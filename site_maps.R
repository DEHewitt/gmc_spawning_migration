library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(ggspatial)
library(patchwork)
library(cowplot)
library(gridExtra)
library(sf)
library(rgdal)

# read in shapefiles for checking
oz <- ne_states(country = "australia", returnclass = "sf")
blr.shape <- st_as_sf(readOGR("data_processed/BLR/blr.shp"))
clr.shape <- st_as_sf(readOGR("data_processed/CLR/clr_clip.shp"))
blr.stations <- read_csv("data_raw/receiver_locations.csv") %>% filter(estuary == "BLR")
clr.stations <- read_csv("data_raw/receiver_locations.csv") %>% filter(estuary == "CLR")
estuaries <- data.frame(estuary = c("clr", "blr"), lat = c(-29.427133, -30.500000), lon = c(153.371989, 153.030000))

blr <- ggplot() +
  geom_sf(data = blr.shape) +
  geom_point(data = blr.stations,
             aes(x = long,
                 y = lat,
                 shape = final.rsp),
             size = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  coord_sf(expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_scale()

clr <- ggplot() +
  geom_rect(aes(xmin = 152.9, xmax = 153.25, ymin = -29.7, ymax = -29.458), fill = "red", alpha = 0.5) +
  geom_rect(aes(xmin = 152.9, xmax = 153.18, ymin = -29.458, ymax = -29.35), fill = "red", alpha = 0.5) +
  geom_sf(data = clr.shape) +
  geom_point(data = clr.stations,
             aes(x = long,
                 y = lat,
                 shape = final.rsp),
             size = 2) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  coord_sf(xlim = c(min(clr.stations$long), max(clr.stations$long)), 
           ylim = c(min(clr.stations$lat), max(clr.stations$lat))) +
  scale_x_continuous(breaks = seq(153, 153.35, by = 0.1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_scale()

inset <- ggplot() +
  geom_sf(data = oz) +
  geom_point(data = estuaries,
             aes(x = lon,
                 y = lat)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  coord_sf(xlim = c(145, 154), ylim = c(-37.5, -15),
           label_axes = "-NE-") +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  scale_y_continuous(position = "right")

site.map <- (clr/blr)|inset

site.map <- site.map + plot_annotation(tag_levels = "a")

ggsave("figures/site_map.png", 
       plot = site.map, 
       device = "png", 
       width = 29, # a4 dimensions
       height = 20, 
       units = "cm", 
       dpi = 600)
