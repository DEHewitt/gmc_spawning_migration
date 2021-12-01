library(tidyverse)
library(lubridate)
library(sf)
library(rgdal)
library(raster)
library(fasterize)
library(patchwork)
library(DHARMa)
library(lunar)
library(scales)
library(ggspatial)
library(mgcv)
library(gratia)
library(zoo)

# load custom functions
source("R/salinity_converter.R")
source("R/distance_to_sea.R")
#source("R/harm.R")

# read in the RSP output
blr.rsp <- readRDS("data_processed/BLR/blr_rsp.rds")
clr.rsp <- readRDS("data_processed/CLR/clr_rsp.rds")

# read in shapefiles for checking
blr.shape <- st_as_sf(readOGR("data_processed/BLR/blr.shp"))
clr.shape <- st_as_sf(readOGR("data_processed/CLR/clr_clip.shp"))

# load water quality data
blr.wq <- read_csv("data_raw/blr_wq_6.csv") %>% 
  dplyr::select(date = Date_text, temp = Temp_CAL, cond = Cond_CAL) %>% # just get the datetime and calibrated measurements
  mutate(date = dmy(date)) %>% # convert date to lubridate format
  mutate(estuary = "Kalang River") %>% # add a column to identify which estuary
  salinity_converter() %>% # convert conductivity to salinity
  mutate(mean.temp_1week = rollmean(temp, k = 7*24, align = "right", fill = NA), # calculate 7-day average
         mean.cond_1week = rollmean(cond, k = 7*24, align = "right", fill = NA),
         mean.sal_1week = rollmean(sal, k = 7*24, align = "right", fill = NA),
         mean.temp_3day = rollmean(temp, k = 3*24, align = "right", fill = NA), # calculate 3-day average
         mean.cond_3day = rollmean(cond, k = 3*24, align = "right", fill = NA),
         mean.sal_3day = rollmean(sal, k = 3*24, align = "right", fill = NA)) %>%
  group_by(date) %>% 
  mutate(mean.temp_1day = mean(temp), # calculate mean temp
         mean.cond_1day = mean(cond), # calculate mean condutivity
         mean.sal_1day = mean(sal)) %>% # calculate mean salinity
  ungroup() %>%
  distinct(date, .keep_all = TRUE) %>%
  arrange(date) %>% # put in date order
  mutate(delta.temp_1day = mean.temp_1day - lag(mean.temp_1day, n = 1), # calculate daily difference in temperature
         delta.cond_1day = mean.cond_1day - lag(mean.cond_1day, n = 1), # calculate daily difference in conductivity
         delta.sal_1day = mean.sal_1day - lag(mean.sal_1day, n = 1), # calculate daily difference in salinity
         delta.temp_1week = mean.temp_1day - mean.temp_1week, # difference between daily mean and past week
         delta.cond_1week = mean.cond_1day - mean.cond_1week,
         delta.sal_1week = mean.sal_1day - mean.sal_1week,
         delta.temp_3day = mean.temp_1day - mean.temp_3day,
         delta.cond_3day = mean.cond_1day - mean.cond_3day,
         delta.sal_3day = mean.sal_1day - mean.sal_3day) %>% 
  ungroup()

# note: the delta.vars have NAs for 4 crabs in cohort 1 as they were detected on the same day
# that the hobo loggers were deployed - these instances can be ommitted when modelling or 
# i could just make it < 7 days (but still doesn't fix the first day). not sure of other options atm

clr.wq <- read_csv("data_raw/clr_wq_3.csv") %>% 
  dplyr::select(date = Date_text, temp = Temp_CAL, cond = Cond_CAL) %>% # just get the datetime and calibrated measurements
  mutate(date = dmy(date)) %>% # convert date to lubridate format
  mutate(estuary = "Clarence River") %>% # add a column to identify which estuary
  salinity_converter() %>% # convert conductivity to salinity
  mutate(mean.temp_1week = rollmean(temp, k = 7*24, align = "right", fill = NA), # calculate 7-day average
         mean.cond_1week = rollmean(cond, k = 7*24, align = "right", fill = NA),
         mean.sal_1week = rollmean(sal, k = 7*24, align = "right", fill = NA),
         mean.temp_3day = rollmean(temp, k = 3*24, align = "right", fill = NA), # calculate 3-day average
         mean.cond_3day = rollmean(cond, k = 3*24, align = "right", fill = NA),
         mean.sal_3day = rollmean(sal, k = 3*24, align = "right", fill = NA)) %>%
  group_by(date) %>% 
  mutate(mean.temp_1day = mean(temp), # calculate mean temp
         mean.cond_1day = mean(cond), # calculate mean condutivity
         mean.sal_1day = mean(sal)) %>% # calculate mean salinity
  ungroup() %>%
  distinct(date, .keep_all = TRUE) %>%
  arrange(date) %>% # put in date order
  mutate(delta.temp_1day = mean.temp_1day - lag(mean.temp_1day, n = 1), # calculate daily difference in temperature
         delta.cond_1day = mean.cond_1day - lag(mean.cond_1day, n = 1), # calculate daily difference in conductivity
         delta.sal_1day = mean.sal_1day - lag(mean.sal_1day, n = 1), # calculate daily difference in salinity
         delta.temp_1week = mean.temp_1day - mean.temp_1week, # difference between daily mean and past week
         delta.cond_1week = mean.cond_1day - mean.cond_1week,
         delta.sal_1week = mean.sal_1day - mean.sal_1week,
         delta.temp_3day = mean.temp_1day - mean.temp_3day,
         delta.cond_3day = mean.cond_1day - mean.cond_3day,
         delta.sal_3day = mean.sal_1day - mean.sal_3day) %>%
  ungroup()

# combine into a single dataframe
wq.data <- bind_rows(blr.wq, clr.wq)

# get out the tracks for kalang
blr.tracks <- blr.rsp$detections %>% 
  bind_rows() %>% 
  mutate(Transmitter = as.character(Transmitter)) %>% 
  mutate(estuary = "Kalang River",
         date = date(Timestamp))

# join with water quality
blr.tracks <- blr.tracks %>% left_join(blr.wq)

# get out the tracks for clarence
clr.tracks <- clr.rsp$detections %>% 
  bind_rows() %>% 
  mutate(Transmitter = as.character(Transmitter)) %>% 
  mutate(estuary = "Clarence River",
         date = date(Timestamp))

# join with water quality
clr.tracks <- clr.tracks %>% left_join(clr.wq)

# get size, release date info from rsp output
meta <- data.frame(crab = c(blr.rsp$bio$Transmitter, clr.rsp$bio$Transmitter),
                   length_mm = c(blr.rsp$bio$Length.mm, clr.rsp$bio$Length.mm),
                   tagging_date = c(blr.rsp$bio$Release.date, clr.rsp$bio$Release.date))

# combine them into one dataframe
tracks <- bind_rows(clr.tracks, blr.tracks) %>%
  dplyr::select(crab = Transmitter, time = Timestamp, lon = Longitude,
                lat = Latitude, estuary = estuary, receiver = Array, date) %>%
  mutate(year = year(time))

tracks <- tracks %>% left_join(meta)

# this bit is just so the indexing is consistent
# later on when doing the distance to sea calcs
# not sure why, but it works
blr <- tracks %>% filter(estuary == "Kalang River")
clr <- tracks %>% filter(estuary == "Clarence River")
tracks <- bind_rows(blr, clr)

# don't need to produce these plots anymore
# idea to make this more automated is to:
# 1. figure out which receivers = mouth
# 2. produce a list of crabs that do make it their (some sort of conditional filter)
# 3. use setdiff() on the list of crabs that do show up at the mouth and the full 
# list of tagged crabs, the result of which (i.e., those not detected at the mouth)
# are the ones to remove (using filter())

if (Sys.Date() == "2020-01-01"){
  ggplot() + 
    geom_sf(data = blr.shape) + 
    geom_path(data = blr, aes(x = lon, y = lat), colour = "red", size = 2) + 
    facet_wrap(vars(crab)) +
    theme_classic()
  
  ggplot() + 
    geom_sf(data = clr.shape) +
    geom_path(data = clr, aes(x = lon, y = lat)) + 
    facet_wrap(vars(crab))
}

# remove the crabs that didn't migrate
# based on visual inspection (above)
tracks <- tracks %>%
  filter(crab != "A69-1602-33406") %>%
  filter(crab != "A69-1602-33453") %>%
  filter(crab != "A69-1602-33407") %>%
  filter(crab != "A69-1602-33455") %>%
  filter(crab != "A69-1602-12309") %>%
  filter(crab != "A69-1602-12320") %>%
  filter(crab != "A69-1602-12316") %>%
  filter(crab != "A69-1602-12327") %>%
  filter(crab != "A69-1602-33415") %>%
  filter(crab != "A69-1602-12344") %>%
  filter(crab != "A69-1602-12353") %>%
  filter(crab != "A69-1602-12356") %>%
  filter(crab != "A69-1602-33402") %>%
  filter(crab != "A69-1602-33404") %>%
  filter(crab != "A69-1602-33444") %>%
  filter(crab != "A69-1602-33446") %>%
  filter(crab != "A69-1602-12322") %>%
  filter(crab != "A69-1602-12341") %>%
  filter(crab != "A69-1602-12318") %>%
  filter(crab != "A69-1602-33408") %>%
  filter(crab != "A69-1602-33447") %>%
  filter(crab != "A69-1602-33398") %>%
  filter(crab != "A69-1602-12354") %>%
  filter(crab != "A69-1602-12324") %>%
  filter(crab != "A69-1602-12346") %>%
  filter(crab != "A69-1602-33448")

# add a column that identifies the cohort (called "year" in the ms)
tracks <- tracks %>%
  mutate(cohort = case_when(estuary == "Kalang River" & year(tagging_date) == "2019" ~ 1,
                            estuary == "Kalang River" & month(tagging_date) == "1" ~ 2,
                            estuary == "Kalang River" & month(tagging_date) == "3" ~ 2,
                            estuary == "Clarence River" & year(tagging_date) == "2019" ~ 1,
                            estuary == "Clarence River" & year(tagging_date) == "2020" ~ 2))

# calculate the total number of crabs that migrated
tracks <- tracks %>%
  group_by(cohort, estuary) %>%
  mutate(total.migrators = length(unique(crab))) %>%
  ungroup()

# calculate distance to sea of each point
blr.dist <- distance_to_sea(data = tracks, 
                            shapefile = blr.shape, 
                            estuary = "Kalang River")

clr.dist <- distance_to_sea(data = tracks, 
                            shapefile = clr.shape, 
                            estuary = "Clarence River")

# add the two data frames together
dist <- bind_rows(blr.dist, clr.dist) %>% dplyr::select(-estuary)

# add back to the interpolated points
tracks <- tracks %>% bind_cols(dist)

# calculate some metrics
# interval and speed metrics are probably not need or useful
tracks <- tracks %>%
  dplyr::group_by(crab) %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(dist.to.sea = round(dist.to.sea, 2)) %>% # in kilometres
  dplyr::mutate(distance.metres = (dist.to.sea - lag(dist.to.sea, 1))*1000) %>% # in metres, negative = downstream
  ungroup()

###### example track and distance to sea plots ######
# extract a crab from each estuary for example plots
clr.crab <- tracks %>% filter(crab == "A69-1602-12321") 
blr.crab <- tracks %>% filter(crab == "A69-1602-33397")

blr.track <- ggplot() +
  geom_sf(data = blr.shape) +
  geom_path(data = blr.crab, aes(x = lon, y = lat), size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()) +
  coord_sf(xlim = c(152.96, 153.053), expand = FALSE) +
  annotation_scale() +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering)

# just get daily positions
blr.dist <- blr.crab %>% distinct(date, .keep_all = TRUE)

blr.dist.plot <- ggplot() +
  geom_path(data = blr.dist, aes(x = date, y = dist.to.sea), size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 12),
        axis.title.y = element_blank()) +
  xlab("Date") +
  scale_y_continuous(breaks = seq(0, 7.5, 1.5), limits = c(0, 8)) +
  scale_x_date(date_labels = "%d %b %Y")

blr.track|blr.dist.plot # put plots together

clr.track <- ggplot() +
  geom_sf(data = clr.shape) +
  geom_path(data = clr.crab, aes(x = lon, y = lat), size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()) +
  coord_sf(ylim = c(-29.5, -29.35), xlim = c(153.15, 153.4), expand = FALSE) +
  annotation_scale()

# just get daily positions
clr.dist <- clr.crab %>% distinct(date, .keep_all = TRUE)

clr.dist.plot <- ggplot() +
  geom_path(data = clr.dist, aes(x = date, y = dist.to.sea), size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12)) +
  xlab("Date") +
  ylab("Distance to sea (km)") +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 30)) +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "2 weeks")

clr.track|clr.dist.plot # put plots together

# set up axis labels
ylab <- "Distance to sea (km)"

ylab <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = ylab, angle = 90) +
  coord_cartesian(clip = "off")+
  theme_void()

example.tracks <- (clr.dist.plot|clr.track)/(blr.dist.plot|blr.track)

example.tracks <- (ylab|example.tracks)  + plot_layout(widths = c(.1, 1))

example.tracks <- example.tracks + plot_annotation(tag_levels = list(c("", "a", "b", "c", "d")))

ggsave(paste0("figures/example_tracks.png"), 
       plot = example.tracks, device = "png", 
       width = 17, height = 13, units = "cm", dpi = 600)
############

###### tracks + distance plots for all crabs ######
blr.tracks <- tracks %>% filter(estuary == "Kalang River")
clr.tracks <- tracks %>% filter(estuary == "Clarence River")

blr.all.tracks <- ggplot() +
  geom_sf(data = blr.shape) +
  geom_path(data = blr.tracks,
            aes(x = lon, y = lat),
            size = 1) +
  facet_wrap(vars(crab)) +
  coord_sf(expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(breaks = seq(152.96, 153.04, 0.04)) +
  ylab("Latitude") +
  xlab("Longitude")

ggsave(paste0("figures/blr_all_tracks.png"), 
       plot = blr.all.tracks, device = "png", 
       width = 29, height = 20, units = "cm", dpi = 600)

clr.all.tracks <- ggplot() +
  geom_sf(data = clr.shape) +
  geom_path(data = clr.tracks,
            aes(x = lon, y = lat),
            size = 1) +
  facet_wrap(vars(crab)) +
  coord_sf(expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(breaks = seq(153.0, 153.3, 0.15)) +
  ylab("Latitude") +
  xlab("Longitude")

ggsave(paste0("figures/clr_all_tracks.png"), 
       plot = clr.all.tracks, device = "png", 
       width = 29, height = 20, units = "cm", dpi = 600)

blr.tracks <- blr.tracks %>% 
  group_by(crab) %>% 
  mutate(day = yday(date)) %>% 
  distinct(day, .keep_all = TRUE) %>% 
  ungroup() 

clr.tracks <- clr.tracks %>% 
  group_by(crab) %>% 
  mutate(day = yday(date)) %>% 
  distinct(day, .keep_all = TRUE) %>% 
  ungroup()

blr.all.dist <- ggplot() +
  geom_path(data = blr.tracks, 
            aes(x = date, y = dist.to.sea), size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 12),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Date") +
  ylab("Distance to sea (km)") +
  facet_wrap(vars(crab), scales = "free_x") +
  scale_x_date(breaks = pretty_breaks())

ggsave(paste0("figures/blr_all_dist.png"), 
       plot = blr.all.dist, device = "png", 
       width = 20, height = 29, units = "cm", dpi = 600)

clr.all.dist <- ggplot() +
  geom_path(data = clr.tracks, 
            aes(x = date, y = dist.to.sea), size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 12),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Date") +
  ylab("Distance to sea (km)") +
  facet_wrap(vars(crab), scales = "free_x") +
  scale_x_date(breaks = pretty_breaks())

ggsave(paste0("figures/clr_all_dist.png"), 
       plot = clr.all.dist, device = "png", 
       width = 20, height = 29, units = "cm", dpi = 600)
############

###### summarise tracks by crab and cohort ######
tracks.summary <- tracks %>%
  group_by(crab, estuary, cohort) %>%
  dplyr::summarise(start = min(time),
                   end = max(time),
                   duration = as.numeric(end-start),
                   total.distance.km = sum(distance.metres, na.rm = TRUE)/1000) %>%
  ungroup()

# get measures of speed + n detections for each crab
tracks.speed <- tracks %>%
  filter(!is.na(receiver)) %>%
  group_by(crab, estuary, cohort) %>%
  arrange(time) %>%
  mutate(distance.metres = (dist.to.sea - lag(dist.to.sea, 1))*1000,
         interval = as.numeric(time - lag(time, 1))*60,
         speed.metres.second = abs(distance.metres/interval)) %>%
  summarise(max.speed = max(speed.metres.second, na.rm = TRUE),
            n = n()) %>%
  ungroup()

# join back to main summary table
tracks.summary <- tracks.summary %>% left_join(tracks.speed)

# get the summary for cohorts (for main ms)
cohort.sum <- tracks.summary %>% 
  group_by(estuary, cohort) %>%
  summarise(tracking.mean.duration = mean(duration),
            tracking.sd.duration = sd(duration),
            mean.distance = mean(total.distance.km),
            sd.distance = sd(total.distance.km)) %>%
  ungroup()

migration.dur <- tracks %>%
  filter(distance.metres < 0) %>%
  group_by(crab) %>%
  mutate(start = min(time),
         end = max(time),
         duration = as.numeric(end-start)) %>%
  ungroup() %>%
  group_by(estuary, cohort) %>%
  summarise(mig.mean.duration = mean(duration),
            mig.sd.duration = sd(duration)) %>%
  ungroup()

# put it all together
cohort.sum <- cohort.sum %>% left_join(migration.dur)

# table summarising things for each crab that migrated
migration.summary <- tracks %>%
  filter(!is.na(receiver)) %>% # remove interpolated points so that this is only based on real detections
  group_by(crab) %>%
  arrange(time) %>%
  mutate(distance.metres = (dist.to.sea - lag(dist.to.sea, 1))*1000,
         interval = as.numeric(time - lag(time, 1))*60,
         speed.metres.second = abs(distance.metres/interval)) %>%
  group_by(crab, estuary, cohort) %>%
  summarise(start = min(time),
            end = max(time),
            mean.speed = mean(speed.metres.second, na.rm = TRUE),
            sd.speed = sd(speed.metres.second, na.rm = TRUE),
            max.speed = max(speed.metres.second, na.rm = TRUE)) %>%
  ungroup()

mig.start <- tracks %>%
  group_by(crab) %>%
  filter(distance.metres < 0) %>%
  summarise(mig.start = min(time),
            tagging_date = tagging_date) %>%
  distinct(crab, .keep_all = TRUE)

migration.summary <- migration.summary %>% 
  left_join(mig.start) %>%
  mutate(track.duration = as.numeric(end - start)/24,
         mig.duration = as.numeric(end - mig.start)/24)

write.table(migration.summary, file = "output/migration_summary.txt", sep = ",")
############

###### create dataframe of migrators ######
migration <- tracks %>%
  group_by(crab, date, estuary, cohort, total.migrators, tagging_date) %>%
  summarise(direction = sum(distance.metres, na.rm = TRUE),
            behaviour = case_when(direction < 0 ~ "migrating", direction >= 0 ~ "not migrating")) %>%
  ungroup() %>% 
  mutate(day = yday(date), tagging_day = yday(tagging_date)) %>%
  group_by(crab) %>%
  mutate(final_day = max(day)) %>%
  ungroup() %>%
  mutate(lunar = lunar.phase(date, shift = 10),
         lunar.scale = scale(lunar)) %>% 
  filter(behaviour == "migrating")

ggplot() + geom_histogram(data = migration, aes(x = day), binwidth = 1) + facet_wrap(vars(estuary, cohort), scales = "free")

###### NO LONGER NEEDED ######
if (Sys.Date() == "2020-01-01"){
  ###### finite mixture modelling ######
  estuaries <- migration$estuary %>% unique()
  cohorts <- migration$cohort %>% unique()
  
  set.seed(42)
  
  mixtures <- list()
  for (i in 1:length(estuaries)){
    for (j in 1:length(cohorts)){
      
      # run one per estuary/cohort combination
      data <- migration %>% filter(estuary == estuaries[i] & cohort == cohorts[j])
      
      # gaussian clustering
      mix <- stepFlexmix(day ~ 1, data = data, k = 1:5, nrep = 500, control = list(verbose = 0))
      
      # aic comparison
      aic <- AIC(mix)
      
      # bic comparison
      bic <- BIC(mix)
      
      # so the list is named
      name <- paste(estuaries[i], cohorts[j], sep = "_")
      
      list <- list(model = mix, AIC = aic, BIC = bic)
      
      mixtures[[name]] <- list
    }
  }
  
  # save output
  saveRDS(mixtures, file = "output/mixtures.rds")
  
  # plot BIC selected distributions
  # clr, cohort 1
  par <- parameters(mixtures$`Clarence River_1`$model@models$`2`)
  
  i = 1
  j = 1
  
  data <- migration %>% filter(estuary == estuaries[i] & cohort == cohorts[j])
  
  plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
  }
  
  lam <- table(clusters(mixtures$`Clarence River_1`$model@models$`2`))
  
  clr.1 <- ggplot(data) +
    geom_histogram(aes(day), binwidth = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[1], par[2], lam[1]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[3], par[4], lam[2]),
                  colour = "red", lwd = 1) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
    geom_vline(xintercept = c(min(data$tagging_day), min(data$final_day)),
               linetype = c("dashed", "dotted")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(colour = "black", size = 12))
  
  # clr, cohort 2
  par <- parameters(mixtures$`Clarence River_2`$model@models$`2`)
  
  i = 1
  j = 2
  
  data <- migration %>% filter(estuary == estuaries[i] & cohort == cohorts[j])
  
  lam <- table(clusters(mixtures$`Clarence River_2`$model@models$`2`))
  
  clr.2 <- ggplot(data) +
    geom_histogram(aes(day), binwidth = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[1], par[2], lam[1]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[3], par[4], lam[2]),
                  colour = "red", lwd = 1) +
    geom_vline(xintercept = c(min(data$tagging_day), min(data$final_day)),
               linetype = c("dashed", "dotted")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(colour = "black", size = 12))
  
  # blr, cohort 1
  par <- parameters(mixtures$`Kalang River_1`$model@models$`4`)
  
  i = 2
  j = 1
  
  data <- migration %>% filter(estuary == estuaries[i] & cohort == cohorts[j])
  
  lam <- table(clusters(mixtures$`Kalang River_1`$model@models$`4`))
  
  blr.1 <- ggplot(data) +
    geom_histogram(aes(day), binwidth = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[1], par[2], lam[1]/4),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[3], par[4], lam[2]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[5], par[6], lam[3]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[7], par[8], lam[4]),
                  colour = "red", lwd = 1) +
    geom_vline(xintercept = c(min(data$tagging_day), min(data$final_day)),
               linetype = c("dashed", "dotted")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(colour = "black", size = 12))
  
  # blr, cohort 2
  par <- parameters(mixtures$`Kalang River_2`$model@models$`3`)
  
  i = 2
  j = 2
  
  data <- migration %>% filter(estuary == estuaries[i] & cohort == cohorts[j])
  
  lam <- table(clusters(mixtures$`Kalang River_2`$model@models$`3`))
  
  blr.2 <- ggplot(data) +
    geom_histogram(aes(day), binwidth = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[1], par[2], lam[1]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[3], par[4], lam[2]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(par[5], par[6], lam[3]),
                  colour = "red", lwd = 1) +
    geom_vline(xintercept = c(min(data$tagging_day), min(data$final_day)),
               linetype = c("dashed", "dotted")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(colour = "black", size = 12))
  
  # set up axis labels
  ylab <- "Count"
  
  ylab <- ggplot() + 
    annotate(geom = "text", x = 1, y = 1, label = ylab, angle = 90) +
    coord_cartesian(clip = "off")+
    theme_void()
  
  xlab <- "Day of the year"
  
  xlab <- ggplot() + 
    annotate(geom = "text", x = 1, y = 1, label = xlab) +
    coord_cartesian(clip = "off")+
    theme_void()
  
  # put all the plots together
  mix.plots <- (clr.1|clr.2)/(blr.1|blr.2)
  
  mix.plots <- (ylab|mix.plots|plot_spacer()) + plot_layout(widths = c(.1, 1, .1))
  
  mix.plots <- mix.plots/xlab + plot_layout(heights = c(1, .1)) 
  
  mix.plots <- mix.plots + plot_annotation(tag_levels = list(c("", "a", "b", "c", "d", "", "")))
  
  # save the plot
  ggsave(paste0("figures/mixtures.png"), 
         plot = mix.plots, device = "png", 
         width = 17, height = 13, units = "cm", dpi = 600)
  
  # create a data frame that summarises all the distributions
  a <- as.data.frame(t(parameters(mixtures$`Clarence River_1`$model@models$`2`))) %>% 
    mutate(estuary = "Clarence River", 
           cohort = 1)
  a <- cbind(mix = rownames(a), a)
  
  b <- as.data.frame(t(parameters(mixtures$`Clarence River_2`$model@models$`2`))) %>% 
    mutate(estuary = "Clarence River", 
           cohort = 2)
  b <- cbind(mix = rownames(b), b)
  
  c <- as.data.frame(t(parameters(mixtures$`Kalang River_1`$model@models$`4`))) %>% 
    mutate(estuary = "Kalang River", 
           cohort = 1)
  c <- cbind(mix = rownames(c), c)
  
  d <- as.data.frame(t(parameters(mixtures$`Kalang River_2`$model@models$`3`))) %>% 
    mutate(estuary = "Kalang River", 
           cohort = 2)
  d <- cbind(mix = rownames(d), d)
  
  mix.pars <- bind_rows(a, b, c, d)
  rownames(mix.pars) <- NULL
  
  saveRDS(mix.pars, "output/mixture_parameters.rds")
  ############
}
##############################

############
migration.binom <- migration %>%
  group_by(crab) %>%
  filter(direction < 0) %>%
  filter(day == min(day)) %>% 
  ungroup() %>%
  dplyr::select(crab, date, day, estuary, cohort, total.migrators, tagging_date, tagging_day, final_day, lunar) %>%
  mutate(tagging_date = date(tagging_date))

# add a column that is a counter for how many crabs have begun migrating - cumsum()
# use this to calculate proportions
# generate zeroes on days when no crabs started
# check that every day in the tracking period is represented
# there is a possibility that days where 
x <- migration.binom %>% filter(estuary == "Clarence River" & cohort == 1)
exists <- x$date %>% unique() # dates when crabs started migrating
should.exist <- seq.Date(from = min(x$tagging_date), to = max(x$date), by = 1)  
clr1.diff <- data.frame(date = as.Date(dplyr::setdiff(should.exist, exists), origin = "1970-01-01"),
                        estuary = "Clarence River", 
                        cohort = 1,
                        mig.n = as.numeric(0))

x <- migration.binom %>% filter(estuary == "Clarence River" & cohort == 2)
exists <- x$date %>% unique()
should.exist <- seq.Date(from = min(x$tagging_date), to = max(x$date), by = 1)
clr2.diff <- data.frame(date = as.Date(dplyr::setdiff(should.exist, exists), origin = "1970-01-01"),
                        estuary = "Clarence River", 
                        cohort = 2,
                        mig.n = as.numeric(0))

x <- migration.binom %>% filter(estuary == "Kalang River" & cohort == 1)
exists <- x$date %>% unique()
should.exist <- seq.Date(from = min(x$tagging_date), to = max(x$date), by = 1)
blr1.diff <- data.frame(date = as.Date(dplyr::setdiff(should.exist, exists), origin = "1970-01-01"),
                        estuary = "Kalang River", 
                        cohort = 1,
                        mig.n = as.numeric(0))

x <- migration.binom %>% filter(estuary == "Kalang River" & cohort == 2)
exists <- x$date %>% unique()
should.exist <- seq.Date(from = min(x$tagging_date), to = max(x$date), by = 1)
blr2.diff <- data.frame(date = as.Date(dplyr::setdiff(should.exist, exists), origin = "1970-01-01"),
                        estuary = "Kalang River", 
                        cohort = 2,
                        mig.n = as.numeric(0))

migration.binom <- migration.binom %>%
  add_column(mig.n = 1) %>%
  bind_rows(clr1.diff, clr2.diff, blr1.diff, blr2.diff) %>%
  left_join(wq.data) %>%
  group_by(estuary, cohort) %>%
  mutate(total.migrators = if_else(is.na(total.migrators), na.locf(total.migrators), total.migrators),
         tagging_date = if_else(is.na(tagging_date), na.locf(tagging_date), tagging_date),
         tagging_day = if_else(is.na(tagging_day), na.locf(tagging_day), tagging_day),
         lunar = lunar.phase(date, shift = 10)) %>%
  ungroup() %>%
  group_by(date, estuary, cohort) %>%
  arrange(date) %>%
  mutate(mig.n = sum(mig.n)) %>%
  ungroup() %>%
  group_by(estuary, cohort) %>%
  distinct(date, .keep_all = TRUE) %>%
  mutate(day = yday(date),
         mig.total = cumsum(mig.n),
         mig.N = if_else(mig.total == 0, total.migrators, as.integer(total.migrators-lag(mig.total))),
         mig.prop = mig.n/mig.N,
         lunar = lunar.phase(date, shift = 10),
         harm.lunar = 2*pi*(lunar/max(lunar, na.rm = TRUE)),
         estuary = factor(estuary),
         cohort = factor(cohort),
         delta.temp_1week = if_else(is.na(delta.temp_1week), delta.temp_1day, delta.temp_1week),
         delta.cond_1week = if_else(is.na(delta.cond_1week), delta.cond_1day, delta.cond_1week),
         delta.temp_3day = if_else(is.na(delta.temp_3day), delta.temp_1day, delta.temp_3day),
         delta.cond_3day = if_else(is.na(delta.cond_3day), delta.cond_1day, delta.cond_3day),
         scale.delta.temp_1week = scale(delta.temp_1week, center = TRUE, scale = TRUE),
         scale.delta.temp_1day = scale(delta.temp_1day, center = TRUE, scale = TRUE),
         scale.mean.temp_1day = scale(mean.temp_1day, center = TRUE, scale = TRUE),
         scale.delta.cond_1week = scale(delta.cond_1week, center = TRUE, scale = TRUE),
         scale.delta.cond_1day = scale(delta.cond_1day, center = TRUE, scale = TRUE),
         scale.mean.cond_1day = scale(mean.cond_1day, center = TRUE, scale = TRUE),
         scale.delta.temp_3day = scale(delta.temp_3day, center = TRUE, scale = TRUE),
         scale.delta.cond_3day = scale(delta.cond_3day, center = TRUE, scale = TRUE)) %>%
  filter(!is.na(delta.temp_1day)) %>%
  ungroup()

###### cumulative migration curve ######
ggplot() + 
  geom_path(data = migration.binom, 
            aes(x = date, y = mig.total/total.migrators)) + 
  geom_path(data = migration.binom,
            aes(x = date, y = mean.temp_1day/max(mean.temp_1day)),
            colour = "red") +
  geom_path(data = migration.binom,
            aes(x = date, y = mean.cond_1day/max(mean.cond_1day)),
            colour = "blue") +
  facet_wrap(vars(estuary, cohort), scales = "free") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(migration.binom$mean.temp_1day)))

a <- migration.binom %>% filter(estuary == "Clarence River" & cohort == 1)
b <- migration.binom %>% filter(estuary == "Clarence River" & cohort == 2)
c <- migration.binom %>% filter(estuary == "Kalang River" & cohort == 1)
d <- migration.binom %>% filter(estuary == "Kalang River" & cohort == 2)

ylab1 <- "Proportion migrated"
ylab2.1 <- bquote(paste("Temperature (", degree, C, ") ")) 
ylab2.2 <- bquote(paste("and"))
ylab2.3 <- bquote(paste("       conductivity (mS cm"^-1, ")"))
xlab <- "Date"

ylab1 <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = ylab1, angle = 90) +
  coord_cartesian(clip = "off")+
  theme_void()

ylab2 <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = ylab2.1, angle = 270, colour = "red", hjust = 1) +
  annotate(geom = "text", x = 1, y = 1, label = ylab2.2, angle = 270, colour = "black", hjust = 0) +
  annotate(geom = "text", x = 1, y = 1, label = ylab2.3, angle = 270, colour = "blue", hjust = 0) +
  coord_cartesian(clip = "off") +
  theme_void()

xlab <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = xlab, angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

a1 <- ggplot() + 
  geom_path(data = a, 
            aes(x = date, y = mig.total/total.migrators),
            size = 1) + 
  geom_path(data = a,
            aes(x = date, y = mean.temp_1day/max(mean.temp_1day)),
            colour = "red",
            size = 1) +
  geom_path(data = a,
            aes(x = date, y = mean.cond_1day/max(mean.cond_1day)),
            colour = "blue",
            size = 1) +
  geom_vline(xintercept = min(a$tagging_date), linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(migration.binom$mean.temp_1day))) +
  scale_x_date(date_breaks = "month", date_labels = "%d %b") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.y.right = element_blank())

b1 <- ggplot() + 
  geom_path(data = b, 
            aes(x = date, y = mig.total/total.migrators),
            size = 1) + 
  geom_path(data = b,
            aes(x = date, y = mean.temp_1day/max(mean.temp_1day)),
            colour = "red",
            size = 1) +
  geom_path(data = b,
            aes(x = date, y = mean.cond_1day/max(mean.cond_1day)),
            colour = "blue",
            size = 1) +
  geom_vline(xintercept = min(b$tagging_date), linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(migration.binom$mean.temp_1day))) +
  scale_x_date(date_labels = "%d %b") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.y.left = element_blank())

c1 <- ggplot() + 
  geom_path(data = c, 
            aes(x = date, y = mig.total/total.migrators),
            size = 1) + 
  geom_path(data = c,
            aes(x = date, y = mean.temp_1day/max(mean.temp_1day)),
            colour = "red",
            size = 1) +
  geom_path(data = c,
            aes(x = date, y = mean.cond_1day/max(mean.cond_1day)),
            colour = "blue",
            size = 1) +
  geom_vline(xintercept = min(c$tagging_date)+1, linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(migration.binom$mean.temp_1day))) +
  scale_x_date(date_breaks = "month", date_labels = "%d %b") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.y.right = element_blank())

d1 <- ggplot() + 
  geom_path(data = d, 
            aes(x = date, y = mig.total/total.migrators),
            size = 1) + 
  geom_path(data = d,
            aes(x = date, y = mean.temp_1day/max(mean.temp_1day)),
            colour = "red",
            size = 1) +
  geom_path(data = d,
            aes(x = date, y = mean.cond_1day/max(mean.cond_1day)),
            colour = "blue",
            size = 1) +
  geom_vline(xintercept = min(d$tagging_date), linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(migration.binom$mean.temp_1day))) +
  scale_x_date(date_breaks = "month", date_labels = "%d %b") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_blank(),
        axis.text.y.left = element_blank())

mig.cum <- (a1|b1)/(c1|d1)

mig.cum <- (ylab1|mig.cum|ylab2) + plot_layout(widths = c(.1, 1, .1))

mig.cum <- mig.cum/xlab + plot_layout(heights = c(1, .1)) 

mig.cum <- mig.cum + plot_annotation(tag_levels = list(c("", "a", "b", "c", "d", "", "")))

ggsave(paste0("figures/cum_migration.png"), 
       plot = mig.cum, device = "png", 
       width = 17, height = 13, units = "cm", dpi = 600)
############

###### gam modelling ######
# hgam for effects of environmental variation
m1 <- gam(mig.prop ~ 
            s(mean.cond_1day, bs = "tp") +
            s(mean.cond_1day, cohort, bs = "fs") +
            s(delta.cond_3day, bs = "tp") + 
            s(delta.cond_3day, cohort, bs = "fs") +
            s(mean.temp_1day, bs = "tp") +
            s(mean.temp_1day, cohort, bs = "fs") +
            s(delta.temp_3day, bs = "tp") + 
            s(delta.temp_3day, cohort, bs = "fs") +
            s(lunar, bs = "cc") + 
            s(lunar, cohort, bs = "fs"),
          knots = list(lunar = c(0.5, 6.5)),
          data = migration.binom,
          family = binomial,
          weights = mig.N,
          method = "REML", 
          select = TRUE)

m2 <- gam(mig.prop ~ 
            s(mean.cond_1day, bs = "tp") +
            s(mean.cond_1day, cohort, bs = "fs") +
            s(delta.cond_1day, bs = "tp") + 
            s(delta.cond_1day, cohort, bs = "fs") +
            s(mean.temp_1day, bs = "tp") +
            s(mean.temp_1day, cohort, bs = "fs") +
            s(delta.temp_1day, bs = "tp") + 
            s(delta.temp_1day, cohort, bs = "fs") +
            s(lunar, bs = "cc") + 
            s(lunar, cohort, bs = "fs"),
          knots = list(lunar = c(0.5, 6.5)),
          data = migration.binom,
          family = binomial,
          weights = mig.N,
          method = "REML", 
          select = TRUE)

m3 <- gam(mig.prop ~ 
            s(mean.cond_1day, bs = "tp") +
            s(mean.cond_1day, cohort, bs = "fs") +
            s(delta.cond_1week, bs = "tp") + 
            s(delta.cond_1week, cohort, bs = "fs") +
            s(mean.temp_1day, bs = "tp") +
            s(mean.temp_1day, cohort, bs = "fs") +
            s(delta.temp_1week, bs = "tp") + 
            s(delta.temp_1week, cohort, bs = "fs") +
            s(lunar, bs = "cc") + 
            s(lunar, cohort, bs = "fs"),
          knots = list(lunar = c(0.5, 6.5)),
          data = migration.binom,
          family = binomial,
          weights = mig.N,
          method = "REML", 
          select = TRUE)

###### model selection ######
# which model is best at explaining the variation in the data
AIC(m1, m2, m3)
############

###### model summary ######
summary(m1)
summary(m2)
summary(m3)

draw(m1, fun = inv_link(m1))
draw(m2, fun = inv_link(m2))
draw(m3, fun = inv_link(m3))

###### model diagnostics ######
# check basis dimension (k) choice is good
gam.check(m1, rep = 1000)

# check model assumptions
gam.diag <- appraise(m1, method = "simulate") 
gam.diag 

m1.res <- simulateResiduals(m1)
plot(m1.res)

#qq plot
gam.qq <- gam.diag[[1]] +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_blank())
  
# resid v. pred
gam.res.v.pred <- gam.diag[[2]] +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_blank())

# resid hist
gam.res.hist <- gam.diag[[3]] +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_blank())

# obs v fit
gam.obs.v.fit <- gam.diag[[4]] +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        title = element_blank())

gam.diagnostics <- (gam.qq|gam.res.v.pred)/(gam.res.hist|gam.obs.v.fit)

ggsave(paste0("figures/gamm_diagnostics.plot.png"), 
       plot = gam.diagnostics, device = "png", 
       width = 17, height = 17, units = "cm", dpi = 600)
############

###### plot model results ######
# look at summary
summary(m1)

# effects
draw(m1, 
     rug = TRUE, 
     fun = inv_link(m1), 
     ci_level = 0.95,
     unconditional = TRUE,
     overall_uncertainty = TRUE, 
     scales = "fixed")

draw(m1, 
     rug = TRUE, 
     fun = inv_link(m1), 
     ci_level = 0.95,
     unconditional = TRUE,
     overall_uncertainty = TRUE, 
     scales = "fixed")

visreg(m1, type = "contrast", scale = "response")

# get the inverse of the logit link function
inverse.logit <- inv_link(m1)

temp.smooth <- bind_rows(evaluate_smooth(m1, smooth = 's(mean.temp_1day)')) %>% 
  mutate(lower = est - (2 * se), 
         upper = est + (2 * se),
         #cohort = if_else(is.na(cohort), "global", as.character(cohort)),
         est.response = inverse.logit(est),
         upper.response = inverse.logit(upper),
         lower.response = inverse.logit(lower))

temp <- ggplot() +
  geom_line(data = temp.smooth, 
            aes(x = mean.temp_1day, 
                y = est),
            size = 1) +
  geom_line(data = temp.smooth, 
            aes(x = mean.temp_1day, 
                y = upper),
            size = 1,
            linetype = "dashed") +
  geom_line(data = temp.smooth, 
            aes(x = mean.temp_1day, 
                y = lower),
            size = 1, 
            linetype = "dashed") +
  geom_rug(data = migration.binom,
           aes(x = mean.temp_1day)) +
  #geom_point(data = res, aes(x = mean.temp_1day, y = `s(mean.temp_1day)`, colour = cohort), alpha = 0.3) +
  theme_bw() +
  #geom_hline(yintercept = max(temp.smooth$lower)) +
  #facet_wrap(vars(cohort), nrow = 1) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  scale_colour_manual(values = c("red", "blue", "black"))  +
  scale_y_continuous(limits = c(-6, 6.5)) +
  xlab(bquote(paste("Mean temperature (", degree, C, ")"))) +
  ylab("Log-odds(migrating)")
  
cond.smooth <- bind_rows(evaluate_smooth(m1, smooth = 's(delta.cond_3day)')) %>% 
  mutate(lower = est - (2 * se), 
         upper = est + (2 * se),
         #cohort = if_else(is.na(cohort), "global", as.character(cohort)),
         est.response = inverse.logit(est),
         upper.response = inverse.logit(upper),
         lower.response = inverse.logit(lower))

cond <- ggplot() +
  geom_line(data = cond.smooth, 
            aes(x = delta.cond_3day, 
                y = est),
            size = 1) +
  geom_line(data = cond.smooth, 
            aes(x = delta.cond_3day, 
                y = upper),
            size = 1, 
            linetype = "dashed") +
  geom_line(data = cond.smooth, 
            aes(x = delta.cond_3day, 
                y = lower),
            size = 1, 
            linetype = "dashed") +
  geom_rug(data = migration.binom,
           aes(x = delta.cond_3day)) +
  #geom_point(data = res, aes(x = mean.cond_3day, y = `s(mean.cond_3day)`, colour = cohort), alpha = 0.3) +
  #geom_hline(yintercept = max(cond.smooth$lower)) +
  theme_bw() +
  #facet_wrap(vars(cohort), nrow = 1) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  scale_colour_manual(values = c("red", "blue", "black"))  +
  scale_y_continuous(limits = c(-6, 6.5)) +
  xlab(bquote(paste(Delta, "-conductivity (mS cm"^-1, ")")))

lunar.smooth <- bind_rows(evaluate_smooth(m1, smooth = 's(lunar)')) %>% 
  mutate(lower = est - (2 * se), 
         upper = est + (2 * se),
         #cohort = if_else(is.na(cohort), "global", as.character(cohort)),
         est.response = inverse.logit(est),
         upper.response = inverse.logit(upper),
         lower.response = inverse.logit(lower))

lunar <- ggplot() +
  geom_line(data = lunar.smooth, 
            aes(x = lunar, 
                y = est),
            size = 1) +
  geom_line(data = lunar.smooth, 
            aes(x = lunar, 
                y = upper),
            size = 1, 
            linetype = "dashed") +
  geom_line(data = lunar.smooth, 
            aes(x = lunar, 
                y = lower),
            size = 1, 
            linetype = "dashed") +
  geom_rug(data = migration.binom,
           aes(x = lunar)) +
  #geom_point(data = res, aes(x = lunar, y = `s(lunar,cohort)`, colour = cohort), alpha = 0.3) +
  theme_bw() +
  #geom_hline(yintercept = 0.550860296, colour = "red") +
  #geom_hline(yintercept = -0.715027772, colour = "blue") +
  #facet_wrap(vars(cohort), nrow = 1) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.825, 0.825),
        #axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  scale_x_continuous(breaks = c(0, pi/2, pi, 3*pi/2), 
                     labels = c("0", 
                                expression(paste(pi,"/2")), 
                                expression(paste(pi)),
                                expression(paste("3", pi, "/2")))) +
  scale_y_continuous(limits = c(-6, 6.5)) +
  scale_colour_manual(values = c("red", "blue", "black"), name = "Cohort", labels = c("1 (2019)", "2 (2020)", "Global")) +
  xlab("Lunar phase (rad)") +
  ylab("Log-odds(migrating)") 

effect.plot <- (temp|cond)/(lunar|plot_spacer())

effect.plot <- effect.plot + plot_annotation(tag_levels = "a")

effect.plot

ggsave(paste0("figures/gamm_effects.plot.png"), 
       plot = effect.plot, device = "png", 
       width = 17, height = 17, units = "cm", dpi = 600)

temp <- ggplot() + 
  geom_ribbon(data = temp.smooth, 
              aes(x = mean.temp_1day, 
                  ymin = lower.response,
                  ymax = upper.response),
              fill = "grey") +
  geom_line(data = temp.smooth, 
            aes(x = mean.temp_1day, 
                y = est.response),
            size = 1) +
  geom_rug(data = migration.binom,
           aes(x = mean.temp_1day)) +
  theme_bw() +
  #geom_hline(yintercept = max(temp.smooth$lower)) +
  #facet_wrap(vars(cohort), nrow = 1) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12),
        axis.title.y = element_blank()) +
  scale_colour_manual(values = c("red", "blue", "black"))  +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(20, 30, 2)) +
  xlab(bquote(paste("Mean temperature (", degree, C, ")")))# +
  #ylab("Pr(migrating)")

cond <- ggplot() +
  geom_ribbon(data = cond.smooth, 
              aes(x = delta.cond_3day, 
                  ymin = lower.response,
                  ymax = upper.response),
              fill = "grey") +
  geom_line(data = cond.smooth, 
            aes(x = delta.cond_3day, 
                y = est.response),
            size = 1) +
  geom_rug(data = migration.binom,
           aes(x = delta.cond_3day)) +
  #geom_point(data = res, aes(x = delta.cond_1day, y = `s(delta.cond_1day)`), alpha = 0.3) +
  #geom_hline(yintercept = max(cond.smooth$lower)) +
  theme_bw() +
  #facet_wrap(vars(cohort), nrow = 1) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  scale_colour_manual(values = c("red", "blue", "black"))  +
  scale_y_continuous(limits = c(0, 1)) +
  xlab(bquote(paste(Delta, "-conductivity (mS cm"^-1, ")")))

lunar <- ggplot() +
  geom_ribbon(data = lunar.smooth, 
              aes(x = lunar, 
                  ymin = lower.response,
                  ymax = upper.response),
              fill = "grey") +
  geom_line(data = lunar.smooth, 
            aes(x = lunar, 
                y = est.response),
            size = 1) +
  geom_rug(data = migration.binom,
           aes(x = lunar)) +
  #geom_point(data = res, aes(x = lunar, y = `s(lunar,cohort)`, colour = cohort), alpha = 0.3) +
  theme_bw() +
  #geom_hline(yintercept = 0.550860296, colour = "red") +
  #geom_hline(yintercept = -0.715027772, colour = "blue") +
  #facet_wrap(vars(cohort), nrow = 1) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.825, 0.825),
        axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  scale_x_continuous(breaks = c(0, pi/2, pi, 3*pi/2), 
                     labels = c("0", 
                                expression(paste(pi,"/2")), 
                                expression(paste(pi)),
                                expression(paste("3", pi, "/2")))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c("red", "blue", "black"), name = "Cohort", labels = c("1 (2019)", "2 (2020)", "Global")) +
  xlab("Lunar phase (rad)") #+
  #ylab("Pr(migrating)") 

effect.plot <- (temp|cond)/(lunar|plot_spacer())

effect.plot <- effect.plot + plot_annotation(tag_levels = "a")

ylab1 <- "Probability of migrating"
                                             
ylab1 <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = ylab1, angle = 90, size = 4) +
  coord_cartesian(clip = "off")+
  theme_void()

effect.plot <- (ylab1|effect.plot) + plot_layout(widths = c(0.1, 1))

effect.plot <- effect.plot + plot_annotation(tag_levels = list(c("", "a", "b", "c")))

ggsave(paste0("figures/gamm_effects_response_plot.png"), 
       plot = effect.plot, device = "png", 
       width = 17, height = 17, units = "cm", dpi = 600)
