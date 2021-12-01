# Script for analysis of Giant Mud Crab telemetry data in Kalang River
# Daniel Hewitt, UNSW

library(RSP)
library(actel)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(sf)
library(rgdal)

# Things to trial:
# 1. Interpolate more points (i.e. lower distance) will need to increase window.size and margin (I think... proportional to the number of points)
# 2. Make time.step lower so that quick moves between K4 and K5 get some interpolation
# 2. Remove multiple detections at intermediate receivers. This artificially slows the crab down

#################### Set up study site and checks ####################
# Study area details
tz <- "Australia/Sydney" # timezone
UTM <- 56 # UTM zone

# Crab speed check
# error prompting intervention if they go faster than this. Value from Alberts-Hubatsch thesis
speed.error <- 0.56 # cm/s
# How many arrays (receivers) will we let a crab 'jump' (i.e. not be detected on)? # let's go with 1 to be conservative
#jump.error <- 1 
# What about inactivity? We know that crabs generally just sit around...
# 3 days as an arbitray choice
#inactive.error <- 3 

################## Load in the spatial data #################
if (Sys.info()[6] == "Dan"){
  setwd("C:/Users/Dan/Documents/PhD/Spawning_migration/data_processed/BLR")
}

blr.shape <- st_as_sf(readOGR("blr.shp"))

# Load the shapefile
if(file.exists("blr_base_raster.rds")){
  blr.base.raster <- readRDS("blr_base_raster.rds")
} else {
  blr.base.raster <- loadShape(shape = "blr.shp", # estuary
                               size = 0.0001, # pixel size (decimal degrees) ~ 11.1 m 
                               spatial = "spatial.csv", # the spatial file with info about receivers/releases
                               coord.x = "Longitude", # names of x/y columns in 'spatial.csv'
                               coord.y = "Latitude",
                               type = "land")#,
                              # buffer = 0.05) # 20 seconds
  
  saveRDS(blr.base.raster, file = "blr_base_raster.rds")
}

# create a raster for katana
#blr.base.raster.k <- loadShape(shape = "blr.shp", # estuary
 #                            size = 0.0001, # pixel size (decimal degrees) ~ 11.1 m 
  #                           spatial = "spatial.csv", # the spatial file with info about receivers/releases
   #                          coord.x = "Longitude", # names of x/y columns in 'spatial.csv'
    #                         coord.y = "Latitude",
     #                        type = "land",
      #                       buffer = 0.05) # 20 seconds

#saveRDS(blr.base.raster.k, file = "blr_base_raster_katana.rds")

# Convert it to a transition layer
# 5 minutes
if(file.exists("blr_transition_layer.rds")){
  blr.transition.layer <- readRDS("blr_transition_layer.rds")
} else {
  blr.transition.layer <- transitionLayer(blr.base.raster, # object from above
                                          directions = 16) # how many directions can an animal move (16 = max)
  
  saveRDS(blr.transition.layer, file = "blr_transition_layer.rds")
}


#################### Filtering with actel ####################
# Create the distances matrix
blr.dist.mat <- distancesMatrix(blr.transition.layer, 
                                coord.x = "Longitude", 
                                coord.y = "Latitude")

y # Do you want to save the distances matrix?

y # Overwrite if it already exists?

# Now we can perform some preliminary checks
blr <- explore(tz = tz,
               speed.error = speed.error,
               save.tables.locally = TRUE,
               report = TRUE) # create a report to send around
  
y # save stray tags?

# event (in)validation
comment
released close to receiver
n

comment
released close to receiver
n

comment
released close to receiver
n

y  # save actel_explore_results.RData? 

#y # Would you like to save a copy of the analysis log?(y/n)

# Check that everything looks right
# e.g. recievers in the water, good pixel size (to resolve channels etc)
plotRaster(input = blr, # the results from actel::explore()
           base.raster = blr.base.raster, # the estuary shapefile
           coord.x = "Longitude",
           coord.y = "Latitude")

##################### Calculate shortest paths #####################
# Open questions for RSP:
# What distance between locations do we want to interpolate?
blr.rsp <- runRSP(input = blr, # results from actel::explore()
                  t.layer = blr.transition.layer, # transition layer
                  coord.x = "Longitude", 
                  coord.y = "Latitude",
                  recaptures = TRUE,
                  time.step = 10, # interval before interpolation starts
                  distance = 100, # How far between locations?
                  min.time = 10, # min. time (default = 10 min) between detections before RSP calc'd (i.e. detections seperated by 5 mins have no interpolated points between them)
                  max.time = 10000, # threshold for splitting into tracks - artificially high so they don't get split
                  er.ad = 0.0001) # super small so that RSP doesn't crash

########## Save the output ##########
file <- "blr_rsp.rds"

if (file.exists(file)){
  file.remove(file)
}

# Save the output
saveRDS(blr.rsp,
     file = file)
#####################################

####### get out tracks for plotting #######
tracks <- blr.rsp$detections %>% bind_rows()

plot_track <- function(data, crab){
  data %>%
    filter(Transmitter == crab) %>%
    ggplot(data = .,) +
    geom_sf(data = blr.shape) +
    geom_point(aes(x = Longitude,
                   y = Latitude,
                   colour = Position,
                   shape = Position)) +
    theme_bw() +
    theme(axis.title = element_text(size = 12, colour = "black"),
          axis.text = element_text(size = 10, colour = "black"),
          panel.grid = element_blank()) +
    coord_sf(xlim = c(min(data$Longitude), max(data$Longitude)), 
             ylim = c(min(data$Latitude), max(data$Latitude)),
             expand = T) +
    ggtitle(crab)
}

plot_track(data = tracks, crab = "A69-1602-12343")

save_plot <- function(data, crab, path = "figures/") {
  dir.create(path, FALSE, TRUE)
  plot_track(data, crab)
  ggsave(paste0(path, crab, ".png"))
}

crabs <- unique(tracks$Transmitter)

for (crab in crabs) {
  save_plot(tracks, crab)
}
######################################################



if(Sys.info()[6] == "Dan"){
  setwd("C:/Users/Dan/Documents/PhD/Spawning_migration/data_processed/CLR")
} else {
  setwd("/home/z5278054/Spawning_migration/CLR")
}

clr.shape <- st_as_sf(readOGR("clr_clip.shp"))

# Create distances matrix
# Load the shapefile
if(file.exists("clr_base_raster.rds")){
  clr.base.raster <- readRDS("clr_base_raster.rds")
} else {
  clr.base.raster <- loadShape(shape = "clr2.shp", # estuary
                               size = 0.0005, # pixel size - this takes about 10 minutes (units - decimal degrees, 0.0001 ~ 11 m)
                               spatial = "spatial.csv", # the spatial file with info about receivers/releases
                               coord.x = "Longitude", # names of x/y columns in 'spatial.csv'
                               coord.y = "Latitude",
                               type = "land")#,
  #buffer = 0.05)
  
  saveRDS(clr.base.raster, file = "clr_base_raster.rds")
}

# create a raster for katana
#clr.base.raster.k <- loadShape(shape = "clr2.shp", # estuary
 #                              size = 0.0001, # pixel size (decimal degrees) ~ 11.1 m 
  #                             spatial = "spatial.csv", # the spatial file with info about receivers/releases
   #                            coord.x = "Longitude", # names of x/y columns in 'spatial.csv'
    #                           coord.y = "Latitude",
     #                          type = "land",
      #                         buffer = 0.05) # 20 seconds

#saveRDS(clr.base.raster.k, file = "clr_base_raster_katana.rds")


# Convert it to a transition layer
if(file.exists("clr_transition_layer.rds")){
  clr.transition.layer <- readRDS("clr_transition_layer.rds")
} else {
  clr.transition.layer <- transitionLayer(clr.base.raster, directions = 16) 
  
  saveRDS(clr.transition.layer, file = "clr_transition_layer.rds")
}

# Create the distances matrix
clr.dist.mat <- distancesMatrix(clr.transition.layer, 
                                coord.x = "Longitude", 
                                coord.y = "Latitude")

y # do you want to save the distances matrix?

y # do you want to overwrite the distances file if it exists?

# When was the first tag deployed in CLR
clr.start.time <- "2019-02-26 15:37:00"

# Now we can perform some preliminary checks
# Run the actel explore() function - when these change from .warning to .error interaction will be required
clr <- explore(tz = tz,
               speed.error = speed.error, # creating errors means I will need to intervene
               #jump.error = jump.error,
               #inactive.error = inactive.error,
               save.tables.locally = T,
               report = T) # create a report to send around

b # discard orphans
b 
b
b
b

y # save stray_tags.csv

# event (in)validation
comment
only just over the speed limit
n

comment
only just over the speed limit
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
only just over (12 -> 10), and getting fast for (4 -> 1) but the water moves quick here
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
only just over (6 -> 4 and 4 -> 1), pretty fast from 11 -> 10 but within the realm of possibility
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

comment
moving fast but believable - the water through here moves quick
n

y # Would you like to save a copy of the results to actel_explore_results.RData?(y/n)

y # Would you like to save a copy of the analysis log to 2021-02-10.14.08.18.actel.log.txt?(y/n)

# Check that everything looks right
# e.g. recievers in the water, good pixel size (to resolve channels etc)
plotRaster(input = clr, # the results from actel::explore()
           base.raster = clr.base.raster, # the estuary shapefile
           coord.x = "Longitude",
           coord.y = "Latitude")

# Calculate shortest paths
clr.rsp <- runRSP(input = clr, # results from actel::explore()
                  t.layer = clr.transition.layer, # transition layer
                  coord.x = "Longitude", 
                  coord.y = "Latitude", 
                  recaptures = TRUE,
                  time.step = 10, # interval before interpolation starts
                  distance = 100, # How far between locations?
                  min.time = 10, # min. time (default = 10 min) between detections before RSP calc'd (i.e. detections seperated by 5 mins have no interpolated points between them)
                  max.time = 10000, # threshold for splitting into tracks - artificially high so they don't get split
                  er.ad = 0.0001) # super small so that RSP doesn't crash

########## Save the output ##########
file <- "clr_rsp.rds"

if (file.exists(file)){
  file.remove(file)
}

# Save the output
saveRDS(clr.rsp,
        file = file)
#####################################

####### get out tracks for plotting #######
tracks <- clr.rsp$detections %>% bind_rows()

plot_track <- function(data, crab){
  data %>%
    filter(Transmitter == crab) %>%
    ggplot(data = .,) +
    geom_sf(data = clr.shape) +
    geom_point(aes(x = Longitude,
                   y = Latitude,
                   colour = Position,
                   shape = Position)) +
    theme_bw() +
    theme(axis.title = element_text(size = 12, colour = "black"),
          axis.text = element_text(size = 10, colour = "black"),
          panel.grid = element_blank()) +
    coord_sf(xlim = c(min(data$Longitude), max(data$Longitude)), 
             ylim = c(min(data$Latitude), max(data$Latitude)),
             expand = T) +
    ggtitle(crab)
}

plot_track(data = tracks, crab =  "A69-1602-12345")

crabs <- unique(tracks$Transmitter)

for (crab in crabs) {
  save_plot(tracks, crab)
}






#################### PLaying with extraction of brownian variance
brownian <- data.frame(brownian = dbbmm$dbbmm$Year_2@DBMvar@means)
coords <- data.frame(lon = dbbmm$dbbmm$Year_2@DBMvar@coords[,1],
                     lat = dbbmm$dbbmm$Year_2@DBMvar@coords[,2])
time <- data.frame(date.time = ymd_hms(dbbmm$dbbmm$Year_2@DBMvar@timestamps)) %>%
  mutate(date = date(date.time))
breaks <- dbbmm$dbbmm$Year_2@DBMvar@break.list
data <- data.frame(bind_cols(brownian, coords, time))


stations <- data.frame(cbind(station = dbbmm$spatial$stations$Station.name,
                             lat = dbbmm$spatial$stations$Latitude,
                             lon = dbbmm$spatial$stations$Longitude)) %>%
  mutate(lat = as.numeric(lat)) %>%
  mutate(lon = as.numeric(lon))

#data <- data %>%
 # mutate(coords.x1 = as.numeric(coords.x1)) %>%
  #mutate(coords.x2 = as.numeric(coords.x2)) %>%
  #mutate(time)

utmcoord <- SpatialPoints(data[, 2:3],
                          proj4string = CRS("+proj=utm +zone=56 ellps=WGS84"))
llcoord <- spTransform(utmcoord,
                       CRS("+proj=longlat +datum=WGS84"))

data$lon <- attr(llcoord, "coords")[,1]
data$lat <- attr(llcoord, "coords")[,2]

ggplot() +
  geom_path(data = data,
             aes(x = lon,
                 y = lat,
                 colour = brownian),
            size = 2) +
 # scale_colour_viridis_b() +
  geom_point(data = stations,
             aes(x = lon,
                 y = lat)) +
  geom_text(data = stations,
            aes(x = lon, 
                y = lat,
                label = station)) 

ggplot() +
  geom_path(data = data,
             aes(x = date,
                 y = brownian,
                 colour = brownian)) +
  scale_x_date(breaks = "1 week") +
  scale_colour_viridis_b() + theme_classic()
