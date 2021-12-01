library(tidyverse)
library(lubridate)

ocean.crabs <- read_csv("data_processed/compiled_offshore_detections.csv")

ocean.crabs <- ocean.crabs %>%
  mutate(time = dmy_hm(time),
         last.estuary.detection = dmy_hm(last.estuary.detection))

ocean.crab.sum <- ocean.crabs %>%
  filter(crab != "A69-1602-12365") %>%
  group_by(crab, station.name) %>%
  mutate(n = n(),
         km.per.day = lat.distance.km/duration) %>%
  distinct(crab, .keep_all = TRUE)
