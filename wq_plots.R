library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)

# load clarence logger datasets
clr.wq.1 <- read_csv("data_raw/clr_wq_1.csv") %>% rename(Temp_CAL = Temp_C, Cond_CAL = `Calibrated_Conductivity_mS/cm`)
clr.wq.2 <- read_csv("data_raw/clr_wq_2.csv") 
clr.wq.3 <- read_csv("data_raw/clr_wq_3.csv") 
clr.wq.4 <- read_csv("data_raw/clr_wq_4.csv")
clr.wq.5 <- read_csv("data_raw/clr_wq_5.csv")

# reformat and summarise to daily mean values
clr.wq.1 <- clr.wq.1 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-28" & date < "2020-12-18")

clr.wq.2 <- clr.wq.2 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-28" & date < "2020-12-18")

clr.wq.3 <- clr.wq.3 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-28" & date < "2020-12-18")

clr.wq.4 <- clr.wq.4 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-28" & date < "2020-12-18")

clr.wq.5 <- clr.wq.5 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-28" & date < "2020-12-18")

clr.temp <- ggplot() +
  geom_line(data = clr.wq.1,
            aes(x = date, y = mean.temp),
            alpha = 0.5) +
  geom_line(data = clr.wq.2,
            aes(x = date,  y = mean.temp),
            alpha = 0.5) +
  geom_line(data = clr.wq.4,
            aes(x = date, y = mean.temp),
            alpha = 0.5) +
  geom_line(data = clr.wq.5,
            aes(x = date, y = mean.temp),
            alpha = 0.5) +
  geom_line(data = clr.wq.3,
            aes(x = date, y = mean.temp),
            colour = "red", size = 1) +
  geom_vline(xintercept = c(date("2019-02-28"), 
                            date("2019-05-20"), 
                            date("2020-11-25"), 
                            date("2020-12-18")), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste("Temperature (",degree,"C)"))) +
  xlab("Date") +
  scale_x_date(date_labels = "%d %b %Y")

clr.cond <- ggplot() +
  geom_line(data = clr.wq.1,
            aes(x = date, y = mean.cond),
            alpha = 0.5) +
  geom_line(data = clr.wq.2,
            aes(x = date, y = mean.cond),
            alpha = 0.5) +
  geom_line(data = clr.wq.3,
            aes(x = date, y = mean.cond),
            alpha = 0.5) +
  geom_line(data = clr.wq.5,
            aes(x = date, y = mean.cond),
            alpha = 0.5) +
  geom_line(data = clr.wq.3,
            aes(x = date, y = mean.cond),
            colour = "blue", size = 1) +
  geom_vline(xintercept = c(date("2019-02-28"), 
                            date("2019-05-20"), 
                            date("2020-11-25"), 
                            date("2020-12-18")), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste("Conductivity (", mS~cm^{-1},")"))) +
  xlab("Date") +
  scale_x_date(date_labels = "%d %b %Y")

# load kalang logger data sets
blr.wq.4 <- read_csv("data_raw/blr_wq_4.csv")
blr.wq.6 <- read_csv("data_raw/blr_wq_6.csv")
blr.wq.9 <- read_csv("data_raw/blr_wq_9.csv")

# reformat and summarise to daily mean values
blr.wq.4 <- blr.wq.4 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-20" & date < "2020-03-08")

blr.wq.6 <- blr.wq.6 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-20" & date < "2020-03-08")

blr.wq.9 <- blr.wq.9 %>%
  mutate(date = dmy(Date_text)) %>%
  group_by(date) %>%
  summarise(mean.temp = mean(Temp_CAL),
            sd.temp = sd(Temp_CAL),
            mean.cond = mean(Cond_CAL),
            sd.cond = sd(Cond_CAL)) %>%
  ungroup() %>%
  filter(date > "2019-02-20" & date < "2020-03-08")

blr.temp <- ggplot() +
  geom_line(data = blr.wq.4,
            aes(x = date, y = mean.temp),
            alpha = 0.5) +
  geom_line(data = blr.wq.9,
            aes(x = date,  y = mean.temp),
            alpha = 0.5) +
  geom_line(data = blr.wq.6,
            aes(x = date, y = mean.temp),
            colour = "red", size = 1) +
  geom_vline(xintercept = c(date("2019-02-20"), 
                            date("2019-06-03"), 
                            date("2020-01-30"), 
                            date("2020-03-08")), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste("Temperature (",degree,"C)"))) +
  xlab("Date") +
  scale_x_date(date_labels = "%d %b %Y")

blr.cond <- ggplot() +
  geom_line(data = blr.wq.4,
            aes(x = date, y = mean.cond),
            alpha = 0.5) +
  geom_line(data = blr.wq.9,
            aes(x = date,  y = mean.cond),
            alpha = 0.5) +
  geom_line(data = blr.wq.6,
            aes(x = date, y = mean.cond),
            colour = "blue", size = 1) +
  geom_vline(xintercept = c(date("2019-02-20"), 
                            date("2019-06-03"), 
                            date("2020-01-30"), 
                            date("2020-03-08")), 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste("Conductivity (", mS~cm^{-1},")"))) +
  xlab("Date") +
  scale_x_date(date_labels = "%d %b %Y")

wq.plot <- (clr.temp|clr.cond)/(blr.temp|blr.cond) + plot_annotation(tag_levels = "a")

ggsave("figures/wq_timeseries.png", 
       plot = wq.plot, 
       device = "png", 
       width = 29, # a4 dimensions
       height = 20, 
       units = "cm", 
       dpi = 600)
