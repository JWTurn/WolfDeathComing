library(lubridate)
library(raster)
library(amt)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)




dat <- read.csv("data/RMNPwolf_rarified.csv") 
  #drop_na() %>%     #####filter !is.na did not work
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")
dat<- dat %>% dplyr::select(x = "X", y = "Y", t = 'datetime', id = "WolfID") %>%
  filter(id == "W06") 



# utm zone 14n
# wgs84 32614
# nad83 26914
# project raster

dat <- amt::make_track(dat, x, y, t, crs = sp::CRS("+init=epsg:32614")) %>%
  amt::transform_coords(sp::CRS("+init=epsg:26914")) # may not need to transform


summarize_sampling_rate(dat)


stps <- amt::track_resample(dat, rate = minutes(120), tolerance = minutes(10)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = T)

land_use <- raster("data/Opendeciduous100m.tif")
wet_dist <- raster("data/Water_Dist.tif")
od <- land_use == 11
names(od) <- "od"

eda1 <- stps %>% extract_covariates(wet_dist, where = "both") %>%
  mutate(lnwetdist_start = log(Water_Dist_start+1)) %>%
  mutate(lnwetdist_end = log(Water_Dist_end+1))

p1 <- eda1 %>% dplyr::select(lnwetdist_end, tod = tod_end_, sl_, ta_) %>%
  
  gather(key, val, -lnwetdist_end, -tod) %>%
  filter(key == "sl_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ lnwetdist_end, nrow = 2) +
  xlab("Step length [m]") + theme_light() +
  ylab("Density") +
  theme(legend.title = element_blank())

p2 <- eda1 %>% dplyr::select(lnwetdist_end, tod = tod_end_, sl_, ta_) %>%
  gather(key, val, -lnwetdist_end, -tod) %>%
  filter(key == "ta_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ lnwetdist_end, nrow = 2) +
  xlab("Turn angle") + theme_light() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank())

library(cowplot)
pg1 <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), rel_widths = c(1, 1)
)

leg <- get_legend(p1)
plot_grid(pg1, leg, rel_widths = c(1, 0.1))
