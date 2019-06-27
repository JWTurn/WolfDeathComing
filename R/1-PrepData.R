### Data ====
# Julie Turner
# Started: June 20 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
libs <- c('data.table', 'spatsoc', 'dplyr', 'amt', 'lubridate', 'raster')
lapply(libs, require, character.only = TRUE)

### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolf_rarified.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")

#### determining nearest neighbor (nn) at orginal point ####
dat.grptimes <- group_times(DT=dat, datetime = "datetime", threshold = "2 hours")

dat.nn <- edge_nn(
  DT = dat.grptimes,
  id = 'WolfID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  # returnDist = TRUE, #(not yet available)
  splitBy = 'PackID'
)

dat.nn <- merge(dat.grptimes,dat.nn, by.x = c('WolfID','PackID', 'timegroup'), 
                by.y = c('ID','PackID', 'timegroup'), all.x = T)

### create steps first, then do nn, but only compare to original nn for each random step, 
### not whoever may be the new closest  

# utm zone 14n
# wgs84 32614
# nad83 26914
# project raster

#### test ####
DT.prep <- dat.nn %>% dplyr::select(x = "X", y = "Y", t = 'datetime', id = "WolfID", nn = 'NN') %>%
  filter(id == "W06") %>%
  filter(t >= as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S") %m-% months(1))

DT <- amt::make_track(DT.prep, x, y, t, crs = sp::CRS("+init=epsg:32614")) # %>% 
# amt::transform_coords(sp::CRS("+init=epsg:26914")) # may not need to transform


summarize_sampling_rate(dat)


stps <- amt::track_resample(DT, rate = minutes(120), tolerance = minutes(10)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = T)

quantile(stps$sl_)


land <- raster(paste0(raw, 'RMNPlandcover.tif'), )
cland <- fread(paste0(raw, 'rcl_cowu.csv'))
land <- raster::reclassify(land, cland)
closed <- land == 1
names(closed) <- "closed"
open <- land == 2
names(open) <- "open"
wet <- land == 3
names(wet) <- "wet"
parkYN <- raster(paste0(raw, 'parkYN.tif'))
parkbound_dist <- raster(paste0(raw, 'boundary_dist.tif'))
roads <- raster(paste0(raw, 'roadall_LinFeat_dist.tif'))

stps.parkdist <- stps %>% extract_covariates(parkbound_dist, where = "both") %>%
  mutate(lnparkdist_start = log(boundary_dist_start + 1)) %>%
  mutate(lnparkdist_end = log(boundary_dist_end + 1))


stps.park <- stps %>% extract_covariates(parkYN, where = "both") %>%
  mutate(parkYN_start = factor(parkYN_start, levels = c(0, 1), labels = c("park", "outside park"))) %>%
  mutate(parkYN_end = factor(parkYN_end, levels = c(0, 1), labels = c("park", "outside park")))

p1.park <- stps.park %>% dplyr::select(parkYN_end, tod = tod_end_, sl_, ta_) %>%
  
  gather(key, val, -parkYN_end, -tod) %>%
  filter(key == "sl_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ parkYN_end, nrow = 2) +
  xlab("Step length [m]") + theme_light() +
  ylab("Density") +
  theme(legend.title = element_blank())

stps.wet <- stps %>% extract_covariates(wet, where = "both") %>%
  mutate(wet_start = factor(wet_start, levels = c(0,1), labels = c('other', 'wet'))) %>%
  mutate(wet_end = factor(wet_end, levels = c(0,1), labels = c('other', 'wet')))


stps.closed <- stps %>% extract_covariates(closed, where = "both") %>%
  mutate(closed_start = factor(closed_start, levels = c(0,1), labels = c('other', 'closed'))) %>%
  mutate(closed_end = factor(closed_end, levels = c(0,1), labels = c('other', 'closed')))

stps.land <- stps %>% extract_covariates(land, where = "both") %>%
  mutate(land_start = factor(RMNPlandcover_start, levels = 1:3, labels = c("closed", "open", 'wet'))) %>%
  mutate(land_end = factor(RMNPlandcover_end, levels = 1:3, labels = c("closed", "open", 'wet')))

p1.land <- stps.land %>% dplyr::select(land_end, tod = tod_end_, sl_, ta_) %>%
  
  gather(key, val, -land_end, -tod) %>%
  filter(key == "sl_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ land_end, nrow = 2) +
  xlab("Step length [m]") + theme_light() +
  ylab("Density") +
  theme(legend.title = element_blank())

#### calculate nn dist -- from Alec ####
# group_times
# DT.prep[,'round.t'] <- as.POSIXct(round(DT.prep$t, units = 'hours'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")
stps.nn <- merge(stps, DT.prep, by.x = c('x1_', 'y1_', 't1_'), by.y = c('x', 'y', 't'))


# measure within some maximum distance
edist <- edge_dist(
  DT = stps,
  threshold = 50000,
  id = 'WolfID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  returnDist = TRUE,
  splitBy = 'PackID'
)

# measure nearest neighbours
enn <- edge_nn(
  DT = DT,
  threshold = 50000,
  id = 'ID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  # returnDist = TRUE, #(not yet available)
  splitBy = 'packID'
)

merge(enn,
      edist,
      by.x = c('ID', 'NN'),
      by.y = c('ID1', 'ID2'),
      all.x = TRUE)


#### extract covariate for SSF ####
ssf1 <- stps %>% random_steps(n = 9) %>%
  amt::extract_covariates(land, where = "both")  %>%
  amt::extract_covariates(parkYN, where = "both") %>%
  amt::extract_covariates(parkbound_dist, where = "both") %>%
  amt::time_of_day(include.crepuscule = T) %>%
  mutate(land_start = factor(RMNPlandcover_start, levels = 1:3, labels = c("closed", "open", 'wet')),
        land_end = factor(RMNPlandcover_end, levels = 1:3, labels = c("closed", "open", 'wet')),
        parkYN_start = factor(parkYN_start, levels = c(0, 1), labels = c("park", "out-park")),
        parkYN_end = factor(parkYN_end, levels = c(0, 1), labels = c("park", "out-park")),
        lnparkdist_start = log(boundary_dist_start + 1),
        lnparkdist_end = log(boundary_dist_end + 1),
        log_sl = log(sl_),
        cos_ta = cos(ta_)) 
ssf1[,'ttd1'] <- as.duration(ssf1$t1_ %--% as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S"))/ddays(1) 
ssf1[,'ttd2'] <- as.duration(ssf1$t2_ %--% as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S"))/ddays(1) 
DT.nn <- dplyr::select(DT.prep, t, nn)
ssf <- merge(ssf1, DT.nn, by.x = 't1_', by.y = 't', all.x = T)
colnames(ssf)[colnames(ssf)=="nn"] <- "nn1"
ssf <- merge(ssf, DT.nn, by.x = 't2_', by.y = 't', all.x = T)
colnames(ssf)[colnames(ssf)=="nn"] <- "nn2"


### calculate all distances ###
remotes::install_github('ropensci/spatsoc')
library(spatsoc)
DT.grp <- select(dat.grptimes, WolfID, PackID, X, Y, datetime, timegroup)
DT.randpts <- merge(ssf, DT.grp, by.y = c('datetime', 'X', 'Y'), by.x = c('t2_', 'x2_', 'y2_'), all.x =T)
edist2 <- edge_dist(
  DT = DT.randpts,
  threshold = 100000, # I don't know how to pick this
  id = 'WolfID',
  coords = c('x2_', 'y2_'),
  timegroup = 't2_',
  returnDist = TRUE,
  splitBy = 'PackID'
)

### didn't work, can't read it for some reason

#######
group_times(dat, datetime = 'datetime', threshold = '20 minutes')

dat_nn <- edge_nn(dat, id = 'WolfID', coords = c('X', 'Y'),
                 timegroup = 'timegroup', splitBy = 'PackID')

