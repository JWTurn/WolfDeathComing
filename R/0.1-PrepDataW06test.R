### Data ====
# Julie Turner
# Started: June 20 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'spatsoc', 'dplyr', 'amt', 'lubridate', 'raster', 'tidyr', 'ggplot')
lapply(libs, require, character.only = TRUE)

### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolf_rarified.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")

dat.meta <- fread(paste0(raw, 'wolf_metadata.csv'))
dat.meta$end_date<- as.POSIXct(dat.meta$end_date, tz = 'UTC', "%Y-%m-%d")

#### determining nearest neighbor (nn) at orginal point ####
dat.grptimes <- group_times(DT=dat, datetime = "datetime", threshold = "2 hours")

enn <- edge_nn(
  DT = dat.grptimes,
  id = 'WolfID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  # returnDist = TRUE, #(not yet available)
  splitBy = 'PackID'
)

### calculating distances and original point
edist1 <- edge_dist(
  DT = dat.grptimes,
  threshold = 10000000, # I don't know how to pick this
  id = 'WolfID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  returnDist = TRUE,
  splitBy = 'PackID'
)

dat.nn <- merge(enn, edist1, by.x = c('ID', 'NN','PackID', 'timegroup'), 
                by.y = c('ID1', 'ID2', 'PackID', 'timegroup'))

#colnames(dat.nn)[colnames(dat.nn)=="distance"] <- "distance1"


dat.nn <- merge(dat.grptimes,dat.nn, by.x = c('WolfID','PackID', 'timegroup'), 
                by.y = c('ID','PackID', 'timegroup'), all.x = T)

### OBJECTIVE: determine nn first from actual data, then create (random) steps
### determine nn dist only from actual data, not who may be nearest at new random step

# utm zone 14n
# wgs84 32614
# nad83 26914
# project raster

#### test w/ W06 ####
DT.prep <- dat.nn %>% dplyr::select(x = "X", y = "Y", t = 'datetime', id = "WolfID", nn = 'NN', distance1 = 'distance') %>%
  filter(id == "W06") %>%
  filter(t >= as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S") %m-% months(1))

# making track
DT <- amt::make_track(DT.prep, x, y, t, crs = sp::CRS("+init=epsg:32614")) # %>% 
# amt::transform_coords(sp::CRS("+init=epsg:26914")) # may not need to transform


summarize_sampling_rate(dat)

# making steps from track
stps <- amt::track_resample(DT, rate = minutes(120), tolerance = minutes(10)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = T) # need to fix this because missing dusk (KK?)

quantile(stps$sl_)

#### layers ####

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


#### exploratory checks ####

stps.parkdist <- stps %>% extract_covariates(parkbound_dist, where = "both") %>%
  mutate(lnparkdist_start = log(boundary_dist_start + 1)) %>%
  mutate(lnparkdist_end = log(boundary_dist_end + 1))


stps.park <- stps %>% extract_covariates(parkYN, where = "both") %>%
  mutate(parkYN_start = factor(parkYN_start, levels = c(0, 1), labels = c("park", "outside park"))) %>%
  mutate(parkYN_end = factor(parkYN_end, levels = c(0, 1), labels = c("park", "outside park")))

p1.park <- stps.park %>% dplyr::select(parkYN_end, tod = tod_end_, sl_, ta_) %>%
  
  tidyr::gather(key, val, -parkYN_end, -tod) %>%
  dplyr::filter(key == "sl_") %>%
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
  
  tidyr::gather(key, val, -land_end, -tod) %>%
  filter(key == "sl_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ land_end, nrow = 2) +
  xlab("Step length [m]") + theme_light() +
  ylab("Density") +
  theme(legend.title = element_blank())





#### extract covariate for SSF ####
ssf1 <- stps %>% random_steps(n = 9) %>%
  amt::extract_covariates(land, where = "both")  %>%
  amt::extract_covariates(parkYN, where = "both") %>%
  amt::extract_covariates(parkbound_dist, where = "both") %>%
  amt::time_of_day(include.crepuscule = T, where = 'start') %>%
  mutate(land_start = factor(RMNPlandcover_start, levels = 1:3, labels = c("closed", "open", 'wet')),
        land_end = factor(RMNPlandcover_end, levels = 1:3, labels = c("closed", "open", 'wet')),
        parkYN_start = factor(parkYN_start, levels = c(0, 1), labels = c("park", "out-park")),
        parkYN_end = factor(parkYN_end, levels = c(0, 1), labels = c("park", "out-park")),
        lnparkdist_start = log(boundary_dist_start + 1),
        lnparkdist_end = log(boundary_dist_end + 1),
        log_sl = log(sl_),
        cos_ta = cos(ta_)) 

# adding time to death at point 1 and 2
ssf1[,'ttd1'] <- as.duration(ssf1$t1_ %--% as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S"))/ddays(1) 
ssf1[,'ttd2'] <- as.duration(ssf1$t2_ %--% as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S"))/ddays(1) 
# adding nn at point 1 and 2
DT.nn.dist <- dplyr::select(DT.prep, t, id, nn, distance1)
ssf <- merge(ssf1, DT.nn.dist, by.x = 't1_', by.y = 't', all.x = T)
colnames(ssf)[colnames(ssf)=="nn"] <- "nn1"
DT.nn <- dplyr::select(DT.prep, t, id, nn)
ssf <- merge(ssf, DT.nn, by.x = c('t2_', 'id'), by.y = c('t', 'id'), all.x = T)
colnames(ssf)[colnames(ssf)=="nn"] <- "nn2"



### calculate all distances for step 2 ###
# need to combine all the distance data with the (random) steps for the focal indiv so neighbor dists can be calculated
dat.grp <- dplyr::select(dat.grptimes, WolfID, PackID, X, Y, datetime, timegroup)
dat.randpts <- merge(ssf, dat.grp, by.y = c('datetime', 'WolfID'), 
                     by.x = c('t2_', 'id'), all =T)
# use the gps pts from the track unless missing (as with all non focal indivs)
dat.randpts$x2_ <- ifelse(!is.na(dat.randpts$x2_), dat.randpts$x2_, dat.randpts$X)
dat.randpts$y2_ <- ifelse(!is.na(dat.randpts$y2_), dat.randpts$y2_, dat.randpts$Y)

# switching to data.table so can run edge_dist() and easier formatting 
DT.randpts <-setDT(dat.randpts)

# limiting it to just time for this test indiv W06
DT.randpts <- DT.randpts[t2_ >= as.POSIXct("2016-06-26 23:59:59", tz = 'UTC', "%Y-%m-%d %H:%M:%S") %m-% months(1)]
# enumerating each random step so can pair distances correctly
DT.randpts <- DT.randpts[, step_id_rank:= 1:.N, by=.(step_id_)]
# making all regular steps (without random steps) have a value of 1
DT.randpts$step_id_rank <- ifelse(!is.na(DT.randpts$step_id_), DT.randpts$step_id_rank, 1)

# merging id with info that will help us pair distances with the correct steps and individuals after they're calculated
# edge_dist() doesn't output a lot of the info needed to do this, so will separated it again when done
DT.randpts[,'id.step'] <- paste(DT.randpts$id, DT.randpts$step_id_rank, sep = '.')
DT.randpts[,'id.step.t'] <- paste(DT.randpts$id.step, DT.randpts$timegroup, sep = '.')

# calculates distances for all individuals within the same pack at each timegroup (this is just for step 2 right now)
edist2 <- edge_dist(
  DT = DT.randpts,
  threshold = 10000000, # I don't know how to pick this
  id = 'id.step.t',
  coords = c('x2_', 'y2_'),
  timegroup = 'timegroup',
  returnDist = TRUE,
  splitBy = 'PackID'
)


# just checking with my test individual of W06, won't need to subdivide later
edist2.ww <- edist2[PackID =='WW']
edist2.ww.sep <- tidyr::separate(edist2.ww, col = ID1, into = c('id1', 'step1', 't1'))
edist2.ww.sep <- tidyr::separate(edist2.ww.sep, col = ID2, into = c('id2', 'step2', 't2'))

#had to make step an integer so that they match between data frames (don't know why step1 and step2 think they're chr)
edist2.ww.sep[,'step'] <- as.integer(edist2.ww.sep$step1)


# Keeping only the necessary columns
edist2.ww.sep<- edist2.ww.sep[,.(timegroup, step, id1, id2, distance)]




# making the times match the ones amt made for the track
DT.grp <- dplyr::select(dat.grp, datetime, timegroup)
DT.grp[,'round.t'] <- as.POSIXct(round(DT.grp$datetime, units = 'mins'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")
DT.grp2 <- setDT(DT.grp)[,.(round.t), by = .(timegroup)] 

# getting a list of all timegroups and associated times to help match points correctly
DT.grp2 <- unique(DT.grp2)

# numbering the steps and random steps to help match the social partners and their distances
ssf.soc <- setDT(ssf)[, step_id_rank:= 1:.N, by=.(step_id_)]
ssf.soc<- merge(ssf.soc, DT.grp2, by.x = 't2_', by.y = 'round.t', all.x = T)
colnames(ssf.soc)[colnames(ssf.soc)=="timegroup"] <- "timegroup2"
ssf.soc<- merge(ssf.soc, DT.grp2, by.x = 't1_', by.y = 'round.t', all.x = T)
colnames(ssf.soc)[colnames(ssf.soc)=="timegroup"] <- "timegroup1"


ssf.soc <- merge(ssf.soc, edist2.ww.sep, by.x = c('id','nn2', 'step_id_rank', 'timegroup2'),
                 by.y = c('id1','id2', 'step', 'timegroup'), all.x = T)
colnames(ssf.soc)[colnames(ssf.soc)=="distance"] <- "distance2"


ssf.W06 <- ssf.soc[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta, tod_start_, 
                      parkYN_start, parkYN_end, lnparkdist_start, lnparkdist_end, land_start, land_end, 
                      id, nn1, nn2, distance1, distance2, timegroup1, timegroup2,
                      ttd1, ttd2)]

#saveRDS(ssf.W06, 'data/derived-data/ssfW06.Rds')

#### fit models ####

### questions
  # if I'm using ttd1, do I need all interactions to be with the other start points?
  # do I need to log transform distance to conspecifics? or ttd? (looks like need to do for ttd)
  # do I need cos_ta or cos_ta:tod

m.core <- ssf.W06 %>% amt::fit_issf(case_ ~ log_sl:tod_start_ + cos_ta:tod_start_ + land_end + log_sl:land_end + strata(step_id_))

m.full <- ssf.W06 %>% amt::fit_issf(case_ ~ log_sl:tod_start_ + cos_ta:tod_start_ + land_end + log_sl:land_end + 
                                      ttd1:log_sl + ttd1:cos_ta +
                                      ttd1:land_start + ttd1:lnparkdist_start:parkYN_start + 
                                      ttd1:distance2 + 
                                      strata(step_id_))
# full not converging

m.move <- ssf.W06 %>% amt::fit_issf(case_ ~ log_sl:tod_start_ + cos_ta:tod_start_ + land_end + log_sl:land_end + 
                                      ttd1:log_sl + ttd1:cos_ta +
                                      strata(step_id_))


m.habitat.parkdist <- ssf.W06 %>% amt::fit_issf(case_ ~ log_sl:tod_start_ + cos_ta:tod_start_ + land_end + log_sl:land_end + 
                                      log(ttd1):land_end + log(ttd1):lnparkdist_start  + 
                                      strata(step_id_))

# ttd1:lnparkdist_start:parkYN_start won't converge

m.habitat.parkYN <- ssf.W06 %>% amt::fit_issf(case_ ~ log_sl:tod_start_ + cos_ta:tod_start_ + land_end + log_sl:land_end + 
                                                  log(ttd1):land_end + log(ttd1):parkYN_start  + 
                                                  strata(step_id_))
## didn't converge


m.social <- ssf.W06 %>% amt::fit_issf(case_ ~ log_sl:tod_start_ + cos_ta:tod_start_ + land_end + log_sl:land_end + 
                                      log(ttd1):log(distance2) + 
                                      strata(step_id_))


bbmle::AICtab(m.core$model, m.move$model, m.habitat.parkdist$model, m.social$model)


