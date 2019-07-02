### Data ====
# Julie Turner
# Started: July 1 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'spatsoc', 'dplyr', 'amt', 'lubridate', 'raster', 'tidyr', 'ggplot2')
lapply(libs, require, character.only = TRUE)

### Functions ----
#function CreateTrack will put dataset into track/step format,
#calculates steps, creates random steps, and extracts stoich layer. 
CreateTrack <- function(x.col, y.col, date.col,  ID, NumbRandSteps, sl_distr, ta_distr) {
  #print(ID)
  #create track from dataset
  trk <- track(x.col, y.col, date.col, ID, crs14) %>% 
    #function turns locs into steps
    steps() 
  #remove any steps that span more than 2hr
  trk <- subset(trk, trk$dt_ > 1.9 & trk$dt_ < 2.1, drop = T)
  #generate random steps
  trk %>%
    random_steps(n = NumbRandSteps, sl_distr, ta_distr) %>%
    amt::time_of_day(include.crepuscule = T, where = 'start') ### check with KK b/c this isn't right
}



### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolf_rarified.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")

dat.meta <- fread(paste0(raw, 'wolf_metadata.csv'))
dat.meta$end_date <- paste(dat.meta$end_date, '23:59:59', sep = ' ')
dat.meta$end_date<- as.POSIXct(dat.meta$end_date, tz = 'UTC', "%Y-%m-%d %H:%M:%S")

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




dat.nn <- merge(dat.grptimes,dat.nn, by.x = c('WolfID','PackID', 'timegroup'), 
                by.y = c('ID','PackID', 'timegroup'), all.x = T)

dat.nn <- merge(dat.nn,dat.meta, by.x = c('WolfID','PackID'), 
                by.y = c('WolfID','PackID'), all.x = T)


### OBJECTIVE: determine nn first from actual data, then create (random) steps
### determine nn dist only from actual data, not who may be nearest at new random step

# utm zone 14n
# wgs84 32614
# nad83 26914
# project raster

utm14N <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
crs14 <- sp::CRS("+init=epsg:32614")

#### filter those that die ####
dat.focal <- setDT(dat.meta)[status == 'dead' & !is.na(PackID)]
focals <- dat.focal$WolfID

DT.prep <- dat.nn %>% dplyr::select(x = "X", y = "Y", t = 'datetime', id = "WolfID", nn = 'NN', distance1 = 'distance',
                                    'status', 'end_date', 'COD') %>%
  filter(id %in% focals) 
  
  

# making track
#not working data.table way



trk <- track(DT.prep$x, DT.prep$y, DT.prep$t, DT.prep$id, crs =  crs) %>% steps() 
  trk <- subset(trk, trk$dt_ > 1.9 & trk$dt_ < 2.1, drop = T)
  trk %>% random_steps(n = NumbRandSteps, sl_distr, ta_distr) %>%
  amt::time_of_day(include.crepuscule = T, where = 'start')


DT <- setDT(DT.prep)[, CreateTrack(x.col=x, y.col=y, date.col=t,  ID=id, NumbRandSteps=9, 
                                   sl_distr = "gamma", ta_distr = "vonmises"), by = 'id']
# can't figure out how to summarize sampling rate this way

DT.trk <- setDT(DT.prep)[,track(x, y, t, crs = sp::CRS("+init=epsg:32614")), by = 'id']


stps <- DT[,amt::track_resample(DT.trk, rate = minutes(120), tolerance = minutes(10)), by = 'id']
#%>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = T) # need to fix this because missing dusk (KK?)


dat_all <- DT.prep %>% group_by(id) %>% nest()

dat_all <- dat_all %>%
  mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:32614")) 
  }))  

dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest

dat.unnest <- dat_all %>% unnest

# stps <- dat_all %>%
#   mutate(steps = map(trk, function(x) {
#     x %>% amt::track_resample(rate = minutes(120), tolerance = minutes(10)) %>%
#       amt::filter_min_n_burst(min_n = 3) %>%
#       amt::steps_by_burst()}))
# 
# stps.unnest<- stps %>% unnest(.key = id)

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


#### making steps from track and extracting covariates ####
ssf <- dat_all %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = minutes(120), tolerance = minutes(10)) %>%
      amt::filter_min_n_burst(min_n = 3) %>%
      amt::steps_by_burst() %>% amt::random_steps(n=9) %>%
      amt::extract_covariates(land, where = "both")  %>%
      amt::extract_covariates(parkYN, where = "both") %>%
      amt::extract_covariates(parkbound_dist, where = "both") %>%
      amt::time_of_day(include.crepuscule = T, where = 'start') %>%  ####check with KK on doing this better
      mutate(land_start = factor(RMNPlandcover_start, levels = 1:3, labels = c("closed", "open", 'wet')),
             land_end = factor(RMNPlandcover_end, levels = 1:3, labels = c("closed", "open", 'wet')),
             parkYN_start = factor(parkYN_start, levels = c(0, 1), labels = c("park", "out-park")),
             parkYN_end = factor(parkYN_end, levels = c(0, 1), labels = c("park", "out-park")),
             lnparkdist_start = log(boundary_dist_start + 1),
             lnparkdist_end = log(boundary_dist_end + 1),
             log_sl = log(sl_),
             cos_ta = cos(ta_))
  }))

#stps.unnest <- stps %>% unnest


ssf.all <- ssf %>% dplyr::select(id, steps) %>% unnest
  








# adding time to death at point 1 and 2
ssf1 <- merge(ssf1, dat.focal, by )
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
