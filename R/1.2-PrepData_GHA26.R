### Data GHA26 ====
# Julie Turner
# Started: July 31 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'spatsoc', 'dplyr', 'amt', 'lubridate', 'raster', 'tidyr', 'ggplot2')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### combining years of data ----
dat2014_15 <- readRDS(paste0(raw, 'GHA26_Full_2014-15_2hr.Rds'))
dat2016 <- readRDS(paste0(raw, 'GHA26_Full_2016_2hr.Rds'))
dat2016_st <- readRDS(paste0(raw, 'GHA26_ST.Full_2016_2hr.Rds'))
dat2017 <- readRDS(paste0(raw, 'GHA26_Full_2017_2hr.Rds'))
dat2018 <- readRDS(paste0(raw, 'GHA26_Full_2018_2hr.Rds'))

dat2014_15 <- dat2014_15[,.(Collar_ID, MBts = ts, Latitude, Longitude)]
dat2014_15$MBts <- as.POSIXct(dat2014_15$MBts, tz = 'MST', "%Y-%m-%d %H:%M:%S")
dat2016 <- dat2016[,.(Collar_ID, MBts, Latitude, Longitude)]
dat2016_st <- dat2016_st[,.(Collar_ID, MBts = ts, Latitude, Longitude)]
dat2017 <- dat2017[,.(Collar_ID, MBts, Latitude, Longitude)]
dat2018 <- dat2018[,.(Collar_ID, MBts, Latitude, Longitude)]

dat <- rbind(dat2014_15, dat2016, dat2016_st, dat2017, dat2018)

# utm zone 14n
# wgs84 32614
# nad83 26914
# project raster

utm14N <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
crs14 <- sp::CRS("+init=epsg:32614")

set.seed(57)

dat[, Longitude := as.numeric(Longitude)]
dat[, Latitude := as.numeric(Latitude)]
dat[, c('X', 'Y') := as.data.table(rgdal::project(cbind(Longitude, Latitude), utm14N))]
dat[,'year'] <- as.numeric(format(dat$MBts,'%Y'))

#write.csv(dat, paste0(derived, 'GHA26_wolf_pts.csv'))

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'), header = T)
dat.meta$death_date <- paste(dat.meta$death_date, '23:59:59', sep = ' ')
dat.meta$death_date<- as.POSIXct(dat.meta$death_date, tz = 'UTC', "%Y-%m-%d %H:%M:%S")
dat.meta[,'year'] <- as.numeric(format(dat.meta$death_date, '%Y'))

#dat.collpack <- dat.meta[,.(WolfID, PackID, year)]

collarID <- fread(paste0(raw, 'GHA26_coll_wolf_id.csv'), header = T)
collarID$Collar_ID <- as.character(collarID$Collar_ID)

#dat.collpack <- dat.collpack[Collar_ID %chin% dat$WolfID]
dat.all <- merge(dat, collarID, by = c('Collar_ID', 'year'), all.x = T)
dat.all[,'WolfID'] <- ifelse(dat.all$Collar_ID=='5148', 'W17', dat.all$WolfID)
dat.all[,'PackID'] <- ifelse(dat.all$Collar_ID=='5148', 'MA', dat.all$PackID)

colnames(dat.all)[colnames(dat.all)=="MBts"] <- "datetime"

dat.all <- dat.all[!is.na(WolfID)]

# range <- dat.all[,.(min=min(datetime), max=max(datetime)), by=.(WolfID)]  
# range <- range[,.(min, max, len=max-min), .(WolfID)]
#### determining nearest neighbor (nn) at orginal point ####
dat.grptimes <- group_times(DT=dat.all, datetime = "datetime", threshold = "2 hours")

enn <- edge_nn(
  DT = dat.grptimes,
  id = 'WolfID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  returnDist = TRUE, 
  splitBy = 'PackID'
)

enn<- unique(enn)




dat.nn <- merge(dat.grptimes,enn, by.x = c('WolfID','PackID', 'timegroup'), 
                by.y = c('ID','PackID', 'timegroup'), all.x = T)
dat.nn <- dat.nn[!is.na(WolfID)]

dat.nn <- merge(dat.nn,dat.meta, by.x = c('WolfID','PackID'), 
                by.y = c('WolfID','PackID'), all.x = T)
dat.nn <- dat.nn[,year.y:=NULL]
colnames(dat.nn)[colnames(dat.nn)=="year.x"] <- "year"


### OBJECTIVE: determine nn first from actual data, then create (random) steps
### determine nn dist only from actual data, not who may be nearest at new random step


#### filter to those that die ####
dat.focal <- setDT(dat.meta)[pop=='GHA26' & use!='n' & !is.na(PackID)]
dat.focal[,'packbound'] <- dat.focal$PackID
focals <- dat.focal$WolfID

#### setting time to death ####
ttd = 61

# DT.prep <- dat.nn %>% dplyr::select(x = "X", y = "Y", t = 'datetime', id = "WolfID", nn = 'NN', distance1 = 'distance',
#                                     'COD', 'death_date') %>%
#   filter(id %in% focals) 
#   
DT.prep <- dat.nn[WolfID %in% focals,.(x = X, y = Y, t = datetime, id = WolfID, nn = NN, distance1 = distance,
                                       status, end_date, COD, death_date)]
DT.prep[, t2d:=(as.duration(t %--% death_date)/ddays(1))]
DT.prep <- DT.prep[t2d >=0 & t2d<=ttd]

DT.prep[,uniqueN(t), by =.(id)]


dat_all <- DT.prep %>% group_by(id) %>% nest()

dat_all <- dat_all %>%
  mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:32614")) 
  }))  

dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest(cols = c(sr))



#### layers ####

land <- raster(paste0(raw, 'GHA26landcover2015_wgs84.tif'), )
#cland <- fread(paste0(raw, 'rcl_cowu.csv'))
cland2 <- fread(paste0(raw, 'rcl_fine.csv'))
land <- raster::reclassify(land, cland2)
# closed <- land == 1
# names(closed) <- "closed"
# open <- land == 2
# names(open) <- "open"
# wet <- land == 3
# names(wet) <- "wet"
coniferous <- land == 1
names(coniferous) <- "coniferous"
deciduous <- land == 2
names(deciduous) <- "deciduous"
mixed <- land == 3
names(mixed) <- "mixed"
shrub <- land == 4
names(shrub) <- "shrub"
open <- land == 5
names(open) <- "open"
wet <- land == 6
names(wet) <- "wet"
urban <- land == 7
names(urban) <- "urban"
## This creates an object which can be used to make a layer of specified diameter
# The d value is what determines the buffer size if you want to change it.
## If you're doing multiple landcover classes, you only need to run this line once, as long as each of the habitat variables has the same resolution
# (e.g., the "Wetland" here could be any cover type)
Buff100 <- focalWeight(land, d=100, type = 'circle')
## This generates a new raster where each cell corresponds to the mean wetland within a 100m buffer.
# Since it's all 1s and 0s, this is the same as the proportion of wetland surrounding the focal variable
propwet <- focal(wet, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
propconif <- focal(coniferous, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
propdecid <- focal(deciduous, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
propmixed <- focal(mixed, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
propopen <- focal(open, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
propshrub <- focal(shrub, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
propurban <- focal(urban, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
roads <- raster(paste0(raw, 'GHA26_roadsPS_dist.tif'))  

FM <- raster(paste0(raw, 'FM_kde.tif')) 
FM <- reclassify(FM, cbind(NA, 0))
FM_dist <- raster(paste0(raw, 'FM_kde_dist.tif'))
GF <- raster(paste0(raw, 'GF_kde.tif')) 
GF <- reclassify(GF, cbind(NA, 0))
GF_dist <- raster(paste0(raw, 'GF_kde_dist.tif'))
GR <- raster(paste0(raw, 'GR_kde.tif')) 
GR <- reclassify(GR, cbind(NA, 0))
GR_dist <- raster(paste0(raw, 'GR_kde_dist.tif'))
HL <- raster(paste0(raw, 'HL_kde.tif')) 
HL <- reclassify(HL, cbind(NA, 0))
HL_dist <- raster(paste0(raw, 'HL_kde_dist.tif'))
MA <- raster(paste0(raw, 'MA_kde.tif')) 
MA <- reclassify(MA, cbind(NA, 0))
MA_dist <- raster(paste0(raw, 'MA_kde_dist.tif')) 
MW <- raster(paste0(raw, 'MW_kde.tif')) 
MW <- reclassify(MW, cbind(NA, 0))
MW_dist <- raster(paste0(raw, 'MW_kde_dist.tif'))
PO <- raster(paste0(raw, 'PO_kde.tif')) 
PO <- reclassify(PO, cbind(NA, 0))
PO_dist <- raster(paste0(raw, 'PO_kde_dist.tif'))
QB <- raster(paste0(raw, 'QB_kde.tif')) 
QB <- reclassify(QB, cbind(NA, 0))
QB_dist <- raster(paste0(raw, 'QB_kde_dist.tif'))
SA <- raster(paste0(raw, 'SA_kde.tif')) 
SA <- reclassify(SA, cbind(NA, 0))
SA_dist <- raster(paste0(raw, 'SA_kde_dist.tif'))
TU <- raster(paste0(raw, 'TU_kde.tif')) 
TU <- reclassify(TU, cbind(NA, 0))
TU_dist <- raster(paste0(raw, 'TU_kde_dist.tif'))






#### making steps from track and extracting covariates ####
ssf <- dat_all %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = minutes(120), tolerance = minutes(10)) %>%
      amt::filter_min_n_burst(min_n = 3) %>%
      amt::steps_by_burst() %>% amt::random_steps(n=10) %>%
      amt::extract_covariates(land, where = "both")  %>%
      amt::extract_covariates(roads, where = "both") %>%
      amt::extract_covariates(FM, where = "both") %>%
      amt::extract_covariates(FM_dist, where = "both") %>%
      amt::extract_covariates(GF, where = "both") %>%
      amt::extract_covariates(GF_dist, where = "both") %>%
      amt::extract_covariates(GR, where = "both") %>%
      amt::extract_covariates(GR_dist, where = "both") %>%
      amt::extract_covariates(HL, where = "both") %>%
      amt::extract_covariates(HL_dist, where = "both") %>%
      amt::extract_covariates(MA, where = "both") %>%
      amt::extract_covariates(MA_dist, where = "both") %>%
      amt::extract_covariates(MW, where = "both") %>%
      amt::extract_covariates(MW_dist, where = "both") %>%
      amt::extract_covariates(PO, where = "both") %>%
      amt::extract_covariates(PO_dist, where = "both") %>%
      amt::extract_covariates(QB, where = "both") %>%
      amt::extract_covariates(QB_dist, where = "both") %>%
      amt::extract_covariates(SA, where = "both") %>%
      amt::extract_covariates(SA_dist, where = "both") %>%
      amt::extract_covariates(TU, where = "both") %>%
      amt::extract_covariates(TU_dist, where = "both") %>%
      mutate(land_start = factor(GHA26landcover2015_wgs84_start, levels = 1:7, labels = c("coniferous", 'deciduous', "mixed", 'shrub', "open", 'wet', 'urban')),
             land_end = factor(GHA26landcover2015_wgs84_end, levels = 1:7, labels = c("coniferous", 'deciduous', "mixed", 'shrub', "open", 'wet', 'urban')),
             FM_kde_start = factor(FM_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             FM_kde_end = factor(FM_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             GF_kde_start = factor(GF_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             GF_kde_end = factor(GF_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             GR_kde_start = factor(GR_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             GR_kde_end = factor(GR_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             HL_kde_start = factor(HL_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             HL_kde_end = factor(HL_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             MA_kde_start = factor(MA_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             MA_kde_end = factor(MA_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             MW_kde_start = factor(MW_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             MW_kde_end = factor(MW_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             PO_kde_start = factor(PO_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             PO_kde_end = factor(PO_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             QB_kde_start = factor(QB_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             QB_kde_end = factor(QB_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             SA_kde_start = factor(SA_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             SA_kde_end = factor(SA_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             TU_kde_start = factor(TU_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             TU_kde_end = factor(TU_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             log_sl = log(sl_ +1),
             cos_ta = cos(ta_))
  }))



ssf.all <- ssf %>% dplyr::select(id, steps) %>% unnest(cols = c(steps))


#### movement parmeters ####
SLdistr <- function(x.col, y.col, date.col, crs, ID, NumbRandSteps, sl_distr, ta_distr) {
  print(ID)
  #create track from dataset
  trk <- track(x.col, y.col, date.col, ID, crs) %>%
    #function turns locs into steps
    steps()
  #remove any steps that span more than 2hr
  trk$dt_ <- difftime(trk$t2_, trk$t1_, unit='hours')
  trk <- subset(trk, trk$dt_ > 1.9 & trk$dt_ < 2.1, drop = T)
  #generate random steps
  trk %>%
    random_steps() %>%
    sl_distr_params()
}

TAdistr <- function(x.col, y.col, date.col, crs, ID, NumbRandSteps, sl_distr, ta_distr) {
  print(ID)
  #create track from dataset
  trk <- track(x.col, y.col, date.col, ID, crs) %>%
    #function turns locs into steps
    steps()
  #remove any steps that span more than 2hr
  trk$dt_ <- difftime(trk$t2_, trk$t1_, unit='hours')
  trk <- subset(trk, trk$dt_ > 1.9 & trk$dt_ < 2.1, drop = T)
  #generate random steps
  trk %>%
    random_steps() %>%
    ta_distr_params()
}


# 
# DT.prep <- merge(dat.all[WolfID %in% focals],dat.meta, by.x = c('WolfID','PackID'), 
#                  by.y = c('WolfID','PackID'), all.x = T)
# DT.prep <- DT.prep[,.(x = X, y = Y, t = datetime, id = WolfID, status, COD, death_date)]
# DT.prep[, ttd:=(as.duration(t %--% death_date)/ddays(1))]


slParams <- DT.prep[, SLdistr(x.col = x, y.col = y, date.col = t, crs = utm14N, ID = id, 
                                     sl_distr = "gamma", ta_distr = "vonmises"),
                    by = id]

taParams <- DT.prep[, TAdistr(x.col = x, y.col = y, date.col = t, crs = utm14N, ID = id, 
                                     sl_distr = "gamma", ta_distr = "vonmises"),
                    by = id]

Params <- merge(slParams, taParams[,.(id,kappa)], by = 'id')



## proportions didn't pull right because of layer name, don't know how to fix
locs_start <- sp::SpatialPoints(data.frame(ssf.all$x1_, ssf.all$y1_))
locs_end <- sp::SpatialPoints(data.frame(ssf.all$x2_, ssf.all$y2_))

# ssf.all[,'propwet_start'] <- raster::extract(propwet, locs_start)
# #ssf.all[,'propclosed_start'] <- raster::extract(propclosed, locs_start)
# ssf.all[,'propconif_start'] <- raster::extract(propconif, locs_start)
# ssf.all[,'propmixed_start'] <- raster::extract(propmixed, locs_start)
# ssf.all[,'propopen_start'] <- raster::extract(propopen, locs_start)

ssf.all[,'propwet_end'] <- raster::extract(propwet, locs_end)
#ssf.all[,'propclosed_end'] <- raster::extract(propclosed, locs_end)
ssf.all[,'propconif_end'] <- raster::extract(propconif, locs_end)
ssf.all[,'propdecid_end'] <- raster::extract(propdecid, locs_end)
ssf.all[,'propmixed_end'] <- raster::extract(propmixed, locs_end)
ssf.all[,'propshrub_end'] <- raster::extract(propshrub, locs_end)
ssf.all[,'propopen_end'] <- raster::extract(propopen, locs_end)
ssf.all[,'propurban_end'] <- raster::extract(propurban, locs_end)

# adding nn and dist for step 1, also adding attributed data about death
ssf.all <- merge(ssf.all, DT.prep, by.x = c('x1_', 'y1_', 't1_', 'id'), by.y = c('x', 'y', 't', 'id'), all.x = T)
colnames(ssf.all)[colnames(ssf.all)=="nn"] <- "nn1"
#colnames(ssf.all)[colnames(ssf.all)=="end_date"] <- "death_date"

# adding nn at step 2
DT.nn <- dplyr::select(DT.prep, t, id, nn)
ssf.all <- merge(ssf.all, DT.nn, by.x = c('t2_', 'id'), by.y = c('t', 'id'), all.x = T)
colnames(ssf.all)[colnames(ssf.all)=="nn"] <- "nn2"

# calc time to death (ttd)
ssf.all[,'ttd1'] <- as.duration(ssf.all$t1_ %--% ssf.all$death_date)/ddays(1) 
ssf.all[,'ttd2'] <- as.duration(ssf.all$t2_  %--% ssf.all$death_date)/ddays(1) 


colnames(ssf.all)[colnames(ssf.all)=="GHA26_roadsPS_dist_start"] <- "roadDist_start"
colnames(ssf.all)[colnames(ssf.all)=="GHA26_roadsPS_dist_end"] <- "roadDist_end"

#saveRDS(ssf.all, 'data/derived-data/ssfRaw_GHA26.Rds')




# Cleaning of grouptime data
dat.grp <- dplyr::select(dat.grptimes, WolfID, PackID, X, Y, datetime, timegroup)
dat.grp[,'round.t'] <- as.POSIXct(round(dat.grp$datetime, units = 'mins'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")

# getting a list of all timegroups and associated times to help match points correctly
DT.grp <- dplyr::select(dat.grp, round.t, timegroup)
DT.grp2 <- unique(DT.grp)

DT.pack <- setDT(dat.focal)[,.(WolfID, packbound)]


#ssf.df = as.data.frame(ssf.all)
#wolf = 'W01'

createSSFnnbyFocal <- function(ssf.df, wolf){
  ssf.sub <- as.data.frame(subset(ssf.df, id == wolf))
  pack <- dat.focal[WolfID == wolf, .(packbound)]$packbound
  pack_kde_start <-paste(pack,"kde_start", sep = '_')
  
  ssf.sub[,'packYN_start'] <- ssf.sub[,(paste(pack,"kde_start", sep = '_'))]
  ssf.sub[,'packYN_end'] <- ssf.sub[,(paste(pack,"kde_end", sep = '_'))]
  
  ssf.sub[,'packDist_start'] <- ssf.sub[,paste(pack,"kde_dist_start", sep = '_')]
  ssf.sub[,'packDist_end'] <- ssf.sub[,paste(pack,"kde_dist_end", sep = '_')]
  
  ssf.sub[,'round.t2'] <- as.POSIXct(round(ssf.sub$t2_, units = 'mins'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")
  dat.rand.sub <- merge(ssf.sub, dat.grp, by.y = c('round.t', 'WolfID'), 
                        by.x = c('round.t2', 'id'), all =T)
  # use the gps pts from the track unless missing (as with all non focal indivs)
  dat.rand.sub$x2_ <- ifelse(!is.na(dat.rand.sub$x2_), dat.rand.sub$x2_, dat.rand.sub$X)
  dat.rand.sub$y2_ <- ifelse(!is.na(dat.rand.sub$y2_), dat.rand.sub$y2_, dat.rand.sub$Y)

  # enumerating each random step so can pair distances correctly
  DT.rand.sub <- setDT(dat.rand.sub)[, step_id_rank:= 1:.N, by=.(step_id_)]
  # making all regular steps (without random steps) have a value of 1
  DT.rand.sub$step_id_rank <- ifelse(!is.na(DT.rand.sub$step_id_), DT.rand.sub$step_id_rank, 1)
  
  # merging id with info that will help us pair distances with the correct steps and individuals after they're calculated
  # edge_dist() doesn't output a lot of the info needed to do this, so will separated it again when done
  DT.rand.sub[,'id.step'] <- paste(DT.rand.sub$id, DT.rand.sub$step_id_rank, sep = '.')
  DT.rand.sub[,'id.step.t'] <- paste(DT.rand.sub$id.step, DT.rand.sub$timegroup, sep = '.')
  
  
  # calculates distances for all individuals within the same pack at each timegroup (this is just for step 2 right now)
  
  #DT.rand.w02.dates <- DT.rand.w02[round.t2>=death_date-ddays(ttd) & round.t2 <= death_date]
  lower.date <- dat.focal[WolfID== wolf, death_date-ddays(ttd)]
  upper.date <- dat.focal[WolfID== wolf, death_date]
  edist2.sub <- DT.rand.sub[round.t2>=lower.date & round.t2 <= upper.date, edge_dist(
    DT = .SD,
    threshold = 1000000, # I don't know how to pick this
    id = 'id.step.t',
    coords = c('x2_', 'y2_'),
    timegroup = 'timegroup',
    returnDist = TRUE,
    splitBy = 'PackID'
  )]
  
  edist2.sub.sep <- tidyr::separate(edist2.sub, col = ID1, into = c('id1', 'step1', 'timegroup1'))
  edist2.sub.sep <- tidyr::separate(edist2.sub.sep, col = ID2, into = c('id2', 'step2', 'timegroup2'))
  
  #had to make step an integer so that they match between data frames (don't know why step1 and step2 think they're chr)
  edist2.sub.sep[,'step'] <- as.integer(edist2.sub.sep$step1)
  
  # Keeping only the necessary columns
  edist2.sub.sep<- edist2.sub.sep[,.(timegroup, step, id1, id2, distance)]
  
  # numbering the steps and random steps to help match the social partners and their distances
  ssf.soc.sub <- setDT(ssf.sub)[, step_id_rank:= 1:.N, by=.(step_id_)]
  ssf.soc.sub<- merge(ssf.soc.sub, DT.grp2, by.x = 'round.t2', by.y = 'round.t', all.x = T)
  colnames(ssf.soc.sub)[colnames(ssf.soc.sub)=="timegroup"] <- "timegroup2"
  ssf.soc.sub[,'round.t1'] <- as.POSIXct(round(ssf.soc.sub$t1_, units = 'mins'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")
  ssf.soc.sub<- merge(ssf.soc.sub, DT.grp2, by.x = 'round.t1', by.y = 'round.t', all.x = T)
  colnames(ssf.soc.sub)[colnames(ssf.soc.sub)=="timegroup"] <- "timegroup1"
  ssf.soc.sub <- ssf.soc.sub[ttd1 <= ttd]
  
  ssf.soc.sub <- merge(ssf.soc.sub, edist2.sub.sep, by.x = c('id','nn2', 'step_id_rank', 'timegroup2'),
                       by.y = c('id1','id2', 'step', 'timegroup'), all.x = T)
  colnames(ssf.soc.sub)[colnames(ssf.soc.sub)=="distance"] <- "distance2"
  
  ssf.soc.sub[,'wet_end'] <- ifelse(ssf.soc.sub$land_end =='wet', 1,0)
  ssf.soc.sub[,'open_end'] <- ifelse(ssf.soc.sub$land_end =='open', 1,0)
  #ssf.soc.sub[,'closed_end'] <- ifelse(ssf.soc.sub$land_end =='closed', 1,0)
  ssf.soc.sub[,'conif_end'] <- ifelse(ssf.soc.sub$land_end =='coniferous', 1,0)
  ssf.soc.sub[,'mixed_end'] <- ifelse(ssf.soc.sub$land_end =='mixed', 1,0)
  ssf.soc.sub[,'decid_end'] <- ifelse(ssf.soc.sub$land_end =='deciduous', 1,0)
  ssf.soc.sub[,'shrub_end'] <- ifelse(ssf.soc.sub$land_end =='shrub', 1,0)
  ssf.soc.sub[,'urban_end'] <- ifelse(ssf.soc.sub$land_end =='urban', 1,0)
  ssf.soc.sub[,'packDistadj_end'] <-ifelse(ssf.soc.sub$packYN_end == 'pack', ssf.soc.sub$packDist_end, (-1)*ssf.soc.sub$packDist_end)
  
  
  
  
  ssf.wolf <- ssf.soc.sub[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta,# tod_start_, 
                            roadDist_start, roadDist_end,
                            land_start, land_end, propwet_end, propopen_end, propconif_end, propmixed_end, propdecid_end, propshrub_end, propurban_end,
                            wet_end, open_end, conif_end, mixed_end, decid_end, shrub_end, urban_end,
                            id, nn1, nn2, distance1, distance2, timegroup1, timegroup2, packYN_start, packYN_end, packDist_start, packDist_end, packDistadj_end,
                            ttd1, ttd2)]
  return(ssf.wolf)
}

unique(sort(dat.focal$WolfID))

unique(setDT(ssf.all)[!is.na(distance1), .(id)])


ssfW01 <- createSSFnnbyFocal(ssf.all, "W01")
ssfW03 <- createSSFnnbyFocal(ssf.all, "W03")
ssfW04 <- createSSFnnbyFocal(ssf.all, "W04")
ssfW05 <- createSSFnnbyFocal(ssf.all, "W05")
ssfW06 <- createSSFnnbyFocal(ssf.all, "W06")
ssfW09 <- createSSFnnbyFocal(ssf.all, "W09")
ssfW10 <- createSSFnnbyFocal(ssf.all, "W10")
ssfW11 <- createSSFnnbyFocal(ssf.all, "W11")
ssfW12 <- createSSFnnbyFocal(ssf.all, "W12")
ssfW13 <- createSSFnnbyFocal(ssf.all, "W13")
ssfW14 <- createSSFnnbyFocal(ssf.all, "W14")
ssfW15 <- createSSFnnbyFocal(ssf.all, "W15")
ssfW16 <- createSSFnnbyFocal(ssf.all, "W16")
#ssfW18 <- createSSFnnbyFocal(ssf.all, "W18")
ssfW22 <- createSSFnnbyFocal(ssf.all, "W22")
ssfW23 <- createSSFnnbyFocal(ssf.all, "W23")
ssfW24 <- createSSFnnbyFocal(ssf.all, "W24")
ssfW25 <- createSSFnnbyFocal(ssf.all, "W25")
ssfW26 <- createSSFnnbyFocal(ssf.all, "W26")
ssfW27 <- createSSFnnbyFocal(ssf.all, "W27")
ssfW29 <- createSSFnnbyFocal(ssf.all, "W29")
ssfW31 <- createSSFnnbyFocal(ssf.all, "W31")
ssfW32 <- createSSFnnbyFocal(ssf.all, "W32")
ssfW34 <- createSSFnnbyFocal(ssf.all, "W34")
ssfW35 <- createSSFnnbyFocal(ssf.all, "W35")
ssfW36 <- createSSFnnbyFocal(ssf.all, "W36")
ssfW37 <- createSSFnnbyFocal(ssf.all, "W37")
ssfW38 <- createSSFnnbyFocal(ssf.all, "W38")
ssfW39 <- createSSFnnbyFocal(ssf.all, "W39")


ssf.soc <- rbind(ssfW01, ssfW03, ssfW04, ssfW05, ssfW06, ssfW09, ssfW10, ssfW11, ssfW12, ssfW13, 
                 ssfW14, ssfW15, ssfW16, ssfW22, ssfW23, ssfW24, ssfW25, ssfW26, ssfW27, ssfW29, ssfW31, 
                 ssfW32, ssfW34, ssfW35, ssfW36, ssfW37, ssfW38, ssfW39)
ssf.soc <- merge(ssf.soc, dat.focal[,.(WolfID, PackID, COD)], by.x = 'id', by.y = 'WolfID', all.x = T)


saveRDS(ssf.soc, 'data/derived-data/ssfAll_2mo_GHA26.Rds')


moveRMNP <- readRDS('data/derived-data/moveParams_2mo_RMNP.Rds')
moveRMNP[,wolfID := paste('RMNP', id, sep = '_')]
Params[,wolfID := paste('GHA26', id, sep = '_')]
Params <- rbind(moveRMNP, Params)
saveRDS(Params, 'data/derived-data/moveParams_2mo_all.Rds')



