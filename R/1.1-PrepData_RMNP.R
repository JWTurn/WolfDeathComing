### Data ====
# Julie Turner
# Started: July 1 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'spatsoc', 'dplyr', 'amt', 'lubridate', 'raster', 'tidyr', 'ggplot2')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolf_rarified.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")
dat[,.(min=min(gmtDate),max=max(gmtDate)), by=.(WolfID)]

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'))
dat.meta$death_date <- paste(dat.meta$death_date, '23:59:59', sep = ' ')
dat.meta$death_date<- as.POSIXct(dat.meta$death_date, tz = 'UTC', "%Y-%m-%d %H:%M:%S")

# range <- dat[,.(min=min(datetime), max=max(datetime)), by=.(WolfID)]  
# range <- range[,.(min, max, len=max-min), .(WolfID)]

#### determining nearest neighbor (nn) at orginal point ####
dat.grptimes <- group_times(DT=dat, datetime = "datetime", threshold = "2 hours")

enn <- edge_nn(
  DT = dat.grptimes,
  id = 'WolfID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  returnDist = TRUE, 
  splitBy = 'PackID'
)



dat.nn <- merge(dat.grptimes,enn, by.x = c('WolfID','PackID', 'timegroup'), 
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

set.seed(53)
#### setting time to death ####
ttd = 61

#### filter to those that die ####
dat.focal <- setDT(dat.meta)[pop=='RMNP' & use!='n' & !is.na(PackID)]
dat.focal <- setnames(dat.focal, old = 'PackID', new = 'packbound')
dat.focal <- dat.focal[!(WolfID %chin% c('W08'))] #doesn't have enough data
focals <- dat.focal$WolfID

# DT.prep <- dat.nn %>% dplyr::select(x = "X", y = "Y", t = 'datetime', id = "WolfID", nn = 'NN', distance1 = 'distance',
#                                     'status', 'end_date', 'COD', 'death_date') %>%
#   filter(id %in% focals) 
#   
  
# DT.prep <- merge(dat[WolfID %in% focals],dat.meta, by.x = c('WolfID','PackID'),
#                  by.y = c('WolfID','PackID'), all.x = T)
DT.prep <- dat.nn[WolfID %in% focals,.(x = X, y = Y, t = datetime, id = WolfID, nn = NN, distance1 = distance,
                     status, end_date, COD, death_date)]
#DT.prep[, t2d:=(as.duration(t %--% death_date)/ddays(1))]
#DT.prep <- DT.prep[t2d >=0 & t2d<=ttd]

DT.prep[,uniqueN(t), by =.(id)]

dat_all <- DT.prep %>% group_by(id) %>% nest()

dat_all <- dat_all %>%
  mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:32614")) 
  }))  



dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest(cols = c(sr))


#### layers ####

land <- raster(paste0(raw, 'RMNPlandcover2015_wgs84.tif'), resolution = c(30, 30))
# land <- raster::projectRaster(from=prepland, crs = crs14)
# res(prepland) <- c(30, 30)
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
parkYN <- raster(paste0(raw, 'parkYN.tif'))
parkbound_dist <- raster(paste0(raw, 'boundary_dist.tif'))
roads <- raster(paste0(raw, 'roadall_LinFeat_dist.tif'))
AD <- raster(paste0(raw, 'AD_kde.tif'))
AD <- reclassify(AD, cbind(NA, 0))
AD_dist <- raster(paste0(raw, 'AD_kde_dist.tif'))
BD <- raster(paste0(raw, 'BD_kde.tif'))
BD <- reclassify(BD, cbind(NA, 0))
BD_dist <- raster(paste0(raw, 'BD_kde_dist.tif'))
BL <- raster(paste0(raw, 'BL_kde.tif'))
BL <- reclassify(BL, cbind(NA, 0))
BL_dist <- raster(paste0(raw, 'BL_kde_dist.tif'))
BT <- raster(paste0(raw, 'BT_kde.tif'))
BT <- reclassify(BT, cbind(NA, 0))
BT_dist <- raster(paste0(raw, 'BT_kde_dist.tif'))
GL <- raster(paste0(raw, 'GL_kde.tif'))
GL <- reclassify(GL, cbind(NA, 0))
GL_dist <- raster(paste0(raw, 'GL_kde_dist.tif'))
RC <- raster(paste0(raw, 'RC_kde.tif'))
RC <- reclassify(RC, cbind(NA, 0))
RC_dist <- raster(paste0(raw, 'RC_kde_dist.tif'))
SL <- raster(paste0(raw, 'SL_kde.tif'))
SL <- reclassify(SL, cbind(NA, 0))
SL_dist <- raster(paste0(raw, 'SL_kde_dist.tif'))
WW <- raster(paste0(raw, 'WW_kde.tif'))
WW <- reclassify(WW, cbind(NA, 0))
WW_dist <- raster(paste0(raw, 'WW_kde_dist.tif'))






#### making steps from track and extracting covariates ####
ssf <- dat_all %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = minutes(120), tolerance = minutes(10)) %>%
      amt::filter_min_n_burst(min_n = 3) %>%
      amt::steps_by_burst() %>% amt::random_steps(n=10) %>%
      amt::extract_covariates(land, where = "both")  %>%
      amt::extract_covariates(parkYN, where = "both") %>%
      amt::extract_covariates(parkbound_dist, where = "both") %>%
      amt::extract_covariates(roads, where = "both") %>%
      amt::extract_covariates(AD, where = "both") %>%
      amt::extract_covariates(AD_dist, where = "both") %>%
      amt::extract_covariates(BD, where = "both") %>%
      amt::extract_covariates(BD_dist, where = "both") %>%
      amt::extract_covariates(BL, where = "both") %>%
      amt::extract_covariates(BL_dist, where = "both") %>%
      amt::extract_covariates(BT, where = "both") %>%
      amt::extract_covariates(BT_dist, where = "both") %>%
      amt::extract_covariates(GL, where = "both") %>%
      amt::extract_covariates(GL_dist, where = "both") %>%
      amt::extract_covariates(RC, where = "both") %>%
      amt::extract_covariates(RC_dist, where = "both") %>%
      amt::extract_covariates(SL, where = "both") %>%
      amt::extract_covariates(SL_dist, where = "both") %>%
      amt::extract_covariates(WW, where = "both") %>%
      amt::extract_covariates(WW_dist, where = "both") %>%
      mutate(land_start = factor(RMNPlandcover2015_wgs84_start, levels = 1:7, labels = c("coniferous", 'deciduous', "mixed", 'shrub', "open", 'wet', 'urban')),
             land_end = factor(RMNPlandcover2015_wgs84_end, levels = 1:7, labels = c("coniferous", 'deciduous', "mixed", 'shrub', "open", 'wet', 'urban')),
             parkYN_start = factor(parkYN_start, levels = c(0, 1), labels = c("park", "out-park")),
             parkYN_end = factor(parkYN_end, levels = c(0, 1), labels = c("park", "out-park")),
             AD_kde_start = factor(AD_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             AD_kde_end = factor(AD_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             BD_kde_start = factor(BD_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             BD_kde_end = factor(BD_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             BL_kde_start = factor(BL_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             BL_kde_end = factor(BL_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             BT_kde_start = factor(BT_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             BT_kde_end = factor(BT_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             GL_kde_start = factor(GL_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             GL_kde_end = factor(GL_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             RC_kde_start = factor(RC_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             RC_kde_end = factor(RC_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             SL_kde_start = factor(SL_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             SL_kde_end = factor(SL_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             WW_kde_start = factor(WW_kde_start, levels = c(1, 0), labels = c("pack", "out-pack")),
             WW_kde_end = factor(WW_kde_end, levels = c(1, 0), labels = c("pack", "out-pack")),
             log_sl = log(sl_ + 1),
             cos_ta = cos(ta_))
  }))



ssf.all <- ssf %>% dplyr::select(id, steps) %>% unnest(cols = c(steps))

# dat_all <- dat_all %>%
#   mutate(steps = map(trk, function(x) {
#     x %>% amt::track_resample(rate = minutes(120), tolerance = minutes(10)) %>%
#       amt::filter_min_n_burst(min_n = 3) %>%
#       amt::steps_by_burst() %>% amt::random_steps(n=10)
#   }))

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






#badwolf <- c('W07')
badwolf <- c('W16')
slParams <- DT.prep[!(id %in% badwolf), SLdistr(x.col = x, y.col = y, date.col = t, crs = utm14N, ID = id, 
                                                sl_distr = "gamma", ta_distr = "vonmises"),
                    by = id]

taParams <- DT.prep[!(id %in% badwolf), TAdistr(x.col = x, y.col = y, date.col = t, crs = utm14N, ID = id, 
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
ssf.all[,'nn2'] <- ifelse(!(is.na(ssf.all$nn2)), ssf.all$nn2, ifelse(!(is.na(ssf.all$nn1)), ssf.all$nn1, NA))

# calc time to death (ttd)
ssf.all[,'ttd1'] <- as.duration(ssf.all$t1_ %--% ssf.all$death_date)/ddays(1) 
ssf.all[,'ttd2'] <- as.duration(ssf.all$t2_  %--% ssf.all$death_date)/ddays(1) 


colnames(ssf.all)[colnames(ssf.all)=="roadall_LinFeat_dist_start"] <- "roadDist_start"
colnames(ssf.all)[colnames(ssf.all)=="roadall_LinFeat_dist_end"] <- "roadDist_end"

#saveRDS(ssf.all, 'data/derived-data/ssfRaw_RMNP.Rds')




# Cleaning of grouptime data
dat.grp <- dplyr::select(dat.grptimes, WolfID, PackID, X, Y, datetime, timegroup)
dat.grp[,'round.t'] <- as.POSIXct(round(dat.grp$datetime, units = 'mins'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")

# getting a list of all timegroups and associated times to help match points correctly
DT.grp <- dplyr::select(dat.grp, round.t, timegroup)
DT.grp2 <- unique(DT.grp)

DT.pack <- setDT(dat.focal)[,.(WolfID, packbound)]


# ssf.df = ssf.all
# wolf = 'W02'

createSSFnnbyFocal <- function(ssf.df, wolf){
  ssf.sub <- subset(ssf.df, id == wolf)
  pack <- dat.focal[WolfID == wolf, .(packbound)]$packbound
  ssf.sub[,'pack'] <- pack
  
  ssf.sub[,'packYN_start'] <- as.data.frame(ssf.sub)[,paste(pack,"kde_start", sep = '_')]
  ssf.sub[,'packYN_end'] <- as.data.frame(ssf.sub)[,paste(pack,"kde_end", sep = '_')]
  
  ssf.sub[,'packDist_start'] <- as.data.frame(ssf.sub)[,paste(pack,"kde_dist_start", sep = '_')]
  ssf.sub[,'packDist_end'] <- as.data.frame(ssf.sub)[,paste(pack,"kde_dist_end", sep = '_')]
  
  ssf.sub[,'round.t2'] <- as.POSIXct(round(ssf.sub$t2_, units = 'mins'), tz = 'UTC', "%Y-%m-%d %H:%M:%S")
  dat.rand.sub <- merge(ssf.sub, dat.grp[PackID == pack], by.y = c('round.t', 'WolfID'), 
                        by.x = c('round.t2', 'id'), all = T)
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
    threshold = 10000000, # I don't know how to pick this
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
  ssf.soc.sub[,'packDistadj_end'] <-ifelse(ssf.soc.sub$packYN_end == 'pack', ssf.soc.sub$packDist_end, ssf.soc.sub$packDist_end*(-1))
  ssf.soc.sub[,'parkDistadj_end'] <-ifelse(ssf.soc.sub$parkYN_end == 'park', ssf.soc.sub$boundary_dist_end, (-1)*ssf.soc.sub$boundary_dist_end)
  
  
  
  ssf.wolf <- ssf.soc.sub[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta, #tod_start_, 
                            parkYN_start, parkYN_end, roadDist_start, roadDist_end, parkDist_end = boundary_dist_end, parkDistadj_end,
                            land_start, land_end, propwet_end, propopen_end, propconif_end, propmixed_end, propdecid_end, propshrub_end, propurban_end,
                            wet_end, open_end, conif_end, mixed_end, decid_end, shrub_end, urban_end,
                            id, nn1, nn2, distance1, distance2, timegroup1, timegroup2, packYN_start, packYN_end, packDist_start, packDist_end, packDistadj_end,
                            ttd1, ttd2)]
  return(ssf.wolf)
}

unique(sort(dat.focal$WolfID))
ssfW02 <- createSSFnnbyFocal(ssf.all, "W02")
ssfW03 <- createSSFnnbyFocal(ssf.all, "W03")
ssfW04 <- createSSFnnbyFocal(ssf.all, "W04")
ssfW05 <- createSSFnnbyFocal(ssf.all, "W05")
ssfW06 <- createSSFnnbyFocal(ssf.all, "W06")
ssfW07 <- createSSFnnbyFocal(ssf.all, "W07") 
ssfW09 <- createSSFnnbyFocal(ssf.all, "W09")
ssfW10 <- createSSFnnbyFocal(ssf.all, "W10")
ssfW11 <- createSSFnnbyFocal(ssf.all, "W11")
ssfW12 <- createSSFnnbyFocal(ssf.all, "W12")
ssfW14 <- createSSFnnbyFocal(ssf.all, "W14")
ssfW15 <- createSSFnnbyFocal(ssf.all, "W15")
ssfW16 <- createSSFnnbyFocal(ssf.all, "W16") 
ssfW19 <- createSSFnnbyFocal(ssf.all, "W19")
ssfW20 <- createSSFnnbyFocal(ssf.all, "W20")
ssfW22 <- createSSFnnbyFocal(ssf.all, "W22")
ssfW25 <- createSSFnnbyFocal(ssf.all, "W25")
ssfW26 <- createSSFnnbyFocal(ssf.all, "W26")
ssfW27 <- createSSFnnbyFocal(ssf.all, "W27")

# , ssfW20 removed from list, didn't run
ssf.soc <- rbind(ssfW02, ssfW03, ssfW04, ssfW05, ssfW06, ssfW07, ssfW09, ssfW10, ssfW11, ssfW12, ssfW14, ssfW15, ssfW16, ssfW19, ssfW20, ssfW22, ssfW25, ssfW26, ssfW27)
ssf.soc <- merge(ssf.soc, dat.focal[,.(WolfID, PackID=packbound, COD)], by.x = 'id', by.y = 'WolfID', all.x = T)



saveRDS(ssf.soc, 'data/derived-data/ssfAll.Rds')

#saveRDS(Params, 'data/derived-data/moveParams_2mo_RMNP.Rds')



