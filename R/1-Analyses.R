### Analyses ====
# Julie Turner
# Started: July 3 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival','forcats')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### SSF data with all covariates ----
set.seed(53)
dat <- readRDS('data/derived-data/ssfAllCov.Rds')

dat <- dat[, stepjum := forcats::fct_shuffle(as.character(step_id_)), by = .(id)]

dat[,.(nstep = uniqueN(as.character(step_id_)), nstepjum = uniqueN(stepjum)), by = .(id)]


#### adding prop landcover, but move to data prep ####
# Generate spatial points object from UTM coordinates
locs <- SpatialPoints(data.frame(dat$x1_, dat$y1_))
locs_end <- SpatialPoints(data.frame(dat$x2_, dat$y2_))

# Read in the landcover map with different habitat types
lc <-  raster(paste0(raw, 'RMNPlandcover.tif'), )
cland <- fread(paste0(raw, 'rcl_cowu.csv'))
lc <- raster::reclassify(lc, cland)
# To get rid of NAs, I make them a new category (10), and then add it to the legend
lc[is.na(lc)] <- 5

closed <- lc == 1
names(closed) <- "closed"
open <- lc == 2
names(open) <- "open"
wet <- lc == 3
names(wet) <- "wet"

## This creates an object which can be used to make a layer of specified diameter
# The d value is what determines the buffer size if you want to change it.
## If you're doing multiple landcover classes, you only need to run this line once, as long as each of the habitat variables has the same resolution
# (e.g., the "Wetland" here could be any cover type)

Buff100 <- focalWeight(closed, d=100, type = 'circle')

## This generates a new raster where each cell corresponds to the mean wetland within a 100m buffer.
# Since it's all 1s and 0s, this is the same as the proportion of wetland surrounding the focal variable
wetBuff100 <- focal(wet, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
closedBuff100 <- focal(closed, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
openBuff100 <- focal(open, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)

## You can then use the "extract" function as you would on any habitat layer and add it to your data:
dat[,'propwet_start'] <- extract(wetBuff100, locs)
dat[,'propclosed_start'] <- extract(closedBuff100, locs)
dat[,'propopen_start'] <- extract(openBuff100, locs)

dat[,'propwet_end'] <- extract(wetBuff100, locs_end)
dat[,'propclosed_end'] <- extract(closedBuff100, locs_end)
dat[,'propopen_end'] <- extract(openBuff100, locs_end)

#saveRDS(dat, 'data/derived-data/ssfAllCov.Rds')

# dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')
# dat[,'wet'] <- ifelse(dat$land_end == 'wet', 'wet', 'not')
# dat.wet<-dat[land_end=='wet', .(id, ttd1, ua, propwet_end)]
# dat.wet[, .(id, .N), by= c('ua', 'id')]
# 
# ggplot(dat, aes(x=factor(ttd1),y=wet, group=(wet)))+
#   stat_summary(aes(color=(wet)),fun.y=length, geom="line")+
#   scale_color_discrete("wet",labels=c("not","wet"))+
#   labs(x="",y="Frequency")

dat.wet<-dat.wet[,dtd := as.integer(ttd1)]
dat.wet<-dat.wet[ua=='used',wetpd:= uniqueN(ttd1), by =.(id,dtd)]
dat.wet.used<-unique(dat.wet[ua=='used', .(id,dtd,wetpd, propwet_end)])

ggplot(dat.wet.used, aes(x=factor(dtd),y=(wetpd), group=(id)))+
     labs(x="",y="Frequency")

lattice::xyplot(wetpd ~ dtd, data = dat.wet.used, groups = id, type ='l')
lattice::xyplot(propwet_end ~ dtd, data = dat.wet.used, groups = id, type ='l')

dat.meta <- fread(paste0(raw, 'wolf_metadata.csv'))

dat.wet[,cor(propwet_end, ttd1), by=.(id)]

#### fit models ####

### questions
# if I'm using ttd1, do I need all interactions to be with the other start points?  --- no
# do I need to log transform distance to conspecifics? or ttd? (looks like need to do for ttd) --- yes
# do I need cos_ta or cos_ta:tod --- not necessarily

# dat.w06 <- dat[id=='W06']
# 
# m.W06 <- clogit(case_ ~ log_sl:ToD_start + land_end + log_sl:land_end + strata(step_id_), dat.w06)
# 
# sum.w06 <- summary(m.W06)$coefficients
# tab.sum.w06 <- sum.w06$coefficients
# tab.sum.w06[,3]
# 
# out <- data.table(t(tab.sum.w06))
# out[,'term']<-c('coef','hr','se','z','p')

#### CORE ####

Core <- function(y, sl, ToD, closed, open, wet, strata1) {
    # Make the model
    model <- clogit(y ~ sl:ToD + closed + open + wet + 
                      #sl:closed + sl:open + sl:wet + 
                      strata(strata1))
    sum.model <- summary(model)$coefficients
    # Transpose the coef of the model and cast as data.table
    term <- c('coef','hr','se','z','p')
    coefOut <- data.table(t(sum.model))

    # Return combined columns
    # print(summary(model))
    print(AIC(model))
    return(data.table(term, coefOut, AIC=AIC(model)))
  }

#run iSSA model by ID
# Core Model
coreOUT<- dat[, Core(case_, log_sl, ToD_start, propclosed_end, propopen_end, propwet_end, stepjum), by = id]
m.core <- merge(coreOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.core.p<-m.core[term=='coef' | term=='p']


m.full <- data.all %>% amt::fit_issf(case_ ~ log_sl:ToD_start + land_end + log_sl:land_end + 
                                      ttd1:log_sl + ttd1:cos_ta +
                                      ttd1:land_start + ttd1:lnparkdist_start:parkYN_start + 
                                      ttd1:distance2 + 
                                      strata(step_id_))
# full not converging based on 1 indiv

#### MOVEMENT ####

Move <- function(y, sl, ToD, land, ttd, ta, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    log(ttd+1):sl + log(ttd+1):ta + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

moveOUT <- dat[, Move(case_, log_sl, ToD_start, land_end, ttd1, cos_ta, stepjum), by = id]
m.move <- merge(moveOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.move.p<-m.move[term=='coef' | term=='p', .(id, `sl:ttd`, `ttd:ta`, COD)]



#### HABITAT ####

Habitat <- function(y, sl, ToD, land, ttd, parkdist, parkYN, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    log(ttd+1):land + log(ttd+1):parkdist:parkYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

habOUT <- dat[, Habitat(case_, log_sl, ToD_start, land_end, ttd1, lnparkdist_end, parkYN_end, stepjum), by = id]
m.habitat <- merge(habOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitat.p <- m.habitat[term=='coef'|term=='p', .(id, term, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
                       `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, COD)]

Habitat.land <- function(y, sl, ToD, land, ttd, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    log(ttd+1):land + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

hablandOUT <- dat[, Habitat.land(case_, log_sl, ToD_start, land_end, ttd1, stepjum), by = id]
m.habitat.land <- merge(hablandOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitat.land.p <- m.habitat.land[term=='p', .(id, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
                                      `log(ttd + 1):parkdist:parkYNpark`,`log(ttd + 1):parkdist:parkYNout-park`, COD)]

# didn't fix the problem with wet


Habitat.park <- function(y, sl, ToD, land, ttd, parkdist, parkYN, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    log(ttd+1):parkdist:parkYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

habparkOUT <- dat[, Habitat.park(case_, log_sl, ToD_start, land_end, ttd1, lnparkdist_end, parkYN_end, stepjum), by = id]
m.habitat.park <- merge(habparkOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitat.park.p <- m.habitat.park[term=='p', .(id, `log(ttd + 1):parkdist:parkYNpark`,`log(ttd + 1):parkdist:parkYNout-park`, COD)]


#### SOCIAL ####

m.social <- data.all %>% amt::fit_issf(case_ ~ log_sl:ToD_start + land_end + log_sl:land_end + 
                                        log(ttd1):log(distance2) + 
                                        strata(step_id_))

Social <- function(y, sl, ToD, land, ttd, nndist, packYN, packdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    log(ttd+1):log(nndist+1) + log(ttd+1):log(packdist+1):packYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

dat[id !='W13' & id != 'W11' & id != 'W19' & id != 'W15',unique(id)]
# W13, W11, W27,  not converge

socOUT <- dat[id !='W13' & id != 'W11' & id != 'W19' & id != 'W15', Social(case_, log_sl, ToD_start, land_end, ttd1, distance2, packYN_end, packDist_end, stepjum), by = id]
m.social <- merge(socOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.social.p <- m.social[term=='coef' |term=='p', .(id, term,`ttd:log(nndist + 1)`, `ttd:log(packdist + 1):packYNpack`, `ttd:log(packdist + 1):packYNout-pack`, COD)]



bbmle::AICtab(m.core$model, m.move$model, m.habitat.parkdist$model, m.social$model)


