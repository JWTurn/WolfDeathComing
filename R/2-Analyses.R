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
dat.RMNP <- readRDS('data/derived-data/ssfAllCov_RMNP.Rds')
dat.GHA26 <- readRDS('data/derived-data/ssfAllCov_GHA26.Rds')


dat.RMNP[,'pop'] <- 'RMNP'
dat.RMNP$wolfID <- paste(dat.RMNP$pop, dat.RMNP$id, sep = '_')
dat.GHA26[,'pop'] <- 'GHA26'
dat.GHA26$wolfID <- paste(dat.GHA26$pop, dat.GHA26$id, sep = '_')
dat<-rbind(dat.RMNP, dat.GHA26, fill=T)

dat <- dat[, stepjum := forcats::fct_shuffle(as.factor(step_id_)), by = .(wolfID)]


dat[,.(nstep = uniqueN(as.character(step_id_)), nstepjum = uniqueN(stepjum)), by = .(wolfID)]

#dat.steps <-dat[id=='W06',.(step_id_, stepjum)]

dat<-dat[ttd1>=0 & ttd2>=0]

#dat[,'wtd1'] <- as.integer(dat$ttd1/7)

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'))
dat.meta[,'wolfpop'] <- paste(dat.meta$pop, dat.meta$WolfID, sep = '_')



dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')

dat<- merge(dat, dat.meta, by.x = c('id', 'pop'), by.y = c('WolfID', 'pop'))


#dat[,'ttd1_adj'] <- ifelse(dat$COD =='none', dat$ttd1*730, dat$ttd1)

dat.meta[,uniqueN(wolfpop), by= .(pop, COD)]

# dat[,'wet'] <- ifelse(dat$land_end == 'wet', 'wet', 'not')
# dat.wet<-dat[land_end=='wet', .(id, ttd1, ua, propwet_end)]
# dat.wet[, .(id, .N), by= c('ua', 'id')]
# 
# ggplot(dat, aes(x=factor(ttd1),y=wet, group=(wet)))+
#   stat_summary(aes(color=(wet)),fun.y=length, geom="line")+
#   scale_color_discrete("wet",labels=c("not","wet"))+
#   labs(x="",y="Frequency")

# dat.wet<-dat.wet[,dtd := as.integer(ttd1)]
# dat.wet<-dat.wet[ua=='used',wetpd:= uniqueN(ttd1), by =.(id,dtd)]
# dat.wet.used<-unique(dat.wet[ua=='used', .(id,dtd,wetpd, propwet_end)])
# 
# ggplot(dat.wet.used, aes(x=factor(dtd),y=(wetpd), group=(id)))+
#      labs(x="",y="Frequency")
# 
# lattice::xyplot(wetpd ~ dtd, data = dat.wet.used, groups = id, type ='l')
# lattice::xyplot(propwet_end ~ dtd, data = dat.wet.used, groups = id, type ='l')
# 
# 
# dat.wet[ua=='used',cor(propwet_end, ttd1), by=.(id)]




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


### dummy vars didn't, Can run prop with this same one
Core.land <- function(y, sl, ToD, closed, open, wet, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + closed + open + wet +
                    sl:closed + sl:open + sl:wet +
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

Core <- function(y, sl, ToD, land, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
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




dat[,mean(propwet_end), by=.(wolfID,ua)]
dat[,mean(propconif_end), by=.(wolfID,ua)]
dat[,mean(propmixed_end), by=.(wolfID,ua)]
dat[,mean(propopen_end), by=.(wolfID,ua)]
    
#run iSSA model by ID
# Core Model

unique(dat$land_end)

dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet', 
                               ifelse(dat$land_end == 'coniferous'|dat$land_end == 'mixed'|dat$land_end == 'deciduous', 'forest','open'))
dat[,'propforest_end_adj'] <- dat$propconif_end+dat$propmixed_end +dat$propdecid_end 
dat[,'propopen_end_adj'] <- dat$propopen_end #+dat$propshrub_end

# dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet', 
#                                ifelse(dat$land_end == 'coniferous', 'coniferous',
#                                       ifelse(dat$land_end == 'deciduous'| dat$land_end == 'mixed'|dat$land_end == 'shrub', 'mixed','open')))
unique(dat$land_end_adj)
dat[,'forest_end_adj'] <- ifelse(dat$land_end_adj == 'forest', 1, 0)
#dat[,'mixed_end_adj'] <- ifelse(dat$land_end_adj == 'mixed', 1, 0)
dat[,'open_end_adj'] <- ifelse(dat$land_end_adj == 'open', 1, 0)
dat[,'wet_end_adj'] <- ifelse(dat$land_end_adj == 'wet', 1, 0)

unique(dat[,.(wolfID)])
#id!='W13', id!='W03' & id!='W14' & id!='W25' 
#wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & wolfID!='RMNP_W12'

# Maybe add RMNPW11 back in? (wet doesn't work right)
#

core2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' # & wolfID!='GHA26_W24', 
                ,Core(case_, log_sl, ToD_start, land_end_adj, stepjum), by = .(wolfID)]


core1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' # & wolfID!='GHA26_W24' 
                 ,Core(case_, log_sl, ToD_start, land_end_adj, stepjum), by = .(wolfID)]



corepropland2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11' 
                    # & wolfID!='GHA26_W35', 
                 ,Core.land(case_, log_sl, ToD_start, propforest_end_adj, propopen_end_adj, propwet_end, stepjum), by = .(wolfID)]


corepropland1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' & wolfID!='GHA26_W24'# & wolfID!='RMNP_W11'
                    # & wolfID!='GHA26_W35'  & wolfID!='GHA26_W25' & wolfID!='GHA26_W26', 
                 ,Core.land(case_, log_sl, ToD_start, propforest_end_adj, propopen_end_adj, propwet_end, stepjum), by = .(wolfID)]

unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32',.(wolfID)])


#### shorter human mods ####

dat.hum <- dat[COD =='human'|COD=='none']

unique(dat.hum[wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
               & wolfID!='GHA26_W34'
               ,.(wolfID)])

corehum2moOUT<- dat.hum[ttd1>31 & ttd1<=38 & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' 
                        & wolfID!='GHA26_W29' 
                        ,Core(case_, log_sl, ToD_start, land_end_adj, stepjum), by = .(wolfID)]

corehum2wkOUT<- dat.hum[ttd1>7 & ttd1<=14 & wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
                         & wolfID!='GHA26_W34' & wolfID!='GHA26_W39'
                        ,Core(case_, log_sl, ToD_start, land_end_adj, stepjum), by = .(wolfID)]

corehum1wkOUT<- dat.hum[ttd1<=7 & wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16'
                         & wolfID!='GHA26_W32' & wolfID!='GHA26_W34'
                        ,Core(case_, log_sl, ToD_start, land_end_adj, stepjum), by = .(wolfID)]

unique(dat.hum[wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16',.(wolfID)])

#### MOVEMENT ####

#### testing probs ####
# dat$wolfID <- as.factor(dat$wolfID)
# dat <- dat.RMNP[, stepjum := forcats::fct_shuffle(as.factor(step_id_)), by = .(id)]
# 
# W06 <- dat[id=='W06', .(case_, log_sl, ToD_start, land_end, ttd1, cos_ta, stepjum, step_id_)]
# W06 <- select(W06, case_, log_sl, ToD_start, land_end, ttd1_adj, cos_ta, stepjum, step_id_)
# write.csv(W06, file = paste0(derived, 'W06_troubleshoot_dat.csv'))
# W06.model <- W06[ttd1>31, clogit(case_ ~ log_sl:ToD_start + land_end +
#                   log_sl:land_end +
#                   log(ttd1+1):log_sl + log(ttd1+1):cos_ta + strata(stepjum))]

Move <- function(y, sl, ToD, land, ttd, ta, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                   # log(ttd+1):sl + log(ttd+1):ta + strata(strata1))
                   ttd:sl + ttd:ta + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}



Move.land <- function(y, sl, ToD, closed, open, wet, ttd, ta, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + closed + open + wet +
                    sl:closed + sl:open + sl:wet +
                    ttd:sl + ttd:ta + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}



move2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                    ,Move(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]


move1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32'# & wolfID!='RMNP_W11'
                   ,Move(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]




movepropland2moOUT <- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11' & wolfID!='GHA26_W35'
                         ,Move.land(case_, log_sl, ToD_start, propforest_end_adj, propopen_end_adj, propwet_end, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]


unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32' & wolfID!='RMNP_W11'
           & wolfID!='GHA26_W35' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25',.(wolfID)])

movepropland1moOUT <- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' & wolfID!='GHA26_W24'# & wolfID!='RMNP_W11' 
                        # & wolfID!='GHA26_W35'  & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
                         ,Move.land(case_, log_sl, ToD_start, propforest_end_adj, propopen_end_adj, propwet_end, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]



#### shorter human mods ####


unique(dat.hum[wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
               & wolfID!='GHA26_W34'
               ,.(wolfID)])

movehum2moOUT<- dat.hum[ttd1>31 & ttd1<=38 & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' 
                        & wolfID!='GHA26_W29' 
                        ,Move(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]

movehum2wkOUT<- dat.hum[ttd1>7 & ttd1<=14 & wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
                        & wolfID!='GHA26_W34' & wolfID!='GHA26_W39'
                        ,Move(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]

movehum1wkOUT<- dat.hum[ttd1<=7 & wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16'
                        & wolfID!='GHA26_W32' & wolfID!='GHA26_W34'
                        ,Move(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), cos_ta, stepjum), by = .(wolfID)]




#### HABITAT ####
Habitat.park <- function(y, sl, ToD, land, ttd, parkdist, rddist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                    ttd:land + ttd:parkdist + ttd:rddist + strata(strata1))
                    #log(ttd+1):land + log(ttd+1):parkdist:parkYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

Habitat<- function(y, sl, ToD, land, ttd, rddist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                    ttd:land + ttd:rddist+ + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}


Habitat.park.land <- function(y, sl, ToD, closed, open, wet, ttd, parkdist, rddist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + closed + open + wet +
                    sl:closed + sl:open + sl:wet +
                    ttd:closed + ttd:open + ttd:wet + ttd:log(parkdist+1) + ttd:log(rddist+1) + strata(strata1))
  #log(ttd+1):land + log(ttd+1):parkdist:parkYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

Habitat.land<- function(y, sl, ToD, closed, open, wet, ttd, rddist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + closed + open + wet +
                    sl:closed + sl:open + sl:wet +
                    ttd:closed + ttd:open + ttd:wet + ttd:rddist + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}




hab2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                 ,Habitat(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+roadDist_end), stepjum), by = .(wolfID)]


unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32' & wolfID!='RMNP_W11',.(wolfID)])

hab1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32'# & wolfID!='RMNP_W11' & wolfID!='GHA26_W24'
                 ,Habitat(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+roadDist_end), stepjum), by = .(wolfID)]







unique(dat[ pop != 'GHA26' & wolfID!='RMNP_W09',.(wolfID)])

habpark2moOUT<- dat[ttd1>31 & pop != 'GHA26' & wolfID!='RMNP_W11', 
                Habitat.park(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+parkDistadj_end), log(1+roadDist_end), stepjum), by = .(wolfID)]

# wolfID!='RMNP_W09' &

habpark1moOUT<- dat[ttd1<=31 & pop != 'GHA26' & wolfID!='RMNP_W11' , 
                Habitat.park(case_, log_sl, ToD_start, land_end_adj,log(1+ttd1), log(1+parkDistadj_end), log(1+roadDist_end), stepjum), by = .(wolfID)]


unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32'# & wolfID!='RMNP_W11' & wolfID!='GHA26_W35'
           ,.(wolfID)])

habpropland2moOUT <- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' & wolfID!='RMNP_W05'# & wolfID!='GHA26_W35'
                         ,Habitat.land(case_, log_sl, ToD_start, scale(propforest_end_adj), scale(propopen_end_adj), scale(propwet_end), log(1+ttd1), scale(roadDist_end), stepjum), by = .(wolfID)]


unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
           #& wolfID!='GHA26_W35' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25'
           ,.(wolfID)])

habpropland1moOUT <- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' & wolfID!='GHA26_W24'# & wolfID!='RMNP_W11' 
                        # & wolfID!='GHA26_W35'  & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
                         ,Habitat.land(case_, log_sl, ToD_start, scale(propforest_end_adj), scale(propopen_end_adj), scale(propwet_end), log(1+ttd1), scale(roadDist_end), stepjum), by = .(wolfID)]



#### shorter human mods ####


unique(dat.hum[wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
               & wolfID!='GHA26_W34'
               ,.(wolfID)])

habhum2moOUT<- dat.hum[ttd1>31 & ttd1<=38 & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' 
                        & wolfID!='GHA26_W29' 
                        ,Habitat(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), scale(roadDist_end), stepjum), by = .(wolfID)]

habhum2wkOUT<- dat.hum[ttd1>7 & ttd1<=14 & wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
                        & wolfID!='GHA26_W34' & wolfID!='GHA26_W39'
                        ,Habitat(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), scale(roadDist_end), stepjum), by = .(wolfID)]

habhum1wkOUT<- dat.hum[ttd1<=7 & wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16'
                        & wolfID!='GHA26_W32' & wolfID!='GHA26_W34'
                       & wolfID!='GHA26_W16'
                        ,Habitat(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), scale(roadDist_end), stepjum), by = .(wolfID)]


unique(dat.hum[wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16'
               & wolfID!='GHA26_W32' & wolfID!='GHA26_W34',.(wolfID)])


#### SOCIAL ####



Social <- function(y, sl, ToD, land, ttd, nndist, packdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                    ttd:nndist+ ttd:packdist + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}


Social.land <- function(y, sl, ToD, closed, open, wet, ttd, nndist, packdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + closed + open + wet + 
                    sl:closed + sl:open + sl:wet +
                    ttd:nndist + ttd:packdist + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]

unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
           & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
           & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
           & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
           & wolfID!='GHA26_W38' & wolfID!='RMNP_W05'
           ,.(wolfID)])

soc2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
                & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                & wolfID!='GHA26_W38' & wolfID!='RMNP_W05' 
                ,Social(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]

socwet2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
                & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                & wolfID!='GHA26_W38' 
                ,Social(case_, log_sl, ToD_start, wet_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]



unique(dat[wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
           & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
           & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
           & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
           & wolfID!='GHA26_W14' & wolfID!='RMNP_W14'
           ,.(wolfID)])

soc1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                   & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                   & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
                   & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                   & wolfID!='GHA26_W14' & wolfID!='RMNP_W14' & wolfID!='GHA26_W31' & wolfID!='GHA26_W38'
                   ,Social(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]


socwet1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
                & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                & wolfID!='GHA26_W14' & wolfID!='RMNP_W14' & wolfID!='GHA26_W31' & wolfID!='GHA26_W38'
                ,Social(case_, log_sl, ToD_start, wet_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]





unique(dat[  wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
             & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
             & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
             & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
             & wolfID!='RMNP_W05' & wolfID!='RMNP_W07' & wolfID!='GHA26_W38'
             ,.(wolfID)])

socpropland2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
                        & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                        & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
                        & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                        & wolfID!='RMNP_W05' & wolfID!='RMNP_W07' & wolfID!='GHA26_W38'
                         ,Social.land(case_, log_sl, ToD_start, propforest_end_adj, propopen_end_adj, propwet_end, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]


socpropland1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
            & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
            & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
            & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
            & wolfID!='RMNP_W10' & wolfID!='GHA26_W14' & wolfID!='RMNP_W14' & wolfID!='GHA26_W24' & wolfID!='GHA26_W31' & wolfID!='GHA26_W38'
            ,Social.land(case_, log_sl, ToD_start, propforest_end_adj, propopen_end_adj, propwet_end, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]

unique(dat[  wolfID!='GHA26_W27' & wolfID!='GHA26_W32' #& wolfID!='RMNP_W11'
             & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
             & wolfID!='GHA26_W15'  & wolfID!='RMNP_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='RMNP_W19'
             & wolfID!='GHA26_W23' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
             & wolfID!='RMNP_W10' & wolfID!='GHA26_W14' & wolfID!='RMNP_W14' & wolfID!='GHA26_W24' & wolfID!='GHA26_W31'
             ,.(wolfID)])




#### shorter human mods ####

 


unique(dat.hum[wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' 
               & wolfID!='GHA26_W29'
               & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
               & wolfID!='GHA26_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
               & wolfID!='GHA26_W14' & wolfID!='RMNP_W14'
               ,.(wolfID)])

sochum2moOUT<- dat.hum[ttd1>31 & ttd1<=38 & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26' 
                       & wolfID!='GHA26_W29' 
                       & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                       & wolfID!='GHA26_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                       & wolfID!='GHA26_W14' & wolfID!='RMNP_W14' & wolfID!='GHA26_W38'
                       ,Social(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]

sochum2wkOUT<- dat.hum[ttd1>7 & ttd1<=14 & wolfID!='GHA26_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W24' & wolfID!='GHA26_W25' & wolfID!='GHA26_W26'
                       & wolfID!='GHA26_W34' & wolfID!='GHA26_W39'
                       & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
                       & wolfID!='GHA26_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
                       & wolfID!='GHA26_W14' & wolfID!='RMNP_W14' & wolfID!='RMNP_W27' & wolfID!='GHA26_W38'
                       ,Social(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]

sochum1wkOUT<- dat.hum[ttd1<=7 & wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16'
               & wolfID!='GHA26_W32' & wolfID!='GHA26_W34'
               & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
               & wolfID!='GHA26_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
               & wolfID!='GHA26_W14' & wolfID!='RMNP_W14'  & wolfID!='GHA26_W22' & wolfID!='GHA26_W24' & wolfID!='GHA26_W26'  
               & wolfID!='RMNP_W27' & wolfID!='GHA26_W38'
               ,Social(case_, log_sl, ToD_start, land_end_adj, log(1+ttd1), log(1+distance2), log(1+packDistadj_end), stepjum), by = .(wolfID)]

unique(dat.hum[wolfID!='RMNP_W03' & wolfID!='GHA26_W11' & wolfID!='RMNP_W16'
               & wolfID!='GHA26_W32' & wolfID!='GHA26_W34'
               & wolfID!='GHA26_W01' & wolfID!='GHA26_W03' & wolfID!='GHA26_W05' & wolfID!='GHA26_W06'
               & wolfID!='GHA26_W15' & wolfID!='GHA26_W16'& wolfID!='RMNP_W16' & wolfID!='GHA26_W23' & wolfID!='GHA26_W36' & wolfID!='GHA26_W39'
               & wolfID!='GHA26_W14' & wolfID!='RMNP_W14'  & wolfID!='GHA26_W22' & wolfID!='GHA26_W24' & wolfID!='GHA26_W26'  
               & wolfID!='RMNP_W27' 
               ,.(wolfID)])

#### AICs ####
m.movepropland2mo.aic <- unique(movepropland2moOUT[,.(AIC_move2mo =AIC), by=.(wolfID)])
m.movepropland1mo.aic <- unique(movepropland1moOUT[,.(AIC_move1mo =AIC), by=.(wolfID)])
move.aic <- merge(m.movepropland2mo.aic, m.movepropland1mo.aic, all = T)
move.aic <- merge(move.aic, dat.meta[,.(COD), by = .(wolfpop)], by.x = 'wolfID', by.y = 'wolfpop')

m.habpropland2mo.aic <- unique(habpropland2moOUT[,.(AIC_hab2mo =AIC), by=.(wolfID)])
m.habpropland1mo.aic <- unique(habpropland1moOUT[,.(AIC_hab1mo =AIC), by=.(wolfID)])
hab.aic <- merge(m.habpropland2mo.aic, m.habpropland1mo.aic, all = T)
hab.aic <- merge(hab.aic, dat.meta[,.(COD), by = .(wolfpop)], by.x = 'wolfID', by.y = 'wolfpop')

m.socpropland2mo.aic <- unique(socpropland2moOUT[,.(AIC_soc2mo =AIC), by=.(wolfID)])
m.socpropland1mo.aic <- unique(socpropland1moOUT[,.(AIC_soc1mo =AIC), by=.(wolfID)])
soc.aic <- merge(m.socpropland2mo.aic, m.socpropland1mo.aic, all = T)
soc.aic <- merge(soc.aic, dat.meta[,.(COD), by = .(wolfpop)], by.x = 'wolfID', by.y = 'wolfpop')

m.all.aic <- merge(move.aic, hab.aic, by=c('wolfID', 'COD'), all = T)
m.all.aic <- merge(m.all.aic, soc.aic, by=c('wolfID', 'COD'), all = T)

#### evidence ratios ####

m.all.aic.2mo<-melt(m.all.aic[,.(wolfID, COD ,move=AIC_move2mo, habitat = AIC_hab2mo, social = AIC_soc2mo)])

m.aic.wts <- m.all.aic.2mo[!is.na(value),.(COD, model=variable, AIC=value, weight=MuMIn::Weights(value)), by=.(wolfID)]

m.aic.wts[, model_rank := frank(AIC, ties.method = 'dense', na.last = 'keep'), by = wolfID]



m.aic.wts.t<-dcast(m.aic.wts,wolfID+COD +model~ model_rank, value.var = c( 'weight'))
m.aic.t<-dcast(m.aic.wts,wolfID+COD +model~ model_rank, value.var = c( 'AIC'))
#m.aic.wts.t[,.(er=)]
m.er <- merge(m.aic.wts[model_rank==1, .(wolfID,COD,model.1=model, AIC.1=AIC)], m.aic.wts[model_rank==2, .(wolfID,COD,model.2=model, AIC.2=AIC)], all =T)
m.er <- merge(m.er, m.aic.wts[model_rank==3, .(wolfID,COD,model.3=model, AIC.3=AIC)], all =T)
m.er[,'er1.2']<- exp((m.er$AIC.2-m.er$AIC.1)/2)
m.er[,'er1.3']<- ifelse(!is.na(m.er$AIC.3),exp((m.er$AIC.3-m.er$AIC.1)/2),NA)
m.er$er1.2<-ifelse(is.infinite(m.er$er1.2), 10^300, m.er$er1.2)

er.mean <-m.er[COD=='cdv'|COD=='human'|COD=='none',.(best=uniqueN(wolfID),mean=mean(er1.2, na.rm=T)), by=.(COD,model.1)]

m.er.cdv <- m.er[COD=='cdv']

#### ICC ####
require(ICC)

m.movepropland2mo.icc <- movepropland2moOUT[term=='coef',-'AIC']
m.movepropland2mo.icc <- m.movepropland2mo.icc[, .( wolfID , lnSL.control = `sl:ttd`, cosTA.control = `ttd:ta`)]

m.movepropland2mo.icc<- merge(m.movepropland2mo.icc, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)



m.movepropland1mo.icc <- movepropland1moOUT[term=='coef',-'AIC']
m.movepropland1mo.icc <- m.movepropland1mo.icc[, .( wolfID , lnSL.case = `sl:ttd`, cosTA.case = `ttd:ta`)]

m.movepropland1mo.icc<- merge(m.movepropland1mo.icc, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)


m.movepropland.icc <- merge(m.movepropland2mo.icc[,.(wolfID, COD, lnSL.control, cosTA.control)], 
                            m.movepropland1mo.icc[,.(wolfID, COD, lnSL.case, cosTA.case)], 
                            by = c('wolfID','COD'), all = T)

#m.movepropland.icc<- melt(m.movepropland.icc)

m.movepropland.icc.cdv <- m.movepropland.icc[COD=='cdv']
m.movepropland.icc.human <- m.movepropland.icc[COD=='human']
m.movepropland.icc.none <- m.movepropland.icc[COD=='none']

cor.test(m.movepropland.icc.cdv$lnSL.control, m.movepropland.icc.cdv$lnSL.case)
cor.test(m.movepropland.icc.human$lnSL.control, m.movepropland.icc.human$lnSL.case)
cor.test(m.movepropland.icc.none$lnSL.control, m.movepropland.icc.none$lnSL.case)

cor.test(m.movepropland.icc.cdv$cosTA.control, m.movepropland.icc.cdv$cosTA.case)
cor.test(m.movepropland.icc.human$cosTA.control, m.movepropland.icc.human$cosTA.case)
cor.test(m.movepropland.icc.none$cosTA.control, m.movepropland.icc.none$cosTA.case)



m.habpropland2mo.icc <- habpropland2moOUT[term=='coef',-'AIC']
m.habpropland2mo.icc <- m.habpropland2mo.icc[, .( wolfID , forest.control = `closed:ttd`,open.control = `open:ttd`, wet.control = `wet:ttd`, road.control = `ttd:rddist`)]

m.habpropland2mo.icc<- merge(m.habpropland2mo.icc, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)



m.habpropland1mo.icc <- habpropland1moOUT[term=='coef',-'AIC']
m.habpropland1mo.icc <- m.habpropland1mo.icc[, .( wolfID , forest.case = `closed:ttd`,open.case = `open:ttd`, wet.case = `wet:ttd`, road.case = `ttd:rddist`)]

m.habpropland1mo.icc<- merge(m.habpropland1mo.icc, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)


m.habpropland.icc <- merge(m.habpropland2mo.icc[,.(wolfID, COD, forest.control, open.control, wet.control, road.control)], 
                            m.habpropland1mo.icc[,.(wolfID, COD, forest.case, open.case, wet.case, road.case)], 
                            by = c('wolfID','COD'), all = T)


m.socpropland2mo.icc <- socpropland2moOUT[term=='coef',-'AIC']
m.socpropland2mo.icc <- m.socpropland2mo.icc[, .( wolfID , nnXttd.control = `ttd:nndist`, packDistXttd.control = `ttd:packdist`)]

m.socpropland2mo.icc<- merge(m.socpropland2mo.icc, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)



m.socpropland1mo.icc <- socpropland1moOUT[term=='coef',-'AIC']
m.socpropland1mo.icc <- m.socpropland1mo.icc[, .( wolfID , nnXttd.case = `ttd:nndist`, packDistXttd.case = `ttd:packdist`)]

m.socpropland1mo.icc<- merge(m.socpropland1mo.icc, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)


m.socpropland.icc <- merge(m.socpropland2mo.icc[,.(wolfID, COD, nnXttd.control, packDistXttd.control)], 
                           m.socpropland1mo.icc[,.(wolfID, COD, nnXttd.case, packDistXttd.case)], 
                           by = c('wolfID','COD'), all = T)



#### delta betas ####

m.movepropland.icc[,'d_lnSL'] <- ifelse(!is.na(m.movepropland.icc$lnSL.control)&!is.na(m.movepropland.icc$lnSL.case), 
                                        m.movepropland.icc$lnSL.control- m.movepropland.icc$lnSL.case, NA)
m.movepropland.icc[,'d_cosTA'] <- ifelse(!is.na(m.movepropland.icc$cosTA.control)&!is.na(m.movepropland.icc$cosTA.case), 
                                         m.movepropland.icc$cosTA.control- m.movepropland.icc$cosTA.case, NA)
m.movepropland.icc[,'model'] <- 'move'

m.habpropland.icc[,'d_forest']<- ifelse(!is.na(m.habpropland.icc$forest.control)&!is.na(m.habpropland.icc$forest.case),
                                        m.habpropland.icc$forest.control -m.habpropland.icc$forest.case, NA)
m.habpropland.icc[,'d_open']<-ifelse(!is.na(m.habpropland.icc$open.control)&!is.na(m.habpropland.icc$open.case),
                                     m.habpropland.icc$open.control -m.habpropland.icc$open.case, NA)
m.habpropland.icc[,'d_wet']<- ifelse(!is.na(m.habpropland.icc$wet.control)&!is.na(m.habpropland.icc$wet.case),
                                     m.habpropland.icc$wet.control -m.habpropland.icc$wet.case, NA)
m.habpropland.icc[,'d_road']<- ifelse(!is.na(m.habpropland.icc$road.control)&!is.na(m.habpropland.icc$road.case),
                                      m.habpropland.icc$road.control -m.habpropland.icc$road.case, NA)
m.habpropland.icc[,'model'] <- 'habitat'


m.socpropland.icc[,'d_nn'] <- ifelse(!is.na(m.socpropland.icc$nnXttd.control)&!is.na(m.socpropland.icc$nnXttd.case),
                                     m.socpropland.icc$nnXttd.control -  m.socpropland.icc$nnXttd.case, NA)
m.socpropland.icc[,'d_pack'] <-ifelse(!is.na(m.socpropland.icc$packDistXttd.control)&!is.na(m.socpropland.icc$packDistXttd.case),
                                      m.socpropland.icc$packDistXttd.control -  m.socpropland.icc$packDistXttd.case, NA)

m.socpropland.icc[,'model'] <- 'social'

dbetas <- merge(m.movepropland.icc[,.(wolfID,COD, d_lnSL, d_cosTA)], 
                m.habpropland.icc[,.(wolfID,COD, d_forest,d_open,d_wet,d_road)], by=c('wolfID','COD'), all=T)

dbetas <- merge(dbetas, m.socpropland.icc[,.(wolfID,COD, d_nn, d_pack)], by=c('wolfID','COD'), all=T)
dbetas.t <- melt(dbetas)
dbetas.t[,'model'] <- ifelse(dbetas.t$variable == 'd_lnSL'|dbetas.t$variable == 'd_cosTA','move',
                             ifelse(dbetas.t$variable == 'd_nn'|dbetas.t$variable == 'd_pack', 'social', 'habitat'))

# dbetas.cdv <- dbetas.t[COD=='cdv']
# dbetas.human <- dbetas.t[COD=='human']
# dbetas.none <- dbetas.t[COD=='none']

dbetas.t <- dbetas.t[COD=='cdv'|COD=='human'|COD=='none'] 
dbetas.t$COD <- factor(dbetas.t$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))


cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")

ggplot(dbetas.t[model=='move'], aes(variable, (value), fill = COD)) +
  geom_boxplot(aes(fill = COD),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = COD),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('delta movement') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-.3,.3)


ggplot(dbetas.t[model=='habitat'], aes(variable, (-1*value), fill = COD)) +
  geom_boxplot(aes(fill = COD),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = COD),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('delta selection') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-5,5)

ggplot(dbetas.t[model=='habitat'& variable=='d_road'], aes(variable, (-1*value), fill = COD)) +
  geom_boxplot(aes(fill = COD),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = COD),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('delta selection') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)# + ylim(-5,5)


ggplot(dbetas.t[model=='social'], aes(variable, (-1*value), fill = COD)) +
  geom_boxplot(aes(fill = COD),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = COD),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('delta selection') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-.8,.8)




#### % in/out park ####
pack2mo <-dat[ttd1>31 & ua=='used' & pop=='RMNP', .(COD, tot= uniqueN(t2_)), by=.(wolfID, packYN_end)]
pack2mo.t<- dcast(pack2mo, wolfID+COD~packYN_end)
pack2mo.t$`out-pack`<- ifelse(is.na(pack2mo.t$`out-pack`), 0, pack2mo.t$`out-pack`)
pack2mo.t[,'total_steps'] <- pack2mo.t$`out-pack`+pack2mo.t$pack
pack2mo.t[,'propout'] <- pack2mo.t$`out-pack`/pack2mo.t$total_steps


pack1mo <-dat[ttd1<=31 & ua=='used' & pop=='RMNP', .(COD, tot= uniqueN(t2_)), by=.(wolfID, packYN_end)]
pack1mo.t<- dcast(pack1mo, wolfID+COD~packYN_end)
pack1mo.t$`out-pack`<- ifelse(is.na(pack1mo.t$`out-pack`), 0, pack1mo.t$`out-pack`)
pack1mo.t[,'total_steps'] <- pack1mo.t$`out-pack`+pack1mo.t$pack
pack1mo.t[,'propout'] <- pack1mo.t$`out-pack`/pack1mo.t$total_steps

pack.propout <- merge(pack2mo.t[,.(wolfID, COD, propout2mo = propout)], pack1mo.t[,.(wolfID, COD, propout1mo = propout)])
pack.propout <- melt(pack.propout)
pack.propout <- pack.propout[COD=='cdv'|COD=='human'|COD=='none']

pack.propmean <- pack.propout[,mean(value), .(COD, variable)]


ggplot(pack.propout, aes(variable, (value), fill = COD)) +
  geom_boxplot(aes(fill = COD),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = COD),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('delta selection') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(0,.25)

#### GRAPHS ####
#### movement graphs ####
m.move2mo.coef <- move2moOUT[term=='coef',-'AIC']
#m.move2mo.coef <- m.move2mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.move2mo.coef <- m.move2mo.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.move2mo.coef<- merge(m.move2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.move2mo.coef <- melt(m.move2mo.coef)
m.move2mo.coef[,'ttd'] <- '2mo'

m.move1mo.coef <- move1moOUT[term=='coef',-'AIC']
#m.move1mo.coef <- m.move1mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.move1mo.coef <- m.move1mo.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.move1mo.coef<- merge(m.move1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.move1mo.coef<- melt(m.move1mo.coef)
m.move1mo.coef[,'ttd'] <- '1mo'

m.move.coef <- rbind(m.move1mo.coef, m.move2mo.coef)
m.move.coef[,'test'] <- ifelse(m.move.coef$ttd =='1mo', 'case', 'control')
m.move.coef$test <- factor(m.move.coef$test, levels = c('control','case'))

m.move.coef.cdv <- m.move.coef[COD=='cdv']
m.move.coef.human <- m.move.coef[COD=='human']
m.move.coef.none <- m.move.coef[COD=='none']

#color = c("#0072B2", "#D55E00", "#009E73")
color = c("darkviolet", "aquamarine4")

ggplot(m.move.coef.cdv, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
        axis.title = element_text(size = 16, color = 'black'),
        axis.text = element_text(size = 14, color = 'black'),
        plot.title=element_text(size = 16, hjust=0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(colour="black", size = 1, fill = "white"),
        strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-.25,.25)


ggplot(m.move.coef.human, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-.5,.5)




ggplot(m.move.coef.none, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-.25,.25)



#### movement propland graphs ####

m.movepropland2mo.coef <- movepropland2moOUT[term=='coef',-'AIC']
#m.movepropland2mo.coef <- m.movepropland2mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movepropland2mo.coef <- m.movepropland2mo.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.movepropland2mo.coef<- merge(m.movepropland2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movepropland2mo.coef <- melt(m.movepropland2mo.coef)
m.movepropland2mo.coef[,'ttd'] <- '2mo'

m.movepropland1mo.coef <- movepropland1moOUT[term=='coef',-'AIC']
#m.movepropland1mo.coef <- m.movepropland1mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movepropland1mo.coef <- m.movepropland1mo.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.movepropland1mo.coef<- merge(m.movepropland1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movepropland1mo.coef<- melt(m.movepropland1mo.coef)
m.movepropland1mo.coef[,'ttd'] <- '1mo'

m.movepropland.coef <- rbind(m.movepropland1mo.coef, m.movepropland2mo.coef)
m.movepropland.coef[,'test'] <- ifelse(m.movepropland.coef$ttd =='1mo', 'case', 'control')
m.movepropland.coef$test <- factor(m.movepropland.coef$test, levels = c('control','case'))

m.movepropland.coef.cdv <- m.movepropland.coef[COD=='cdv']
m.movepropland.coef.human <- m.movepropland.coef[COD=='human']
m.movepropland.coef.none <- m.movepropland.coef[COD=='none']

#color = c("#0072B2", "#D55E00", "#009E73")
color = c("darkviolet", "aquamarine4")

ggplot(m.movepropland.coef.cdv, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-.25,.25)


ggplot(m.movepropland.coef.human, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-.2,.4)


test<-m.movepropland.coef.none[,head(.SD,10),by=.(variable,test)]

ggplot(m.movepropland.coef.none[,head(.SD,10),by=.(variable,test)], aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-.25,.25)


#### movement human mod graphs ####

m.movehum2mo.coef <- movehum2moOUT[term=='coef',-'AIC']
#m.movehum2mo.coef <- m.movehum2mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movehum2mo.coef <- m.movehum2mo.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.movehum2mo.coef<- merge(m.movehum2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movehum2mo.coef <- melt(m.movehum2mo.coef)
m.movehum2mo.coef[,'ttd'] <- '2mo'

m.movehum2wk.coef <- movehum2wkOUT[term=='coef',-'AIC']
#m.movehum2wk.coef <- m.movehum2wk.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movehum2wk.coef <- m.movehum2wk.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.movehum2wk.coef<- merge(m.movehum2wk.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movehum2wk.coef<- melt(m.movehum2wk.coef)
m.movehum2wk.coef[,'ttd'] <- '2wk'

m.movehum1wk.coef <- movehum1wkOUT[term=='coef',-'AIC']
#m.movehum1wk.coef <- m.movehum1wk.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movehum1wk.coef <- m.movehum1wk.coef[, .( wolfID , lnSL = `sl:ttd`, cosTA = `ttd:ta`)]

m.movehum1wk.coef<- merge(m.movehum1wk.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movehum1wk.coef<- melt(m.movehum1wk.coef)
m.movehum1wk.coef[,'ttd'] <- '1wk'

m.movehum.coef <- rbind(m.movehum1wk.coef, m.movehum2wk.coef, m.movehum2mo.coef)
#m.movehum.coef[,'test'] <- ifelse(m.movehum.coef$ttd =='1mo', 'case', 'control')
m.movehum.coef$ttd <- factor(m.movehum.coef$ttd, levels = c('2mo','2wk', '1wk'))

m.movehum.coef.human <- m.movehum.coef[COD=='human']
m.movehum.coef.none <- m.movehum.coef[COD=='none']

#color = c("#0072B2", "#D55E00", "#009E73")
color = c("darkviolet", "violet", "aquamarine4")

ggplot(m.movehum.coef.human, aes(variable, (-1*value), fill = ttd)) +
  geom_boxplot(aes(fill = ttd),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = ttd),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) #+ ylim(-.25,.25)


ggplot(m.movehum.coef.none, aes(variable, (-1*value), fill = ttd)) +
  geom_boxplot(aes(fill = ttd),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = ttd),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('movement') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-1,1)




#### habitat graphs ####

m.hab2mo.coef <- hab2moOUT[term=='coef',-'AIC']
m.hab2mo.coef <- m.hab2mo.coef[, .( wolfID, forest = `landforest:ttd`, open = `landopen:ttd`, wet=`landwet:ttd`, roads = `ttd:rddist`)]
m.hab2mo.coef<- merge(m.hab2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.hab2mo.coef <- melt(m.hab2mo.coef)
m.hab2mo.coef[,'ttd'] <- '2mo'

m.hab1mo.coef <- hab1moOUT[term=='coef',-'AIC']
m.hab1mo.coef <- m.hab1mo.coef[, .(  wolfID, forest = `landforest:ttd`, open = `landopen:ttd`, wet=`landwet:ttd`, roads = `ttd:rddist`)]
m.hab1mo.coef<- merge(m.hab1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.hab1mo.coef<- melt(m.hab1mo.coef)
m.hab1mo.coef[,'ttd'] <- '1mo'

m.hab.coef <- rbind(m.hab1mo.coef, m.hab2mo.coef)
m.hab.coef[,'test'] <- ifelse(m.hab.coef$ttd =='1mo', 'case', 'control')
m.hab.coef$test <- factor(m.hab.coef$test, levels = c('control','case'))

m.hab.coef.cdv <- m.hab.coef[COD == 'cdv']
m.hab.coef.human <- m.hab.coef[COD == 'human']
m.hab.coef.none <- m.hab.coef[COD == 'none']
#m.hab.coef.land <- m.hab.coef[variable=='closed' | variable=='open' | variable=='wet' & COD == 'cdv']
#m.hab.coef.rd <- m.hab.coef[variable=='roads' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.hab.coef.cdv, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-8,4)


# ggplot(m.hab.coef.rd, aes(variable, value, fill = test)) +
#   geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
#                outlier.color = NA, lwd = 0.6,
#                alpha = 0.25) +
#   geom_jitter(aes(color = test),
#               position = position_jitterdodge(.35),
#               size = 2, alpha = 0.4) +
#   #ggtitle('Interaction with community identity') +
#   geom_hline(aes(yintercept = 0), lty = 2) +
#   theme(#legend.position = 'none',
#     axis.title = element_text(size = 16, color = 'black'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title=element_text(size = 16, hjust=0),
#     axis.line = element_line(colour = "black"),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     strip.background = element_rect(colour="black", size = 1, fill = "white"),
#     strip.text = element_text(size = 14)) +
#   xlab('') +
#   ylab('Selection') +
#   scale_fill_manual(values = color) +
#   scale_color_manual(values = color) +
#   ylim(-.1,.15)
# 

ggplot(m.hab.coef.human, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-5,5)


ggplot(m.hab.coef.none, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-5,5)


#### habitat propland graphs ####

m.habpropland2mo.coef <- habpropland2moOUT[term=='coef',-'AIC']
m.habpropland2mo.coef <- m.habpropland2mo.coef[, .( wolfID, forest = `closed:ttd`, open = `open:ttd`, wet=`wet:ttd`, roads = `ttd:rddist`)]
m.habpropland2mo.coef<- merge(m.habpropland2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habpropland2mo.coef <- melt(m.habpropland2mo.coef)
m.habpropland2mo.coef[,'ttd'] <- '2mo'

m.habpropland1mo.coef <- habpropland1moOUT[term=='coef',-'AIC']
m.habpropland1mo.coef <- m.habpropland1mo.coef[, .(  wolfID, forest = `closed:ttd`, open = `open:ttd`, wet=`wet:ttd`, roads = `ttd:rddist`)]
m.habpropland1mo.coef<- merge(m.habpropland1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habpropland1mo.coef<- melt(m.habpropland1mo.coef)
m.habpropland1mo.coef[,'ttd'] <- '1mo'

m.habpropland.coef <- rbind(m.habpropland1mo.coef, m.habpropland2mo.coef)
m.habpropland.coef[,'test'] <- ifelse(m.habpropland.coef$ttd =='1mo', 'case', 'control')
m.habpropland.coef$test <- factor(m.habpropland.coef$test, levels = c('control','case'))

m.habpropland.coef.cdv <- m.habpropland.coef[COD == 'cdv']
m.habpropland.coef.human <- m.habpropland.coef[COD == 'human']
m.habpropland.coef.none <- m.habpropland.coef[COD == 'none']
#m.habpropland.coef.land <- m.habpropland.coef[variable=='closed' | variable=='open' | variable=='wet' & COD == 'cdv']
#m.habpropland.coef.rd <- m.habpropland.coef[variable=='roads' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.habpropland.coef.cdv, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-50,50)

ggplot(m.habpropland.coef.cdv[variable=='roads'], aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-5,4)



ggplot(m.habpropland.coef.human, aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-20,30)


ggplot(m.habpropland.coef.human[variable=='roads'], aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-3,5)

ggplot(m.habpropland.coef.none[,head(.SD,10),by=.(variable,test)], aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-5,5)


ggplot(m.habpropland.coef.none[variable=='roads',head(.SD,10),by=.(variable,test)], aes(variable, (-1*value), fill = test)) +
  geom_boxplot(aes(fill = test),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-5,5)
### MAKE PARK GRAPHS

#### habitat human mod graphs ####

m.habhum2mo.coef <- habhum2moOUT[term=='coef',-'AIC']
#m.habhum2mo.coef <- m.habhum2mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.habhum2mo.coef <- m.habhum2mo.coef[, .( wolfID , forest = `landforest:ttd`, open = `landopen:ttd`, wet=`landwet:ttd`, roads = `ttd:rddist`)]

m.habhum2mo.coef<- merge(m.habhum2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habhum2mo.coef <- melt(m.habhum2mo.coef)
m.habhum2mo.coef[,'ttd'] <- '2mo'

m.habhum2wk.coef <- habhum2wkOUT[term=='coef',-'AIC']
#m.habhum2wk.coef <- m.habhum2wk.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.habhum2wk.coef <- m.habhum2wk.coef[, .( wolfID ,forest = `landforest:ttd`, open = `landopen:ttd`, wet=`landwet:ttd`, roads = `ttd:rddist`)]

m.habhum2wk.coef<- merge(m.habhum2wk.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habhum2wk.coef<- melt(m.habhum2wk.coef)
m.habhum2wk.coef[,'ttd'] <- '2wk'

m.habhum1wk.coef <- habhum1wkOUT[term=='coef',-'AIC']
#m.habhum1wk.coef <- m.habhum1wk.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.habhum1wk.coef <- m.habhum1wk.coef[, .( wolfID , forest = `landforest:ttd`, open = `landopen:ttd`, wet=`landwet:ttd`, roads = `ttd:rddist`)]

m.habhum1wk.coef<- merge(m.habhum1wk.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habhum1wk.coef<- melt(m.habhum1wk.coef)
m.habhum1wk.coef[,'ttd'] <- '1wk'

m.habhum.coef <- rbind(m.habhum1wk.coef, m.habhum2wk.coef, m.habhum2mo.coef)
#m.habhum.coef[,'test'] <- ifelse(m.habhum.coef$ttd =='1mo', 'case', 'control')
m.habhum.coef$ttd <- factor(m.habhum.coef$ttd, levels = c('2mo','2wk', '1wk'))

m.habhum.coef.human <- m.habhum.coef[COD=='human']
m.habhum.coef.none <- m.habhum.coef[COD=='none']

#color = c("#0072B2", "#D55E00", "#009E73")
#color = c("darkviolet", "violet", "aquamarine4")

ggplot(m.habhum.coef.human, aes(variable, (-1*value), fill = ttd)) +
  geom_boxplot(aes(fill = ttd),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = ttd),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-5,5)


ggplot(m.habhum.coef.none, aes(variable, (-1*value), fill = ttd)) +
  geom_boxplot(aes(fill = ttd),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = ttd),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-1,1)



#### social graphs ####

m.soc2mo.coef <- soc2moOUT[term=='coef',-'AIC']
m.soc2mo.coef <- m.soc2mo.coef[, .( wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]
m.soc2mo.coef<- merge(m.soc2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.soc2mo.coef <- melt(m.soc2mo.coef)
m.soc2mo.coef[,'ttd'] <- '2mo'

m.soc1mo.coef <- soc1moOUT[term=='coef',-'AIC']
m.soc1mo.coef <- m.soc1mo.coef[, .(wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]
m.soc1mo.coef<- merge(m.soc1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.soc1mo.coef<- melt(m.soc1mo.coef)
m.soc1mo.coef[,'ttd'] <- '1mo'

m.soc.coef <- rbind(m.soc1mo.coef, m.soc2mo.coef)
m.soc.coef[,'test'] <- ifelse(m.soc.coef$ttd =='1mo', 'case', 'control')
m.soc.coef$test <- factor(m.soc.coef$test, levels = c('control','case'))

m.soc.coef.cdv <- m.soc.coef[COD == 'cdv']
m.soc.coef.human <- m.soc.coef[COD == 'human']
m.soc.coef.none <- m.soc.coef[COD == 'none']

# m.soc.coef.nn <- m.soc.coef[variable =='nnXttd' & COD == 'cdv']
# m.soc.coef.pack <- m.soc.coef[variable =='packDistXttd' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.soc.coef.cdv, aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-.25,.05)


ggplot(m.soc.coef.human, aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-.25,.05)



ggplot(m.soc.coef.none[,head(.SD,10),by=.(variable,test)], aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-1,2.5)




#### social wet graphs ####

m.socwet2mo.coef <- socwet2moOUT[term=='coef',-'AIC']
m.socwet2mo.coef <- m.socwet2mo.coef[, .( wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]
m.socwet2mo.coef<- merge(m.socwet2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.socwet2mo.coef <- melt(m.socwet2mo.coef)
m.socwet2mo.coef[,'ttd'] <- '2mo'

m.socwet1mo.coef <- socwet1moOUT[term=='coef',-'AIC']
m.socwet1mo.coef <- m.socwet1mo.coef[, .(wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]
m.socwet1mo.coef<- merge(m.socwet1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.socwet1mo.coef<- melt(m.socwet1mo.coef)
m.socwet1mo.coef[,'ttd'] <- '1mo'

m.socwet.coef <- rbind(m.socwet1mo.coef, m.socwet2mo.coef)
m.socwet.coef[,'test'] <- ifelse(m.socwet.coef$ttd =='1mo', 'case', 'control')
m.socwet.coef$test <- factor(m.socwet.coef$test, levels = c('control','case'))

m.socwet.coef.cdv <- m.socwet.coef[COD == 'cdv']
m.socwet.coef.human <- m.socwet.coef[COD == 'human']
m.socwet.coef.none <- m.socwet.coef[COD == 'none']

# m.socwet.coef.nn <- m.socwet.coef[variable =='nnXttd' & COD == 'cdv']
# m.socwet.coef.pack <- m.socwet.coef[variable =='packDistXttd' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.socwet.coef.cdv, aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-.25,.05)


ggplot(m.socwet.coef.human, aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-.25,.05)



ggplot(m.socwet.coef.none[,head(.SD,10),by=.(variable,test)], aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-1,2.5)




#### social propland graphs ####


m.socpropland2mo.coef <- socpropland2moOUT[term=='coef',-'AIC']
m.socpropland2mo.coef <- m.socpropland2mo.coef[, .( wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]
m.socpropland2mo.coef<- merge(m.socpropland2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.socpropland2mo.coef <- melt(m.socpropland2mo.coef)
m.socpropland2mo.coef[,'ttd'] <- '2mo'

m.socpropland1mo.coef <- socpropland1moOUT[term=='coef',-'AIC']
m.socpropland1mo.coef <- m.socpropland1mo.coef[, .(wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]
m.socpropland1mo.coef<- merge(m.socpropland1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.socpropland1mo.coef<- melt(m.socpropland1mo.coef)
m.socpropland1mo.coef[,'ttd'] <- '1mo'

m.socpropland.coef <- rbind(m.socpropland1mo.coef, m.socpropland2mo.coef)
m.socpropland.coef[,'test'] <- ifelse(m.socpropland.coef$ttd =='1mo', 'case', 'control')
m.socpropland.coef$test <- factor(m.socpropland.coef$test, levels = c('control','case'))

m.socpropland.coef.cdv <- m.socpropland.coef[COD == 'cdv']
m.socpropland.coef.human <- m.socpropland.coef[COD == 'human']
m.socpropland.coef.none <- m.socpropland.coef[COD == 'none']

# m.socpropland.coef.nn <- m.socpropland.coef[variable =='nnXttd' & COD == 'cdv']
# m.socpropland.coef.pack <- m.socpropland.coef[variable =='packDistXttd' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.socpropland.coef.cdv, aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-.25,.05)


ggplot(m.socpropland.coef.human, aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)# + ylim(-.25,.05)



ggplot(m.socpropland.coef.none[,head(.SD,10),by=.(variable,test)], aes(variable, value, fill = test)) +
  geom_boxplot(aes(fill = test), #notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = test),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('Selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-1,1)




#### social human mod graphs ####

m.sochum2mo.coef <- sochum2moOUT[term=='coef',-'AIC']
#m.sochum2mo.coef <- m.sochum2mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.sochum2mo.coef <- m.sochum2mo.coef[, .( wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]

m.sochum2mo.coef<- merge(m.sochum2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.sochum2mo.coef <- melt(m.sochum2mo.coef)
m.sochum2mo.coef[,'ttd'] <- '2mo'

m.sochum2wk.coef <- sochum2wkOUT[term=='coef',-'AIC']
#m.sochum2wk.coef <- m.sochum2wk.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.sochum2wk.coef <- m.sochum2wk.coef[, .( wolfID ,nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]

m.sochum2wk.coef<- merge(m.sochum2wk.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.sochum2wk.coef<- melt(m.sochum2wk.coef)
m.sochum2wk.coef[,'ttd'] <- '2wk'

m.sochum1wk.coef <- sochum1wkOUT[term=='coef',-'AIC']
#m.sochum1wk.coef <- m.sochum1wk.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.sochum1wk.coef <- m.sochum1wk.coef[, .( wolfID , nnXttd = `ttd:nndist`, packDistXttd = `ttd:packdist`)]

m.sochum1wk.coef<- merge(m.sochum1wk.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.sochum1wk.coef<- melt(m.sochum1wk.coef)
m.sochum1wk.coef[,'ttd'] <- '1wk'

m.sochum.coef <- rbind(m.sochum1wk.coef, m.sochum2wk.coef, m.sochum2mo.coef)
#m.sochum.coef[,'test'] <- ifelse(m.sochum.coef$ttd =='1mo', 'case', 'control')
m.sochum.coef$ttd <- factor(m.sochum.coef$ttd, levels = c('2mo','2wk', '1wk'))

m.sochum.coef.human <- m.sochum.coef[COD=='human']
m.sochum.coef.none <- m.sochum.coef[COD=='none']

#color = c("#0072B2", "#D55E00", "#009E73")
#color = c("darkviolet", "violet", "aquamarine4")

ggplot(m.sochum.coef.human, aes(variable, (value), fill = ttd)) +
  geom_boxplot(aes(fill = ttd),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = ttd),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-5,5)


ggplot(m.sochum.coef.none, aes(variable, (value), fill = ttd)) +
  geom_boxplot(aes(fill = ttd),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = ttd),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('') +
  ylab('selection') +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) + ylim(-1,1)



#### moveVis ####
require(moveVis)
require(move)

move_df <- unique(dat[wolfID == 'RMNP_W06'| wolfID == 'RMNP_W02'| wolfID == 'RMNP_W12',
                .(wolfID, x1_, y1_, t1_)])

utm14N <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

move_df <- move_df[, c('Long', 'Lat') := as.data.table(proj4::project(cbind(x1_, y1_), utm14N, inverse= T))]


move_dat <- df2move(move_df, proj = utm14N, x = 'Long', y = 'Lat', time = 't1_', track_id = 'wolfID')
m <- align_move(move_dat, unit = 'hours')

frames <- frames_spatial(m, path_colours = c("purple", "green", "blue"),
                         map_service = "osm",
                         map_type = "watercolor")
                         #map_service = "mapbox",
                         #map_token = 'pk.eyJ1Ijoiand0dXJuZXIiLCJhIjoiY2p6NXJtbDNzMGFpMTNjb3V0eG00MzNuaCJ9.gYv2XPxEGbcFGO3OMu948Q', 
                         #map_type = "satellite")# %>% 
  # add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  # add_northarrow() %>% 
  # add_scalebar() %>% 
  # add_timestamps(m, type = "label") %>% 
  # add_progress()





