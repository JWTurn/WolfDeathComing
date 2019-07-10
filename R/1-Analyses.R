### Analyses ====
# Julie Turner
# Started: July 3 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### SSF data with all covariates ----
dat <- readRDS('data/derived-data/ssfAllCov.Rds')
dat[,'wtd1'] <- round(dat$ttd1/7)
dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')
dat[,'wet'] <- ifelse(dat$land_end == 'wet', 'wet', 'not')
dat.wet<-dat[land_end=='wet', .(id, .N, wtd1, ttd1, ua)]
dat.wet[, .(id, .N), by= c('ua', 'id')]

ggplot(dat, aes(x=factor(wtd1),y=wet, group=(wet)))+
  stat_summary(aes(color=(wet)),fun.y=length, geom="line")+
  scale_color_discrete("wet",labels=c("not","wet"))+
  labs(x="",y="Frequency")

dat.meta <- fread(paste0(raw, 'wolf_metadata.csv'))


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

Core <- function(y, sl, ToD, land, strata1) {
    # Make the model
    model <- clogit(y ~ sl:ToD + land + sl:land + strata(strata1))
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
coreOUT<- dat[, Core(case_, log_sl, ToD_start, land_end, step_id_), by = id]
m.core <- merge(coreOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)


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

moveOUT <- dat[, Move(case_, log_sl, ToD_start, land_end, wtd1, cos_ta, step_id_), by = id]
m.move <- merge(moveOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.move[term=='p', .(id, `sl:ttd`, `ttd:ta`, COD)]



#### HABITAT ####

Habitat <- function(y, sl, ToD, land, ttd, parkdist, parkYN, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    ttd:land + ttd:parkdist:parkYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

habOUT <- dat[, Habitat(case_, log_sl, ToD_start, land_end, wtd1, lnparkdist_end, parkYN_end, step_id_), by = id]
m.habitat <- merge(habOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitat.p <- m.habitat[term=='p', .(id, `landclosed:ttd`, `landopen:ttd`, `landwet:ttd`, 
                       `ttd:parkdist:parkYNpark`,`ttd:parkdist:parkYNout-park`, COD)]

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

hablandOUT <- dat[, Habitat.land(case_, log_sl, ToD_start, land_end, ttd1, step_id_), by = id]
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

habparkOUT <- dat[, Habitat.park(case_, log_sl, ToD_start, land_end, ttd1, lnparkdist_end, parkYN_end, step_id_), by = id]
m.habitat.park <- merge(habparkOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitat.park.p <- m.habitat.park[term=='p', .(id, `log(ttd + 1):parkdist:parkYNpark`,`log(ttd + 1):parkdist:parkYNout-park`, COD)]


#### SOCIAL ####

m.social <- data.all %>% amt::fit_issf(case_ ~ log_sl:ToD_start + land_end + log_sl:land_end + 
                                        log(ttd1):log(distance2) + 
                                        strata(step_id_))

Social <- function(y, sl, ToD, land, ttd, nndist, packYN, packdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land + sl:land + 
                    ttd:log(nndist+1) + ttd:log(packdist+1):packYN + strata(strata1))
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

socOUT <- dat[id !='W13' & id != 'W11' & id != 'W19' & id != 'W15', Social(case_, log_sl, ToD_start, land_end, wtd1, distance2, packYN_end, packDist_end, step_id_), by = id]
m.social <- merge(socOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.social.p <- m.social[term=='coef' |term=='p', .(id, term,`ttd:log(nndist + 1)`, `ttd:log(packdist + 1):packYNpack`, `ttd:log(packdist + 1):packYNout-pack`, COD)]



bbmle::AICtab(m.core$model, m.move$model, m.habitat.parkdist$model, m.social$model)


