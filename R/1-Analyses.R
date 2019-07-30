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
set.seed(37)
dat <- readRDS('data/derived-data/ssfAllCov.Rds')

dat <- dat[, stepjum := forcats::fct_shuffle(as.factor(step_id_)), by = .(id)]

dat[,.(nstep = uniqueN(as.character(step_id_)), nstepjum = uniqueN(stepjum)), by = .(id)]

dat.steps <-dat[id=='W06',.(step_id_, stepjum)]

dat<-dat[ttd1>=0 & ttd2>=0]

dat[,'wtd1'] <- as.integer(dat$ttd1/7)

dat.meta <- fread(paste0(raw, 'wolf_metadata.csv'))



dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')

dat<- merge(dat, dat.meta, by.x = 'id', by.y = 'WolfID')


dat[,'ttd1_adj'] <- ifelse(dat$COD =='none', dat$ttd1*730, dat$ttd1)

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


dat.wet[ua=='used',cor(propwet_end, ttd1), by=.(id)]




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
model <- clogit(case_ ~ log_sl:ToD_start + land_end + 
                  #sl:closed + sl:open + sl:wet + 
                  log_sl:land_end +
                  strata(stepjum), data = dat[id=='W06'])
sum.mod <- summary(model)
sum.mod$conf.int

### dummy vars didn't, Can run prop with this same one
Core.land <- function(y, sl, ToD, conif, mixed, open, wet, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + conif + mixed + open + wet +
                    sl:conif + sl:mixed + sl:open + sl:wet +
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




dat[,mean(propwet_end), by=.(id,ua)]
dat[,mean(propconif_end), by=.(id,ua)]
dat[,mean(propmixed_end), by=.(id,ua)]
dat[,mean(propopen_end), by=.(id,ua)]
    
#run iSSA model by ID
# Core Model

#id!='W13'
coreOUT<- dat[id!='W03' & id!='W14' & id!='W25', Core(case_, log_sl, ToD_start, land_end, stepjum), by = id]
m.core <- merge(coreOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)


m.core.co <- coreOUT[term=='coef',-'AIC']
m.core.co <-merge(m.core.co, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

ggplot(melt(m.core.co)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter() +
  ylim(-1,1)

ggplot(melt(coreOUT)) +
  geom_density(aes(value), fill = 'dodgerblue', alpha =0.5) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  facet_wrap(~variable, scale = "free")


corewetOUT<- dat[id!='W03' & id!='W14' & id!='W25', Core(case_, log_sl, ToD_start, wet_end, stepjum), by = id]
m.corewet <- merge(corewetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

corepropwetOUT<- dat[id!='W03' & id!='W14' & id!='W25', Core(case_, log_sl, ToD_start, propwet_end, stepjum), by = id]
m.corepropwet <- merge(corepropwetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)



# wet still not work, many fitter error
corelandOUT<- dat[id!='W03' & id!='W14' & id!='W25', Core.land(case_, log_sl, ToD_start, conif_end, mixed_end, open_end, wet_end, stepjum), by = id]
m.coreland <- merge(corelandOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.coreland.p<-m.coreland[term=='coef' | term=='p']
m.coreland.coef<- m.coreland[term=='coef' ]



corepropOUT<- dat[id!='W03' & id!='W14' & id!='W25', Core.land(case_, log_sl, ToD_start, propconif_end, propmixed_end, propopen_end, propwet_end, stepjum), by = id]
m.coreprop <- merge(corepropOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.coreprop.p<-m.coreprop[term=='coef' | term=='p']




#### MOVEMENT ####

Move <- function(y, sl, ToD, land, ttd, ta, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                   # log(ttd+1):sl + log(ttd+1):ta + strata(strata1))
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

Move.td <- function(y, sl, ToD, land, ttd, ta, strata1) {
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


Move.prop <- function(y, sl, ToD, wet, open, ttd, ta, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + wet + open +
                    sl:wet + sl:open +
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





unique(dat[id!='W03' & id!='W14' & id!='W25',.(id)])
# still prob w/ var on both sides
movewetOUT <- dat[id!='W03'  & id!='W14' & id!='W25',Move(case_, log_sl, ToD_start, wet_end, ttd1_adj, cos_ta, stepjum), by=.(id)]
m.movewet <- merge(movewetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

m.movewet.coef <- movewetOUT[term=='coef',-'AIC']
m.movewet.coef <- merge(m.movewet.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

ggplot(melt(m.movewet.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) # +
  #ylim(-1,1)


# not working
movewet.td.OUT <- dat[id!='W03'  & id!='W14' & id!='W25',Move.td(case_, log_sl, ToD_start, wet_end, as.factor(wtd1), cos_ta, stepjum), by=.(id)]
m.movewet <- merge(movewetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

movepropwetOUT <- dat[id!='W03'  & id!='W14' & id!='W25',Move(case_, log_sl, ToD_start, propwet_end, ttd1_adj, cos_ta, stepjum), by=.(id)]
m.movepropwet <- merge(movepropwetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

moveOUT <- dat[wtd1==0 | wtd1==4, Move(case_, log_sl, ToD_start, land_end, as.factor(wtd1), cos_ta, stepjum), by = id]
m.move <- merge(moveOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.move.p<-m.move[term=='coef' | term=='p', .(id, term, `sl:log(ttd + 1)`, `log(ttd + 1):ta`, COD)]
m.move.aic<-m.move[term=='coef' | term=='se']
m.move.coef<-m.move[term=='coef']

m.move.coef <- moveOUT[term=='coef',-'AIC']
m.move.coef <- merge(m.move.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

ggplot(melt(m.move.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) +
  ylim(-1,1)

ggplot(melt(moveOUT)) +
  geom_density(aes(value), fill = 'dodgerblue', alpha =0.5) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  facet_wrap(~variable, scale = "free")


movepropOUT <- dat[, Move.prop(case_, log_sl, ToD_start, propwet_end, propopen_end, ttd1, cos_ta, as.integer(stepjum)), by = id]
m.moveprop <- merge(movepropOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.moveprop.p<-m.moveprop[term=='coef' | term=='p', .(id, term, `sl:log(ttd + 1)`, `log(ttd + 1):ta`, COD)]



#### HABITAT ####
Habitat <- function(y, sl, ToD, land, ttd, parkdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                    log(ttd+1):land + log(ttd+1):log(parkdist+1) + strata(strata1))
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

Habitat.prop <- function(y, sl, ToD, wet, open, ttd, parkdist, parkYN, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + wet + open +
                    sl:wet + sl:open +
                    log(ttd+1):wet + log(ttd+1):open + log(ttd+1):parkdist:parkYN + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}

habwetOUT <- dat[id!='W03' & id!='W14' & id!='W25', Habitat(case_, log_sl, ToD_start, wet_end, ttd1_adj, parkDistadj_end, stepjum), by = id]

m.habwet <- merge(habwetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

m.habwet.coef <- habwetOUT[term=='coef',-'AIC']
m.habwet.coef <- merge(m.habwet.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

ggplot(melt(m.habwet.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter(aes(color=COD))  +
  ylim(-1,1)



habOUT <- dat[, Habitat(case_, log_sl, ToD_start, land_end, as.factor(wtd1), parkDistAdj_end, parkYN_end, stepjum), by = id]
m.habitat <- merge(habOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitat.p <- m.habitat[term=='coef'|term=='p', .(id, term, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
                       `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, COD)]
m.habitat.aic <- m.habitat[term=='coef'|term=='se', .(id, term, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
                                                   `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, AIC, COD)]
m.habitat.coef <- m.habitat[term=='coef']
m.hab.coef <- habOUT[term=='coef',-'AIC']
m.hab.coef<- merge(m.hab.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)


ggplot(melt(m.hab.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) +
  ylim(-1,1)

ggplot(melt(m.hab.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter() +
  ylim(-10,5)

ggplot(melt(habOUT)) +
  geom_density(aes(value), fill = 'dodgerblue', alpha =0.5) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  facet_wrap(~variable, scale = "free")


habpropOUT <- dat[, Habitat.prop(case_, log_sl, ToD_start, propwet_end, propopen_end, ttd1, lnparkdist_end, parkYN_end, stepjum), by = id]
m.habitatprop <- merge(habpropOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.habitatprop.p <- m.habitatprop[term=='coef'|term=='p', .(id, term, `wet:log(ttd + 1)`, `open:log(ttd + 1)`, 
                                                  `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, COD)]
m.habitatprop.coef <- m.habitatprop[term=='coef']
m.habprop.coef <- habpropOUT[term=='coef',-'AIC']


ggplot(melt(m.habprop.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter() +
  ylim(-10,5)


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



Social <- function(y, sl, ToD, land, ttd, nndist, packdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + land +
                    sl:land +
                    log(ttd+1):log(nndist+1) + log(ttd+1):log(packdist+1) + strata(strata1))
  sum.model <- summary(model)$coefficients
  # Transpose the coef of the model and cast as data.table
  term <- c('coef','hr','se','z','p')
  coefOut <- data.table(t(sum.model))
  
  # Return combined columns
  # print(summary(model))
  print(AIC(model))
  return(data.table(term, coefOut, AIC=AIC(model)))
}


Social.prop <- function(y, sl, ToD, wet, open, ttd, nndist, packYN, packdist, strata1) {
  # Make the model
  model <- clogit(y ~ sl:ToD + wet + open +
                    sl:wet + sl:open +
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

dat[,unique(id)] 
dat[id!='W03' & id!='W14' & id!='W25' & id!='W15' & id!='W19',unique(id)]
# W13, W11, W27,  not converge

### only wet hab now
socOUT <- dat[id!='W03' & id!='W14' & id!='W25' & id!='W15'  & id!='W19', Social(case_, log_sl, ToD_start, wet_end, ttd1, distance2, packDistadj_end, stepjum), by = id]
m.social <- merge(socOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
#m.social.p <- m.social[term=='coef' |term=='p', .(id, term, `log(ttd + 1):log(nndist + 1)`, `log(ttd + 1):log(packdist + 1):packYNpack`, `log(ttd + 1):log(packdist + 1):packYNout-pack`, COD)]
m.social.aic <- m.social[term=='coef' |term=='se']
m.social.coef <- m.social[term=='coef']

m.soc.coef <- socOUT[term=='coef',-'AIC']
m.soc.coef<- merge(m.soc.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

ggplot(melt(m.soc.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) #+ ylim(-1,1)


plot(term)
ggplot(melt(m.soc.coef)) + aes(variable, value) +
       geom_boxplot() +
       geom_jitter() +
      ylim(-10,5)

socpropOUT <- dat[id !='W05' & id != 'W03' & id != 'W07' & id != 'W19' & id != 'W15' & id != 'W14', Social.prop(case_, log_sl, ToD_start, propwet_end, propopen_end, ttd1, distance2, packYN_end, packDist_end, as.integer(stepjum)), by = id]
m.socialprop <- merge(socpropOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
m.socialprop.p <- m.socialprop[term=='coef' |term=='p', .(id, term, `log(ttd + 1):log(nndist + 1)`, `log(ttd + 1):log(packdist + 1):packYNpack`, `log(ttd + 1):log(packdist + 1):packYNout-pack`, COD)]

ggplot(melt(socOUT)) +
  geom_density(aes(value), fill = 'dodgerblue', alpha =0.5) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  facet_wrap(~variable, scale = "free")

dat.focal<- dat.meta[status=='dead',.(WolfID, death=as.Date(end_date))]
dat.focal <- dat.focal[WolfID %chin% unique(dat$id)] 
range(dat.focal$death)
plot(dat.focal$death)

bbmle::AICtab(m.core$model, m.move$model, m.habitat.parkdist$model, m.social$model)



range(dat$step_id_)
range(as.integer(dat$stepjum))

dat.steps <- dat[,.(id, step_id_, stepjum, stepjumnum = as.integer(stepjum))]






