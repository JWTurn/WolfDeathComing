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
dat.RMNP <- readRDS('data/derived-data/ssfAllCov.Rds')
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

dat[,'wtd1'] <- as.integer(dat$ttd1/7)

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'))
dat.meta[,'wolfpop'] <- paste(dat.meta$pop, dat.meta$WolfID, sep = '_')



dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')

dat<- merge(dat, dat.meta, by.x = c('id', 'pop'), by.y = c('WolfID', 'pop'))


dat[,'ttd1_adj'] <- ifelse(dat$COD =='none', dat$ttd1*730, dat$ttd1)

dat.meta[,uniqueN(wolfpop), by= .(status, COD)]

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
# model <- clogit(case_ ~ log_sl:ToD_start + land_end + 
#                   #sl:closed + sl:open + sl:wet + 
#                   log_sl:land_end +
#                   strata(stepjum), data = dat[id=='W06'])
# sum.mod <- summary(model)
# sum.mod$conf.int

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




dat[,mean(propwet_end), by=.(wolfID,ua)]
dat[,mean(propconif_end), by=.(wolfID,ua)]
dat[,mean(propmixed_end), by=.(wolfID,ua)]
dat[,mean(propopen_end), by=.(wolfID,ua)]
    
#run iSSA model by ID
# Core Model


unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='RMNP_W04' & wolfID!='RMNP_W05' & wolfID!='GHA26_W09',.(wolfID)])
#id!='W13', id!='W03' & id!='W14' & id!='W25' 
#wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & wolfID!='RMNP_W12'
unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37',.(wolfID)])
core2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & wolfID!='RMNP_W14' & 
                wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                Core(case_, log_sl, ToD_start, land_end, stepjum), by = .(wolfID)]

unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & 
             wolfID!='GHA26_W31' & wolfID!='GHA26_W37',.(wolfID)])
# & wolfID!='RMNP_W27' 
core1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & wolfID!='RMNP_W14' &
                   wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                 Core(case_, log_sl, ToD_start, land_end, stepjum), by = .(wolfID)]



coreOUT<- dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='RMNP_W04' & wolfID!='RMNP_W05' & wolfID!='GHA26_W09', Core(case_, log_sl, ToD_start, land_end, stepjum), by = .(wolfID)]
m.core <- merge(coreOUT, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)


m.core.co <- coreOUT[term=='coef',-'AIC']
m.core.co <-merge(m.core.co, dat.meta, by.x = 'wolfID', by.y = 'WolfID', all.x = T)

ggplot(melt(m.core.co)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter() +
  ylim(-1,1)

ggplot(melt(coreOUT)) +
  geom_density(aes(value), fill = 'dodgerblue', alpha =0.5) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  facet_wrap(~variable, scale = "free")

# 
corewet2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                Core(case_, log_sl, ToD_start, wet_end, stepjum), by = .(wolfID)]


corewet1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Core(case_, log_sl, ToD_start, wet_end, stepjum), by = .(wolfID)]




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



move2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Move(case_, log_sl, ToD_start, land_end, ttd1_adj, cos_ta, stepjum), by = .(wolfID)]


move1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Move(case_, log_sl, ToD_start, land_end, ttd1_adj, cos_ta, stepjum), by = .(wolfID)]




movewet2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Move(case_, log_sl, ToD_start, wet_end, ttd1_adj, cos_ta, stepjum), by = .(wolfID)]


movewet1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Move(case_, log_sl, ToD_start, wet_end, ttd1_adj, cos_ta, stepjum), by = .(wolfID)]




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


unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
             wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & 
             wolfID!='GHA26_W37' & wolfID!='RMNP_W05' & wolfID!='RMNP_W09' & wolfID!='RMNP_W11' ,.(wolfID)])
hab2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                   wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & 
                  wolfID!='GHA26_W37'  & wolfID!='RMNP_W05' & wolfID!='RMNP_W09' & wolfID!='RMNP_W11'  , 
                 Habitat(case_, log_sl, ToD_start, land_end, ttd1_adj, parkDistadj_end, stepjum), by = .(wolfID)]


   
unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
             wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' & 
             wolfID!='RMNP_W15' & wolfID!='RMNP_W09' & wolfID!='RMNP_W11' & wolfID!='RMNP_W13'
           ,.(wolfID)])
hab1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                   wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' & 
                  wolfID!='RMNP_W15' & wolfID!='RMNP_W09' & wolfID!='RMNP_W11' & wolfID!='RMNP_W13', 
                 Habitat(case_, log_sl, ToD_start, land_end, ttd1_adj, parkDistadj_end, stepjum), by = .(wolfID)]




habwet2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Habitat(case_, log_sl, ToD_start, wet_end, ttd1_adj, parkDistadj_end, stepjum), by = .(wolfID)]


habwet1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37', 
                    Habitat(case_, log_sl, ToD_start, wet_end, ttd1_adj, parkDistadj_end, stepjum), by = .(wolfID)]





# habwetOUT <- dat[id!='W03' & id!='W14' & id!='W25', Habitat(case_, log_sl, ToD_start, wet_end, ttd1_adj, parkDistadj_end, stepjum), by = id]
# 
# m.habwet <- merge(habwetOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
# 
# m.habwet.coef <- habwetOUT[term=='coef',-'AIC']
# m.habwet.coef <- merge(m.habwet.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)

ggplot(melt(m.habwet.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter(aes(color=COD))  +
  ylim(-1,1)



# habOUT <- dat[, Habitat(case_, log_sl, ToD_start, land_end, as.factor(wtd1), parkDistAdj_end, parkYN_end, stepjum), by = id]
# m.habitat <- merge(habOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
# m.habitat.p <- m.habitat[term=='coef'|term=='p', .(id, term, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
#                        `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, COD)]
# m.habitat.aic <- m.habitat[term=='coef'|term=='se', .(id, term, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
#                                                    `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, AIC, COD)]
# m.habitat.coef <- m.habitat[term=='coef']
# m.hab.coef <- habOUT[term=='coef',-'AIC']
# m.hab.coef<- merge(m.hab.coef, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)


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


# habpropOUT <- dat[, Habitat.prop(case_, log_sl, ToD_start, propwet_end, propopen_end, ttd1, lnparkdist_end, parkYN_end, stepjum), by = id]
# m.habitatprop <- merge(habpropOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
# m.habitatprop.p <- m.habitatprop[term=='coef'|term=='p', .(id, term, `wet:log(ttd + 1)`, `open:log(ttd + 1)`, 
#                                                   `log(ttd + 1):parkdist:parkYNpark`, `log(ttd + 1):parkdist:parkYNout-park`, COD)]
# m.habitatprop.coef <- m.habitatprop[term=='coef']
# m.habprop.coef <- habpropOUT[term=='coef',-'AIC']
# 

ggplot(melt(m.habprop.coef)) + aes(variable, value) +
  geom_boxplot() +
  geom_jitter() +
  ylim(-10,5)


# Habitat.land <- function(y, sl, ToD, land, ttd, strata1) {
#   # Make the model
#   model <- clogit(y ~ sl:ToD + land + sl:land + 
#                     log(ttd+1):land + strata(strata1))
#   sum.model <- summary(model)$coefficients
#   # Transpose the coef of the model and cast as data.table
#   term <- c('coef','hr','se','z','p')
#   coefOut <- data.table(t(sum.model))
#   
#   # Return combined columns
#   # print(summary(model))
#   print(AIC(model))
#   return(data.table(term, coefOut, AIC=AIC(model)))
# }
# 
# hablandOUT <- dat[, Habitat.land(case_, log_sl, ToD_start, land_end, ttd1, stepjum), by = id]
# m.habitat.land <- merge(hablandOUT, dat.meta, by.x = 'id', by.y = 'WolfID', all.x = T)
# m.habitat.land.p <- m.habitat.land[term=='p', .(id, `landclosed:log(ttd + 1)`, `landopen:log(ttd + 1)`, `landwet:log(ttd + 1)`, 
#                                       `log(ttd + 1):parkdist:parkYNpark`,`log(ttd + 1):parkdist:parkYNout-park`, COD)]

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


dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
     wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' & 
      wolfID!='RMNP_W05' & wolfID!='RMNP_W13' & wolfID!='RMNP_W15' & wolfID!='RMNP_W19'
   ,unique(wolfID)]

soc2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                   wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' & 
                  wolfID!='RMNP_W05' & wolfID!='RMNP_W13' & wolfID!='RMNP_W15' & wolfID!='RMNP_W19', 
                 Social(case_, log_sl, ToD_start, land_end, ttd1_adj, distance1, packDistadj_end, stepjum), by = .(wolfID)]

unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
             wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' & 
             wolfID!='RMNP_W15' & wolfID!='RMNP_W19'
           ,.(wolfID)])


soc1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                   wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' &
                  wolfID!='RMNP_W15' & wolfID!='RMNP_W19', 
                 Social(case_, log_sl, ToD_start, land_end, ttd1_adj, distance1, packDistadj_end, stepjum), by = .(wolfID)]


unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
             wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' & 
             wolfID!='RMNP_W13' & wolfID!='RMNP_W15' & wolfID!='RMNP_W19'
           ,.(wolfID)])


socwet2moOUT<- dat[ttd1>31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='RMNP_W27' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' 
                     & wolfID!='RMNP_W13' & wolfID!='RMNP_W15' & wolfID!='RMNP_W19', 
                    Social(case_, log_sl, ToD_start, wet_end, ttd1_adj, distance1, packDistadj_end, stepjum), by = .(wolfID)]


unique(dat[wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
             wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37' 
             & wolfID!='RMNP_W15' & wolfID!='RMNP_W19'
           ,.(wolfID)])


socwet1moOUT<- dat[ttd1<=31 & wolfID!='GHA26_W03' & wolfID!='RMNP_W03' & wolfID!='GHA26_W09' & 
                      wolfID!='RMNP_W14' & wolfID!='RMNP_W25' & wolfID!='GHA26_W31' & wolfID!='GHA26_W37'
                   & wolfID!='RMNP_W15' & wolfID!='RMNP_W19', 
                    Social(case_, log_sl, ToD_start, wet_end, ttd1_adj, distance1, packDistadj_end, stepjum), by = .(wolfID)]



#### AICs ####
m.movewet2mo.aic <- unique(movewet2moOUT[,.(AIC_move2mo =AIC), by=.(wolfID)])
m.movewet1mo.aic <- unique(movewet1moOUT[,.(AIC_move1mo =AIC), by=.(wolfID)])
move.aic <- merge(m.movewet2mo.aic, m.movewet1mo.aic, all = T)
move.aic <- merge(move.aic, dat.meta[,.(COD), by = .(wolfpop)], by.x = 'wolfID', by.y = 'wolfpop')

m.habwet2mo.aic <- unique(habwet2moOUT[,.(AIC_hab2mo =AIC), by=.(wolfID)])
m.habwet1mo.aic <- unique(habwet1moOUT[,.(AIC_hab1mo =AIC), by=.(wolfID)])
hab.aic <- merge(m.habwet2mo.aic, m.habwet1mo.aic, all = T)
hab.aic <- merge(hab.aic, dat.meta[,.(COD), by = .(wolfpop)], by.x = 'wolfID', by.y = 'wolfpop')

m.socwet2mo.aic <- unique(socwet2moOUT[,.(AIC_soc2mo =AIC), by=.(wolfID)])
m.socwet1mo.aic <- unique(socwet1moOUT[,.(AIC_soc1mo =AIC), by=.(wolfID)])
soc.aic <- merge(m.socwet2mo.aic, m.socwet1mo.aic, all = T)
soc.aic <- merge(soc.aic, dat.meta[,.(COD), by = .(wolfpop)], by.x = 'wolfID', by.y = 'wolfpop')

m.all.aic <- merge(move.aic, hab.aic, by=c('wolfID', 'COD'), all = T)
m.all.aic <- merge(m.all.aic, soc.aic, by=c('wolfID', 'COD'), all = T)




#### GRAPHS ####
m.movewet2mo.coef <- movewet2moOUT[term=='coef',-'AIC']
m.movewet2mo.coef <- m.movewet2mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movewet2mo.coef<- merge(m.movewet2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movewet2mo.coef <- melt(m.movewet2mo.coef)
m.movewet2mo.coef[,'ttd'] <- '2mo'

m.movewet1mo.coef <- movewet1moOUT[term=='coef',-'AIC']
m.movewet1mo.coef <- m.movewet1mo.coef[, .( wolfID , lnSL = `sl:log(ttd + 1)`, cosTA = `log(ttd + 1):ta`)]
m.movewet1mo.coef<- merge(m.movewet1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.movewet1mo.coef<- melt(m.movewet1mo.coef)
m.movewet1mo.coef[,'ttd'] <- '1mo'

m.movewet.coef <- rbind(m.movewet1mo.coef, m.movewet2mo.coef)
m.movewet.coef[,'test'] <- ifelse(m.movewet.coef$ttd =='1mo', 'case', 'control')
m.movewet.coef$test <- factor(m.movewet.coef$test, levels = c('control','case'))

m.movewet.coef.cdv <- m.movewet.coef[COD=='cdv']

#color = c("#0072B2", "#D55E00", "#009E73")
color = c("darkviolet", "aquamarine4")

ggplot(m.movewet.coef.cdv, aes(variable, value, fill = test)) +
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
  scale_color_manual(values = color) +
  ylim(-.25,.25)

ggplot(m.movewet.coef) + aes(variable, value, fill=ttd) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) + ylim(-.3,.3)


m.movewet.coef.cdv <- m.movewet.coef[COD=='cdv']

ggplot(m.movewet.coef.cdv) + aes(variable, value, fill=ttd) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= ttd, color = ttd)) + ylim(-.25,.2)


ggplot(m.movewet1mo.coef, aes(variable, value)) +
  geom_boxplot(aes( fill=COD)) +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= COD, color = COD)) +
  ylim(-.25,.2)

ggplot(m.movewet2mo.coef, aes(variable, value)) +
  geom_boxplot(aes( fill=COD)) +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= COD, color = COD))+
  ylim(-.25,.2)


#### habitat graphs ####

m.habwet2mo.coef <- habwet2moOUT[term=='coef',-'AIC']
m.habwet2mo.coef <- m.habwet2mo.coef[, .( wolfID , wet = `land:log(ttd + 1)`, parkDist = `log(ttd + 1):log(parkdist + 1)`)]
m.habwet2mo.coef<- merge(m.habwet2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habwet2mo.coef <- melt(m.habwet2mo.coef)
m.habwet2mo.coef[,'ttd'] <- '2mo'

m.habwet1mo.coef <- habwet1moOUT[term=='coef',-'AIC']
m.habwet1mo.coef <- m.habwet1mo.coef[, .( wolfID , wet = `land:log(ttd + 1)`, parkDist = `log(ttd + 1):log(parkdist + 1)`)]
m.habwet1mo.coef<- merge(m.habwet1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.habwet1mo.coef<- melt(m.habwet1mo.coef)
m.habwet1mo.coef[,'ttd'] <- '1mo'

m.habwet.coef <- rbind(m.habwet1mo.coef, m.habwet2mo.coef)
m.habwet.coef[,'test'] <- ifelse(m.habwet.coef$ttd =='1mo', 'case', 'control')
m.habwet.coef$test <- factor(m.habwet.coef$test, levels = c('control','case'))

m.habwet.coef.wet <- m.habwet.coef[variable=='wet' & COD == 'cdv']
m.habwet.coef.park <- m.habwet.coef[variable=='parkDist' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.habwet.coef.wet, aes(variable, value, fill = test)) +
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
  scale_color_manual(values = color) +
  ylim(-1,1)


ggplot(m.habwet.coef.park, aes(variable, value, fill = test)) +
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
  scale_color_manual(values = color) +
  ylim(-.1,.15)


ggplot(m.habwet.coef) + aes(variable, value, fill=ttd) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) + ylim(-1,1)


m.habwet.coef.cdv <- m.habwet.coef[COD=='cdv']

ggplot(m.habwet.coef.cdv) + aes(variable, value, fill=ttd) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= ttd, color = ttd)) + ylim(-.75,.85)




ggplot(m.habwet1mo.coef, aes(variable, value)) +
  geom_boxplot(aes( fill=COD)) +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= COD, color = COD)) +
  ylim(-.7,.5)

m.habwet2mo.coef.g <- m.habwet2mo.coef[COD!='human']
ggplot(m.habwet2mo.coef.g, aes(variable, value)) +
  geom_boxplot(aes( fill=COD)) +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= COD, color = COD))+
  ylim(-.7,1)


#### social graphs ####

m.socwet2mo.coef <- socwet2moOUT[term=='coef',-'AIC']
m.socwet2mo.coef <- m.socwet2mo.coef[, .( wolfID , nnXttd = `log(ttd + 1):log(nndist + 1)`, packDistXttd = `log(ttd + 1):log(packdist + 1)`)]
m.socwet2mo.coef<- merge(m.socwet2mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.socwet2mo.coef <- melt(m.socwet2mo.coef)
m.socwet2mo.coef[,'ttd'] <- '2mo'

m.socwet1mo.coef <- socwet1moOUT[term=='coef',-'AIC']
m.socwet1mo.coef <- m.socwet1mo.coef[, .(wolfID , nnXttd = `log(ttd + 1):log(nndist + 1)`, packDistXttd = `log(ttd + 1):log(packdist + 1)`)]
m.socwet1mo.coef<- merge(m.socwet1mo.coef, dat.meta, by.x = 'wolfID', by.y = 'wolfpop', all.x = T)
m.socwet1mo.coef<- melt(m.socwet1mo.coef)
m.socwet1mo.coef[,'ttd'] <- '1mo'

m.socwet.coef <- rbind(m.socwet1mo.coef, m.socwet2mo.coef)
m.socwet.coef[,'test'] <- ifelse(m.socwet.coef$ttd =='1mo', 'case', 'control')
m.socwet.coef$test <- factor(m.socwet.coef$test, levels = c('control','case'))

m.socwet.coef.nn <- m.socwet.coef[variable =='nnXttd' & COD == 'cdv']
m.socwet.coef.pack <- m.socwet.coef[variable =='packDistXttd' & COD == 'cdv']


#color = c("#0072B2", "#D55E00", "#009E73")

ggplot(m.socwet.coef.nn, aes(variable, value, fill = test)) +
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
  scale_color_manual(values = color) +
  ylim(-.25,.05)


ggplot(m.socwet.coef.pack, aes(variable, value, fill = test)) +
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
  scale_color_manual(values = color) +
  ylim(-.3,.75)


ggplot(m.socwet.coef) + aes(variable, value, fill=ttd) +
  geom_boxplot() +
  geom_jitter(aes(color=COD)) + ylim(-.5,.5)


m.socwet.coef.cdv <- m.socwet.coef[COD=='cdv']

ggplot(m.socwet.coef.cdv) + aes(variable, value, fill=ttd) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= ttd, color = ttd)) + ylim(-.4,.75)


m.socwet1mo.coef.g <- m.socwet1mo.coef[COD!='wolf']
ggplot(m.socwet1mo.coef.g, aes(variable, value)) +
  geom_boxplot(aes( fill=COD)) +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= COD, color = COD)) +
  ylim(-.3,.75)

ggplot(m.socwet2mo.coef, aes(variable, value)) +
  geom_boxplot(aes( fill=COD)) +
  geom_dotplot(binaxis = 'y', stackdir='center',
               position=position_dodge(1) , aes( fill= COD, color = COD))+
  ylim(-.3,.2)



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





