### Analyses -- Mixed effects ====
# Julie Turner
# Started: November 18 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival','forcats', 'glmmTMB', 'tibble', 'bbmle')
lapply(libs, require, character.only = TRUE)

### Function ----
se <- function(x){
  sd(x, na.rm = T)/ sqrt(length(na.omit(x)))
}


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


dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet', 
                               ifelse(dat$land_end == 'mixed'|dat$land_end == 'deciduous'|dat$land_end == 'coniferous', 'forest','open'))

dat[,'forest']
# dat[,'propforest_end_adj'] <- dat$propconif_end+dat$propmixed_end +dat$propdecid_end 
# dat[,'propopen_end_adj'] <- dat$propopen_end +dat$propshrub_end + dat$propurban_end
# 
# 
# dat[,mean(propwet_end), by=.(wolfID,ua)]
# dat[,mean(propforest_end_adj), by=.(wolfID,ua)]
# dat[,mean(propopen_end_adj), by=.(wolfID,ua)]
# dat<-dat[COD == 'cdv'|COD == 'human'|COD == 'none']

dat[,'wolf_step_id'] <- paste(dat$wolfID, dat$step_id_, sep = '_')


#### behav predict COD ####

total <- mclogit::mclogit(COD ~ log_sl:ToD_start +
                  log_sl:land_end_adj +
                  log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                  land_end_adj + log(1+roadDist_end) +
                  log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                  log(1+distance2) + log(1+packDistadj_end) +
                  log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                  strata(wolf_step_id), dat )


#### full model #### 
#### all time ####
#### human ####
dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
#dat.soc <- dat[wolfID %chin% dat.wnn$wolfID]
full.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                      log_sl:land_end_adj +
                      log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                      (1|wolf_step_id) +
                      (0 + (log_sl)|wolfID) +
                      (0 + (cos_ta)|wolfID) +
                      (0 + (log(ttd1+1):log_sl)|wolfID) +
                      (0 + (log(ttd1+1):cos_ta)|wolfID) +
                      land_end_adj + #log(1+roadDist_end) +
                      log(ttd1+1):land_end_adj + # log(ttd1+1):log(1+roadDist_end) +
                      
                      (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                      #(0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                      
                      log(1+distance2) + log(1+packDistadj_end) +
                      log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                      
                      (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                      (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                    , family=poisson(),
                    data = dat[wolfID %chin% dat.wnn$wolfID & COD == 'human'], doFit=FALSE)

full.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.human$parameters$theta)
full.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.human <- glmmTMB:::fitTMB(full.human)
pophuman <- summary(full.human)

summary(full.human)$coef$cond[-1, "Estimate"]
saveRDS(pophuman, 'data/derived-data/pophuman_noroad.Rds')
summary(full.human)$varcor


full.all.indiv.human <- coef(full.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")

full.human.both <- glmmTMB(case_ ~ log_sl:ToD_start +
                        log_sl:land_end_adj +
                        log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                        (1|wolf_step_id) +
                        (0 + (log_sl)|wolfID) +
                        (0 + (cos_ta)|wolfID) +
                        (0 + (log(ttd1+1):log_sl)|wolfID) +
                        (0 + (log(ttd1+1):cos_ta)|wolfID) +
                        land_end_adj + log(1+roadDist_end) +
                        log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                        
                        (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                        (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                        
                        log(1+packDistadj_end) +
                        log(ttd1+1):log(1+packDistadj_end) +
                        (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                      , family=poisson(),
                      data = dat[COD == 'human'], doFit=FALSE)

full.human.both$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.human.both$parameters$theta)
full.human.both$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.human.both <- glmmTMB:::fitTMB(full.human.both)
pophuman.both <- summary(full.human.both)

summary(full.human.both)$coef$cond[-1, "Estimate"]
saveRDS(pophuman.both, 'data/derived-data/pophumanboth.Rds')
summary(full.human.both)$varcor


full.all.indiv.human.both <- coef(full.human.both)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")

#### cdv ####
#dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
#dat.soc <- dat[wolfID %chin% dat.wnn$wolfID]
full.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                      log_sl:land_end_adj +
                      log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                      (1|wolf_step_id) +
                      (0 + (log_sl)|wolfID) +
                      (0 + (cos_ta)|wolfID) +
                      (0 + (log(ttd1+1):log_sl)|wolfID) +
                      (0 + (log(ttd1+1):cos_ta)|wolfID) +
                      land_end_adj  + log(1+roadDist_end) +
                      log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                      
                      (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                      (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                      
                      log(1+distance2) + log(1+packDistadj_end) +
                      log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                      
                      (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                      (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                    , family=poisson(),
                    data = dat[wolfID %chin% dat.wnn$wolfID & COD == 'cdv'], doFit=FALSE)

full.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.cdv$parameters$theta)
full.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.cdv <- glmmTMB:::fitTMB(full.cdv)
summary(full.cdv)
popcdv<- summary(full.cdv)$coef$cond[-1, 1:2]
saveRDS(popcdv, 'data/derived-data/popcdv.Rds')

summary(full.cdv)$coef$cond[-1, "Estimate"]
summary(full.cdv)$varcor


full.all.indiv.cdv <- coef(full.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



#### none ####

full.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                      log_sl:land_end_adj +
                      log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                      (1|wolf_step_id) +
                      (0 + (log_sl)|wolfID) +
                      (0 + (cos_ta)|wolfID) +
                      (0 + (log(ttd1+1):log_sl)|wolfID) +
                      (0 + (log(ttd1+1):cos_ta)|wolfID) +
                      land_end_adj + log(1+roadDist_end) +
                      log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                      
                      (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                      (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                      
                      log(1+distance2) + log(1+packDistadj_end) +
                      log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                      
                      (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                      (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                    , family=poisson(),
                    data = dat[wolfID %chin% dat.wnn$wolfID & COD == 'none'], doFit=FALSE)

full.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.none$parameters$theta)
full.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.none <- glmmTMB:::fitTMB(full.none)
summary(full.none)

summary(full.none)$coef$cond[-1, "Estimate"]
summary(full.none)$varcor


full.all.indiv.none <- coef(full.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### Last month ####
dat.wnn.lastmo <- dat[ttd1<=31 & !is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
### human ####
full.lastmo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                             #log_sl:land_end_adj +
                             log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                             (1|wolf_step_id) +
                             (0 + (log_sl)|wolfID) +
                             (0 + (cos_ta)|wolfID) +
                             (0 + (log(ttd1+1):log_sl)|wolfID) +
                             (0 + (log(ttd1+1):cos_ta)|wolfID) +
                             land_end_adj + log(1+roadDist_end) +
                             log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                             
                             (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                             (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                             
                             log(1+distance2) + log(1+packDistadj_end) +
                             log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                             
                             (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                             (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                           , family=poisson(),
                           data = dat[ttd1<=31 &wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'human'], doFit=FALSE)

full.lastmo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.lastmo.human$parameters$theta)
full.lastmo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.lastmo.human <- glmmTMB:::fitTMB(full.lastmo.human)
summary(full.lastmo.human)

summary(full.lastmo.human)$coef$cond[-1, "Estimate"]
summary(full.lastmo.human)$varcor


full.lastmo.all.indiv.human <- coef(full.lastmo.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



#### cdv ####
full.lastmo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                             #log_sl:land_end_adj +
                             log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                             (1|wolf_step_id) +
                             (0 + (log_sl)|wolfID) +
                             (0 + (cos_ta)|wolfID) +
                             (0 + (log(ttd1+1):log_sl)|wolfID) +
                             (0 + (log(ttd1+1):cos_ta)|wolfID) +
                             land_end_adj + log(1+roadDist_end) +
                             log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                             
                             (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                             (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                             
                             log(1+distance2) + log(1+packDistadj_end) +
                             log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                             
                             (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                             (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                           , family=poisson(),
                           data = dat[ttd1<=31 &wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'cdv'], doFit=FALSE)

full.lastmo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.lastmo.cdv$parameters$theta)
full.lastmo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.lastmo.cdv <- glmmTMB:::fitTMB(full.lastmo.cdv)
summary(full.lastmo.cdv)

popcdv.lastmo <- summary(full.lastmo.cdv)$coef$cond[-1, 1:2]
saveRDS(popcdv.lastmo, 'data/derived-data/popcdv_lastmo.Rds')
summary(full.lastmo.cdv)$coef$cond[-1, "Estimate"]
summary(full.lastmo.cdv)$varcor


full.lastmo.all.indiv.cdv <- coef(full.lastmo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### none ####
full.lastmo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                             #log_sl:land_end_adj +
                             log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                             (1|wolf_step_id) +
                             (0 + (log_sl)|wolfID) +
                             (0 + (cos_ta)|wolfID) +
                             (0 + (log(ttd1+1):log_sl)|wolfID) +
                             (0 + (log(ttd1+1):cos_ta)|wolfID) +
                             land_end_adj + log(1+roadDist_end) +
                             log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                             
                             (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                             (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                             
                             log(1+distance2) + log(1+packDistadj_end) +
                             log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                             
                             (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                             (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                           , family=poisson(),
                           data = dat[ttd1<=31 &wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'none'], doFit=FALSE)

full.lastmo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.lastmo.none$parameters$theta)
full.lastmo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.lastmo.none <- glmmTMB:::fitTMB(full.lastmo.none)
summary(full.lastmo.none)

summary(full.lastmo.none)$coef$cond[-1, "Estimate"]
summary(full.lastmo.none)$varcor


full.lastmo.all.indiv.none <- coef(full.lastmo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



#### Time split ####
#### intx COD ####

full.2mo.COD <- glmmTMB(case_ ~ log_sl:ToD_start +
                            log_sl:land_end_adj +
                            cos_ta + 
                            (1|wolf_step_id) + 
                            (0 + log_sl|wolfID) + (0 + (log_sl:COD)|wolfID) +
                            (0 + cos_ta|wolfID) + (0 + (cos_ta:COD)|wolfID) +
                            land_end_adj + log(1+roadDist_end) +
                            
                            
                            (0 + land_end_adj|wolfID) + (0 + (land_end_adj:COD)|wolfID) + 
                            (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(1+roadDist_end):COD)|wolfID) +
                            
                            log(1+distance2) + log(1+packDistadj_end) +
                            (0 + (log(1+distance2))|wolfID) + (0 + (log(1+distance2):COD)|wolfID) + 
                            (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(1+packDistadj_end):COD)|wolfID)
                          , family=poisson(), 
                          data = dat[ttd1>31], doFit=FALSE)



full.2mo.COD$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.2mo.COD$parameters$theta)
full.2mo.COD$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.2mo.COD <- glmmTMB:::fitTMB(full.2mo.COD)
summary(full.2mo.COD)

summary(full.2mo.COD)$coef$cond[-1, "Estimate"]
summary(full.2mo.COD)$varcor


full.indiv.2mo.COD <- coef(full.2mo.COD)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



full.1mo.COD <- glmmTMB(case_ ~ log_sl:ToD_start +
                          log_sl:land_end_adj +
                          cos_ta + 
                          (1|wolf_step_id) + 
                          (0 + log_sl|wolfID) + (0 + (log_sl:COD)|wolfID) +
                          (0 + cos_ta|wolfID) + (0 + (cos_ta:COD)|wolfID) +
                          land_end_adj + log(1+roadDist_end) +
                          
                          
                          (0 + land_end_adj|wolfID) + (0 + (land_end_adj:COD)|wolfID) + 
                          (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(1+roadDist_end):COD)|wolfID) +
                          
                          log(1+distance2) + log(1+packDistadj_end) +
                          (0 + (log(1+distance2))|wolfID) + (0 + (log(1+distance2):COD)|wolfID) + 
                          (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(1+packDistadj_end):COD)|wolfID)
                        , family=poisson(), 
                        data = dat[ttd1<=31], doFit=FALSE)



full.1mo.COD$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.1mo.COD$parameters$theta)
full.1mo.COD$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.1mo.COD <- glmmTMB:::fitTMB(full.1mo.COD)
summary(full.1mo.COD)

summary(full.1mo.COD)$coef$cond[-1, "Estimate"]
summary(full.1mo.COD)$varcor


full.indiv.1mo.COD <- coef(full.1mo.COD)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



#### human ####
full.2mo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                          log_sl:land_end_adj +
                          cos_ta + 
                          (1|wolf_step_id) + 
                          (0 + log_sl|wolfID) + 
                          (0 + cos_ta|wolfID) +
                          land_end_adj + log(1+roadDist_end) +
                          
                          
                          (0 + land_end_adj|wolfID) + 
                          (0 + (log(1+roadDist_end))|wolfID) +
                          
                          log(1+distance2) + log(1+packDistadj_end) +
                          
                          
                          (0 + (log(1+distance2))|wolfID) + 
                          (0 + (log(1+packDistadj_end))|wolfID) 
                        , family=poisson(), 
                        data = dat[ttd1>31 & COD == 'human'], doFit=FALSE)



full.2mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.2mo.human$parameters$theta)
full.2mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.2mo.human <- glmmTMB:::fitTMB(full.2mo.human)
summary(full.2mo.human)

summary(full.2mo.human)$coef$cond[-1, "Estimate"]
summary(full.2mo.human)$varcor


full.indiv.2mo.human <- coef(full.2mo.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



full.1mo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                            log_sl:land_end_adj +
                            cos_ta + 
                            (1|wolf_step_id) + 
                            (0 + log_sl|wolfID) + 
                            (0 + cos_ta|wolfID) +
                            land_end_adj + log(1+roadDist_end) +
                            
                            
                            (0 + land_end_adj|wolfID) + 
                            (0 + (log(1+roadDist_end))|wolfID) +
                            
                            log(1+distance2) + log(1+packDistadj_end) +
                            
                            
                            (0 + (log(1+distance2))|wolfID) + 
                            (0 + (log(1+packDistadj_end))|wolfID) 
                          , family=poisson(), 
                          data = dat[ttd1<=31 & COD == 'human'], doFit=FALSE)



full.1mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.1mo.human$parameters$theta)
full.1mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.1mo.human <- glmmTMB:::fitTMB(full.1mo.human)
summary(full.1mo.human)

summary(full.1mo.human)$coef$cond[-1, "Estimate"]
summary(full.1mo.human)$varcor


full.indiv.1mo.human <- coef(full.1mo.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### cdv ####

full.2mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                          log_sl:land_end_adj +
                          cos_ta + 
                          (1|wolf_step_id) + 
                          (0 + log_sl|wolfID) + 
                          (0 + cos_ta|wolfID) +
                          land_end_adj + log(1+roadDist_end) +
                          
                          
                          (0 + land_end_adj|wolfID) + 
                          (0 + (log(1+roadDist_end))|wolfID) +
                          
                          log(1+distance2) + log(1+packDistadj_end) +
                          
                          
                          (0 + (log(1+distance2))|wolfID) + 
                          (0 + (log(1+packDistadj_end))|wolfID) 
                        , family=poisson(), 
                        data = dat[ttd1>31 & COD == 'cdv'], doFit=FALSE)



full.2mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.2mo.cdv$parameters$theta)
full.2mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.2mo.cdv <- glmmTMB:::fitTMB(full.2mo.cdv)
summary(full.2mo.cdv)

summary(full.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(full.2mo.cdv)$varcor


full.indiv.2mo.cdv <- coef(full.2mo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


full.1mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                          log_sl:land_end_adj +
                          log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                          (1|wolf_step_id) + 
                          (0 + (log(ttd1+1):log_sl)|wolfID) + 
                          (0 + (log(ttd1+1):cos_ta)|wolfID) +
                          land_end_adj + log(1+roadDist_end) +
                          log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                          
                          (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                          (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                          
                          log(1+distance2) + log(1+packDistadj_end) +
                          log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                          
                          (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                          (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                        , family=poisson(), 
                        data = dat[ttd1<=31 & COD == 'cdv'], doFit=FALSE)

full.1mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                          log_sl:land_end_adj +
                          cos_ta + 
                          (1|wolf_step_id) + 
                          (0 + log_sl|wolfID) + 
                          (0 + cos_ta|wolfID) +
                          land_end_adj + log(1+roadDist_end) +
                          
                          
                          (0 + land_end_adj|wolfID) + 
                          (0 + (log(1+roadDist_end))|wolfID) +
                          
                          log(1+distance2) + log(1+packDistadj_end) +
                          
                          
                          (0 + (log(1+distance2))|wolfID) + 
                          (0 + (log(1+packDistadj_end))|wolfID) 
                        , family=poisson(), 
                        data = dat[ttd1<=31 & COD == 'cdv'], doFit=FALSE)


full.1mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.1mo.cdv$parameters$theta)
full.1mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.1mo.cdv <- glmmTMB:::fitTMB(full.1mo.cdv)
summary(full.1mo.cdv)
confint(summary(full.2mo.cdv))

summary(full.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(full.1mo.cdv)$varcor


full.indiv.1mo.cdv <- coef(full.1mo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")





#### none ####
full.2mo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                            log_sl:land_end_adj +
                            cos_ta + 
                            (1|wolf_step_id) + 
                            (0 + log_sl|wolfID) + 
                            (0 + cos_ta|wolfID) +
                            land_end_adj + log(1+roadDist_end) +
                            
                            
                            (0 + land_end_adj|wolfID) + 
                            (0 + (log(1+roadDist_end))|wolfID) +
                            
                            log(1+distance2) + log(1+packDistadj_end) +
                            
                            
                            (0 + (log(1+distance2))|wolfID) + 
                            (0 + (log(1+packDistadj_end))|wolfID) 
                          , family=poisson(), 
                          data = dat[ttd1>31 & COD == 'none'], doFit=FALSE)



full.2mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.2mo.none$parameters$theta)
full.2mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.2mo.none <- glmmTMB:::fitTMB(full.2mo.none)
summary(full.2mo.none)

summary(full.2mo.none)$coef$cond[-1, "Estimate"]
summary(full.2mo.none)$varcor


full.indiv.2mo.none <- coef(full.2mo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



full.1mo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                            log_sl:land_end_adj +
                            cos_ta + 
                            (1|wolf_step_id) + 
                            (0 + log_sl|wolfID) + 
                            (0 + cos_ta|wolfID) +
                            land_end_adj + log(1+roadDist_end) +
                            
                            
                            (0 + land_end_adj|wolfID) + 
                            (0 + (log(1+roadDist_end))|wolfID) +
                            
                            log(1+distance2) + log(1+packDistadj_end) +
                            
                            
                            (0 + (log(1+distance2))|wolfID) + 
                            (0 + (log(1+packDistadj_end))|wolfID) 
                          , family=poisson(), 
                          data = dat[ttd1<=31 & COD == 'none'], doFit=FALSE)



full.1mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(full.1mo.none$parameters$theta)
full.1mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
full.1mo.none <- glmmTMB:::fitTMB(full.1mo.none)
summary(full.1mo.none)

summary(full.1mo.none)$coef$cond[-1, "Estimate"]
summary(full.1mo.none)$varcor


full.indiv.1mo.none <- coef(full.1mo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



#### MOMVEMENT ####
#### human ####
move.human <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                            #log_sl:land_end_adj +
                            log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                            #(0 + (log_sl:land_end_adj)|wolfID) +
                            (1|wolf_step_id) +
                            (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                            (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                          , family=poisson(), 
                          data = dat[COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
# TMBStruc$parameters$theta[4] = log(1e3) 
# TMBStruc$mapArg = list(theta=factor(c(1:3, NA)))
nvar_parm <- length(move.human$parameters$theta)
move.human$parameters$theta[1] <- log(1e3)

move.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.human.fit <- glmmTMB:::fitTMB(move.human)
summary(move.human.fit)

move.indiv.human <- coef(move.human.fit)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")

move.2mo.human <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                           #log_sl:land_end_adj +
                            log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                            #(0 + (log_sl:land_end_adj)|wolfID) +
                            (1|wolf_step_id) +
                            (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                            (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                         , family=poisson(), 
                         data = dat[ttd1>31 & COD == 'human'], doFit=FALSE)


move.2mo.human <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                            #log_sl:land_end_adj +
                             cos_ta + 
                            #(0 + (log_sl:land_end_adj)|wolfID) +
                            (1|wolf_step_id) +
                            (0 + log_sl|wolfID)  + 
                            (0 + cos_ta|wolfID) 
                          , family=poisson(), 
                          data = dat[ttd1>31 & COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
# TMBStruc$parameters$theta[4] = log(1e3) 
# TMBStruc$mapArg = list(theta=factor(c(1:3, NA)))
nvar_parm <- length(move.2mo.human$parameters$theta)
move.2mo.human$parameters$theta[1] <- log(1e3)
#nvar_parm <- length(move.2mo.human$parameters$theta)
move.2mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.2mo.human.fit <- glmmTMB:::fitTMB(move.2mo.human)
summary(move.2mo.human.fit)

move.indiv.2mo.human <- coef(move.2mo.human.fit)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



move.1mo.human <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                            #log_sl:land_end_adj +
                            log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                            #(0 + (log_sl:land_end_adj)|wolfID) +
                            (1|wolf_step_id) +
                            (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                            (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                         , family=poisson(), 
                         data = dat[ttd1<=31 & COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
move.1mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(move.1mo.human$parameters$theta)
move.1mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.1mo.human <- glmmTMB:::fitTMB(move.1mo.human)
summary(move.1mo.human)

move.indiv.1mo.human <- coef(move.1mo.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### cdv ####
move.2mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                          #log_sl:land_end_adj +
                          log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                          #(0 + (log_sl:land_end_adj)|wolfID) +
                          (1|wolf_step_id) +
                          (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                          (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                       , family=poisson(), 
                     data = dat[ttd1>31 & COD == 'cdv'], doFit=FALSE)


#' Set variance of random intercept to 10^6
move.2mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(move.2mo.cdv$parameters$theta)
move.2mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.2mo.cdv <- glmmTMB:::fitTMB(move.2mo.cdv)
summary(move.2mo.cdv)

summary(move.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(move.2mo.cdv)$varcor


move.indiv.2mo.cdv <- coef(move.2mo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#dat.1mo.cdv <- dat[ttd1<=31 & COD == 'cdv']

move.1mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                          #log_sl:land_end_adj +
                          log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                          #(0 + (log_sl:land_end_adj)|wolfID) +
                          (1|wolf_step_id) +
                          (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                          (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                        , family=poisson(), 
                        data = dat[ttd1<=31 & COD == 'cdv'], doFit=FALSE)

#' Set variance of random intercept to 10^6
move.1mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(move.1mo.cdv$parameters$theta)
move.1mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.1mo.cdv <- glmmTMB:::fitTMB(move.1mo.cdv)
summary(move.1mo.cdv)

summary(move.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(move.1mo.cdv)$varcor


move.indiv.1mo.cdv <- coef(move.1mo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")

#### none ####
move.2mo.none <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                           #log_sl:land_end_adj +
                           log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                           #(0 + (log_sl:land_end_adj)|wolfID) +
                           (1|wolf_step_id) +
                           (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                           (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                        , family=poisson(), 
                        data = dat[ttd1>31 & COD == 'none'], doFit=FALSE)

#' Set variance of random intercept to 10^6
move.2mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(move.2mo.none$parameters$theta)
move.2mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.2mo.none <- glmmTMB:::fitTMB(move.2mo.none)
summary(move.2mo.none)

summary(move.2mo.none)$coef$cond[-1, "Estimate"]
summary(move.2mo.none)$varcor

move.indiv.2mo.none <- coef(move.2mo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



move.1mo.none <- glmmTMB(case_ ~ log_sl:ToD_start + #log_sl +
                           #log_sl:land_end_adj +
                           log(1+ttd1):log_sl + cos_ta + log(1+ttd1):cos_ta +
                           #(0 + (log_sl:land_end_adj)|wolfID) +
                           (1|wolf_step_id) +
                           (0 + log_sl|wolfID) + (0 + log(1+ttd1):log_sl|wolfID) + 
                           (0 + cos_ta|wolfID) + (0 + log(1+ttd1):cos_ta|wolfID)
                         , family=poisson(), 
                         data = dat[ttd1<=31 & COD == 'none'], doFit=FALSE)

#' Set variance of random intercept to 10^6
move.1mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(move.1mo.none$parameters$theta)
move.1mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
move.1mo.none <- glmmTMB:::fitTMB(move.1mo.none)
summary(move.1mo.none)

summary(move.1mo.none)$coef$cond[-1, "Estimate"]
summary(move.1mo.none)$varcor

move.indiv.1mo.none <- coef(move.1mo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")





#### HABITAT ####
#### human ####
habitat.2mo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                            land_end_adj + log(1+roadDist_end) +
                            log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                            (1|wolf_step_id) + 
                            (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                            (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID)
                           ,family=poisson(), 
                          data = dat[ttd1>31 & COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
habitat.2mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(habitat.2mo.human$parameters$theta)
habitat.2mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
habitat.2mo.human <- glmmTMB:::fitTMB(habitat.2mo.human)
summary(habitat.2mo.human)

summary(habitat.2mo.human)$coef$cond[-1, "Estimate"]
summary(habitat.2mo.human)$varcor

habitat.indiv.2mo.human <- coef(habitat.2mo.human)$cond$wolfID[ , -1] %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



habitat.1mo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                               land_end_adj + log(1+roadDist_end) +
                               log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                               (1|wolf_step_id) + 
                               (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                               (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID)
                             , family=poisson(), 
                          data = dat[ttd1<=31 & COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
habitat.1mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(habitat.1mo.human$parameters$theta)
habitat.1mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
habitat.1mo.human <- glmmTMB:::fitTMB(habitat.1mo.human)
summary(habitat.1mo.human)

summary(habitat.1mo.human)$coef$cond[-1, "Estimate"]
summary(habitat.1mo.human)$varcor

habitat.indiv.1mo.human <- coef(habitat.1mo.human)$cond$wolfID[ , -1] %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### cdv ####
habitat.2mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                             land_end_adj + log(1+roadDist_end) +
                             log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                             (1|wolf_step_id) + 
                             (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                             (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID)
                        , family=poisson(), 
                        data = dat[ttd1>31 & COD == 'cdv'], doFit=FALSE)


#' Set variance of random intercept to 10^6
habitat.2mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(habitat.2mo.cdv$parameters$theta)
habitat.2mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
habitat.2mo.cdv <- glmmTMB:::fitTMB(habitat.2mo.cdv)
summary(habitat.2mo.cdv)

summary(habitat.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(habitat.2mo.cdv)$varcor


habitat.indiv.2mo.cdv <- coef(habitat.2mo.cdv)$cond$wolfID[ , -1] %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


dat.1mo.cdv <- dat[ttd1<=31 & COD == 'cdv']

habitat.1mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                             land_end_adj + log(1+roadDist_end) +
                             log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                             (1|wolf_step_id) + 
                             (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                             (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID)
                           , family=poisson(), 
                        data = dat[ttd1<=31 & COD == 'cdv'], doFit=FALSE)

#' Set variance of random intercept to 10^6
habitat.1mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(habitat.1mo.cdv$parameters$theta)
habitat.1mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
habitat.1mo.cdv <- glmmTMB:::fitTMB(habitat.1mo.cdv)
summary(habitat.1mo.cdv)

summary(habitat.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(habitat.1mo.cdv)$varcor


habitat.indiv.1mo.cdv <- coef(habitat.1mo.cdv)$cond$wolfID[ , -1] %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")

#### none ####
habitat.2mo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                              land_end_adj + log(1+roadDist_end) +
                              log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                              (1|wolf_step_id) + 
                              (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                              (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID)
                            , family=poisson(), 
                         data = dat[ttd1>31 & COD == 'none'], doFit=FALSE)

#' Set variance of random intercept to 10^6
habitat.2mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(habitat.2mo.none$parameters$theta)
habitat.2mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
habitat.2mo.none <- glmmTMB:::fitTMB(habitat.2mo.none)
summary(habitat.2mo.none)

summary(habitat.2mo.none)$coef$cond[-1, "Estimate"]
summary(habitat.2mo.none)$varcor

habitat.indiv.2mo.none <- coef(habitat.2mo.none)$cond$wolfID[ , -1] %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



habitat.1mo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                              log_sl:land_end_adj +
                              land_end_adj + log(1+roadDist_end) +
                              log(ttd1+1):land_end_adj + log(ttd1+1):log(1+roadDist_end) +
                              (1|wolf_step_id) + 
                              (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) + 
                              (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID)
                            , family=poisson(), 
                         data = dat[ttd1<=31 & COD == 'none'], doFit=FALSE)

#' Set variance of random intercept to 10^6
habitat.1mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(habitat.1mo.none$parameters$theta)
habitat.1mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
habitat.1mo.none <- glmmTMB:::fitTMB(habitat.1mo.none)
summary(habitat.1mo.none)

summary(habitat.1mo.none)$coef$cond[-1, "Estimate"]
summary(habitat.1mo.none)$varcor

habitat.indiv.1mo.none <- coef(habitat.1mo.none)$cond$wolfID[ , -1] %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")





#### SOCIAL ####

# remove those that don't have packmates at time of analysis
dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.soc <- dat[wolfID %chin% dat.wnn$wolfID]

#### human ####
social.2mo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                               log(1+distance2) + log(1+packDistadj_end) +
                               log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                               (1|wolf_step_id) + 
                               (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                               (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                             ,family=poisson(), 
                             data = dat.soc[ttd1>31 & COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
social.2mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(social.2mo.human$parameters$theta)
social.2mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
social.2mo.human <- glmmTMB:::fitTMB(social.2mo.human)
summary(social.2mo.human)

summary(social.2mo.human)$coef$cond[-1, "Estimate"]
summary(social.2mo.human)$varcor

social.indiv.2mo.human <- coef(social.2mo.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



social.1mo.human <- glmmTMB(case_ ~ log_sl:ToD_start +
                              log(1+distance2) + log(1+packDistadj_end) +
                              log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                              (1|wolf_step_id) + 
                              (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                              (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                            ,family=poisson(), 
                            data = dat.soc[ttd1<=31 & COD == 'human'], doFit=FALSE)

#' Set variance of random intercept to 10^6
social.1mo.human$parameters$theta[1] <- log(1e3)
nvar_parm <- length(social.1mo.human$parameters$theta)
social.1mo.human$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
social.1mo.human <- glmmTMB:::fitTMB(social.1mo.human)
summary(social.1mo.human)

summary(social.1mo.human)$coef$cond[-1, "Estimate"]
summary(social.1mo.human)$varcor

social.indiv.1mo.human <- coef(social.1mo.human)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### cdv ####
social.2mo.cdv <- glmmTMB(case_ ~log_sl:ToD_start +
                            log(1+distance2) + log(1+packDistadj_end) +
                            log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                            (1|wolf_step_id) + 
                            (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                            (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                          ,family=poisson(), 
                          data = dat.soc[ttd1>31 & COD == 'cdv'], doFit=FALSE)


#' Set variance of random intercept to 10^6
social.2mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(social.2mo.cdv$parameters$theta)
social.2mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
social.2mo.cdv <- glmmTMB:::fitTMB(social.2mo.cdv)
summary(social.2mo.cdv)

summary(social.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(social.2mo.cdv)$varcor


social.indiv.2mo.cdv <- coef(social.2mo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#dat.1mo.cdv <- dat[ttd1<=31 & COD == 'cdv']

social.1mo.cdv <- glmmTMB(case_ ~ log_sl:ToD_start +
                            log(1+distance2) + log(1+packDistadj_end) +
                            log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                            (1|wolf_step_id) + 
                            (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                            (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                          ,family=poisson(), 
                          data = dat.soc[ttd1<=31 & COD == 'cdv'], doFit=FALSE)

#' Set variance of random intercept to 10^6
social.1mo.cdv$parameters$theta[1] <- log(1e3)
nvar_parm <- length(social.1mo.cdv$parameters$theta)
social.1mo.cdv$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
social.1mo.cdv <- glmmTMB:::fitTMB(social.1mo.cdv)
summary(social.1mo.cdv)

summary(social.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(social.1mo.cdv)$varcor


social.indiv.1mo.cdv <- coef(social.1mo.cdv)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")

#### none ####
social.2mo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                             log(1+distance2) + log(1+packDistadj_end) +
                             log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                             (1|wolf_step_id) + 
                             (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                             (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                           ,family=poisson(), 
                           data = dat.soc[ttd1>31 & COD == 'none'], doFit=FALSE)

#' Set variance of random intercept to 10^6
social.2mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(social.2mo.none$parameters$theta)
social.2mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
social.2mo.none <- glmmTMB:::fitTMB(social.2mo.none)
summary(social.2mo.none)

summary(social.2mo.none)$coef$cond[-1, "Estimate"]
summary(social.2mo.none)$varcor

social.indiv.2mo.none <- coef(social.2mo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")



social.1mo.none <- glmmTMB(case_ ~ log_sl:ToD_start +
                             log(1+distance2) + log(1+packDistadj_end) +
                             log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                             (1|wolf_step_id) + 
                             (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) + 
                             (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                           ,family=poisson(), 
                           data = dat.soc[ttd1<=31 & COD == 'none'], doFit=FALSE)

#' Set variance of random intercept to 10^6
social.1mo.none$parameters$theta[1] <- log(1e3)
nvar_parm <- length(social.1mo.none$parameters$theta)
social.1mo.none$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
social.1mo.none <- glmmTMB:::fitTMB(social.1mo.none)
summary(social.1mo.none)

summary(social.1mo.none)$coef$cond[-1, "Estimate"]
summary(social.1mo.none)$varcor

social.indiv.1mo.none <- coef(social.1mo.none)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
  pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
  mutate(method = "ME")


#### GATHERING RESULTS ####
#### none ####
#move
summary(move.2mo.none)$coef$cond[-1, "Estimate"]
summary(move.2mo.none)$varcor
move.indiv.2mo.none

summary(move.1mo.none)$coef$cond[-1, "Estimate"]
summary(move.1mo.none)$varcor
move.indiv.1mo.none

move.indiv.none <- merge(setDT(move.indiv.2mo.none)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                         setDT(move.indiv.1mo.none)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

#habitat 
summary(habitat.2mo.none)$coef$cond[-1, "Estimate"]
summary(habitat.2mo.none)$varcor
habitat.indiv.2mo.none

summary(habitat.1mo.none)$coef$cond[-1, "Estimate"]
summary(habitat.1mo.none)$varcor
habitat.indiv.1mo.none

habitat.indiv.none <- merge(setDT(habitat.indiv.2mo.none)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                            setDT(habitat.indiv.1mo.none)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))


#social
summary(social.2mo.none)$coef$cond[-1, "Estimate"]
summary(social.2mo.none)$varcor
social.indiv.2mo.none

summary(social.1mo.none)$coef$cond[-1, "Estimate"]
summary(social.1mo.none)$varcor
social.indiv.1mo.human

social.indiv.none <- merge(setDT(social.indiv.2mo.none)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                           setDT(social.indiv.1mo.none)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))



#### humans ####
#move
summary(move.2mo.human)$coef$cond[-1, "Estimate"]
summary(move.2mo.human)$varcor
move.indiv.2mo.human

summary(move.1mo.human)$coef$cond[-1, "Estimate"]
#move.indiv.1mo.human[,se(estimate), by = .(term)]
summary(move.1mo.human)$varcor
move.indiv.1mo.human

move.indiv.human <- merge(setDT(move.indiv.2mo.human)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                         setDT(move.indiv.1mo.human)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))



#habitat 
summary(habitat.2mo.human)$coef$cond[-1, "Estimate"]
summary(habitat.2mo.human)$varcor
habitat.indiv.2mo.human

summary(habitat.1mo.human)$coef$cond[-1, "Estimate"]
summary(habitat.1mo.human)$varcor
habitat.indiv.1mo.human

habitat.indiv.human <- merge(setDT(habitat.indiv.2mo.human)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                          setDT(habitat.indiv.1mo.human)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

#social
summary(social.2mo.human)$coef$cond[-1, "Estimate"]
summary(social.2mo.human)$varcor
social.indiv.2mo.human

summary(social.1mo.human)$coef$cond[-1, "Estimate"]
summary(social.1mo.human)$varcor
social.indiv.1mo.human

social.indiv.human <- merge(setDT(social.indiv.2mo.human)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                          setDT(social.indiv.1mo.human)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))



#### cdv ####
#move
summary(move.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(move.2mo.cdv)$varcor
move.indiv.2mo.cdv

summary(move.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(move.1mo.cdv)$varcor
move.indiv.1mo.cdv

move.indiv.cdv <- merge(setDT(move.indiv.2mo.cdv)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                        setDT(move.indiv.1mo.cdv)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

#habitat 
summary(habitat.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(habitat.2mo.cdv)$varcor
habitat.indiv.2mo.cdv

summary(habitat.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(habitat.1mo.cdv)$varcor
habitat.indiv.1mo.cdv

habitat.indiv.cdv <- merge(setDT(habitat.indiv.2mo.cdv)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                        setDT(habitat.indiv.1mo.cdv)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

#social
summary(social.2mo.cdv)$coef$cond[-1, "Estimate"]
summary(social.2mo.cdv)$varcor
social.indiv.2mo.cdv

summary(social.1mo.cdv)$coef$cond[-1, "Estimate"]
summary(social.1mo.cdv)$varcor
social.indiv.1mo.cdv

social.indiv.cdv <- merge(setDT(social.indiv.2mo.cdv)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                        setDT(social.indiv.1mo.cdv)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))





#### combine ####
move.indiv.none[,'COD'] <- 'none'
move.indiv.human[,'COD'] <- 'human'
move.indiv.cdv[,'COD'] <- 'cdv'

move.indiv <- rbind(move.indiv.none, move.indiv.human, move.indiv.cdv)
move.indiv$COD <- factor(move.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
move.indiv[,'delta'] <- move.indiv$control - move.indiv$case


cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")

ggplot(move.indiv[term=="log(1 + ttd1):log_sl"|term=="log(1 + ttd1):cos_ta"], aes(term, (delta), fill = COD)) +
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
  scale_color_manual(values = cbPalette) #+ ylim(-.3,.3)

melt(move.indiv)

ggplot(melt(move.indiv)[term=="log(1 + ttd1):log_sl"|term=="log(1 + ttd1):cos_ta"&COD=='human' &variable != 'delta'], aes(term, (value), fill = variable)) +
  geom_boxplot(aes(fill = variable),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = variable),
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
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-.3,.3)



social.indiv.none[,'COD'] <- 'none'
social.indiv.human[,'COD'] <- 'human'
social.indiv.cdv[,'COD'] <- 'cdv'

social.indiv <- rbind(social.indiv.none, social.indiv.human, social.indiv.cdv)
social.indiv$COD <- factor(social.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
social.indiv[,'delta'] <- social.indiv$control - social.indiv$case


cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")

ggplot(social.indiv[term=="log(ttd1 + 1):log(1 + distance2)"|term=="log(ttd1 + 1):log(1 + packDistadj_end)"], aes(term, (delta), fill = COD)) +
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
  ylab('delta socialment') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-.05,.05)



#### full ####
#### intx COD ####

full.indiv.COD <- merge(setDT(full.indiv.2mo.COD)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                          setDT(full.indiv.1mo.COD)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

full.COD <- merge(full.indiv.COD[case!=0|control!=0], dat.meta[,.(wolfpop, COD)], by.x = 'wolfID', by.y = 'wolfpop', all.x = T)


full.COD$COD <- factor(full.COD$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
full.COD[,'delta_abs'] <- abs(full.COD$control) - full.COD$case
full.COD[,'delta'] <- full.COD$control - full.COD$case

# full.COD <- full.COD[term=='cos_ta'|term=='land_end_adjforest'|term=='land_end_adjopen'|term=='land_end_adjwet'|
#                            term=='log(1 + distance2)'|term=='log(1 + packDistadj_end)'|
#                            term=='log(1 + roadDist_end)'|term=='log_sl']
full.terms <- (unique(full.COD$term))
full.terms <- c("log_sl", "cos_ta", "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
                "log(1 + distance2)", "log(1 + packDistadj_end)")
full.COD$term <- factor(full.COD$term, levels = full.terms, labels = c("log_sl", "cos_ta", "forest", "open", "wet", "roadDist",
                                                                           "nnDist", "boundaryDist"))
minmax <- setDT(full.COD)[,.(min= min(case), max=max(case)), by = .(term, COD)]
full.COD.betas <- minmax[min!= max]
full.COD.betas.names <- unique(full.COD.betas$term)


full.COD.betas <- full.COD[term %like% "COD", ]

ggplot(full.COD.betas, aes(term, (control), fill = COD)) +
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
  ylab('beta') +
  ggtitle("a) control") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)# + ylim(-2,2)

ggplot(full.COD.betas, aes(term, (case), fill = COD)) +
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
  ylab('beta') +
  ggtitle("b) case") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)# + ylim(-2,2)


ggplot(full.COD.betas, aes(term, (delta), fill = COD)) +
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
  ylab('delta beta') +
  ggtitle("c) change between control and case") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)# + ylim(-2,2)


full.COD.betas.names <- c("cos_ta:CODnone", "cos_ta", 
                          "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
                          "log(1 + distance2)", "log(1 + packDistadj_end)",
                          "log(ttd1 + 1):log_sl", "log(ttd1 + 1):cos_ta", 
                          "log(ttd1 + 1):land_end_adjforest", "log(ttd1 + 1):land_end_adjopen", "log(ttd1 + 1):land_end_adjwet", "log(ttd1 + 1):log(1 + roadDist_end)", 
                          "log(ttd1 + 1):log(1 + distance2)", "log(ttd1 + 1):log(1 + packDistadj_end)")

#[term %chin% full.all.betas.names]

# full.betas <- c("log_sl", "cos_ta", "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
#                 "log(1 + distance2)", "log(1 + packDistadj_end)")

full.all.indiv.betas$term <- factor(full.all.indiv.betas$term, levels = full.all.betas.names, labels = c("log_sl", "cos_ta", "open", "wet", "roadDist",
                                                                                                         "nnDist", "boundaryDist",
                                                                                                         "log_sl-ttd", "cos_ta-ttd", "forest-ttd", "open-ttd", "wet-ttd", "roadDist-ttd",
                                                                                                         "nnDist-ttd", "boundaryDist-ttd"))




#### COD separate ####
confint(full.2mo.human$fit)
full.indiv.human <- merge(setDT(full.indiv.2mo.human)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                        setDT(full.indiv.1mo.human)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

full.indiv.cdv <- merge(setDT(full.indiv.2mo.cdv)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                        setDT(full.indiv.1mo.cdv)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))

full.indiv.none <- merge(setDT(full.indiv.2mo.none)[term != '(Intercept)',.(wolfID, term, control = estimate)],
                        setDT(full.indiv.1mo.none)[term != '(Intercept)',.(wolfID, term, case = estimate)], by= c('wolfID', 'term'))


full.indiv.none[,'COD'] <- 'none'
full.indiv.human[,'COD'] <- 'human'
full.indiv.cdv[,'COD'] <- 'cdv'

full.indiv <- rbind(full.indiv.none, full.indiv.human, full.indiv.cdv)
full.indiv$COD <- factor(full.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
full.indiv[,'delta_abs'] <- abs(full.indiv$control) - full.indiv$case
full.indiv[,'delta'] <- full.indiv$control - full.indiv$case

full.indiv <- full.indiv[term=='cos_ta'|term=='land_end_adjforest'|term=='land_end_adjopen'|term=='land_end_adjwet'|
                           term=='log(1 + distance2)'|term=='log(1 + packDistadj_end)'|
                           term=='log(1 + roadDist_end)'|term=='log_sl']
full.terms <- (unique(full.indiv$term))
full.terms <- c("log_sl", "cos_ta", "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
                "log(1 + distance2)", "log(1 + packDistadj_end)")
full.indiv$term <- factor(full.indiv$term, levels = full.terms, labels = c("log_sl", "cos_ta", "forest", "open", "wet", "roadDist",
                                                                           "nnDist", "boundaryDist"))


#saveRDS(full.indiv, 'data/derived-data/full_cc_betas.Rds')

cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")

p.delta <- ggplot(full.indiv, aes(term, (delta), fill = COD)) +
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
  ylab('delta beta') +
  ggtitle("c) change between control and case") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)# + ylim(-2,2)


p.control <- ggplot(full.indiv, aes(term, (control), fill = COD)) +
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
  ylab('beta') +
  ggtitle("a) control") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-2,2)



p.case <- ggplot(full.indiv, aes(term, (case), fill = COD)) +
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
  ylab('beta') +
  ggtitle("b) case") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-2,2)

require(patchwork)
p.control/p.case/p.delta


#### full all ####
summary(full.cdv)
full.all.indiv.none[,'COD'] <- 'none'
full.all.indiv.human[,'COD'] <- 'human'
full.all.indiv.cdv[,'COD'] <- 'cdv'

full.all.indiv <- rbind(full.all.indiv.none, full.all.indiv.human, full.all.indiv.cdv)
full.all.indiv$COD <- factor(full.all.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))

minmax <- setDT(full.all.indiv)[,.(min= min(estimate), max=max(estimate)), by = .(term, COD)]
full.all.betas <- minmax[min!= max]
full.all.betas.names <- unique(full.all.betas$term)

full.all.betas.names <- c("log_sl", "cos_ta", 
                           "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
                "log(1 + distance2)", "log(1 + packDistadj_end)",
                "log(ttd1 + 1):log_sl", "log(ttd1 + 1):cos_ta", 
                "log(ttd1 + 1):land_end_adjforest", "log(ttd1 + 1):land_end_adjopen", "log(ttd1 + 1):land_end_adjwet", "log(ttd1 + 1):log(1 + roadDist_end)", 
                "log(ttd1 + 1):log(1 + distance2)", "log(ttd1 + 1):log(1 + packDistadj_end)")

full.all.indiv.betas <- full.all.indiv[term %chin% full.all.betas.names]
# full.betas <- c("log_sl", "cos_ta", "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
#                 "log(1 + distance2)", "log(1 + packDistadj_end)")

full.all.indiv.betas$term <- factor(full.all.indiv.betas$term, levels = full.all.betas.names, labels = c("log_sl", "cos_ta", "open", "wet", "roadDist",
                                                                           "nnDist", "boundaryDist",
                                                                           "log_sl-ttd", "cos_ta-ttd", "forest-ttd", "open-ttd", "wet-ttd", "roadDist-ttd",
                                                                           "nnDist-ttd", "boundaryDist-ttd"))

#saveRDS(full.all.indiv.betas, 'data/derived-data/full_betas.Rds')
full.ttd <- full.all.indiv.betas[term %like% "ttd", ]


ggplot(full.ttd, aes(term, (estimate), fill = COD)) +
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
  ylab('beta') +
 # ggtitle("b) case") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-2,2)

#### full last mo ####

full.lastmo.all.indiv.none[,'COD'] <- 'none'
full.lastmo.all.indiv.human[,'COD'] <- 'human'
full.lastmo.all.indiv.cdv[,'COD'] <- 'cdv'

full.lastmo.indiv <- rbind(full.lastmo.all.indiv.none, full.lastmo.all.indiv.human, full.lastmo.all.indiv.cdv)
full.lastmo.indiv$COD <- factor(full.lastmo.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))

minmax <- setDT(full.lastmo.indiv)[,.(min= min(estimate), max=max(estimate)), by = .(term, COD)]
full.lastmo.betas <- minmax[min!= max]
full.lastmo.betas.names <- unique(full.lastmo.betas$term)

full.lastmo.betas.names <- c("log_sl", "cos_ta", 
                          "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
                          "log(1 + distance2)", "log(1 + packDistadj_end)",
                          "log(ttd1 + 1):log_sl", "log(ttd1 + 1):cos_ta", 
                          "log(ttd1 + 1):land_end_adjforest", "log(ttd1 + 1):land_end_adjopen", "log(ttd1 + 1):land_end_adjwet", "log(ttd1 + 1):log(1 + roadDist_end)", 
                          "log(ttd1 + 1):log(1 + distance2)", "log(ttd1 + 1):log(1 + packDistadj_end)")

full.lastmo.indiv.betas <- full.lastmo.indiv[term %chin% full.lastmo.betas.names]
# full.betas <- c("log_sl", "cos_ta", "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
#                 "log(1 + distance2)", "log(1 + packDistadj_end)")

full.lastmo.indiv.betas$term <- factor(full.lastmo.indiv.betas$term, levels = full.lastmo.betas.names, labels = c("log_sl", "cos_ta", "open", "wet", "roadDist",
                                                                                                         "nnDist", "boundaryDist",
                                                                                                         "log_sl-ttd", "cos_ta-ttd", "forest-ttd", "open-ttd", "wet-ttd", "roadDist-ttd",
                                                                                                         "nnDist-ttd", "boundaryDist-ttd"))

#saveRDS(full.lastmo.indiv.betas, 'data/derived-data/full_lastmo_betas.Rds')

full.lastmo.ttd <- full.lastmo.indiv.betas[term %like% "ttd", ]


ggplot(full.lastmo.ttd, aes(term, (estimate), fill = COD)) +
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
  ylab('beta') +
  # ggtitle("b) case") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-.15,.15)




