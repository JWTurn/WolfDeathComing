### Analyses -- Mixed effects ====
# Julie Turner
# Started: November 18 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival','forcats', 'glmmTMB', 'tibble', 'bbmle', 'patchwork', 'broom.mixed')
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
dat.RMNP <- readRDS('data/derived-data/ssfAllCov_2mo_RMNP.Rds')
dat.GHA26 <- readRDS('data/derived-data/ssfAllCov_2mo_GHA26.Rds')


dat.RMNP[,'pop'] <- 'RMNP'
dat.RMNP$wolfID <- paste(dat.RMNP$pop, dat.RMNP$id, sep = '_')
dat.GHA26[,'pop'] <- 'GHA26'
dat.GHA26$wolfID <- paste(dat.GHA26$pop, dat.GHA26$id, sep = '_')
dat<-rbind(dat.RMNP, dat.GHA26, fill=T)

dat[ttd1>=0 & ttd2>=0, unique(wolfID)]

#dat<-dat[ttd1>=0 & ttd2>=0]

#dat[,'wtd1'] <- as.integer(dat$ttd1/7)

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'))
dat.meta[,'wolfpop'] <- paste(dat.meta$pop, dat.meta$WolfID, sep = '_')

params <- readRDS('data/derived-data/moveParams_2mo_all.Rds')

dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')

dat<- merge(dat, dat.meta, by.x = c('id', 'pop'), by.y = c('WolfID', 'pop'))

dat$packDist_end <- ifelse(dat$packDistadj_end >=0, dat$packDistadj_end, 0)

summary(dat[,.(land_end)])
summary(dat[case_==TRUE,.(land_end)])

# dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet', 
#                                ifelse(dat$land_end == 'mixed'|dat$land_end == 'deciduous', 'forest',
#                                       ifelse(dat$land_end == 'coniferous', 'coniferous', 'open')))
dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet',
                               ifelse(dat$land_end == 'mixed'|dat$land_end == 'deciduous'|dat$land_end == 'coniferous'|dat$land_end == 'shrub',
                                      'forest','open'))

#dat[,'forest']
dat[,'propforest_end_adj'] <- dat$propconif_end+dat$propmixed_end +dat$propdecid_end +dat$propshrub_end # remove shrub?
dat[,'propopen_end_adj'] <- dat$propopen_end + dat$propurban_end 

# 
# dat[,mean(propwet_end), by=.(wolfID,ua)]
# dat[,mean(propforest_end_adj), by=.(wolfID,ua)]
# dat[,mean(propopen_end_adj), by=.(wolfID,ua)]
# dat<-dat[COD == 'cdv'|COD == 'human'|COD == 'none']

dat[,'wolf_step_id'] <- paste(dat$wolfID, dat$step_id_, sep = '_')

dat[ua=='used',.(road=mean(roadDist_end, na.rm = T), nn= mean(distance2, na.rm = T), pack= mean(packDistadj_end, na.rm = T)), by= COD]
dat[!is.na(distance2), .N, by=.(ua)]
dat[is.na(distance2), .N, by=.(ua)]

dat$COD <- factor(dat$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
dat$COD[is.na(dat$COD)] <- "control"

dat$ToD_start <- as.factor(dat$ToD_start)
dat$land_end_adj <- as.factor(dat$land_end_adj)

#### only wolves with packmates ####
dat[,uniqueN(step_id_), by=.(wolfID)]
dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo <- dat[ttd1<=31 & !is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo.cod <- merge(dat.wnn.lastmo, dat.meta[,.(wolfpop, pop, COD)], by.x = 'wolfID', by.y = 'wolfpop', all.x = T)

dat.wnn.lastmo.cod[,.(N=uniqueN(wolfID)), by=.(pop, COD)]
dat.meta[,.(N=uniqueN(wolfpop)), by=.(pop, COD)]

#####

dat[case_ == FALSE & wolfID %chin% dat.wnn.lastmo$wolfID, mean(propwet_end), by=.(wolfID)]
dat[case_ == FALSE & wolfID %chin% dat.wnn.lastmo$wolfID, mean(propforest_end_adj), by=.(wolfID)]
dat[case_ == FALSE & wolfID %chin% dat.wnn.lastmo$wolfID, mean(propopen_end_adj), by=.(wolfID)]

######
ggplot(dat[ua == 'used'], aes(ttd1, sl_,colour = COD)) +
 # geom_point() +
  geom_smooth() 

sl.wolf <- ggplot(dat[ua == 'used'], aes(sl_, colour = wolfID)) +
  geom_density()

ggplot(dat, aes(ta_, colour = wolfID)) +
  geom_density()

ta.wolf <- ggplot(dat[ua == 'used'], aes(ta_, colour = wolfID)) +
  geom_density()
sl.wolf|ta.wolf

sl <- ggplot(dat[ua == 'used'], aes(sl_)) +
  geom_density()
ta <- ggplot(dat[ua == 'used'], aes(ta_)) +
  geom_density()
sl|ta


# gam <- dat[ua == 'used',.(vals = MASS::fitdistr(sl_, "gamma", lower = c(0,0))[[1]],
#                          param = names(MASS::fitdistr(sl_, "gamma", lower = c(0,0))[[1]])), by = .(wolfID)]
# gam.wide <- dcast(gam, wolfID ~ param, value.var = "vals")
# 
# gam <- dat[ua == 'used',.(vals = MASS::fitdistr(sl_, "gamma", lower = c(0,0))[[1]],
#                           param = names(MASS::fitdistr(sl_, "gamma", lower = c(0,0))[[1]])), by = .(wolfID)]

gam.wolves <- merge(dat[case_==TRUE,.(wolfID, sl_)], params, by = c('wolfID'), all.x = T)


ggplot(gam.wolves[wolfID %chin% dat.wnn.lastmo$wolfID], aes(sl_, ..density..)) +
  geom_histogram(binwidth = 50) +
  geom_line(aes( y=dgamma(sl_, shape[1], scale[1])), color="blue", size = 1) +
  facet_wrap(vars(wolfID))


ggplot(dat[ua == 'used' & pop == 'RMNP' & wolfID %chin% dat.wnn.lastmo$wolfID], aes(sl_)) +
  geom_density(color='blue') + #geom_histogram(bins = 500) +
  facet_wrap(vars(wolfID))



ggplot(dat[ ua =='used' & COD=='cdv'], aes(ttd1, (distance2), color = wolfID)) +
  geom_point() + geom_smooth(method = lm, se=F)
 
ggplot(dat[ttd1<=31 & ua =='used' & COD=='cdv'], aes(ttd1, (roadDist_end), color = wolfID)) +
  geom_point() + geom_smooth(method = lm, se=F)

rd.u <-ggplot(dat[ttd1<=31 & ua =='used' & COD=='cdv'], aes(roadDist_end)) +
  geom_histogram(aes(y=..density..), binwidth = 50) + geom_density() +
  ggtitle("used - CDV") 
rd.a<-ggplot(dat[ttd1<=31 & ua =='avail' & COD=='cdv'], aes(roadDist_end)) +
  geom_histogram(aes(y=..density..), binwidth = 50) + geom_density() + 
  ggtitle("avaliable - CDV") 

rd.u.c <-ggplot(dat[ttd1<=31 & ua =='used' & COD=='control'], aes(roadDist_end)) +
  geom_histogram(aes(y=..density..), binwidth = 50) + geom_density() +
  ggtitle("used - none") 
rd.a.c <-ggplot(dat[ttd1<=31 & ua =='avail' & COD=='control'], aes(roadDist_end)) +
  geom_histogram(aes(y=..density..), binwidth = 50) + geom_density() + 
  ggtitle("avaliable - none") 
(rd.u.c|rd.u)/(rd.a.c|rd.a)


nn.u <-ggplot(dat[ttd1<=31 & ua =='used' & COD=='cdv'], aes(distance2)) +
  geom_histogram(aes(y=..density..), binwidth = 100) + geom_density() +
  ggtitle("used - CDV") + xlab('Dist to NN')
nn.a<-ggplot(dat[ttd1<=31 & ua =='avail' & COD=='cdv'], aes(distance2)) +
  geom_histogram(aes(y=..density..), binwidth = 100) + geom_density() + 
  ggtitle("avaliable - CDV") + xlab('Dist to NN')

nn.u.c <-ggplot(dat[ttd1<=31 & ua =='used' & COD=='control'], aes(distance2)) +
  geom_histogram(aes(y=..density..), binwidth = 100) + geom_density() +
  ggtitle("used - none") + xlab('Dist to NN')
nn.a.c <-ggplot(dat[ttd1<=31 & ua =='avail' & COD=='control'], aes(distance2)) +
  geom_histogram(aes(y=..density..), binwidth = 100) + geom_density() + 
  ggtitle("avaliable - none") + xlab('Dist to NN')
(nn.u.c|nn.u)/(nn.a.c|nn.a)




ggplot(dat[ttd1<=31 & ua =='used' & COD=='human'], aes(roadDist_end)) +
  geom_histogram() 
ggplot(dat[ttd1<=31 & ua =='avail' & COD=='human'], aes(roadDist_end)) +
  geom_histogram() 

ggplot(dat[ttd1<=31 & ua =='used' & COD=='cdv'], aes(land_end_adj)) +
  geom_histogram(stat = 'count') 
ggplot(dat[ttd1<=31 & ua =='avail' & COD=='cdv'], aes(land_end_adj)) +
  geom_histogram(stat = 'count') 


dat[COD=='cdv',.(last(packDistadj_end), min(packDistadj_end, na.rm = T), max(packDistadj_end, na.rm = T)), by=.(wolfID)]


#### full model #### 
#### all time ####

quantile(dat$distance2, probs = c(0, 0.05, .95, 1), na.rm = T)
quantile(dat$distance2, na.rm = T)
quantile(dat$packDist_end, probs = c(0, 0.05, .95, 1), na.rm = T)
quantile(dat$packDist_end, na.rm = T)
dat[distance2 >=50000,.(unique(wolfID), .N), by=.(COD)]
dat[,.(unique(wolfID), .N), by=.(COD)]

dat[,'nnDist_end'] <- ifelse(dat$distance2<=30000, dat$distance2, NA)
dat[,'packDist_end_5'] <- ifelse(dat$packDist_end<=50000, dat$packDist_end, NA)

#### everyone ####

everyone <- glmmTMB(case_ ~# pop + 
                      #log_sl:ToD_start +
                      log_sl:cos_ta +
                       # log_sl:land_end_adj +
                      log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                        log_sl:COD + cos_ta:COD + 
                        I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                        (1|wolf_step_id) +
                        (0 + (log_sl)|wolfID) +
                        (0 + (log_sl:cos_ta)|wolfID) +
                      (0 + (cos_ta)|wolfID) +
                        (0 + (I(log(ttd1 + 1)):log_sl)|wolfID) +
                        (0 + (I(log(ttd1 + 1)):cos_ta)|wolfID) +
                        # COD:land_end_adj + I(log(1+roadDist_end)):COD +
                        # COD:I(log(ttd1 + 1)):land_end_adj +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +
                      COD:propforest_end_adj + COD:propopen_end_adj + COD:propwet_end + I(log(1+roadDist_end)):COD +
                      COD:I(log(ttd1 + 1)):propforest_end_adj +  COD:I(log(ttd1 + 1)):propopen_end_adj +  COD:I(log(ttd1 + 1)):propwet_end +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +
                      
                        # (0 + land_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):land_end_adj)|wolfID) +
                      (0 + propforest_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):propforest_end_adj)|wolfID) +
                      (0 + propopen_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):propopen_end_adj)|wolfID) +
                      (0 + propwet_end|wolfID) + (0 + (I(log(ttd1 + 1)):propwet_end)|wolfID) +
                        (0 + I(log(1+roadDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+roadDist_end)))|wolfID) +
                        
                        I(log(1+distance2)):COD + I(log(1+packDist_end)):COD +
                        I(log(ttd1 + 1)):I(log(1+distance2)):COD + I(log(ttd1 + 1)):I(log(1+packDist_end)):COD +
                        
                        (0 + I(log(1+distance2))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+distance2)))|wolfID) +
                        (0 + I(log(1+packDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+packDist_end)))|wolfID)
                      , family=poisson(),
                      data = dat[wolfID %chin% dat.wnn.lastmo$wolfID], 
                    map = list(theta=factor(c(NA,1:17))), start = list(theta=c(log(1000),seq(0,0, length.out = 17))))

# everyone$parameters$theta[1] <- log(1e3)
# nvar_parm <- length(everyone$parameters$theta)
# everyone$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
# everyone <- glmmTMB:::fitTMB(everyone)
summary(everyone)

summary(everyone)$coef$cond[-1, "Estimate"]
popeveryone<- summary(everyone)$coef$cond[-1, 1:2]
#saveRDS(popeveryone, 'data/derived-data/popeveryone_COD.Rds')
sum.everyone<- tidy(everyone)
#saveRDS(sum.everyone, 'data/derived-data/summarypopeveryone_COD.Rds')

everyone.ran_vals <-tidy(everyone, effect= 'ran_vals')
# everyone.ran_pars <-tidy(everyone, effect= 'ran_pars')
everyone.se <-setDT(everyone.ran_vals)[group=='wolfID']


# everyone.all.indiv <- coef(everyone)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
#   pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
#   mutate(method = "ME")




#### by model ####
everyone.move <- glmmTMB(case_ ~ log_sl:ToD_start +
                           log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                              #log_sl:land_end_adj +
                              log_sl:COD + cos_ta:COD +
                              I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                              (1|wolf_step_id) +
                              (0 + (log_sl)|wolfID) +
                              (0 + (cos_ta)|wolfID) +
                              (0 + (I(log(ttd1 + 1)):log_sl)|wolfID) +
                              (0 + (I(log(ttd1 + 1)):cos_ta)|wolfID) # +
                              # COD:land_end_adj + I(log(1+roadDist_end)):COD +
                              # COD:I(log(ttd1 + 1)):land_end_adj +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +
                              # 
                              # (0 + land_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):land_end_adj)|wolfID) +
                              # (0 + I(log(1+roadDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+roadDist_end)))|wolfID)# +
                              # 
                            # I(log(1+distance2)):COD + I(log(1+packDist_end)):COD +
                            # I(log(ttd1 + 1)):I(log(1+distance2)):COD + I(log(ttd1 + 1)):I(log(1+packDist_end)):COD +
                            # 
                            # (0 + I(log(1+distance2))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+distance2)))|wolfID) +
                            # (0 + I(log(1+packDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+packDist_end)))|wolfID)
                            , family=poisson(),
                            data = dat, #[wolfID %chin% dat.wnn.lastmo$wolfID], 
                            map = list(theta=factor(c(NA,1:4))), start = list(theta=c(log(1000),seq(0,0, length.out = 4))))
summary(everyone.move)
popeveryoneMove<- summary(everyone.move)$coef$cond[-1, 1:2]
#saveRDS(popeveryoneMove, 'data/derived-data/popeveryoneMove_COD.Rds')
sum.everyoneMove<- summary(everyone.move)$coef$cond
#saveRDS(sum.everyoneMove, 'data/derived-data/summarypopeveryoneMove_COD.Rds')



everyone.habitat <- glmmTMB(case_ ~ log_sl:ToD_start +
                             log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                             # log_sl:COD + cos_ta:COD + 
                             # I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                             (1|wolf_step_id) +
                             # (0 + (log_sl)|wolfID) +
                             # (0 + (cos_ta)|wolfID) +
                             # (0 + (I(log(ttd1 + 1)):log_sl)|wolfID) +
                             # (0 + (I(log(ttd1 + 1)):cos_ta)|wolfID) +
                             COD:propforest_end_adj + COD:propopen_end_adj + COD:propwet_end + I(log(1+roadDist_end)):COD +
                             COD:I(log(ttd1 + 1)):propforest_end_adj +  COD:I(log(ttd1 + 1)):propopen_end_adj +  COD:I(log(ttd1 + 1)):propwet_end +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +

                             (0 + propforest_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):propforest_end_adj)|wolfID) +
                              (0 + propopen_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):propopen_end_adj)|wolfID) +
                              (0 + propwet_end|wolfID) + (0 + (I(log(ttd1 + 1)):propwet_end)|wolfID) +
                             (0 + I(log(1+roadDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+roadDist_end)))|wolfID)# +

                             # I(log(1+distance2)):COD + I(log(1+packDist_end)):COD +
                             # I(log(ttd1 + 1)):I(log(1+distance2)):COD + I(log(ttd1 + 1)):I(log(1+packDist_end)):COD +
                             # 
                             # (0 + I(log(1+distance2))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+distance2)))|wolfID) +
                             # (0 + I(log(1+packDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+packDist_end)))|wolfID)
                           , family=poisson(),
                           data = dat[wolfID %chin% dat.wnn.lastmo$wolfID], 
                           map = list(theta=factor(c(NA,1:8))), start = list(theta=c(log(1000),seq(0,0, length.out = 8))))
summary(everyone.habitat)
popeveryoneHab<- summary(everyone.habitat)$coef$cond[-1, 1:2]
#saveRDS(popeveryoneHab, 'data/derived-data/popeveryoneHab_COD.Rds')
sum.everyoneHab<- summary(everyone.habitat)$coef$cond
#saveRDS(sum.everyoneHab, 'data/derived-data/summarypopeveryoneHab_COD.Rds')

everyone.social <- glmmTMB(case_ ~ log_sl:ToD_start +
                             log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                      # log_sl:COD + cos_ta:COD + 
                      # I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                      (1|wolf_step_id) +
                      # (0 + (log_sl)|wolfID) +
                      # (0 + (cos_ta)|wolfID) +
                      # (0 + (I(log(ttd1 + 1)):log_sl)|wolfID) +
                      # (0 + (I(log(ttd1 + 1)):cos_ta)|wolfID) +
                      # COD:land_end_adj + I(log(1+roadDist_end)):COD +
                      # COD:I(log(ttd1 + 1)):land_end_adj +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +
                      # 
                      # (0 + land_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):land_end_adj)|wolfID) +
                      # (0 + I(log(1+roadDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+roadDist_end)))|wolfID) +
                      # 
                      I(log(1+distance2)):COD + I(log(1+packDist_end)):COD +
                      I(log(ttd1 + 1)):I(log(1+distance2)):COD + I(log(ttd1 + 1)):I(log(1+packDist_end)):COD +
                      
                      (0 + I(log(1+distance2))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+distance2)))|wolfID) +
                      (0 + I(log(1+packDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+packDist_end)))|wolfID)
                    , family=poisson(),
                    data = dat[wolfID %chin% dat.wnn.lastmo$wolfID], 
                    map = list(theta=factor(c(NA,1:4))), start = list(theta=c(log(1000),seq(0,0, length.out = 4))))
summary(everyone.social)
popeveryoneSoc<- summary(everyone.social)$coef$cond[-1, 1:2]
#saveRDS(popeveryoneSoc, 'data/derived-data/popeveryoneSoc_COD.Rds')
sum.everyoneSoc<- summary(everyone.social)$coef$cond
#saveRDS(sum.everyoneSoc, 'data/derived-data/summarypopeveryoneSoc_COD.Rds')


everyone.pack <- glmmTMB(case_ ~ log_sl:ToD_start +
                             log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                             # log_sl:COD + cos_ta:COD + 
                             # I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                             (1|wolf_step_id) +
                             # (0 + (log_sl)|wolfID) +
                             # (0 + (cos_ta)|wolfID) +
                             # (0 + (I(log(ttd1 + 1)):log_sl)|wolfID) +
                             # (0 + (I(log(ttd1 + 1)):cos_ta)|wolfID) +
                             # COD:land_end_adj + I(log(1+roadDist_end)):COD +
                             # COD:I(log(ttd1 + 1)):land_end_adj +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +
                             # 
                             # (0 + land_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):land_end_adj)|wolfID) +
                             # (0 + I(log(1+roadDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+roadDist_end)))|wolfID) +
                             # 
                            # I(log(1+distance2)):COD + 
                           I(log(1+packDist_end)):COD +
                           #  I(log(ttd1 + 1)):I(log(1+distance2)):COD + 
                           I(log(ttd1 + 1)):I(log(1+packDist_end)):COD +
                             
                            # (0 + I(log(1+distance2))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+distance2)))|wolfID) +
                             (0 + I(log(1+packDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+packDist_end)))|wolfID)
                           , family=poisson(),
                           data = dat[wolfID %chin% dat.wnn.lastmo$wolfID], 
                           map = list(theta=factor(c(NA,1:2))), start = list(theta=c(log(1000),seq(0,0, length.out = 2))))
summary(everyone.pack)
popeveryonePack<- summary(everyone.pack)$coef$cond[-1, 1:2]
#saveRDS(popeveryonePack, 'data/derived-data/popeveryonePack_COD.Rds')
sum.everyonePack<- summary(everyone.pack)$coef$cond
#saveRDS(sum.everyonePack, 'data/derived-data/summarypopeveryonePack_COD.Rds')


AICtab(everyone.move, everyone.habitat, everyone, everyone.social, everyone.pack)


#### GATHERING RESULTS ####

cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")
gcolors <- c("deepskyblue", "purple", "dark green")

#### everyone BETA graphs ####
everyone.indiv<- merge(everyone.se[,.(wolfID= level, term, estimate, se=std.error)], dat.meta[,.(wolfpop, COD)], by.x ='wolfID', by.y= 'wolfpop', all.x=T)
everyone.indiv$COD <- factor(everyone.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
everyone.indiv$COD[is.na(everyone.indiv$COD)] <- "control"

unique(everyone.indiv$term)
everyone.all.betas.names <- c("log_sl", "cos_ta",  'log_sl:cos_ta',
                             "propforest_end_adj", "propopen_end_adj", "propwet_end", "I(log(1 + roadDist_end))",
                             "I(log(1 + distance2))", "I(log(1 + packDist_end))",
                             "I(log(ttd1 + 1)):log_sl", "I(log(ttd1 + 1)):cos_ta", 
                             "I(log(ttd1 + 1)):propforest_end_adj", "I(log(ttd1 + 1)):propopen_end_adj", "I(log(ttd1 + 1)):propwet_end", "I(log(ttd1 + 1)):I(log(1 + roadDist_end))",
                             "I(log(ttd1 + 1)):I(log(1 + distance2))", "I(log(ttd1 + 1)):I(log(1 + packDist_end))")
everyone.indiv$term <- factor(everyone.indiv$term, levels = everyone.all.betas.names, labels = c("log_sl", "cos_ta", 'sl_ta','forest', "open", "wet", "roadDist",
                                                                                                                     "nnDist", "boundaryDist",
                                                                                                                     "log_sl-ttd", "cos_ta-ttd", "forest-ttd", "open-ttd", "wet-ttd", "roadDist-ttd",
                                                                                                                     "nnDist-ttd", "boundaryDist-ttd"))


#saveRDS(everyone.indiv, 'data/derived-data/everyone_betas_COD.Rds')
everyone.indiv<-readRDS('data/derived-data/everyone_betas_COD.Rds')
everyone.ttd <- everyone.indiv[term %like% "ttd", ]

everyone.main <- everyone.indiv[!(term %like% "ttd"), ]




#### ALL BETAS ####

ttd.vars <- ggplot(everyone.ttd, aes(term, (estimate), fill = COD)) +
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


# ggplot(everyone.ttd[term== 'wet-ttd' | term== 'boundaryDist-ttd'], aes(term, (estimate), fill = COD)) +
#   geom_boxplot(aes(fill = COD),# notch = TRUE, notchwidth = 0.7,
#                outlier.color = NA, lwd = 0.6,
#                alpha = 0.25) +
#   geom_jitter(aes(color = COD),
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
#   ylab('beta') +
#   # ggtitle("b) case") +
#   scale_fill_manual(values = cbPalette) +
#   scale_color_manual(values = cbPalette) #+ ylim(-2,2)

main.vars <- ggplot(everyone.main, aes(term, (estimate), fill = COD)) +
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

main.vars/ttd.vars

#### MOVE BETAS ####
unique(sum.everyone$term)
everyone.pop.betas.move <- c("log_sl:CODhuman", "log_sl:CODCDV", "CODcontrol:cos_ta", "CODhuman:cos_ta", "CODCDV:cos_ta",
                             "log_sl:CODcontrol:I(log(ttd1 + 1))",  "log_sl:CODhuman:I(log(ttd1 + 1))", "log_sl:CODCDV:I(log(ttd1 + 1))",
                             "CODcontrol:cos_ta:I(log(ttd1 + 1))", "CODhuman:cos_ta:I(log(ttd1 + 1))", "CODCDV:cos_ta:I(log(ttd1 + 1))")

everyone.pop.move <- setDT(sum.everyone)[term %in% everyone.pop.betas.move]
everyone.pop.move[,'COD'] <- ifelse(everyone.pop.move$term %like% 'CDV', 'CDV', ifelse(everyone.pop.move$term %like% 'human', 'human', 'control'))
everyone.pop.move$term <- factor(everyone.pop.move$term, levels = everyone.pop.betas.move, labels = c("log_sl", "log_sl", "cos_ta", "cos_ta", "cos_ta",
                                                                                                      "log_sl-ttd",  "log_sl-ttd", "log_sl-ttd",
                                                                                                      "cos_ta-ttd", "cos_ta-ttd", "cos_ta-ttd"))
everyone.pop.move$COD <- factor(everyone.pop.move$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))

everyone.move.ttd <- everyone.pop.move[term %like% "ttd", ]

everyone.move.main <- everyone.pop.move[!(term %like% "ttd"), ]

pd <- position_dodge(0.5) # move them .05 to the left and right

move.ttd <-ggplot(everyone.move.ttd, aes(x=term, y=estimate, colour=COD)) + 
  geom_errorbar(aes(ymin=estimate - 1.96*std.error, ymax=estimate + 1.96*std.error), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  scale_x_discrete(breaks=c("log_sl-ttd", "cos_ta-ttd"),
                   labels=c("lnSL x TTD", "cosTA x TTD")) +
  xlab('') +
  ylab('beta') +
  #ggtitle("Movement") +
  scale_fill_manual(values = gcolors) +
  scale_color_manual(values = gcolors) #+ ylim(-.35,.3)

move.ttd 
move.ttd + geom_point(data = everyone.ttd[term=='log_sl-ttd'|term=='cos_ta-ttd'], aes(term, (estimate), fill = COD), position = pd)


ttd.move <- ggplot(everyone.ttd[term=='log_sl-ttd'|term=='cos_ta-ttd'], aes(term, (estimate), fill = COD)) +
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
  ggtitle("Movement") +
  scale_fill_manual(values = gcolors) +
  scale_color_manual(values = gcolors) #+ ylim(-.35,.3)



main.move <- ggplot(everyone.main[term=='log_sl'|term=='cos_ta'], aes(term, (estimate), fill = COD)) +
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
  ggtitle("Movement") +
  scale_fill_manual(values = gcolors) +
  scale_color_manual(values = gcolors)# + ylim(-.35,.3)

ttd.move

#### HABITAT BETAS ####
ttd.hab <- ggplot(everyone.ttd[term=='forest-ttd'|term=='open-ttd'|term=='wet-ttd'], aes(term, (estimate), fill = COD)) +
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
  ggtitle("Habitat") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-8,12)


main.hab <- ggplot(everyone.main[term=='forest'|term=='open'|term=='wet'], aes(term, (estimate), fill = COD)) +
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
  ggtitle("Habitat") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-8,12)


main.hab/ttd.hab

#### DISTANCE BETAS ####
ttd.dist <- ggplot(everyone.ttd[term=='roadDist-ttd'|term=='nnDist-ttd'|term=='boundaryDist-ttd'], aes(term, (estimate), fill = COD)) +
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
  ggtitle("Distance to Rd and Social") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-1,2.5)

main.dist <- ggplot(everyone.main[term=='roadDist'|term=='nnDist'|term=='boundaryDist'], aes(term, (estimate), fill = COD)) +
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
  ggtitle("Distance to Rd and Social") +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) #+ ylim(-1,2.5)

main.dist/ttd.dist



#### RSS ####
#### h2 RSS ####
### CDV ###
CDV.1day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = median(distance2, na.rm = T),
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('CDV', levels = levels(COD)),
                                                           ttd1 = 1,
                                                           wolf_step_id = NA, wolfID= NA)]
p.CDV.1day.2 <- predict(everyone, newdata = CDV.1day.2, type='link', re.form = NA)
p.CDV.1day.2.habitat <- predict(everyone.habitat, newdata = CDV.1day.2, type='link', re.form = NA)
p.CDV.1day.2.social <- predict(everyone.social, newdata = CDV.1day.2, type='link', re.form = NA)
p.CDV.1day.2.pack <- predict(everyone.pack, newdata = CDV.1day.2, type='link', re.form = NA)



CDV.60day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                        log_sl = mean(log_sl),
                                                        cos_ta = mean(cos_ta),
                                                        # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                        propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                        propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                        propwet_end= mean(propwet_end, na.rm = T),      
                                                        roadDist_end = median(roadDist_end, na.rm = T),
                                                        distance2 = median(distance2, na.rm = T),
                                                        nnDist_end = median(nnDist_end, na.rm = T),
                                                        packDist_end = median(packDist_end, na.rm = T),
                                                        COD = factor('CDV', levels = levels(COD)),
                                                        ttd1 = 60,
                                                        wolf_step_id = NA, wolfID= NA)]
p.CDV.60day.2 <- predict(everyone, newdata = CDV.60day.2, type='link', re.form = NA)
p.CDV.60day.2.habitat <- predict(everyone.habitat, newdata = CDV.60day.2, type='link', re.form = NA)
p.CDV.60day.2.social <- predict(everyone.social, newdata = CDV.60day.2, type='link', re.form = NA)
p.CDV.60day.2.pack <- predict(everyone.pack, newdata = CDV.60day.2, type='link', re.form = NA)


### human ###
human.1day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                        log_sl = mean(log_sl),
                                                        cos_ta = mean(cos_ta),
                                                        # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                        propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                        propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                        propwet_end= mean(propwet_end, na.rm = T),      
                                                        roadDist_end = median(roadDist_end, na.rm = T),
                                                        distance2 = median(distance2, na.rm = T),
                                                        nnDist_end = median(nnDist_end, na.rm = T),
                                                        packDist_end = median(packDist_end, na.rm = T),
                                                        COD = factor('human', levels = levels(COD)),
                                                        ttd1 = 1,
                                                        wolf_step_id = NA, wolfID= NA)]
p.human.1day.2 <- predict(everyone, newdata = human.1day.2, type='link', re.form = NA)
p.human.1day.2.habitat <- predict(everyone.habitat, newdata = human.1day.2, type='link', re.form = NA)
p.human.1day.2.social <- predict(everyone.social, newdata = human.1day.2, type='link', re.form = NA)
p.human.1day.2.pack <- predict(everyone.pack, newdata = human.1day.2, type='link', re.form = NA)



human.60day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= mean(propwet_end, na.rm = T),      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('human', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.human.60day.2 <- predict(everyone, newdata = human.60day.2, type='link', re.form = NA)
p.human.60day.2.habitat <- predict(everyone.habitat, newdata = human.60day.2, type='link', re.form = NA)
p.human.60day.2.social <- predict(everyone.social, newdata = human.60day.2, type='link', re.form = NA)
p.human.60day.2.pack <- predict(everyone.pack, newdata = human.60day.2, type='link', re.form = NA)


### control ###
control.1day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                        log_sl = mean(log_sl),
                                                        cos_ta = mean(cos_ta),
                                                        # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                        propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                        propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                        propwet_end= mean(propwet_end, na.rm = T),      
                                                        roadDist_end = median(roadDist_end, na.rm = T),
                                                        distance2 = median(distance2, na.rm = T),
                                                        nnDist_end = median(nnDist_end, na.rm = T),
                                                        packDist_end = median(packDist_end, na.rm = T),
                                                        COD = factor('control', levels = levels(COD)),
                                                        ttd1 = 1,
                                                        wolf_step_id = NA, wolfID= NA)]
p.control.1day.2 <- predict(everyone, newdata = control.1day.2, type='link', re.form = NA)
p.control.1day.2.habitat <- predict(everyone.habitat, newdata = control.1day.2, type='link', re.form = NA)
p.control.1day.2.social <- predict(everyone.social, newdata = control.1day.2, type='link', re.form = NA)
p.control.1day.2.pack <- predict(everyone.pack, newdata = control.1day.2, type='link', re.form = NA)


control.60day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= mean(propwet_end, na.rm = T),      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('control', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.control.60day.2 <- predict(everyone, newdata = control.60day.2, type='link', re.form = NA)
p.control.60day.2.habitat <- predict(everyone.habitat, newdata = control.60day.2, type='link', re.form = NA)
p.control.60day.2.social <- predict(everyone.social, newdata = control.60day.2, type='link', re.form = NA)
p.control.60day.2.pack <- predict(everyone.pack, newdata = control.60day.2, type='link', re.form = NA)


### INDIVs ###
p.h2.indiv <- function(ids, DT, mod, death, t2death){
  lapply(ids, function(i) {
    unique(
      DT[#wolfID == i,
        ,.(h2 = predict(
          mod,
          newdata = .SD[, .(
            ToD_start = factor('day', levels = levels(ToD_start)),
            log_sl = mean(log_sl),
            cos_ta = mean(cos_ta),
            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
            propforest_end_adj = mean(propforest_end_adj, na.rm = T),
            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
            propwet_end= mean(propwet_end, na.rm = T), 
            roadDist_end = median(roadDist_end, na.rm = T),
            distance2 = median(distance2, na.rm = T),
            nnDist_end = median(nnDist_end, na.rm = T),
            packDist_end = median(packDist_end, na.rm = T),
            COD = factor(death, levels = levels(COD)),
            ttd1 = t2death,
            wolf_step_id = NA,
            wolfID = i
          )],
          type = "link",
          re.form = NULL
        ), wolfID = i)]
    )
  })}

### CDV ###
CDV.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='CDV', wolfID])

p.CDV.1day.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)
p.CDV.1day.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 1)
p.CDV.1day.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 1)
p.CDV.1day.2.indiv.pack <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', t2death = 1)

p.CDV.60day.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)
p.CDV.60day.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 60)
p.CDV.60day.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 60)
p.CDV.60day.2.indiv.pack <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', t2death = 60)


### human ###
human.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='human', wolfID])

p.human.1day.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)
p.human.1day.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 1)
p.human.1day.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 1)
p.human.1day.2.indiv.pack <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', t2death = 1)

p.human.60day.2.indiv <- p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)
p.human.60day.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 60)
p.human.60day.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 60)
p.human.60day.2.indiv.pack <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', t2death = 60)


### control ###
control.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='control', wolfID])

p.control.1day.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)
p.control.1day.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 1)
p.control.1day.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 1)
p.control.1day.2.indiv.pack <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', t2death = 1)

p.control.60day.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)
p.control.60day.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 60)
p.control.60day.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 60)
p.control.60day.2.indiv.pack <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', t2death = 60)




#### FOREST ####
#### CDV ####
forest.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),
                                                             roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = median(distance2, na.rm = T),
                                                             nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA,
                                                             wolfID = NA)]

p.forest.CDV.1day.1 <- predict(everyone, newdata = forest.CDV.1day.1, type='link', re.form = NA)
p.forest.CDV.1day.1.habitat <- predict(everyone.habitat, newdata = forest.CDV.1day.1, type='link', re.form = NA)


logRSS.forest.CDV.1day<- p.forest.CDV.1day.1 - p.CDV.1day.2
logRSS.forest.CDV.1day.habitat<- p.forest.CDV.1day.1.habitat - p.CDV.1day.2.habitat

forest.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                              log_sl = mean(log_sl),
                                                              cos_ta = mean(cos_ta),
                                                              #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                              propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                              propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                              propwet_end= mean(propwet_end, na.rm = T),
                                                              roadDist_end = median(roadDist_end, na.rm = T),
                                                              distance2 = median(distance2, na.rm = T),
                                                              nnDist_end = median(nnDist_end, na.rm = T),
                                                              packDist_end = median(packDist_end, na.rm = T),
                                                              COD = factor('CDV', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]
p.forest.CDV.60day.1 <- predict(everyone, newdata = forest.CDV.60day.1, type='link', re.form =NA)
p.forest.CDV.60day.1.habitat <- predict(everyone.habitat, newdata = forest.CDV.60day.1, type='link', re.form =NA)


logRSS.forest.CDV.60day<- p.forest.CDV.60day.1 - p.CDV.60day.2
logRSS.forest.CDV.60day.habitat<- p.forest.CDV.60day.1.habitat - p.CDV.60day.2.habitat

logRSS.forest <- data.frame(COD= 'CDV', ttd = '1 day', var = 'forest', rss = logRSS.forest.CDV.1day, x = seq(from = 0, to = 1, length.out = 100))
logRSS.forest <- rbind(logRSS.forest, data.frame(COD= 'CDV', ttd = '60 days', var = 'forest', rss = logRSS.forest.CDV.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.forest.habitat <- data.frame(COD= 'CDV', ttd = '1 day', var = 'forest', rss = logRSS.forest.CDV.1day.habitat, x = seq(from = 0, to = 1, length.out = 100))
logRSS.forest.habitat <- rbind(logRSS.forest.habitat, data.frame(COD= 'CDV', ttd = '60 days', var = 'forest', rss = logRSS.forest.CDV.60day.habitat, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####

p.forest.h1.indiv <- function(ids, DT, mod, death, t2death){
    lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
          propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
          propwet_end= mean(propwet_end, na.rm = T),
          roadDist_end = median(roadDist_end, na.rm = T),
          distance2 = median(distance2, na.rm = T),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = median(packDist_end, na.rm = T),
          COD = factor(death, levels = levels(COD)),
          ttd1 = t2death,
          wolf_step_id = NA,
          wolfID = i
        )],
        type = "link",
        re.form = NULL
      ), wolfID = i)]
    # )
  })
}

p.forest.CDV.1day.1.indiv <- p.forest.h1.indiv(ids = CDV.wolfID, DT=dat, mod=everyone, death = 'CDV', t2death = 1)
  
h1.forest.CDV.1day.indiv <- data.table(rbindlist(p.forest.CDV.1day.1.indiv),
                                     ttd = '1 day', COD = 'CDV', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                     ttd = '1 day', COD = 'CDV', var = 'forest')



logRSS.forest.CDV.1day.indiv <- merge(h1.forest.CDV.1day.indiv, h2.forest.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.CDV.1day.indiv[,'rss'] <- logRSS.forest.CDV.1day.indiv$h1 - logRSS.forest.CDV.1day.indiv$h2

#### habitat model

p.forest.CDV.1day.1.indiv.habitat <-p.forest.h1.indiv(ids = CDV.wolfID, DT=dat, mod=everyone.habitat, death = 'CDV', t2death = 1)


h1.forest.CDV.1day.indiv.habitat <- data.table(rbindlist(p.forest.CDV.1day.1.indiv.habitat),
                                             ttd = '1 day', COD = 'CDV', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.CDV.1day.indiv.habitat <- data.table(rbindlist(p.CDV.1day.2.indiv.habitat),
                                             ttd = '1 day', COD = 'CDV', var = 'forest')


logRSS.forest.CDV.1day.indiv.habitat <- merge(h1.forest.CDV.1day.indiv.habitat, h2.forest.CDV.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.CDV.1day.indiv.habitat[,'rss'] <- logRSS.forest.CDV.1day.indiv.habitat$h1 - logRSS.forest.CDV.1day.indiv.habitat$h2


### 60 days

p.forest.CDV.60day.1.indiv <- p.forest.h1.indiv(ids = CDV.wolfID, DT=dat, mod=everyone, death = 'CDV', t2death = 60)

h1.forest.CDV.60day.indiv <- data.table(rbindlist(p.forest.CDV.60day.1.indiv),
                                      ttd = '60 days', COD = 'CDV', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.CDV.60day.indiv <- data.table(rbindlist(p.CDV.60day.2.indiv),
                                      ttd = '60 days', COD = 'CDV', var = 'forest')



logRSS.forest.CDV.60day.indiv <- merge(h1.forest.CDV.60day.indiv, h2.forest.CDV.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.CDV.60day.indiv[,'rss'] <- logRSS.forest.CDV.60day.indiv$h1 - logRSS.forest.CDV.60day.indiv$h2

#### habitat model

p.forest.CDV.60day.1.indiv.habitat <- p.forest.h1.indiv(ids = CDV.wolfID, DT=dat, mod=everyone.habitat, death = 'CDV', t2death = 60)


h1.forest.CDV.60day.indiv.habitat <- data.table(rbindlist(p.forest.CDV.60day.1.indiv.habitat),
                                              ttd = '60 days', COD = 'CDV', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))



h2.forest.CDV.60day.indiv.habitat <- data.table(rbindlist(p.CDV.60day.2.indiv.habitat),
                                              ttd = '60 days', COD = 'CDV', var = 'forest')


logRSS.forest.CDV.60day.indiv.habitat <- merge(h1.forest.CDV.60day.indiv.habitat, h2.forest.CDV.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.CDV.60day.indiv.habitat[,'rss'] <- logRSS.forest.CDV.60day.indiv.habitat$h1 - logRSS.forest.CDV.60day.indiv.habitat$h2



#### human ####

forest.human.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('forest', levels = levels(land_end_adj)),
    propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
    propopen_end_adj = mean(propopen_end_adj, na.rm = T),
    propwet_end= mean(propwet_end, na.rm = T),
    roadDist_end = median(roadDist_end, na.rm = T),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('human', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.forest.human.1day.1 <- predict(everyone, newdata = forest.human.1day.1, type='link', re.form = NA)
p.forest.human.1day.1.habitat <- predict(everyone.habitat, newdata = forest.human.1day.1, type='link', re.form = NA)


logRSS.forest.human.1day<- p.forest.human.1day.1 - p.human.1day.2
logRSS.forest.human.1day.habitat<- p.forest.human.1day.1.habitat - p.human.1day.2.habitat


forest.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                                propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                                propwet_end= mean(propwet_end, na.rm = T),
                                                                roadDist_end = median(roadDist_end, na.rm = T),
                                                                distance2 = median(distance2, na.rm = T),
                                                                nnDist_end = median(nnDist_end, na.rm = T),
                                                                packDist_end = median(packDist_end, na.rm = T),
                                                                COD = factor('human', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]

p.forest.human.60day.1 <- predict(everyone, newdata = forest.human.60day.1, type='link', re.form = NA)
p.forest.human.60day.1.habitat <- predict(everyone.habitat, newdata = forest.human.60day.1, type='link', re.form = NA)


logRSS.forest.human.60day<- p.forest.human.60day.1 - p.human.60day.2
logRSS.forest.human.60day.habitat<- p.forest.human.60day.1.habitat - p.human.60day.2.habitat

logRSS.forest <- rbind(logRSS.forest, data.frame(COD= 'human', ttd = '1 day', var = 'forest', rss = logRSS.forest.human.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.forest <- rbind(logRSS.forest, data.frame(COD= 'human', ttd = '60 days', var = 'forest', rss = logRSS.forest.human.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.forest.habitat <- rbind(logRSS.forest.habitat, data.frame(COD= 'human', ttd = '1 day', var = 'forest', rss = logRSS.forest.human.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.forest.habitat <- rbind(logRSS.forest.habitat, data.frame(COD= 'human', ttd = '60 days', var = 'forest', rss = logRSS.forest.human.60day, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####
p.forest.human.1day.1.indiv <- p.forest.h1.indiv(ids = human.wolfID, DT=dat, mod=everyone, death = 'human', t2death = 1)


h1.forest.human.1day.indiv <- data.table(rbindlist(p.forest.human.1day.1.indiv),
                                       ttd = '1 day', COD = 'human', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.human.1day.indiv <- data.table(rbindlist(p.human.1day.2.indiv),
                                       ttd = '1 day', COD = 'human', var = 'forest')



logRSS.forest.human.1day.indiv <- merge(h1.forest.human.1day.indiv, h2.forest.human.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.human.1day.indiv[,'rss'] <- logRSS.forest.human.1day.indiv$h1 - logRSS.forest.human.1day.indiv$h2

#### habitat model

p.forest.human.1day.1.indiv.habitat <- p.forest.h1.indiv(ids = human.wolfID, DT=dat, mod=everyone.habitat, death = 'human', t2death = 1)

h1.forest.human.1day.indiv.habitat <- data.table(rbindlist(p.forest.human.1day.1.indiv.habitat),
                                               ttd = '1 day', COD = 'human', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.human.1day.indiv.habitat <- data.table(rbindlist(p.human.1day.2.indiv.habitat),
                                               ttd = '1 day', COD = 'human', var = 'forest')


logRSS.forest.human.1day.indiv.habitat <- merge(h1.forest.human.1day.indiv.habitat, h2.forest.human.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.human.1day.indiv.habitat[,'rss'] <- logRSS.forest.human.1day.indiv.habitat$h1 - logRSS.forest.human.1day.indiv.habitat$h2


### 60 days

p.forest.human.60day.1.indiv <- p.forest.h1.indiv(ids = human.wolfID, DT=dat, mod=everyone, death = 'human', t2death = 60)


h1.forest.human.60day.indiv <- data.table(rbindlist(p.forest.human.60day.1.indiv),
                                        ttd = '60 days', COD = 'human', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.human.60day.indiv <- data.table(rbindlist(p.human.60day.2.indiv),
                                        ttd = '60 days', COD = 'human', var = 'forest')



logRSS.forest.human.60day.indiv <- merge(h1.forest.human.60day.indiv, h2.forest.human.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.human.60day.indiv[,'rss'] <- logRSS.forest.human.60day.indiv$h1 - logRSS.forest.human.60day.indiv$h2

#### habitat model

p.forest.human.60day.1.indiv.habitat <- p.forest.h1.indiv(ids = human.wolfID, DT=dat, mod=everyone.habitat, death = 'human', t2death = 60)


h1.forest.human.60day.indiv.habitat <- data.table(rbindlist(p.forest.human.60day.1.indiv.habitat),
                                                ttd = '60 days', COD = 'human', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))



h2.forest.human.60day.indiv.habitat <- data.table(rbindlist(p.human.60day.2.indiv.habitat),
                                                ttd = '60 days', COD = 'human', var = 'forest')


logRSS.forest.human.60day.indiv.habitat <- merge(h1.forest.human.60day.indiv.habitat, h2.forest.human.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.human.60day.indiv.habitat[,'rss'] <- logRSS.forest.human.60day.indiv.habitat$h1 - logRSS.forest.human.60day.indiv.habitat$h2




#### control ####

forest.control.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('forest', levels = levels(land_end_adj)),
    propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
    propopen_end_adj = mean(propopen_end_adj, na.rm = T),
    propwet_end= mean(propwet_end, na.rm = T),
    roadDist_end = median(roadDist_end, na.rm = T),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('control', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.forest.control.1day.1 <- predict(everyone, newdata = forest.control.1day.1, type='link', re.form = NA)
p.forest.control.1day.1.habitat <- predict(everyone.habitat, newdata = forest.control.1day.1, type='link', re.form = NA)


logRSS.forest.control.1day<- p.forest.control.1day.1 - p.control.1day.2
logRSS.forest.control.1day.habitat<- p.forest.control.1day.1.habitat - p.control.1day.2.habitat


forest.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                  log_sl = mean(log_sl),
                                                                  cos_ta = mean(cos_ta),
                                                                  #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                  propforest_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                                  propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                                  propwet_end= mean(propwet_end, na.rm = T),
                                                                  roadDist_end = median(roadDist_end, na.rm = T),
                                                                  distance2 = median(distance2, na.rm = T),
                                                                  nnDist_end = median(nnDist_end, na.rm = T),
                                                                  packDist_end = median(packDist_end, na.rm = T),
                                                                  COD = factor('control', levels = levels(COD)),
                                                                  ttd1 = 60,
                                                                  wolf_step_id = NA, wolfID = NA)]

p.forest.control.60day.1 <- predict(everyone, newdata = forest.control.60day.1, type='link', re.form = NA)
p.forest.control.60day.1.habitat <- predict(everyone.habitat, newdata = forest.control.60day.1, type='link', re.form = NA)


logRSS.forest.control.60day<- p.forest.control.60day.1 - p.control.60day.2
logRSS.forest.control.60day.habitat<- p.forest.control.60day.1.habitat - p.control.60day.2.habitat

logRSS.forest <- rbind(logRSS.forest, data.frame(COD= 'control', ttd = '1 day', var = 'forest', rss = logRSS.forest.control.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.forest <- rbind(logRSS.forest, data.frame(COD= 'control', ttd = '60 days', var = 'forest', rss = logRSS.forest.control.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.forest.habitat <- rbind(logRSS.forest.habitat, data.frame(COD= 'control', ttd = '1 day', var = 'forest', rss = logRSS.forest.control.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.forest.habitat <- rbind(logRSS.forest.habitat, data.frame(COD= 'control', ttd = '60 days', var = 'forest', rss = logRSS.forest.control.60day, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####
p.forest.control.1day.1.indiv <- p.forest.h1.indiv(ids = control.wolfID, DT=dat, mod=everyone, death = 'control', t2death = 1)

h1.forest.control.1day.indiv <- data.table(rbindlist(p.forest.control.1day.1.indiv),
                                         ttd = '1 day', COD = 'control', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.control.1day.indiv <- data.table(rbindlist(p.control.1day.2.indiv),
                                         ttd = '1 day', COD = 'control', var = 'forest')



logRSS.forest.control.1day.indiv <- merge(h1.forest.control.1day.indiv, h2.forest.control.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.control.1day.indiv[,'rss'] <- logRSS.forest.control.1day.indiv$h1 - logRSS.forest.control.1day.indiv$h2

#### habitat model

p.forest.control.1day.1.indiv.habitat <- p.forest.h1.indiv(ids = control.wolfID, DT=dat, mod=everyone.habitat, death = 'control', t2death = 1)

h1.forest.control.1day.indiv.habitat <- data.table(rbindlist(p.forest.control.1day.1.indiv.habitat),
                                                 ttd = '1 day', COD = 'control', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.control.1day.indiv.habitat <- data.table(rbindlist(p.control.1day.2.indiv.habitat),
                                                 ttd = '1 day', COD = 'control', var = 'forest')


logRSS.forest.control.1day.indiv.habitat <- merge(h1.forest.control.1day.indiv.habitat, h2.forest.control.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.control.1day.indiv.habitat[,'rss'] <- logRSS.forest.control.1day.indiv.habitat$h1 - logRSS.forest.control.1day.indiv.habitat$h2


### 60 days

p.forest.control.60day.1.indiv <-p.forest.h1.indiv(ids = control.wolfID, DT=dat, mod=everyone, death = 'control', t2death = 60)

h1.forest.control.60day.indiv <- data.table(rbindlist(p.forest.control.60day.1.indiv),
                                          ttd = '60 days', COD = 'control', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))


h2.forest.control.60day.indiv <- data.table(rbindlist(p.control.60day.2.indiv),
                                          ttd = '60 days', COD = 'control', var = 'forest')



logRSS.forest.control.60day.indiv <- merge(h1.forest.control.60day.indiv, h2.forest.control.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.control.60day.indiv[,'rss'] <- logRSS.forest.control.60day.indiv$h1 - logRSS.forest.control.60day.indiv$h2

#### habitat model

p.forest.control.60day.1.indiv.habitat <- p.forest.h1.indiv(ids = control.wolfID, DT=dat, mod=everyone.habitat, death = 'control', t2death = 60)


h1.forest.control.60day.indiv.habitat <- data.table(rbindlist(p.forest.control.60day.1.indiv.habitat),
                                                  ttd = '60 days', COD = 'control', var = 'forest', x = seq(from = 0, to = 1, length.out = 100))



h2.forest.control.60day.indiv.habitat <- data.table(rbindlist(p.control.60day.2.indiv.habitat),
                                                  ttd = '60 days', COD = 'control', var = 'forest')


logRSS.forest.control.60day.indiv.habitat <- merge(h1.forest.control.60day.indiv.habitat, h2.forest.control.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.forest.control.60day.indiv.habitat[,'rss'] <- logRSS.forest.control.60day.indiv.habitat$h1 - logRSS.forest.control.60day.indiv.habitat$h2





logRSS.forest.indiv <- rbind(logRSS.forest.control.1day.indiv, logRSS.forest.control.60day.indiv, 
                           logRSS.forest.human.1day.indiv, logRSS.forest.human.60day.indiv, 
                           logRSS.forest.CDV.1day.indiv, logRSS.forest.CDV.60day.indiv)

logRSS.forest.indiv.habitat <- rbind(logRSS.forest.control.1day.indiv.habitat, logRSS.forest.control.60day.indiv.habitat, 
                                   logRSS.forest.human.1day.indiv.habitat, logRSS.forest.human.60day.indiv.habitat, 
                                   logRSS.forest.CDV.1day.indiv.habitat, logRSS.forest.CDV.60day.indiv.habitat)





#### OPEN ####
#### CDV ####
open.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                               log_sl = mean(log_sl),
                                                               cos_ta = mean(cos_ta),
                                                               #land_end_adj = factor('open', levels = levels(land_end_adj)),
                                                               propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                               propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                               propwet_end= mean(propwet_end, na.rm = T),
                                                               roadDist_end = median(roadDist_end, na.rm = T),
                                                               distance2 = median(distance2, na.rm = T),
                                                               nnDist_end = median(nnDist_end, na.rm = T),
                                                               packDist_end = median(packDist_end, na.rm = T),
                                                               COD = factor('CDV', levels = levels(COD)),
                                                               ttd1 = 1,
                                                               wolf_step_id = NA,
                                                               wolfID = NA)]

p.open.CDV.1day.1 <- predict(everyone, newdata = open.CDV.1day.1, type='link', re.form = NA)
p.open.CDV.1day.1.habitat <- predict(everyone.habitat, newdata = open.CDV.1day.1, type='link', re.form = NA)


logRSS.open.CDV.1day<- p.open.CDV.1day.1 - p.CDV.1day.2
logRSS.open.CDV.1day.habitat<- p.open.CDV.1day.1.habitat - p.CDV.1day.2.habitat

open.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                #land_end_adj = factor('open', levels = levels(land_end_adj)),
                                                              propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                              propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                                propwet_end= mean(propwet_end, na.rm = T),
                                                                roadDist_end = median(roadDist_end, na.rm = T),
                                                                distance2 = median(distance2, na.rm = T),
                                                                nnDist_end = median(nnDist_end, na.rm = T),
                                                                packDist_end = median(packDist_end, na.rm = T),
                                                                COD = factor('CDV', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]
p.open.CDV.60day.1 <- predict(everyone, newdata = open.CDV.60day.1, type='link', re.form =NA)
p.open.CDV.60day.1.habitat <- predict(everyone.habitat, newdata = open.CDV.60day.1, type='link', re.form =NA)


logRSS.open.CDV.60day<- p.open.CDV.60day.1 - p.CDV.60day.2
logRSS.open.CDV.60day.habitat<- p.open.CDV.60day.1.habitat - p.CDV.60day.2.habitat

logRSS.open <- data.frame(COD= 'CDV', ttd = '1 day', var = 'open', rss = logRSS.open.CDV.1day, x = seq(from = 0, to = 1, length.out = 100))
logRSS.open <- rbind(logRSS.open, data.frame(COD= 'CDV', ttd = '60 days', var = 'open', rss = logRSS.open.CDV.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.open.habitat <- data.frame(COD= 'CDV', ttd = '1 day', var = 'open', rss = logRSS.open.CDV.1day.habitat, x = seq(from = 0, to = 1, length.out = 100))
logRSS.open.habitat <- rbind(logRSS.open.habitat, data.frame(COD= 'CDV', ttd = '60 days', var = 'open', rss = logRSS.open.CDV.60day.habitat, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####

p.open.h1.indiv <- function(ids, DT, mod, death, t2death){
    lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          # land_end_adj = factor('open', levels = levels(land_end_adj)),
          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
          propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
          propwet_end= mean(propwet_end, na.rm = T),
          roadDist_end = median(roadDist_end, na.rm = T),
          distance2 = median(distance2, na.rm = T),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = median(packDist_end, na.rm = T),
          COD = factor(death, levels = levels(COD)),
          ttd1 = t2death,
          wolf_step_id = NA,
          wolfID = i
        )],
        type = "link",
        re.form = NULL
      ), wolfID = i)]
    # )
  })
}

p.open.CDV.1day.1.indiv <- p.open.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)

h1.open.CDV.1day.indiv <- data.table(rbindlist(p.open.CDV.1day.1.indiv),
                                       ttd = '1 day', COD = 'CDV', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                       ttd = '1 day', COD = 'CDV', var = 'open')



logRSS.open.CDV.1day.indiv <- merge(h1.open.CDV.1day.indiv, h2.open.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.CDV.1day.indiv[,'rss'] <- logRSS.open.CDV.1day.indiv$h1 - logRSS.open.CDV.1day.indiv$h2

#### habitat model

p.open.CDV.1day.1.indiv.habitat <- p.open.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 1)


h1.open.CDV.1day.indiv.habitat <- data.table(rbindlist(p.open.CDV.1day.1.indiv.habitat),
                                               ttd = '1 day', COD = 'CDV', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.CDV.1day.indiv.habitat <- data.table(rbindlist(p.CDV.1day.2.indiv.habitat),
                                               ttd = '1 day', COD = 'CDV', var = 'open')


logRSS.open.CDV.1day.indiv.habitat <- merge(h1.open.CDV.1day.indiv.habitat, h2.open.CDV.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.CDV.1day.indiv.habitat[,'rss'] <- logRSS.open.CDV.1day.indiv.habitat$h1 - logRSS.open.CDV.1day.indiv.habitat$h2


### 60 days

p.open.CDV.60day.1.indiv <- p.open.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)

h1.open.CDV.60day.indiv <- data.table(rbindlist(p.open.CDV.60day.1.indiv),
                                        ttd = '60 days', COD = 'CDV', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.CDV.60day.indiv <- data.table(rbindlist(p.CDV.60day.2.indiv),
                                        ttd = '60 days', COD = 'CDV', var = 'open')



logRSS.open.CDV.60day.indiv <- merge(h1.open.CDV.60day.indiv, h2.open.CDV.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.CDV.60day.indiv[,'rss'] <- logRSS.open.CDV.60day.indiv$h1 - logRSS.open.CDV.60day.indiv$h2

#### habitat model

p.open.CDV.60day.1.indiv.habitat <- p.open.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 60)


h1.open.CDV.60day.indiv.habitat <- data.table(rbindlist(p.open.CDV.60day.1.indiv.habitat),
                                                ttd = '60 days', COD = 'CDV', var = 'open', x = seq(from = 0, to = 1, length.out = 100))



h2.open.CDV.60day.indiv.habitat <- data.table(rbindlist(p.CDV.60day.2.indiv.habitat),
                                                ttd = '60 days', COD = 'CDV', var = 'open')


logRSS.open.CDV.60day.indiv.habitat <- merge(h1.open.CDV.60day.indiv.habitat, h2.open.CDV.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.CDV.60day.indiv.habitat[,'rss'] <- logRSS.open.CDV.60day.indiv.habitat$h1 - logRSS.open.CDV.60day.indiv.habitat$h2



#### human ####

open.human.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('open', levels = levels(land_end_adj)),
    propforest_end_adj = mean(propforest_end_adj, na.rm = T),
    propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
    propwet_end= mean(propwet_end, na.rm = T),
    roadDist_end = median(roadDist_end, na.rm = T),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('human', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.open.human.1day.1 <- predict(everyone, newdata = open.human.1day.1, type='link', re.form = NA)
p.open.human.1day.1.habitat <- predict(everyone.habitat, newdata = open.human.1day.1, type='link', re.form = NA)


logRSS.open.human.1day<- p.open.human.1day.1 - p.human.1day.2
logRSS.open.human.1day.habitat<- p.open.human.1day.1.habitat - p.human.1day.2.habitat


open.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                  log_sl = mean(log_sl),
                                                                  cos_ta = mean(cos_ta),
                                                                  #land_end_adj = factor('open', levels = levels(land_end_adj)),
                                                                propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                                propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                                  propwet_end= mean(propwet_end, na.rm = T),
                                                                  roadDist_end = median(roadDist_end, na.rm = T),
                                                                  distance2 = median(distance2, na.rm = T),
                                                                  nnDist_end = median(nnDist_end, na.rm = T),
                                                                  packDist_end = median(packDist_end, na.rm = T),
                                                                  COD = factor('human', levels = levels(COD)),
                                                                  ttd1 = 60,
                                                                  wolf_step_id = NA, wolfID = NA)]

p.open.human.60day.1 <- predict(everyone, newdata = open.human.60day.1, type='link', re.form = NA)
p.open.human.60day.1.habitat <- predict(everyone.habitat, newdata = open.human.60day.1, type='link', re.form = NA)


logRSS.open.human.60day<- p.open.human.60day.1 - p.human.60day.2
logRSS.open.human.60day.habitat<- p.open.human.60day.1.habitat - p.human.60day.2.habitat

logRSS.open <- rbind(logRSS.open, data.frame(COD= 'human', ttd = '1 day', var = 'open', rss = logRSS.open.human.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.open <- rbind(logRSS.open, data.frame(COD= 'human', ttd = '60 days', var = 'open', rss = logRSS.open.human.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.open.habitat <- rbind(logRSS.open.habitat, data.frame(COD= 'human', ttd = '1 day', var = 'open', rss = logRSS.open.human.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.open.habitat <- rbind(logRSS.open.habitat, data.frame(COD= 'human', ttd = '60 days', var = 'open', rss = logRSS.open.human.60day, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####
p.open.human.1day.1.indiv <- p.open.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)

h1.open.human.1day.indiv <- data.table(rbindlist(p.open.human.1day.1.indiv),
                                         ttd = '1 day', COD = 'human', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.human.1day.indiv <- data.table(rbindlist(p.human.1day.2.indiv),
                                         ttd = '1 day', COD = 'human', var = 'open')



logRSS.open.human.1day.indiv <- merge(h1.open.human.1day.indiv, h2.open.human.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.human.1day.indiv[,'rss'] <- logRSS.open.human.1day.indiv$h1 - logRSS.open.human.1day.indiv$h2

#### habitat model

p.open.human.1day.1.indiv.habitat <- p.open.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 1)

h1.open.human.1day.indiv.habitat <- data.table(rbindlist(p.open.human.1day.1.indiv.habitat),
                                                 ttd = '1 day', COD = 'human', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.human.1day.indiv.habitat <- data.table(rbindlist(p.human.1day.2.indiv.habitat),
                                                 ttd = '1 day', COD = 'human', var = 'open')


logRSS.open.human.1day.indiv.habitat <- merge(h1.open.human.1day.indiv.habitat, h2.open.human.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.human.1day.indiv.habitat[,'rss'] <- logRSS.open.human.1day.indiv.habitat$h1 - logRSS.open.human.1day.indiv.habitat$h2


### 60 days

p.open.human.60day.1.indiv <- p.open.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)


h1.open.human.60day.indiv <- data.table(rbindlist(p.open.human.60day.1.indiv),
                                          ttd = '60 days', COD = 'human', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.human.60day.indiv <- data.table(rbindlist(p.human.60day.2.indiv),
                                          ttd = '60 days', COD = 'human', var = 'open')



logRSS.open.human.60day.indiv <- merge(h1.open.human.60day.indiv, h2.open.human.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.human.60day.indiv[,'rss'] <- logRSS.open.human.60day.indiv$h1 - logRSS.open.human.60day.indiv$h2

#### habitat model

p.open.human.60day.1.indiv.habitat <- p.open.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 60)


h1.open.human.60day.indiv.habitat <- data.table(rbindlist(p.open.human.60day.1.indiv.habitat),
                                                  ttd = '60 days', COD = 'human', var = 'open', x = seq(from = 0, to = 1, length.out = 100))



h2.open.human.60day.indiv.habitat <- data.table(rbindlist(p.human.60day.2.indiv.habitat),
                                                  ttd = '60 days', COD = 'human', var = 'open')


logRSS.open.human.60day.indiv.habitat <- merge(h1.open.human.60day.indiv.habitat, h2.open.human.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.human.60day.indiv.habitat[,'rss'] <- logRSS.open.human.60day.indiv.habitat$h1 - logRSS.open.human.60day.indiv.habitat$h2




#### control ####

open.control.1day.1 <-dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('open', levels = levels(land_end_adj)),
    propforest_end_adj = mean(propforest_end_adj, na.rm = T),
    propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
    propwet_end= mean(propwet_end, na.rm = T),
    roadDist_end = median(roadDist_end, na.rm = T),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('control', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.open.control.1day.1 <- predict(everyone, newdata = open.control.1day.1, type='link', re.form = NA)
p.open.control.1day.1.habitat <- predict(everyone.habitat, newdata = open.control.1day.1, type='link', re.form = NA)


logRSS.open.control.1day<- p.open.control.1day.1 - p.control.1day.2
logRSS.open.control.1day.habitat<- p.open.control.1day.1.habitat - p.control.1day.2.habitat


open.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                    log_sl = mean(log_sl),
                                                                    cos_ta = mean(cos_ta),
                                                                    #land_end_adj = factor('open', levels = levels(land_end_adj)),
                                                                  propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                                  propopen_end_adj = seq(from = 0, to = 1, length.out = 100),
                                                                    propwet_end= mean(propwet_end, na.rm = T),
                                                                    roadDist_end = median(roadDist_end, na.rm = T),
                                                                    distance2 = median(distance2, na.rm = T),
                                                                    nnDist_end = median(nnDist_end, na.rm = T),
                                                                    packDist_end = median(packDist_end, na.rm = T),
                                                                    COD = factor('control', levels = levels(COD)),
                                                                    ttd1 = 60,
                                                                    wolf_step_id = NA, wolfID = NA)]

p.open.control.60day.1 <- predict(everyone, newdata = open.control.60day.1, type='link', re.form = NA)
p.open.control.60day.1.habitat <- predict(everyone.habitat, newdata = open.control.60day.1, type='link', re.form = NA)


logRSS.open.control.60day<- p.open.control.60day.1 - p.control.60day.2
logRSS.open.control.60day.habitat<- p.open.control.60day.1.habitat - p.control.60day.2.habitat

logRSS.open <- rbind(logRSS.open, data.frame(COD= 'control', ttd = '1 day', var = 'open', rss = logRSS.open.control.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.open <- rbind(logRSS.open, data.frame(COD= 'control', ttd = '60 days', var = 'open', rss = logRSS.open.control.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.open.habitat <- rbind(logRSS.open.habitat, data.frame(COD= 'control', ttd = '1 day', var = 'open', rss = logRSS.open.control.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.open.habitat <- rbind(logRSS.open.habitat, data.frame(COD= 'control', ttd = '60 days', var = 'open', rss = logRSS.open.control.60day, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####
p.open.control.1day.1.indiv <- p.open.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)


h1.open.control.1day.indiv <- data.table(rbindlist(p.open.control.1day.1.indiv),
                                           ttd = '1 day', COD = 'control', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.control.1day.indiv <- data.table(rbindlist(p.control.1day.2.indiv),
                                           ttd = '1 day', COD = 'control', var = 'open')



logRSS.open.control.1day.indiv <- merge(h1.open.control.1day.indiv, h2.open.control.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.control.1day.indiv[,'rss'] <- logRSS.open.control.1day.indiv$h1 - logRSS.open.control.1day.indiv$h2

#### habitat model

p.open.control.1day.1.indiv.habitat <- p.open.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 1)

h1.open.control.1day.indiv.habitat <- data.table(rbindlist(p.open.control.1day.1.indiv.habitat),
                                                   ttd = '1 day', COD = 'control', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.control.1day.indiv.habitat <- data.table(rbindlist(p.control.1day.2.indiv.habitat),
                                                   ttd = '1 day', COD = 'control', var = 'open')


logRSS.open.control.1day.indiv.habitat <- merge(h1.open.control.1day.indiv.habitat, h2.open.control.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.control.1day.indiv.habitat[,'rss'] <- logRSS.open.control.1day.indiv.habitat$h1 - logRSS.open.control.1day.indiv.habitat$h2


### 60 days

p.open.control.60day.1.indiv <- p.open.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)


h1.open.control.60day.indiv <- data.table(rbindlist(p.open.control.60day.1.indiv),
                                            ttd = '60 days', COD = 'control', var = 'open', x = seq(from = 0, to = 1, length.out = 100))


h2.open.control.60day.indiv <- data.table(rbindlist(p.control.60day.2.indiv),
                                            ttd = '60 days', COD = 'control', var = 'open')



logRSS.open.control.60day.indiv <- merge(h1.open.control.60day.indiv, h2.open.control.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.control.60day.indiv[,'rss'] <- logRSS.open.control.60day.indiv$h1 - logRSS.open.control.60day.indiv$h2

#### habitat model

p.open.control.60day.1.indiv.habitat <- p.open.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 60)

h1.open.control.60day.indiv.habitat <- data.table(rbindlist(p.open.control.60day.1.indiv.habitat),
                                                    ttd = '60 days', COD = 'control', var = 'open', x = seq(from = 0, to = 1, length.out = 100))



h2.open.control.60day.indiv.habitat <- data.table(rbindlist(p.control.60day.2.indiv.habitat),
                                                    ttd = '60 days', COD = 'control', var = 'open')


logRSS.open.control.60day.indiv.habitat <- merge(h1.open.control.60day.indiv.habitat, h2.open.control.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.open.control.60day.indiv.habitat[,'rss'] <- logRSS.open.control.60day.indiv.habitat$h1 - logRSS.open.control.60day.indiv.habitat$h2





logRSS.open.indiv <- rbind(logRSS.open.control.1day.indiv, logRSS.open.control.60day.indiv, 
                             logRSS.open.human.1day.indiv, logRSS.open.human.60day.indiv, 
                             logRSS.open.CDV.1day.indiv, logRSS.open.CDV.60day.indiv)

logRSS.open.indiv.habitat <- rbind(logRSS.open.control.1day.indiv.habitat, logRSS.open.control.60day.indiv.habitat, 
                                     logRSS.open.human.1day.indiv.habitat, logRSS.open.human.60day.indiv.habitat, 
                                     logRSS.open.CDV.1day.indiv.habitat, logRSS.open.CDV.60day.indiv.habitat)





#### WET ####
#### CDV ####
wet.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             #land_end_adj = factor('wet', levels = levels(land_end_adj)),
                                                             propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= seq(from = 0, to = 1, length.out = 100),
                                                             roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = median(distance2, na.rm = T),
                                                             nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA,
                                                             wolfID = NA)]

p.wet.CDV.1day.1 <- predict(everyone, newdata = wet.CDV.1day.1, type='link', re.form = NA)
p.wet.CDV.1day.1.habitat <- predict(everyone.habitat, newdata = wet.CDV.1day.1, type='link', re.form = NA)


logRSS.wet.CDV.1day<- p.wet.CDV.1day.1 - p.CDV.1day.2
logRSS.wet.CDV.1day.habitat<- p.wet.CDV.1day.1.habitat - p.CDV.1day.2.habitat

wet.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                              log_sl = mean(log_sl),
                                                              cos_ta = mean(cos_ta),
                                                              #land_end_adj = factor('wet', levels = levels(land_end_adj)),
                                                              propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= seq(from = 0, to = 1, length.out = 100),
                                                              roadDist_end = median(roadDist_end, na.rm = T),
                                                              distance2 = median(distance2, na.rm = T),
                                                              nnDist_end = median(nnDist_end, na.rm = T),
                                                              packDist_end = median(packDist_end, na.rm = T),
                                                              COD = factor('CDV', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]
p.wet.CDV.60day.1 <- predict(everyone, newdata = wet.CDV.60day.1, type='link', re.form =NA)
p.wet.CDV.60day.1.habitat <- predict(everyone.habitat, newdata = wet.CDV.60day.1, type='link', re.form =NA)


logRSS.wet.CDV.60day<- p.wet.CDV.60day.1 - p.CDV.60day.2
logRSS.wet.CDV.60day.habitat<- p.wet.CDV.60day.1.habitat - p.CDV.60day.2.habitat

logRSS.wet <- data.frame(COD= 'CDV', ttd = '1 day', var = 'wet', rss = logRSS.wet.CDV.1day, x = seq(from = 0, to = 1, length.out = 100))
logRSS.wet <- rbind(logRSS.wet, data.frame(COD= 'CDV', ttd = '60 days', var = 'wet', rss = logRSS.wet.CDV.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.wet.habitat <- data.frame(COD= 'CDV', ttd = '1 day', var = 'wet', rss = logRSS.wet.CDV.1day.habitat, x = seq(from = 0, to = 1, length.out = 100))
logRSS.wet.habitat <- rbind(logRSS.wet.habitat, data.frame(COD= 'CDV', ttd = '60 days', var = 'wet', rss = logRSS.wet.CDV.60day.habitat, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####

p.wet.h1.indiv <- function(ids, DT, mod, death, t2death){
    lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          # land_end_adj = factor('wet', levels = levels(land_end_adj)),
          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
          propwet_end= seq(from = 0, to = 1, length.out = 100),
          roadDist_end = median(roadDist_end, na.rm = T),
          distance2 = median(distance2, na.rm = T),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = median(packDist_end, na.rm = T),
          COD = factor(death, levels = levels(COD)),
          ttd1 = t2death,
          wolf_step_id = NA,
          wolfID = i
        )],
        type = "link",
        re.form = NULL
      ), wolfID = i)]
    # )
  })
}

p.wet.CDV.1day.1.indiv <- p.wet.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)

h1.wet.CDV.1day.indiv <- data.table(rbindlist(p.wet.CDV.1day.1.indiv),
                                     ttd = '1 day', COD = 'CDV', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                     ttd = '1 day', COD = 'CDV', var = 'wet')



logRSS.wet.CDV.1day.indiv <- merge(h1.wet.CDV.1day.indiv, h2.wet.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.CDV.1day.indiv[,'rss'] <- logRSS.wet.CDV.1day.indiv$h1 - logRSS.wet.CDV.1day.indiv$h2

#### habitat model

p.wet.CDV.1day.1.indiv.habitat <- p.wet.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 1)


h1.wet.CDV.1day.indiv.habitat <- data.table(rbindlist(p.wet.CDV.1day.1.indiv.habitat),
                                             ttd = '1 day', COD = 'CDV', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.CDV.1day.indiv.habitat <- data.table(rbindlist(p.CDV.1day.2.indiv.habitat),
                                             ttd = '1 day', COD = 'CDV', var = 'wet')


logRSS.wet.CDV.1day.indiv.habitat <- merge(h1.wet.CDV.1day.indiv.habitat, h2.wet.CDV.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.CDV.1day.indiv.habitat[,'rss'] <- logRSS.wet.CDV.1day.indiv.habitat$h1 - logRSS.wet.CDV.1day.indiv.habitat$h2


### 60 days

p.wet.CDV.60day.1.indiv <- p.wet.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)

h1.wet.CDV.60day.indiv <- data.table(rbindlist(p.wet.CDV.60day.1.indiv),
                                      ttd = '60 days', COD = 'CDV', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.CDV.60day.indiv <- data.table(rbindlist(p.CDV.60day.2.indiv),
                                      ttd = '60 days', COD = 'CDV', var = 'wet')



logRSS.wet.CDV.60day.indiv <- merge(h1.wet.CDV.60day.indiv, h2.wet.CDV.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.CDV.60day.indiv[,'rss'] <- logRSS.wet.CDV.60day.indiv$h1 - logRSS.wet.CDV.60day.indiv$h2

#### habitat model

p.wet.CDV.60day.1.indiv.habitat <- p.wet.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 60)


h1.wet.CDV.60day.indiv.habitat <- data.table(rbindlist(p.wet.CDV.60day.1.indiv.habitat),
                                              ttd = '60 days', COD = 'CDV', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))



h2.wet.CDV.60day.indiv.habitat <- data.table(rbindlist(p.CDV.60day.2.indiv.habitat),
                                              ttd = '60 days', COD = 'CDV', var = 'wet')


logRSS.wet.CDV.60day.indiv.habitat <- merge(h1.wet.CDV.60day.indiv.habitat, h2.wet.CDV.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.CDV.60day.indiv.habitat[,'rss'] <- logRSS.wet.CDV.60day.indiv.habitat$h1 - logRSS.wet.CDV.60day.indiv.habitat$h2



#### human ####

wet.human.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('wet', levels = levels(land_end_adj)),
    propforest_end_adj = mean(propforest_end_adj, na.rm = T),
    propopen_end_adj = mean(propopen_end_adj, na.rm = T),
    propwet_end= seq(from = 0, to = 1, length.out = 100),
    roadDist_end = median(roadDist_end, na.rm = T),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('human', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.wet.human.1day.1 <- predict(everyone, newdata = wet.human.1day.1, type='link', re.form = NA)
p.wet.human.1day.1.habitat <- predict(everyone.habitat, newdata = wet.human.1day.1, type='link', re.form = NA)


logRSS.wet.human.1day<- p.wet.human.1day.1 - p.human.1day.2
logRSS.wet.human.1day.habitat<- p.wet.human.1day.1.habitat - p.human.1day.2.habitat


wet.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                #land_end_adj = factor('wet', levels = levels(land_end_adj)),
                                                                propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                               propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                               propwet_end= seq(from = 0, to = 1, length.out = 100),
                                                                roadDist_end = median(roadDist_end, na.rm = T),
                                                                distance2 = median(distance2, na.rm = T),
                                                                nnDist_end = median(nnDist_end, na.rm = T),
                                                                packDist_end = median(packDist_end, na.rm = T),
                                                                COD = factor('human', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]

p.wet.human.60day.1 <- predict(everyone, newdata = wet.human.60day.1, type='link', re.form = NA)
p.wet.human.60day.1.habitat <- predict(everyone.habitat, newdata = wet.human.60day.1, type='link', re.form = NA)


logRSS.wet.human.60day<- p.wet.human.60day.1 - p.human.60day.2
logRSS.wet.human.60day.habitat<- p.wet.human.60day.1.habitat - p.human.60day.2.habitat

logRSS.wet <- rbind(logRSS.wet, data.frame(COD= 'human', ttd = '1 day', var = 'wet', rss = logRSS.wet.human.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.wet <- rbind(logRSS.wet, data.frame(COD= 'human', ttd = '60 days', var = 'wet', rss = logRSS.wet.human.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.wet.habitat <- rbind(logRSS.wet.habitat, data.frame(COD= 'human', ttd = '1 day', var = 'wet', rss = logRSS.wet.human.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.wet.habitat <- rbind(logRSS.wet.habitat, data.frame(COD= 'human', ttd = '60 days', var = 'wet', rss = logRSS.wet.human.60day, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####
p.wet.human.1day.1.indiv <- p.wet.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)

h1.wet.human.1day.indiv <- data.table(rbindlist(p.wet.human.1day.1.indiv),
                                       ttd = '1 day', COD = 'human', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.human.1day.indiv <- data.table(rbindlist(p.human.1day.2.indiv),
                                       ttd = '1 day', COD = 'human', var = 'wet')



logRSS.wet.human.1day.indiv <- merge(h1.wet.human.1day.indiv, h2.wet.human.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.human.1day.indiv[,'rss'] <- logRSS.wet.human.1day.indiv$h1 - logRSS.wet.human.1day.indiv$h2

#### habitat model

p.wet.human.1day.1.indiv.habitat <- p.wet.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 1)


h1.wet.human.1day.indiv.habitat <- data.table(rbindlist(p.wet.human.1day.1.indiv.habitat),
                                               ttd = '1 day', COD = 'human', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.human.1day.indiv.habitat <- data.table(rbindlist(p.human.1day.2.indiv.habitat),
                                               ttd = '1 day', COD = 'human', var = 'wet')


logRSS.wet.human.1day.indiv.habitat <- merge(h1.wet.human.1day.indiv.habitat, h2.wet.human.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.human.1day.indiv.habitat[,'rss'] <- logRSS.wet.human.1day.indiv.habitat$h1 - logRSS.wet.human.1day.indiv.habitat$h2


### 60 days

p.wet.human.60day.1.indiv <-  p.wet.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)

h1.wet.human.60day.indiv <- data.table(rbindlist(p.wet.human.60day.1.indiv),
                                        ttd = '60 days', COD = 'human', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.human.60day.indiv <- data.table(rbindlist(p.human.60day.2.indiv),
                                        ttd = '60 days', COD = 'human', var = 'wet')



logRSS.wet.human.60day.indiv <- merge(h1.wet.human.60day.indiv, h2.wet.human.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.human.60day.indiv[,'rss'] <- logRSS.wet.human.60day.indiv$h1 - logRSS.wet.human.60day.indiv$h2

#### habitat model

p.wet.human.60day.1.indiv.habitat <- p.wet.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 60)


h1.wet.human.60day.indiv.habitat <- data.table(rbindlist(p.wet.human.60day.1.indiv.habitat),
                                                ttd = '60 days', COD = 'human', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))



h2.wet.human.60day.indiv.habitat <- data.table(rbindlist(p.human.60day.2.indiv.habitat),
                                                ttd = '60 days', COD = 'human', var = 'wet')


logRSS.wet.human.60day.indiv.habitat <- merge(h1.wet.human.60day.indiv.habitat, h2.wet.human.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.human.60day.indiv.habitat[,'rss'] <- logRSS.wet.human.60day.indiv.habitat$h1 - logRSS.wet.human.60day.indiv.habitat$h2




#### control ####

wet.control.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('wet', levels = levels(land_end_adj)),
    propforest_end_adj = mean(propforest_end_adj, na.rm = T),
    propopen_end_adj = mean(propopen_end_adj, na.rm = T),
    propwet_end= seq(from = 0, to = 1, length.out = 100),
    roadDist_end = median(roadDist_end, na.rm = T),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('control', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.wet.control.1day.1 <- predict(everyone, newdata = wet.control.1day.1, type='link', re.form = NA)
p.wet.control.1day.1.habitat <- predict(everyone.habitat, newdata = wet.control.1day.1, type='link', re.form = NA)


logRSS.wet.control.1day<- p.wet.control.1day.1 - p.control.1day.2
logRSS.wet.control.1day.habitat<- p.wet.control.1day.1.habitat - p.control.1day.2.habitat


wet.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                  log_sl = mean(log_sl),
                                                                  cos_ta = mean(cos_ta),
                                                                  #land_end_adj = factor('wet', levels = levels(land_end_adj)),
                                                                  propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                                 propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                                 propwet_end= seq(from = 0, to = 1, length.out = 100),
                                                                  roadDist_end = median(roadDist_end, na.rm = T),
                                                                  distance2 = median(distance2, na.rm = T),
                                                                  nnDist_end = median(nnDist_end, na.rm = T),
                                                                  packDist_end = median(packDist_end, na.rm = T),
                                                                  COD = factor('control', levels = levels(COD)),
                                                                  ttd1 = 60,
                                                                  wolf_step_id = NA, wolfID = NA)]

p.wet.control.60day.1 <- predict(everyone, newdata = wet.control.60day.1, type='link', re.form = NA)
p.wet.control.60day.1.habitat <- predict(everyone.habitat, newdata = wet.control.60day.1, type='link', re.form = NA)


logRSS.wet.control.60day<- p.wet.control.60day.1 - p.control.60day.2
logRSS.wet.control.60day.habitat<- p.wet.control.60day.1.habitat - p.control.60day.2.habitat

logRSS.wet <- rbind(logRSS.wet, data.frame(COD= 'control', ttd = '1 day', var = 'wet', rss = logRSS.wet.control.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.wet <- rbind(logRSS.wet, data.frame(COD= 'control', ttd = '60 days', var = 'wet', rss = logRSS.wet.control.60day, x = seq(from = 0, to = 1, length.out = 100)))

logRSS.wet.habitat <- rbind(logRSS.wet.habitat, data.frame(COD= 'control', ttd = '1 day', var = 'wet', rss = logRSS.wet.control.1day, x = seq(from = 0, to = 1, length.out = 100)))
logRSS.wet.habitat <- rbind(logRSS.wet.habitat, data.frame(COD= 'control', ttd = '60 days', var = 'wet', rss = logRSS.wet.control.60day, x = seq(from = 0, to = 1, length.out = 100)))

#### indivs ####
p.wet.control.1day.1.indiv <- p.wet.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)

h1.wet.control.1day.indiv <- data.table(rbindlist(p.wet.control.1day.1.indiv),
                                         ttd = '1 day', COD = 'control', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.control.1day.indiv <- data.table(rbindlist(p.control.1day.2.indiv),
                                         ttd = '1 day', COD = 'control', var = 'wet')



logRSS.wet.control.1day.indiv <- merge(h1.wet.control.1day.indiv, h2.wet.control.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.control.1day.indiv[,'rss'] <- logRSS.wet.control.1day.indiv$h1 - logRSS.wet.control.1day.indiv$h2

#### habitat model

p.wet.control.1day.1.indiv.habitat <- p.wet.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 1)

h1.wet.control.1day.indiv.habitat <- data.table(rbindlist(p.wet.control.1day.1.indiv.habitat),
                                                 ttd = '1 day', COD = 'control', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.control.1day.indiv.habitat <- data.table(rbindlist(p.control.1day.2.indiv.habitat),
                                                 ttd = '1 day', COD = 'control', var = 'wet')


logRSS.wet.control.1day.indiv.habitat <- merge(h1.wet.control.1day.indiv.habitat, h2.wet.control.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.control.1day.indiv.habitat[,'rss'] <- logRSS.wet.control.1day.indiv.habitat$h1 - logRSS.wet.control.1day.indiv.habitat$h2


### 60 days

p.wet.control.60day.1.indiv <- p.wet.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)


h1.wet.control.60day.indiv <- data.table(rbindlist(p.wet.control.60day.1.indiv),
                                          ttd = '60 days', COD = 'control', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))


h2.wet.control.60day.indiv <- data.table(rbindlist(p.control.60day.2.indiv),
                                          ttd = '60 days', COD = 'control', var = 'wet')



logRSS.wet.control.60day.indiv <- merge(h1.wet.control.60day.indiv, h2.wet.control.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.control.60day.indiv[,'rss'] <- logRSS.wet.control.60day.indiv$h1 - logRSS.wet.control.60day.indiv$h2

#### habitat model

p.wet.control.60day.1.indiv.habitat <-  p.wet.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 60)

h1.wet.control.60day.indiv.habitat <- data.table(rbindlist(p.wet.control.60day.1.indiv.habitat),
                                                  ttd = '60 days', COD = 'control', var = 'wet', x = seq(from = 0, to = 1, length.out = 100))



h2.wet.control.60day.indiv.habitat <- data.table(rbindlist(p.control.60day.2.indiv.habitat),
                                                  ttd = '60 days', COD = 'control', var = 'wet')


logRSS.wet.control.60day.indiv.habitat <- merge(h1.wet.control.60day.indiv.habitat, h2.wet.control.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.wet.control.60day.indiv.habitat[,'rss'] <- logRSS.wet.control.60day.indiv.habitat$h1 - logRSS.wet.control.60day.indiv.habitat$h2





logRSS.wet.indiv <- rbind(logRSS.wet.control.1day.indiv, logRSS.wet.control.60day.indiv, 
                           logRSS.wet.human.1day.indiv, logRSS.wet.human.60day.indiv, 
                           logRSS.wet.CDV.1day.indiv, logRSS.wet.CDV.60day.indiv)

logRSS.wet.indiv.habitat <- rbind(logRSS.wet.control.1day.indiv.habitat, logRSS.wet.control.60day.indiv.habitat, 
                                   logRSS.wet.human.1day.indiv.habitat, logRSS.wet.human.60day.indiv.habitat, 
                                   logRSS.wet.CDV.1day.indiv.habitat, logRSS.wet.CDV.60day.indiv.habitat)





#### ROAD DIST ####
#### CDV ####
road.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                          log_sl = mean(log_sl),
                          cos_ta = mean(cos_ta),
                          #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                          propwet_end= mean(propwet_end, na.rm = T),
                          roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                          distance2 = median(distance2, na.rm = T),
                          nnDist_end = median(nnDist_end, na.rm = T),
                          packDist_end = median(packDist_end, na.rm = T),
                          COD = factor('CDV', levels = levels(COD)),
                          ttd1 = 1,
                          wolf_step_id = NA,
                          wolfID = NA)]

p.road.CDV.1day.1 <- predict(everyone, newdata = road.CDV.1day.1, type='link', re.form = NA)
p.road.CDV.1day.1.habitat <- predict(everyone.habitat, newdata = road.CDV.1day.1, type='link', re.form = NA)


logRSS.road.CDV.1day<- p.road.CDV.1day.1 - p.CDV.1day.2
logRSS.road.CDV.1day.habitat<- p.road.CDV.1day.1.habitat - p.CDV.1day.2.habitat

road.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),
                                                             roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                             distance2 = median(distance2, na.rm = T),
                                                             nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID = NA)]
p.road.CDV.60day.1 <- predict(everyone, newdata = road.CDV.60day.1, type='link', re.form =NA)
p.road.CDV.60day.1.habitat <- predict(everyone.habitat, newdata = road.CDV.60day.1, type='link', re.form =NA)


logRSS.road.CDV.60day<- p.road.CDV.60day.1 - p.CDV.60day.2
logRSS.road.CDV.60day.habitat<- p.road.CDV.60day.1.habitat - p.CDV.60day.2.habitat

logRSS.road <- data.frame(COD= 'CDV', ttd = '1 day', var = 'road', rss = logRSS.road.CDV.1day, x = seq(from = 0, to = 3000, length.out = 100))
logRSS.road <- rbind(logRSS.road, data.frame(COD= 'CDV', ttd = '60 days', var = 'road', rss = logRSS.road.CDV.60day, x = seq(from = 0, to = 3000, length.out = 100)))

logRSS.road.habitat <- data.frame(COD= 'CDV', ttd = '1 day', var = 'road', rss = logRSS.road.CDV.1day.habitat, x = seq(from = 0, to = 3000, length.out = 100))
logRSS.road.habitat <- rbind(logRSS.road.habitat, data.frame(COD= 'CDV', ttd = '60 days', var = 'road', rss = logRSS.road.CDV.60day.habitat, x = seq(from = 0, to = 3000, length.out = 100)))

#### indivs ####

p.road.h1.indiv <- function(ids, DT, mod, death, t2death){
  lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
          propwet_end= mean(propwet_end, na.rm = T), 
          roadDist_end = seq(from = 0, to = 3000, length.out = 100),
          distance2 = median(distance2, na.rm = T),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = median(packDist_end, na.rm = T),
          COD = factor(death, levels = levels(COD)),
          ttd1 = t2death,
          wolf_step_id = NA,
          wolfID = i
        )],
        type = "link",
        re.form = NULL
      ), wolfID = i)]
    # )
  })
}

p.road.CDV.1day.1.indiv <- p.road.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)
  
h1.road.CDV.1day.indiv <- data.table(rbindlist(p.road.CDV.1day.1.indiv),
                                   ttd = '1 day', COD = 'CDV', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                   ttd = '1 day', COD = 'CDV', var = 'road')



logRSS.road.CDV.1day.indiv <- merge(h1.road.CDV.1day.indiv, h2.road.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.CDV.1day.indiv[,'rss'] <- logRSS.road.CDV.1day.indiv$h1 - logRSS.road.CDV.1day.indiv$h2

#### habitat model

p.road.CDV.1day.1.indiv.habitat <-  p.road.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 1)

h1.road.CDV.1day.indiv.habitat <- data.table(rbindlist(p.road.CDV.1day.1.indiv.habitat),
                                          ttd = '1 day', COD = 'CDV', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.CDV.1day.indiv.habitat <- data.table(rbindlist(p.CDV.1day.2.indiv.habitat),
                                          ttd = '1 day', COD = 'CDV', var = 'road')


logRSS.road.CDV.1day.indiv.habitat <- merge(h1.road.CDV.1day.indiv.habitat, h2.road.CDV.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.CDV.1day.indiv.habitat[,'rss'] <- logRSS.road.CDV.1day.indiv.habitat$h1 - logRSS.road.CDV.1day.indiv.habitat$h2


### 60 days

p.road.CDV.60day.1.indiv <-  p.road.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)


h1.road.CDV.60day.indiv <- data.table(rbindlist(p.road.CDV.60day.1.indiv),
                                    ttd = '60 days', COD = 'CDV', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.CDV.60day.indiv <- data.table(rbindlist(p.CDV.60day.2.indiv),
                                    ttd = '60 days', COD = 'CDV', var = 'road')



logRSS.road.CDV.60day.indiv <- merge(h1.road.CDV.60day.indiv, h2.road.CDV.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.CDV.60day.indiv[,'rss'] <- logRSS.road.CDV.60day.indiv$h1 - logRSS.road.CDV.60day.indiv$h2

#### habitat model

p.road.CDV.60day.1.indiv.habitat <- p.road.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 60)

h1.road.CDV.60day.indiv.habitat <- data.table(rbindlist(p.road.CDV.60day.1.indiv.habitat),
                                           ttd = '60 days', COD = 'CDV', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))



h2.road.CDV.60day.indiv.habitat <- data.table(rbindlist(p.CDV.60day.2.indiv.habitat),
                                           ttd = '60 days', COD = 'CDV', var = 'road')


logRSS.road.CDV.60day.indiv.habitat <- merge(h1.road.CDV.60day.indiv.habitat, h2.road.CDV.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.CDV.60day.indiv.habitat[,'rss'] <- logRSS.road.CDV.60day.indiv.habitat$h1 - logRSS.road.CDV.60day.indiv.habitat$h2



#### human ####

road.human.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('forest', levels = levels(land_end_adj)),
    propforest_end_adj = mean(propforest_end_adj, na.rm = T),
    propopen_end_adj = mean(propopen_end_adj, na.rm = T),
    propwet_end = mean(propwet_end, na.rm = T),
    roadDist_end = seq(
      from = 0,
      to = 3000,
      length.out = 100
    ),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('human', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.road.human.1day.1 <- predict(everyone, newdata = road.human.1day.1, type='link', re.form = NA)
p.road.human.1day.1.habitat <- predict(everyone.habitat, newdata = road.human.1day.1, type='link', re.form = NA)


logRSS.road.human.1day<- p.road.human.1day.1 - p.human.1day.2
logRSS.road.human.1day.habitat<- p.road.human.1day.1.habitat - p.human.1day.2.habitat


road.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                              log_sl = mean(log_sl),
                                                              cos_ta = mean(cos_ta),
                                                              #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                              propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                              propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                              propwet_end= mean(propwet_end, na.rm = T),
                                                              roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                              distance2 = median(distance2, na.rm = T),
                                                              nnDist_end = median(nnDist_end, na.rm = T),
                                                              packDist_end = median(packDist_end, na.rm = T),
                                                              COD = factor('human', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]

p.road.human.60day.1 <- predict(everyone, newdata = road.human.60day.1, type='link', re.form = NA)
p.road.human.60day.1.habitat <- predict(everyone.habitat, newdata = road.human.60day.1, type='link', re.form = NA)


logRSS.road.human.60day<- p.road.human.60day.1 - p.human.60day.2
logRSS.road.human.60day.habitat<- p.road.human.60day.1.habitat - p.human.60day.2.habitat

logRSS.road <- rbind(logRSS.road, data.frame(COD= 'human', ttd = '1 day', var = 'road', rss = logRSS.road.human.1day, x = seq(from = 0, to = 3000, length.out = 100)))
logRSS.road <- rbind(logRSS.road, data.frame(COD= 'human', ttd = '60 days', var = 'road', rss = logRSS.road.human.60day, x = seq(from = 0, to = 3000, length.out = 100)))

logRSS.road.habitat <- rbind(logRSS.road.habitat, data.frame(COD= 'human', ttd = '1 day', var = 'road', rss = logRSS.road.human.1day, x = seq(from = 0, to = 3000, length.out = 100)))
logRSS.road.habitat <- rbind(logRSS.road.habitat, data.frame(COD= 'human', ttd = '60 days', var = 'road', rss = logRSS.road.human.60day, x = seq(from = 0, to = 3000, length.out = 100)))

#### indivs ####
p.road.human.1day.1.indiv <- p.road.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)

h1.road.human.1day.indiv <- data.table(rbindlist(p.road.human.1day.1.indiv),
                                     ttd = '1 day', COD = 'human', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.human.1day.indiv <- data.table(rbindlist(p.human.1day.2.indiv),
                                     ttd = '1 day', COD = 'human', var = 'road')



logRSS.road.human.1day.indiv <- merge(h1.road.human.1day.indiv, h2.road.human.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.human.1day.indiv[,'rss'] <- logRSS.road.human.1day.indiv$h1 - logRSS.road.human.1day.indiv$h2

#### habitat model

p.road.human.1day.1.indiv.habitat <- p.road.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 1)

h1.road.human.1day.indiv.habitat <- data.table(rbindlist(p.road.human.1day.1.indiv.habitat),
                                             ttd = '1 day', COD = 'human', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.human.1day.indiv.habitat <- data.table(rbindlist(p.human.1day.2.indiv.habitat),
                                             ttd = '1 day', COD = 'human', var = 'road')


logRSS.road.human.1day.indiv.habitat <- merge(h1.road.human.1day.indiv.habitat, h2.road.human.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.human.1day.indiv.habitat[,'rss'] <- logRSS.road.human.1day.indiv.habitat$h1 - logRSS.road.human.1day.indiv.habitat$h2


### 60 days

p.road.human.60day.1.indiv <- p.road.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)


h1.road.human.60day.indiv <- data.table(rbindlist(p.road.human.60day.1.indiv),
                                      ttd = '60 days', COD = 'human', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.human.60day.indiv <- data.table(rbindlist(p.human.60day.2.indiv),
                                      ttd = '60 days', COD = 'human', var = 'road')



logRSS.road.human.60day.indiv <- merge(h1.road.human.60day.indiv, h2.road.human.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.human.60day.indiv[,'rss'] <- logRSS.road.human.60day.indiv$h1 - logRSS.road.human.60day.indiv$h2

#### habitat model

p.road.human.60day.1.indiv.habitat <- p.road.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 60)

h1.road.human.60day.indiv.habitat <- data.table(rbindlist(p.road.human.60day.1.indiv.habitat),
                                              ttd = '60 days', COD = 'human', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))



h2.road.human.60day.indiv.habitat <- data.table(rbindlist(p.human.60day.2.indiv.habitat),
                                              ttd = '60 days', COD = 'human', var = 'road')


logRSS.road.human.60day.indiv.habitat <- merge(h1.road.human.60day.indiv.habitat, h2.road.human.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.human.60day.indiv.habitat[,'rss'] <- logRSS.road.human.60day.indiv.habitat$h1 - logRSS.road.human.60day.indiv.habitat$h2




#### control ####

road.control.1day.1 <-
  dat[wolfID %chin% dat.wnn.lastmo$wolfID, .(
    ToD_start = factor('day', levels = levels(ToD_start)),
    log_sl = mean(log_sl),
    cos_ta = mean(cos_ta),
    # land_end_adj = factor('forest', levels = levels(land_end_adj)),
    propforest_end_adj = mean(propforest_end_adj, na.rm = T),
    propopen_end_adj = mean(propopen_end_adj, na.rm = T),
    propwet_end = mean(propwet_end, na.rm = T),
    roadDist_end = seq(
      from = 0,
      to = 3000,
      length.out = 100
    ),
    distance2 = median(distance2, na.rm = T),
    nnDist_end = median(nnDist_end, na.rm = T),
    packDist_end = median(packDist_end, na.rm = T),
    COD = factor('control', levels = levels(COD)),
    ttd1 = 1,
    wolf_step_id = NA,
    wolfID = NA
  )]

p.road.control.1day.1 <- predict(everyone, newdata = road.control.1day.1, type='link', re.form = NA)
p.road.control.1day.1.habitat <- predict(everyone.habitat, newdata = road.control.1day.1, type='link', re.form = NA)


logRSS.road.control.1day<- p.road.control.1day.1 - p.control.1day.2
logRSS.road.control.1day.habitat<- p.road.control.1day.1.habitat - p.control.1day.2.habitat


road.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                                propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                                propwet_end= mean(propwet_end, na.rm = T),
                                                                roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                distance2 = median(distance2, na.rm = T),
                                                                nnDist_end = median(nnDist_end, na.rm = T),
                                                                packDist_end = median(packDist_end, na.rm = T),
                                                                COD = factor('control', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]

p.road.control.60day.1 <- predict(everyone, newdata = road.control.60day.1, type='link', re.form = NA)
p.road.control.60day.1.habitat <- predict(everyone.habitat, newdata = road.control.60day.1, type='link', re.form = NA)


logRSS.road.control.60day<- p.road.control.60day.1 - p.control.60day.2
logRSS.road.control.60day.habitat<- p.road.control.60day.1.habitat - p.control.60day.2.habitat

logRSS.road <- rbind(logRSS.road, data.frame(COD= 'control', ttd = '1 day', var = 'road', rss = logRSS.road.control.1day, x = seq(from = 0, to = 3000, length.out = 100)))
logRSS.road <- rbind(logRSS.road, data.frame(COD= 'control', ttd = '60 days', var = 'road', rss = logRSS.road.control.60day, x = seq(from = 0, to = 3000, length.out = 100)))

logRSS.road.habitat <- rbind(logRSS.road.habitat, data.frame(COD= 'control', ttd = '1 day', var = 'road', rss = logRSS.road.control.1day, x = seq(from = 0, to = 3000, length.out = 100)))
logRSS.road.habitat <- rbind(logRSS.road.habitat, data.frame(COD= 'control', ttd = '60 days', var = 'road', rss = logRSS.road.control.60day, x = seq(from = 0, to = 3000, length.out = 100)))

#### indivs ####
p.road.control.1day.1.indiv <- p.road.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)

h1.road.control.1day.indiv <- data.table(rbindlist(p.road.control.1day.1.indiv),
                                       ttd = '1 day', COD = 'control', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.control.1day.indiv <- data.table(rbindlist(p.control.1day.2.indiv),
                                       ttd = '1 day', COD = 'control', var = 'road')



logRSS.road.control.1day.indiv <- merge(h1.road.control.1day.indiv, h2.road.control.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.control.1day.indiv[,'rss'] <- logRSS.road.control.1day.indiv$h1 - logRSS.road.control.1day.indiv$h2

#### habitat model

p.road.control.1day.1.indiv.habitat <-  p.road.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 1)


h1.road.control.1day.indiv.habitat <- data.table(rbindlist(p.road.control.1day.1.indiv.habitat),
                                               ttd = '1 day', COD = 'control', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.control.1day.indiv.habitat <- data.table(rbindlist(p.control.1day.2.indiv.habitat),
                                               ttd = '1 day', COD = 'control', var = 'road')


logRSS.road.control.1day.indiv.habitat <- merge(h1.road.control.1day.indiv.habitat, h2.road.control.1day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.control.1day.indiv.habitat[,'rss'] <- logRSS.road.control.1day.indiv.habitat$h1 - logRSS.road.control.1day.indiv.habitat$h2


### 60 days

p.road.control.60day.1.indiv <- p.road.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)


h1.road.control.60day.indiv <- data.table(rbindlist(p.road.control.60day.1.indiv),
                                        ttd = '60 days', COD = 'control', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


h2.road.control.60day.indiv <- data.table(rbindlist(p.control.60day.2.indiv),
                                        ttd = '60 days', COD = 'control', var = 'road')



logRSS.road.control.60day.indiv <- merge(h1.road.control.60day.indiv, h2.road.control.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.control.60day.indiv[,'rss'] <- logRSS.road.control.60day.indiv$h1 - logRSS.road.control.60day.indiv$h2

#### habitat model

p.road.control.60day.1.indiv.habitat <- p.road.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 60)


h1.road.control.60day.indiv.habitat <- data.table(rbindlist(p.road.control.60day.1.indiv.habitat),
                                                ttd = '60 days', COD = 'control', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))



h2.road.control.60day.indiv.habitat <- data.table(rbindlist(p.control.60day.2.indiv.habitat),
                                                ttd = '60 days', COD = 'control', var = 'road')


logRSS.road.control.60day.indiv.habitat <- merge(h1.road.control.60day.indiv.habitat, h2.road.control.60day.indiv.habitat, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.control.60day.indiv.habitat[,'rss'] <- logRSS.road.control.60day.indiv.habitat$h1 - logRSS.road.control.60day.indiv.habitat$h2





logRSS.road.indiv <- rbind(logRSS.road.control.1day.indiv, logRSS.road.control.60day.indiv, 
                           logRSS.road.human.1day.indiv, logRSS.road.human.60day.indiv, 
                           logRSS.road.CDV.1day.indiv, logRSS.road.CDV.60day.indiv)

logRSS.road.indiv.habitat <- rbind(logRSS.road.control.1day.indiv.habitat, logRSS.road.control.60day.indiv.habitat, 
                           logRSS.road.human.1day.indiv.habitat, logRSS.road.human.60day.indiv.habitat, 
                           logRSS.road.CDV.1day.indiv.habitat, logRSS.road.CDV.60day.indiv.habitat)





#### NN DIST ####
#### CDV ####
nn.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(pop = factor('GHA26', levels = levels(pop)), ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),  
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = seq(from = 0, to = 7500, length.out = 150),
                                                           nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA, wolfID = NA)]

p.nn.CDV.1day.1 <- predict(everyone, newdata = nn.CDV.1day.1, type='link', re.form = NA)
p.nn.CDV.1day.1.social <- predict(everyone.social, newdata = nn.CDV.1day.1, type='link', re.form = NA)


logRSS.nn.CDV.1day<- p.nn.CDV.1day.1 - p.CDV.1day.2
logRSS.nn.CDV.1day.social<- p.nn.CDV.1day.1.social - p.CDV.1day.2.social


nn.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(pop = factor('GHA26', levels = levels(pop)), ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T), 
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = seq(from = 0, to = 7500, length.out = 150),
                                                           nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('CDV', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]
p.nn.CDV.60day.1 <- predict(everyone, newdata = nn.CDV.60day.1, type='link', re.form=NA)
p.nn.CDV.60day.1.social <- predict(everyone.social, newdata = nn.CDV.60day.1, type='link', re.form=NA)


logRSS.nn.CDV.60day<- p.nn.CDV.60day.1 - p.CDV.60day.2
logRSS.nn.CDV.60day.social<- p.nn.CDV.60day.1.social - p.CDV.60day.2.social


logRSS.nn <- data.frame(COD= 'CDV', ttd = '1 day', var = 'nn', rss = logRSS.nn.CDV.1day, x = seq(from = 0, to = 7500, length.out = 150))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'CDV', ttd = '60 days', var = 'nn', rss = logRSS.nn.CDV.60day, x = seq(from = 0, to = 7500, length.out = 150)))

logRSS.social <- data.frame(COD= 'CDV', ttd = '1 day', var = 'nn', rss = logRSS.nn.CDV.1day.social, x = seq(from = 0, to = 7500, length.out = 150))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '60 days', var = 'nn', rss = logRSS.nn.CDV.60day.social, x = seq(from = 0, to = 7500, length.out = 150)))

#### indivs ####

p.nn.h1.indiv <- function(ids, DT, mod, death, t2death){
    lapply(ids, function(i) {
    #unique(
      DT[#wolfID == i,
                    ,.(h1 = predict(
                      mod,
                      newdata = .SD[, .(
                        ToD_start = factor('day', levels = levels(ToD_start)),
                        log_sl = mean(log_sl),
                        cos_ta = mean(cos_ta),
                        # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                        propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                        propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                        propwet_end= mean(propwet_end, na.rm = T), 
                        roadDist_end = median(roadDist_end, na.rm = T),
                        distance2 = seq(from = 0, to = 7500, length.out = 150),
                        nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                        packDist_end = median(packDist_end, na.rm = T),
                        COD = factor(death, levels = levels(COD)),
                        ttd1 = t2death,
                        wolf_step_id = NA,
                        wolfID = i
                      )],
                      type = "link",
                      re.form = NULL
                    ), wolfID = i)]
   # )
  })
}

p.nn.CDV.1day.1.indiv <- p.nn.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)
  
h1.nn.CDV.1day.indiv <- data.table(rbindlist(p.nn.CDV.1day.1.indiv),
            ttd = '1 day', COD = 'CDV', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                   ttd = '1 day', COD = 'CDV', var = 'nn')



logRSS.nn.CDV.1day.indiv <- merge(h1.nn.CDV.1day.indiv, h2.nn.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.CDV.1day.indiv[,'rss'] <- logRSS.nn.CDV.1day.indiv$h1 - logRSS.nn.CDV.1day.indiv$h2

#### social model

p.nn.CDV.1day.1.indiv.social <- p.nn.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 1)

h1.nn.CDV.1day.indiv.social <- data.table(rbindlist(p.nn.CDV.1day.1.indiv.social),
                                   ttd = '1 day', COD = 'CDV', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.CDV.1day.indiv.social <- data.table(rbindlist(p.CDV.1day.2.indiv.social),
                                   ttd = '1 day', COD = 'CDV', var = 'nn')


logRSS.nn.CDV.1day.indiv.social <- merge(h1.nn.CDV.1day.indiv.social, h2.nn.CDV.1day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.CDV.1day.indiv.social[,'rss'] <- logRSS.nn.CDV.1day.indiv.social$h1 - logRSS.nn.CDV.1day.indiv.social$h2


### 60 days

p.nn.CDV.60day.1.indiv <- p.nn.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)

h1.nn.CDV.60day.indiv <- data.table(rbindlist(p.nn.CDV.60day.1.indiv),
                                   ttd = '60 days', COD = 'CDV', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.CDV.60day.indiv <- data.table(rbindlist(p.CDV.60day.2.indiv),
                                   ttd = '60 days', COD = 'CDV', var = 'nn')



logRSS.nn.CDV.60day.indiv <- merge(h1.nn.CDV.60day.indiv, h2.nn.CDV.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.CDV.60day.indiv[,'rss'] <- logRSS.nn.CDV.60day.indiv$h1 - logRSS.nn.CDV.60day.indiv$h2

#### social model

p.nn.CDV.60day.1.indiv.social <- p.nn.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 60)

h1.nn.CDV.60day.indiv.social <- data.table(rbindlist(p.nn.CDV.60day.1.indiv.social),
                                          ttd = '60 days', COD = 'CDV', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))



h2.nn.CDV.60day.indiv.social <- data.table(rbindlist(p.CDV.60day.2.indiv.social),
                                          ttd = '60 days', COD = 'CDV', var = 'nn')


logRSS.nn.CDV.60day.indiv.social <- merge(h1.nn.CDV.60day.indiv.social, h2.nn.CDV.60day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.CDV.60day.indiv.social[,'rss'] <- logRSS.nn.CDV.60day.indiv.social$h1 - logRSS.nn.CDV.60day.indiv.social$h2


#### human ####

nn.human.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),  
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = seq(from = 0, to = 7500, length.out = 150),
                                                           nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('human', levels = levels(COD)),
                                                           ttd1 = 1,
                                                           wolf_step_id = NA, wolfID = NA)]

p.nn.human.1day.1 <- predict(everyone, newdata = nn.human.1day.1, type='link', re.form = NA)
p.nn.human.1day.1.social <- predict(everyone.social, newdata = nn.human.1day.1, type='link', re.form = NA)


logRSS.nn.human.1day<- p.nn.human.1day.1 - p.human.1day.2
logRSS.nn.human.1day.social<- p.nn.human.1day.1.social - p.human.1day.2.social


nn.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T), 
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = seq(from = 0, to = 7500, length.out = 150),
                                                            nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('human', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID = NA)]
p.nn.human.60day.1 <- predict(everyone, newdata = nn.human.60day.1, type='link', re.form=NA)
p.nn.human.60day.1.social <- predict(everyone.social, newdata = nn.human.60day.1, type='link', re.form=NA)



logRSS.nn.human.60day<- p.nn.human.60day.1 - p.human.60day.2
logRSS.nn.human.60day.social<- p.nn.human.60day.1.social - p.human.60day.2.social

logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'human', ttd = '1 day', var = 'nn', rss = logRSS.nn.human.1day, x = seq(from = 0, to = 7500, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'human', ttd = '60 days', var = 'nn', rss = logRSS.nn.human.60day, x = seq(from = 0, to = 7500, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '1 day', var = 'nn', rss = logRSS.nn.human.1day.social, x = seq(from = 0, to = 7500, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '60 days', var = 'nn', rss = logRSS.nn.human.60day.social, x = seq(from = 0, to = 7500, length.out = 150)))

#### indivs ####

p.nn.human.1day.1.indiv <- p.nn.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)

h1.nn.human.1day.indiv <- data.table(rbindlist(p.nn.human.1day.1.indiv),
                                   ttd = '1 day', COD = 'human', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))

h2.nn.human.1day.indiv <- data.table(rbindlist(p.human.1day.2.indiv),
                                   ttd = '1 day', COD = 'human', var = 'nn')



logRSS.nn.human.1day.indiv <- merge(h1.nn.human.1day.indiv, h2.nn.human.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.human.1day.indiv[,'rss'] <- logRSS.nn.human.1day.indiv$h1 - logRSS.nn.human.1day.indiv$h2

#### social model

p.nn.human.1day.1.indiv.social <- p.nn.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 1)

h1.nn.human.1day.indiv.social <- data.table(rbindlist(p.nn.human.1day.1.indiv.social),
                                          ttd = '1 day', COD = 'human', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.human.1day.indiv.social <- data.table(rbindlist(p.human.1day.2.indiv.social),
                                          ttd = '1 day', COD = 'human', var = 'nn')


logRSS.nn.human.1day.indiv.social <- merge(h1.nn.human.1day.indiv.social, h2.nn.human.1day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.human.1day.indiv.social[,'rss'] <- logRSS.nn.human.1day.indiv.social$h1 - logRSS.nn.human.1day.indiv.social$h2


### 60 days

p.nn.human.60day.1.indiv <- p.nn.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)

h1.nn.human.60day.indiv <- data.table(rbindlist(p.nn.human.60day.1.indiv),
                                    ttd = '60 days', COD = 'human', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))

h2.nn.human.60day.indiv <- data.table(rbindlist(p.human.60day.2.indiv),
                                    ttd = '60 days', COD = 'human', var = 'nn')



logRSS.nn.human.60day.indiv <- merge(h1.nn.human.60day.indiv, h2.nn.human.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.human.60day.indiv[,'rss'] <- logRSS.nn.human.60day.indiv$h1 - logRSS.nn.human.60day.indiv$h2

#### social model

p.nn.human.60day.1.indiv.social <- p.nn.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 60)

h1.nn.human.60day.indiv.social <- data.table(rbindlist(p.nn.human.60day.1.indiv.social),
                                           ttd = '60 days', COD = 'human', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.human.60day.indiv.social <- data.table(rbindlist(p.human.60day.2.indiv.social),
                                           ttd = '60 days', COD = 'human', var = 'nn')


logRSS.nn.human.60day.indiv.social <- merge(h1.nn.human.60day.indiv.social, h2.nn.human.60day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.human.60day.indiv.social[,'rss'] <- logRSS.nn.human.60day.indiv.social$h1 - logRSS.nn.human.60day.indiv.social$h2


#### control ####

nn.control.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),  
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = seq(from = 0, to = 7500, length.out = 150),
                                                           nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('control', levels = levels(COD)),
                                                           ttd1 = 1,
                                                           wolf_step_id = NA, wolfID = NA)]

p.nn.control.1day.1 <- predict(everyone, newdata = nn.control.1day.1, type='link', re.form = NA)
p.nn.control.1day.1.social <- predict(everyone.social, newdata = nn.control.1day.1, type='link', re.form = NA)


logRSS.nn.control.1day<- p.nn.control.1day.1 - p.control.1day.2
logRSS.nn.control.1day.social<- p.nn.control.1day.1.social - p.control.1day.2.social


nn.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T), 
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = seq(from = 0, to = 7500, length.out = 150),
                                                            nnDist_end = seq(from = 0, to = 7500, length.out = 150),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('control', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID = NA)]
p.nn.control.60day.1 <- predict(everyone, newdata = nn.control.60day.1, type='link', re.form=NA)
p.nn.control.60day.1.social <- predict(everyone.social, newdata = nn.control.60day.1, type='link', re.form=NA)



logRSS.nn.control.60day<- p.nn.control.60day.1 - p.control.60day.2
logRSS.nn.control.60day.social<- p.nn.control.60day.1.social - p.control.60day.2.social

logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'control', ttd = '1 day', var = 'nn', rss = logRSS.nn.control.1day, x = seq(from = 0, to = 7500, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'control', ttd = '60 days', var = 'nn', rss = logRSS.nn.control.60day, x = seq(from = 0, to = 7500, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '1 day', var = 'nn', rss = logRSS.nn.control.1day.social, x = seq(from = 0, to = 7500, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '60 days', var = 'nn', rss = logRSS.nn.control.60day.social, x = seq(from = 0, to = 7500, length.out = 150)))


#### indivs ####

p.nn.control.1day.1.indiv <- p.nn.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)


h1.nn.control.1day.indiv <- data.table(rbindlist(p.nn.control.1day.1.indiv),
                                   ttd = '1 day', COD = 'control', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.control.1day.indiv <- data.table(rbindlist(p.control.1day.2.indiv),
                                   ttd = '1 day', COD = 'control', var = 'nn')



logRSS.nn.control.1day.indiv <- merge(h1.nn.control.1day.indiv, h2.nn.control.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.control.1day.indiv[,'rss'] <- logRSS.nn.control.1day.indiv$h1 - logRSS.nn.control.1day.indiv$h2

#### social model

p.nn.control.1day.1.indiv.social <- p.nn.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 1)

h1.nn.control.1day.indiv.social <- data.table(rbindlist(p.nn.control.1day.1.indiv.social),
                                          ttd = '1 day', COD = 'control', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))


h2.nn.control.1day.indiv.social <- data.table(rbindlist(p.control.1day.2.indiv.social),
                                          ttd = '1 day', COD = 'control', var = 'nn')


logRSS.nn.control.1day.indiv.social <- merge(h1.nn.control.1day.indiv.social, h2.nn.control.1day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.control.1day.indiv.social[,'rss'] <- logRSS.nn.control.1day.indiv.social$h1 - logRSS.nn.control.1day.indiv.social$h2


### 60 days

p.nn.control.60day.1.indiv <- p.nn.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)

h1.nn.control.60day.indiv <- data.table(rbindlist(p.nn.control.60day.1.indiv),
                                    ttd = '60 days', COD = 'control', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))



h2.nn.control.60day.indiv <- data.table(rbindlist(p.control.60day.2.indiv),
                                    ttd = '60 days', COD = 'control', var = 'nn')



logRSS.nn.control.60day.indiv <- merge(h1.nn.control.60day.indiv, h2.nn.control.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.control.60day.indiv[,'rss'] <- logRSS.nn.control.60day.indiv$h1 - logRSS.nn.control.60day.indiv$h2

#### social model

p.nn.control.60day.1.indiv.social <- p.nn.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 60)

h1.nn.control.60day.indiv.social <- data.table(rbindlist(p.nn.control.60day.1.indiv.social),
                                           ttd = '60 days', COD = 'control', var = 'nn', x = seq(from = 0, to = 7500, length.out = 150))



h2.nn.control.60day.indiv.social <- data.table(rbindlist(p.control.60day.2.indiv.social),
                                           ttd = '60 days', COD = 'control', var = 'nn')


logRSS.nn.control.60day.indiv.social <- merge(h1.nn.control.60day.indiv.social, h2.nn.control.60day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.nn.control.60day.indiv.social[,'rss'] <- logRSS.nn.control.60day.indiv.social$h1 - logRSS.nn.control.60day.indiv.social$h2




##### gathering RSS
logRSS.nn.indiv <- rbind(logRSS.nn.control.1day.indiv, logRSS.nn.control.60day.indiv, 
                           logRSS.nn.human.1day.indiv, logRSS.nn.human.60day.indiv, 
                           logRSS.nn.CDV.1day.indiv, logRSS.nn.CDV.60day.indiv)

logRSS.nn.indiv.social <- rbind(logRSS.nn.control.1day.indiv.social, logRSS.nn.control.60day.indiv.social, 
                           logRSS.nn.human.1day.indiv.social, logRSS.nn.human.60day.indiv.social, 
                           logRSS.nn.CDV.1day.indiv.social, logRSS.nn.CDV.60day.indiv.social)

#saveRDS(logRSS.nn.indiv, 'data/derived-data/logRSS_nn_indiv.Rds')
#saveRDS(logRSS.nn.indiv.social, 'data/derived-data/logRSS_nn_indiv_social.Rds')


#saveRDS(logRSS.nn, 'data/derived-data/logRSS_nn.Rds')
#saveRDS(logRSS.social, 'data/derived-data/logRSS_social.Rds')


#### PACK EDGE ####
#### CDV ####
pack.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),  
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = median(distance2, na.rm = T),
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = seq(0, 15000, length.out = 150),
                                                           COD = factor('CDV', levels = levels(COD)),
                                                           ttd1 = 1,
                                                           wolf_step_id = NA, wolfID = NA)]

p.pack.CDV.1day.1 <- predict(everyone, newdata = pack.CDV.1day.1, type='link', re.form = NA)
p.pack.CDV.1day.1.social <- predict(everyone.social, newdata = pack.CDV.1day.1, type='link', re.form = NA)
p.pack.CDV.1day.1.pack <- predict(everyone.pack, newdata = pack.CDV.1day.1, type='link', re.form = NA)


logRSS.pack.CDV.1day<- p.pack.CDV.1day.1 - p.CDV.1day.2
logRSS.pack.CDV.1day.social<- p.pack.CDV.1day.1.social - p.CDV.1day.2.social
logRSS.pack.CDV.1day.pack<- p.pack.CDV.1day.1.pack - p.CDV.1day.2.pack


pack.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T), 
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = seq(0, 15000, length.out = 150),
                                                            COD = factor('CDV', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID = NA)]
p.pack.CDV.60day.1 <- predict(everyone, newdata = pack.CDV.60day.1, type='link', re.form=NA)
p.pack.CDV.60day.1.social <- predict(everyone.social, newdata = pack.CDV.60day.1, type='link', re.form=NA)
p.pack.CDV.60day.1.pack <- predict(everyone.pack, newdata = pack.CDV.60day.1, type='link', re.form=NA)


logRSS.pack.CDV.60day<- p.pack.CDV.60day.1 - p.CDV.60day.2
logRSS.pack.CDV.60day.social<- p.pack.CDV.60day.1.social - p.CDV.60day.2.social
logRSS.pack.CDV.60day.pack<- p.pack.CDV.60day.1.social - p.CDV.60day.2.pack


logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'CDV', ttd = '1 day', var = 'boundary', rss = logRSS.pack.CDV.1day, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'CDV', ttd = '60 days', var = 'boundary', rss = logRSS.pack.CDV.60day, x = seq(from = 0, to = 15000, length.out = 150)))

# logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '1 day', var = 'boundary', rss = logRSS.pack.CDV.1day.social, x = seq(from = 0, to = 15000, length.out = 150)))
# logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '60 days', var = 'boundary', rss = logRSS.pack.CDV.60day.social, x = seq(from = 0, to = 15000, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '1 day', var = 'boundary', rss = logRSS.pack.CDV.1day.pack, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '60 days', var = 'boundary', rss = logRSS.pack.CDV.60day.pack, x = seq(from = 0, to = 15000, length.out = 150)))

#### indivs ####


p.pack.h1.indiv <- function(ids, DT, mod, death, t2death){
    lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
          propwet_end= mean(propwet_end, na.rm = T), 
          roadDist_end = median(roadDist_end, na.rm = T),
          distance2 = median(distance2, na.rm = T),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = seq(0, 15000, length.out = 150),
          COD = factor(death, levels = levels(COD)),
          ttd1 = t2death,
          wolf_step_id = NA,
          wolfID = i
        )],
        type = "link",
        re.form = NULL
      ), wolfID = i)]
    # )
  })
  }

p.pack.CDV.1day.1.indiv <- p.pack.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)

h1.pack.CDV.1day.indiv <- data.table(rbindlist(p.pack.CDV.1day.1.indiv),
                                   ttd = '1 day', COD = 'CDV', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                   ttd = '1 day', COD = 'CDV', var = 'boundary')



logRSS.pack.CDV.1day.indiv <- merge(h1.pack.CDV.1day.indiv, h2.pack.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.CDV.1day.indiv[,'rss'] <- logRSS.pack.CDV.1day.indiv$h1 - logRSS.pack.CDV.1day.indiv$h2

#### social model

p.pack.CDV.1day.1.indiv.social <- p.pack.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 1)

h1.pack.CDV.1day.indiv.social <- data.table(rbindlist(p.pack.CDV.1day.1.indiv.social),
                                          ttd = '1 day', COD = 'CDV', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.CDV.1day.indiv.social <- data.table(rbindlist(p.CDV.1day.2.indiv.social),
                                          ttd = '1 day', COD = 'CDV', var = 'boundary')


logRSS.pack.CDV.1day.indiv.social <- merge(h1.pack.CDV.1day.indiv.social, h2.pack.CDV.1day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.CDV.1day.indiv.social[,'rss'] <- logRSS.pack.CDV.1day.indiv.social$h1 - logRSS.pack.CDV.1day.indiv.social$h2


#### pack model
p.pack.CDV.1day.1.indiv.pack <- p.pack.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', t2death = 1)

h1.pack.CDV.1day.indiv.pack <- data.table(rbindlist(p.pack.CDV.1day.1.indiv.pack),
                                            ttd = '1 day', COD = 'CDV', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.CDV.1day.indiv.pack <- data.table(rbindlist(p.CDV.1day.2.indiv.pack),
                                            ttd = '1 day', COD = 'CDV', var = 'boundary')


logRSS.pack.CDV.1day.indiv.pack <- merge(h1.pack.CDV.1day.indiv.pack, h2.pack.CDV.1day.indiv.pack, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.CDV.1day.indiv.pack[,'rss'] <- logRSS.pack.CDV.1day.indiv.pack$h1 - logRSS.pack.CDV.1day.indiv.pack$h2


### 60 days

p.pack.CDV.60day.1.indiv <- p.pack.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)


h1.pack.CDV.60day.indiv <- data.table(rbindlist(p.pack.CDV.60day.1.indiv),
                                    ttd = '60 days', COD = 'CDV', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.CDV.60day.indiv <- data.table(rbindlist(p.CDV.60day.2.indiv),
                                    ttd = '60 days', COD = 'CDV', var = 'boundary')



logRSS.pack.CDV.60day.indiv <- merge(h1.pack.CDV.60day.indiv, h2.pack.CDV.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.CDV.60day.indiv[,'rss'] <- logRSS.pack.CDV.60day.indiv$h1 - logRSS.pack.CDV.60day.indiv$h2

#### social model

p.pack.CDV.60day.1.indiv.social <- p.pack.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 60)


h1.pack.CDV.60day.indiv.social <- data.table(rbindlist(p.pack.CDV.60day.1.indiv.social),
                                           ttd = '60 days', COD = 'CDV', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))



h2.pack.CDV.60day.indiv.social <- data.table(rbindlist(p.CDV.60day.2.indiv.social),
                                           ttd = '60 days', COD = 'CDV', var = 'boundary')


logRSS.pack.CDV.60day.indiv.social <- merge(h1.pack.CDV.60day.indiv.social, h2.pack.CDV.60day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.CDV.60day.indiv.social[,'rss'] <- logRSS.pack.CDV.60day.indiv.social$h1 - logRSS.pack.CDV.60day.indiv.social$h2

#### pack model

p.pack.CDV.60day.1.indiv.pack <- p.pack.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', t2death = 60)


h1.pack.CDV.60day.indiv.pack <- data.table(rbindlist(p.pack.CDV.60day.1.indiv.pack),
                                             ttd = '60 days', COD = 'CDV', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))



h2.pack.CDV.60day.indiv.pack <- data.table(rbindlist(p.CDV.60day.2.indiv.pack),
                                             ttd = '60 days', COD = 'CDV', var = 'boundary')


logRSS.pack.CDV.60day.indiv.pack <- merge(h1.pack.CDV.60day.indiv.pack, h2.pack.CDV.60day.indiv.pack, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.CDV.60day.indiv.pack[,'rss'] <- logRSS.pack.CDV.60day.indiv.pack$h1 - logRSS.pack.CDV.60day.indiv.pack$h2

#### human ####

pack.human.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),  
                                                             roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = median(distance2, na.rm = T),
                                                             nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = seq(0, 15000, length.out = 150),
                                                             COD = factor('human', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA, wolfID = NA)]

p.pack.human.1day.1 <- predict(everyone, newdata = pack.human.1day.1, type='link', re.form = NA)
p.pack.human.1day.1.social <- predict(everyone.social, newdata = pack.human.1day.1, type='link', re.form = NA)
p.pack.human.1day.1.pack <- predict(everyone.pack, newdata = pack.human.1day.1, type='link', re.form = NA)


logRSS.pack.human.1day<- p.pack.human.1day.1 - p.human.1day.2
logRSS.pack.human.1day.social<- p.pack.human.1day.1.social - p.human.1day.2.social
logRSS.pack.human.1day.pack<- p.pack.human.1day.1.pack - p.human.1day.2.pack


pack.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                              log_sl = mean(log_sl),
                                                              cos_ta = mean(cos_ta),
                                                              # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                              propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                              propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                              propwet_end= mean(propwet_end, na.rm = T), 
                                                              roadDist_end = median(roadDist_end, na.rm = T),
                                                              distance2 = median(distance2, na.rm = T),
                                                              nnDist_end = median(nnDist_end, na.rm = T),
                                                              packDist_end = seq(0, 15000, length.out = 150),
                                                              COD = factor('human', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]
p.pack.human.60day.1 <- predict(everyone, newdata = pack.human.60day.1, type='link', re.form=NA)
p.pack.human.60day.1.social <- predict(everyone.social, newdata = pack.human.60day.1, type='link', re.form=NA)
p.pack.human.60day.1.pack <- predict(everyone.pack, newdata = pack.human.60day.1, type='link', re.form=NA)



logRSS.pack.human.60day<- p.pack.human.60day.1 - p.human.60day.2
logRSS.pack.human.60day.social<- p.pack.human.60day.1.social - p.human.60day.2.social
logRSS.pack.human.60day.pack<- p.pack.human.60day.1.pack - p.human.60day.2.pack

logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'human', ttd = '1 day', var = 'boundary', rss = logRSS.pack.human.1day, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'human', ttd = '60 days', var = 'boundary', rss = logRSS.pack.human.60day, x = seq(from = 0, to = 15000, length.out = 150)))

# logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '1 day', var = 'boundary', rss = logRSS.pack.human.1day.social, x = seq(from = 0, to = 15000, length.out = 150)))
# logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '60 days', var = 'boundary', rss = logRSS.pack.human.60day.social, x = seq(from = 0, to = 15000, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '1 day', var = 'boundary', rss = logRSS.pack.human.1day.pack, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '60 days', var = 'boundary', rss = logRSS.pack.human.60day.pack, x = seq(from = 0, to = 15000, length.out = 150)))

#### indivs ####

p.pack.human.1day.1.indiv <- p.pack.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)

h1.pack.human.1day.indiv <- data.table(rbindlist(p.pack.human.1day.1.indiv),
                                     ttd = '1 day', COD = 'human', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))

h2.pack.human.1day.indiv <- data.table(rbindlist(p.human.1day.2.indiv),
                                     ttd = '1 day', COD = 'human', var = 'boundary')



logRSS.pack.human.1day.indiv <- merge(h1.pack.human.1day.indiv, h2.pack.human.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.human.1day.indiv[,'rss'] <- logRSS.pack.human.1day.indiv$h1 - logRSS.pack.human.1day.indiv$h2

#### social model

p.pack.human.1day.1.indiv.social <- p.pack.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 1)


h1.pack.human.1day.indiv.social <- data.table(rbindlist(p.pack.human.1day.1.indiv.social),
                                            ttd = '1 day', COD = 'human', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.human.1day.indiv.social <- data.table(rbindlist(p.human.1day.2.indiv.social),
                                            ttd = '1 day', COD = 'human', var = 'boundary')


logRSS.pack.human.1day.indiv.social <- merge(h1.pack.human.1day.indiv.social, h2.pack.human.1day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.human.1day.indiv.social[,'rss'] <- logRSS.pack.human.1day.indiv.social$h1 - logRSS.pack.human.1day.indiv.social$h2


#### pack model

p.pack.human.1day.1.indiv.pack <- p.pack.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', t2death = 1)


h1.pack.human.1day.indiv.pack <- data.table(rbindlist(p.pack.human.1day.1.indiv.pack),
                                              ttd = '1 day', COD = 'human', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.human.1day.indiv.pack <- data.table(rbindlist(p.human.1day.2.indiv.pack),
                                              ttd = '1 day', COD = 'human', var = 'boundary')


logRSS.pack.human.1day.indiv.pack <- merge(h1.pack.human.1day.indiv.pack, h2.pack.human.1day.indiv.pack, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.human.1day.indiv.pack[,'rss'] <- logRSS.pack.human.1day.indiv.pack$h1 - logRSS.pack.human.1day.indiv.pack$h2


### 60 days

p.pack.human.60day.1.indiv <- p.pack.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)


h1.pack.human.60day.indiv <- data.table(rbindlist(p.pack.human.60day.1.indiv),
                                      ttd = '60 days', COD = 'human', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))

h2.pack.human.60day.indiv <- data.table(rbindlist(p.human.60day.2.indiv),
                                      ttd = '60 days', COD = 'human', var = 'boundary')



logRSS.pack.human.60day.indiv <- merge(h1.pack.human.60day.indiv, h2.pack.human.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.human.60day.indiv[,'rss'] <- logRSS.pack.human.60day.indiv$h1 - logRSS.pack.human.60day.indiv$h2

#### social model

p.pack.human.60day.1.indiv.social <- p.pack.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 60)

h1.pack.human.60day.indiv.social <- data.table(rbindlist(p.pack.human.60day.1.indiv.social),
                                             ttd = '60 days', COD = 'human', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.human.60day.indiv.social <- data.table(rbindlist(p.human.60day.2.indiv.social),
                                             ttd = '60 days', COD = 'human', var = 'boundary')


logRSS.pack.human.60day.indiv.social <- merge(h1.pack.human.60day.indiv.social, h2.pack.human.60day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.human.60day.indiv.social[,'rss'] <- logRSS.pack.human.60day.indiv.social$h1 - logRSS.pack.human.60day.indiv.social$h2

#### pack model

p.pack.human.60day.1.indiv.pack <- p.pack.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', t2death = 60)

h1.pack.human.60day.indiv.pack <- data.table(rbindlist(p.pack.human.60day.1.indiv.pack),
                                               ttd = '60 days', COD = 'human', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.human.60day.indiv.pack <- data.table(rbindlist(p.human.60day.2.indiv.pack),
                                               ttd = '60 days', COD = 'human', var = 'boundary')


logRSS.pack.human.60day.indiv.pack <- merge(h1.pack.human.60day.indiv.pack, h2.pack.human.60day.indiv.pack, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.human.60day.indiv.pack[,'rss'] <- logRSS.pack.human.60day.indiv.pack$h1 - logRSS.pack.human.60day.indiv.pack$h2

#### control ####

pack.control.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                               log_sl = mean(log_sl),
                                                               cos_ta = mean(cos_ta),
                                                               #land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                               propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                               propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                               propwet_end= mean(propwet_end, na.rm = T),  
                                                               roadDist_end = median(roadDist_end, na.rm = T),
                                                               distance2 = median(distance2, na.rm = T),
                                                               nnDist_end = median(nnDist_end, na.rm = T),
                                                               packDist_end = seq(0, 15000, length.out = 150),
                                                               COD = factor('control', levels = levels(COD)),
                                                               ttd1 = 1,
                                                               wolf_step_id = NA, wolfID = NA)]

p.pack.control.1day.1 <- predict(everyone, newdata = pack.control.1day.1, type='link', re.form = NA)
p.pack.control.1day.1.social <- predict(everyone.social, newdata = pack.control.1day.1, type='link', re.form = NA)
p.pack.control.1day.1.pack <- predict(everyone.pack, newdata = pack.control.1day.1, type='link', re.form = NA)


logRSS.pack.control.1day<- p.pack.control.1day.1 - p.control.1day.2
logRSS.pack.control.1day.social<- p.pack.control.1day.1.social - p.control.1day.2.social
logRSS.pack.control.1day.pack<- p.pack.control.1day.1.pack - p.control.1day.2.pack


pack.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                                propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                                propwet_end= mean(propwet_end, na.rm = T), 
                                                                roadDist_end = median(roadDist_end, na.rm = T),
                                                                distance2 = median(distance2, na.rm = T),
                                                                nnDist_end = median(nnDist_end, na.rm = T),
                                                                packDist_end = seq(0, 15000, length.out = 150),
                                                                COD = factor('control', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]
p.pack.control.60day.1 <- predict(everyone, newdata = pack.control.60day.1, type='link', re.form=NA)
p.pack.control.60day.1.social <- predict(everyone.social, newdata = pack.control.60day.1, type='link', re.form=NA)
p.pack.control.60day.1.pack <- predict(everyone.pack, newdata = pack.control.60day.1, type='link', re.form=NA)



logRSS.pack.control.60day<- p.pack.control.60day.1 - p.control.60day.2
logRSS.pack.control.60day.social<- p.pack.control.60day.1.social - p.control.60day.2.social
logRSS.pack.control.60day.pack<- p.pack.control.60day.1.pack - p.control.60day.2.pack

logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'control', ttd = '1 day', var = 'boundary', rss = logRSS.pack.control.1day, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'control', ttd = '60 days', var = 'boundary', rss = logRSS.pack.control.60day, x = seq(from = 0, to = 15000, length.out = 150)))

# logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '1 day', var = 'boundary', rss = logRSS.pack.control.1day.social, x = seq(from = 0, to = 15000, length.out = 150)))
# logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '60 days', var = 'boundary', rss = logRSS.pack.control.60day.social, x = seq(from = 0, to = 15000, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '1 day', var = 'boundary', rss = logRSS.pack.control.1day.pack, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '60 days', var = 'boundary', rss = logRSS.pack.control.60day.pack, x = seq(from = 0, to = 15000, length.out = 150)))

#### indivs ####

p.pack.control.1day.1.indiv <- p.pack.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)

h1.pack.control.1day.indiv <- data.table(rbindlist(p.pack.control.1day.1.indiv),
                                       ttd = '1 day', COD = 'control', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.control.1day.indiv <- data.table(rbindlist(p.control.1day.2.indiv),
                                       ttd = '1 day', COD = 'control', var = 'boundary')



logRSS.pack.control.1day.indiv <- merge(h1.pack.control.1day.indiv, h2.pack.control.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.control.1day.indiv[,'rss'] <- logRSS.pack.control.1day.indiv$h1 - logRSS.pack.control.1day.indiv$h2

#### social model

p.pack.control.1day.1.indiv.social <- p.pack.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 1)


h1.pack.control.1day.indiv.social <- data.table(rbindlist(p.pack.control.1day.1.indiv.social),
                                              ttd = '1 day', COD = 'control', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.control.1day.indiv.social <- data.table(rbindlist(p.control.1day.2.indiv.social),
                                              ttd = '1 day', COD = 'control', var = 'boundary')


logRSS.pack.control.1day.indiv.social <- merge(h1.pack.control.1day.indiv.social, h2.pack.control.1day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.control.1day.indiv.social[,'rss'] <- logRSS.pack.control.1day.indiv.social$h1 - logRSS.pack.control.1day.indiv.social$h2


#### pack model

p.pack.control.1day.1.indiv.pack <- p.pack.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', t2death = 1)


h1.pack.control.1day.indiv.pack <- data.table(rbindlist(p.pack.control.1day.1.indiv.pack),
                                                ttd = '1 day', COD = 'control', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))


h2.pack.control.1day.indiv.pack <- data.table(rbindlist(p.control.1day.2.indiv.pack),
                                                ttd = '1 day', COD = 'control', var = 'boundary')


logRSS.pack.control.1day.indiv.pack <- merge(h1.pack.control.1day.indiv.pack, h2.pack.control.1day.indiv.pack, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.control.1day.indiv.pack[,'rss'] <- logRSS.pack.control.1day.indiv.pack$h1 - logRSS.pack.control.1day.indiv.pack$h2


### 60 days

p.pack.control.60day.1.indiv <- p.pack.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)


h1.pack.control.60day.indiv <- data.table(rbindlist(p.pack.control.60day.1.indiv),
                                        ttd = '60 days', COD = 'control', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))



h2.pack.control.60day.indiv <- data.table(rbindlist(p.control.60day.2.indiv),
                                        ttd = '60 days', COD = 'control', var = 'boundary')



logRSS.pack.control.60day.indiv <- merge(h1.pack.control.60day.indiv, h2.pack.control.60day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.control.60day.indiv[,'rss'] <- logRSS.pack.control.60day.indiv$h1 - logRSS.pack.control.60day.indiv$h2

#### social model

p.pack.control.60day.1.indiv.social <- p.pack.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 60)

h1.pack.control.60day.indiv.social <- data.table(rbindlist(p.pack.control.60day.1.indiv.social),
                                               ttd = '60 days', COD = 'control', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))



h2.pack.control.60day.indiv.social <- data.table(rbindlist(p.control.60day.2.indiv.social),
                                               ttd = '60 days', COD = 'control', var = 'boundary')


logRSS.pack.control.60day.indiv.social <- merge(h1.pack.control.60day.indiv.social, h2.pack.control.60day.indiv.social, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.control.60day.indiv.social[,'rss'] <- logRSS.pack.control.60day.indiv.social$h1 - logRSS.pack.control.60day.indiv.social$h2


#### pack model

p.pack.control.60day.1.indiv.pack <- p.pack.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', t2death = 60)

h1.pack.control.60day.indiv.pack <- data.table(rbindlist(p.pack.control.60day.1.indiv.pack),
                                                 ttd = '60 days', COD = 'control', var = 'boundary', x = seq(from = 0, to = 15000, length.out = 150))



h2.pack.control.60day.indiv.pack <- data.table(rbindlist(p.control.60day.2.indiv.pack),
                                                 ttd = '60 days', COD = 'control', var = 'boundary')


logRSS.pack.control.60day.indiv.pack <- merge(h1.pack.control.60day.indiv.pack, h2.pack.control.60day.indiv.pack, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.pack.control.60day.indiv.pack[,'rss'] <- logRSS.pack.control.60day.indiv.pack$h1 - logRSS.pack.control.60day.indiv.pack$h2


##### gathering RSS
logRSS.pack.indiv <- rbind(logRSS.pack.control.1day.indiv, logRSS.pack.control.60day.indiv, 
                         logRSS.pack.human.1day.indiv, logRSS.pack.human.60day.indiv, 
                         logRSS.pack.CDV.1day.indiv, logRSS.pack.CDV.60day.indiv)

logRSS.pack.indiv.social <- rbind(logRSS.pack.control.1day.indiv.social, logRSS.pack.control.60day.indiv.social, 
                                logRSS.pack.human.1day.indiv.social, logRSS.pack.human.60day.indiv.social, 
                                logRSS.pack.CDV.1day.indiv.social, logRSS.pack.CDV.60day.indiv.social)

logRSS.pack.indiv.pack <- rbind(logRSS.pack.control.1day.indiv.pack, logRSS.pack.control.60day.indiv.pack, 
                                  logRSS.pack.human.1day.indiv.pack, logRSS.pack.human.60day.indiv.pack, 
                                  logRSS.pack.CDV.1day.indiv.pack, logRSS.pack.CDV.60day.indiv.pack)


#### all ###

logRSS.indiv <- rbind(logRSS.forest.indiv, logRSS.open.indiv, logRSS.wet.indiv, logRSS.road.indiv, 
                      logRSS.nn.indiv, logRSS.pack.indiv)
# logRSS.indiv.model <- rbind(logRSS.forest.indiv.habitat, logRSS.open.indiv.habitat, logRSS.wet.indiv.habitat, logRSS.road.indiv.habitat, 
#                       logRSS.nn.indiv.social, logRSS.pack.indiv.social)
logRSS.indiv.model <- rbind(logRSS.forest.indiv.habitat, logRSS.open.indiv.habitat, logRSS.wet.indiv.habitat, logRSS.road.indiv.habitat, 
                            logRSS.nn.indiv.social, logRSS.pack.indiv.pack)

logRSS <- rbind(logRSS.forest, logRSS.open, logRSS.wet, logRSS.road, logRSS.nn)
logRSS.model <- rbind(logRSS.forest.habitat, logRSS.open.habitat, logRSS.wet.habitat, logRSS.road.habitat, logRSS.social)

#saveRDS(logRSS.indiv, 'data/derived-data/logRSS_indiv2.Rds')
#saveRDS(logRSS.indiv.model, 'data/derived-data/logRSS_indiv_model2.Rds')


#saveRDS(logRSS, 'data/derived-data/logRSS2.Rds')
#saveRDS(logRSS.model, 'data/derived-data/logRSS_model2.Rds')


logRSS.indiv.cor <- dcast(logRSS.indiv.model, wolfID + COD + var + x ~ ttd, value.var = 'rss')
logRSS.indiv.cor[,cor.test(`1 day`,`60 days`), by= .(COD, var)]

#### GRAPHS ####
logRSS.indiv <- readRDS('data/derived-data/logRSS_indiv2.Rds')
logRSS.indiv.model <- readRDS('data/derived-data/logRSS_indiv_model2.Rds')


logRSS <- readRDS('data/derived-data/logRSS2.Rds')
logRSS.model <- readRDS('data/derived-data/logRSS_model2.Rds')

logRSS.indiv <- setDT(logRSS.indiv)
logRSS.indiv[,'COD'] <- factor(logRSS.indiv$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))
logRSS.indiv[,'ttd'] <- as.factor(logRSS.indiv$ttd)

logRSS.indiv.se <- unique(logRSS.indiv[,.(se=se(rss)), by = .(x, var, COD, ttd)])

logRSS.pop <- merge(setDT(logRSS), logRSS.indiv.se, by = c('x', 'COD', 'ttd','var'))
logRSS.pop[,'COD'] <- factor(logRSS.pop$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))


logRSS.indiv.model <- setDT(logRSS.indiv.model)
logRSS.indiv.model[,'COD'] <- factor(logRSS.indiv.model$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))
logRSS.indiv.model[,'ttd'] <- as.factor(logRSS.indiv.model$ttd)

logRSS.indiv.model.se <- unique(logRSS.indiv.model[,.(se=se(rss)), by = .(x, var, COD, ttd)])

logRSS.pop.model <- merge(setDT(logRSS.model), logRSS.indiv.model.se, by = c('x', 'COD', 'ttd','var'))
logRSS.pop.model[,'COD'] <- factor(logRSS.pop.model$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))
logRSS.pop.model[,'ttd'] <- as.factor(logRSS.pop.model$ttd)


#### TIMELINE ####

ttd.days <-data.frame(time = c(-60, -1,0), days = c("-60 days", "-1 day",""), end = c('','','0 days'))
ttd <- ggplot(ttd.days, aes(time, 0))+
  geom_point() +
  geom_line(size=1) +
  geom_text(aes(x = time, y = 0.01, label = days))+
  geom_text(aes(x = time, y = -0.01, label = end))+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
    axis.line.x = element_blank()) +
  theme(plot.title=element_text(size=15,hjust = 0.05),axis.text.x = element_blank(), axis.title = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())+ theme(axis.ticks = element_blank()) +
   ylim(-0.015, 0.015)+
 ggtitle("Time to death") + theme(plot.title = element_text(hjust = 0.5))


#### FOREST not sig ####
forest.1 <- ggplot(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Proportion of forest") +
  ggtitle("1 day")  +
 # ylim(-0.7,1.3) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


forest.60 <- ggplot(data=logRSS.pop[var == 'forest'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Proportion of forest") +
  ggtitle("60 days")  +
 # ylim(-0.7,1.3) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

forest.1|forest.60


 ggplot(data=setDT(logRSS.indiv)[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_smooth() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion of forest") +
  ggtitle("1 day") 


forest.1.b <- ggplot(data=setDT(logRSS.indiv)[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion forest") +
  ggtitle("1 day") + theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

forest.60.b <- ggplot(data=setDT(logRSS.indiv)[var == 'forest'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm', show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion forest") +
  ggtitle("60 days") + theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

forest.60.b|forest.1.b

#### OPEN sig ####
open.1 <- ggplot(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .1), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Proportion of open") +
  ggtitle("1 day")  +
  ylim(-3,1) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


open.60 <- ggplot(data=logRSS.pop[var == 'open'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .1), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Proportion of open") +
  ggtitle("60 days")  +
  ylim(-3,1) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

open.1|open.60

open.1.b <- ggplot(data=setDT(logRSS.indiv)[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion open") +
 # ggtitle("b) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

open.60.b <- ggplot(data=setDT(logRSS.indiv)[var == 'open'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm', show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion open") +
  #ggtitle("a) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

open.60.b|open.1.b


#### WET not sig ####
wet.1 <- ggplot(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Proportion of wet") +
  ggtitle("1 day")  +
  ylim(-1.5,1.7) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


wet.60 <- ggplot(data=logRSS.pop[var == 'wet'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Proportion of wet") +
  ggtitle("60 days")  +
  ylim(-1.5,1.7) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

wet.1|wet.60

wet.1.b <- ggplot(data=setDT(logRSS.indiv)[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion wet") +
  #ggtitle("b) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

wet.60.b <- ggplot(data=setDT(logRSS.indiv)[var == 'wet'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm', show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion wet") +
  #ggtitle("a) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10)+
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

#### ALL HAB GRAPH ####
(forest.60.b|forest.1.b)/(open.60.b|open.1.b)/(wet.60.b|wet.1.b)


# ggplot(data=setDT(logRSS.indiv)[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
#   geom_smooth() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
#   #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
#   ylab("logRSS") + xlab("Proportion of forest") +
#   ggtitle("1 day") 



#### ROAD control sig ####
road.1 <- ggplot(data=logRSS.pop[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("1 day")  +
  ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


road.60 <- ggplot(data=logRSS.pop[var == 'road'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("60 days")  +
  ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

road.1|road.60

ggplot(data=setDT(logRSS.indiv)[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_smooth(method = 'loess') +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Dist to road") +
  ggtitle("1 day") 

ggplot(data=setDT(logRSS.indiv)[var == 'road'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_smooth(method = 'loess') +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Dist to road") +
  ggtitle("60 days") 

road.1.model <- ggplot(data=logRSS.pop.model[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("1 day")  +
  ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


road.60.model <- ggplot(data=logRSS.pop.model[var == 'road'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("60 days")  +
  ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

road.1.model|road.60.model


road.1.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'road'& ttd=='1 day'], aes(x, -rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='dashed', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("b) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-2,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

road.60.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'road'& ttd=='60 days'], aes(x, -rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='dashed', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("a) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-2,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

road.60.b|road.1.b

#### NN  sig ####
nn.1 <- ggplot(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to nn (m)") +
  ggtitle("1 day")  +
 # ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


nn.60 <- ggplot(data=logRSS.pop[var == 'nn'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to nn (m)") +
  ggtitle("60 days")  +
 # ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

nn.1|nn.60


ggplot(data=setDT(logRSS.indiv)[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_smooth(method = 'loess') +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Dist to nn") +
  ggtitle("1 day") 

ggplot(data=setDT(logRSS.indiv)[var == 'nn'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_smooth(method = 'loess') +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Dist to nn") +
  ggtitle("60 days") 

nn.1.model <- ggplot(data=logRSS.pop.model[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to nn (m)") +
  ggtitle("1 day")  +
  ylim(-1,8) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


nn.60.model <- ggplot(data=logRSS.pop.model[var == 'nn'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to nn (m)") +
  ggtitle("60 days")  +
  ylim(-1,8) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

nn.60.model|nn.1.model


nn.1.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
 # geom_line(data=logRSS.pop.model[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to nearest neighbor (m)") +
 ggtitle("1 day") +
  theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-.5,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

nn.60.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'nn'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess', show.legend = F) +
  # geom_line(data=logRSS.pop.model[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to nearest neighbor (m)") +
 ggtitle("60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-.5,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))



(nn.60.model.b|nn.1.model.b)





#### PACK  sig ####
pack.1 <- ggplot(data=logRSS.pop[var == 'boundary'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to pack (m)") +
  ggtitle("1 day")  +
  # ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


pack.60 <- ggplot(data=logRSS.pop[var == 'boundary'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to pack (m)") +
  ggtitle("60 days")  +
  # ylim(-2,12) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

pack.1|pack.60


ggplot(data=setDT(logRSS.indiv)[var == 'boundary'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_smooth(method = 'loess') +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Dist to pack") +
  ggtitle("1 day") 

ggplot(data=setDT(logRSS.indiv)[var == 'boundary'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_smooth(method = 'loess') +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Dist to pack") +
  ggtitle("60 days") 

pack.1.model <- ggplot(data=logRSS.pop.model[var == 'boundary'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to pack (m)") +
  ggtitle("1 day")  +
  #ylim(-1,8) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))


pack.60.model <- ggplot(data=logRSS.pop.model[var == 'boundary'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2), show.legend = F)+
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("logRSS") + xlab("Distance to pack (m)") +
  ggtitle("60 days")  +
 # ylim(-1,8) +
  # scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

pack.1.model|pack.60.model


pack.1.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'boundary'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='dashed', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'boundary'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to pack boundary (m)") +
  ggtitle("b) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
 
  # ylim(-.5,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

pack.60.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'boundary'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='dashed', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'boundary'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to pack boundary (m)") +
  ggtitle("a) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-.5,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

pack.60.model.b|pack.1.model.b

