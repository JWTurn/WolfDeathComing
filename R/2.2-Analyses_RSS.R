### RSS ====
# Julie Turner
# Revised July 6 2020


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
dat.RMNP <- readRDS('data/derived-data/ssfAllCov_RMNP.Rds')
dat.GHA26 <- readRDS('data/derived-data/ssfAllCov_GHA26.Rds')


dat.RMNP[,'pop'] <- 'RMNP'
dat.RMNP$wolfID <- paste(dat.RMNP$pop, dat.RMNP$id, sep = '_')
dat.GHA26[,'pop'] <- 'GHA26'
dat.GHA26$wolfID <- paste(dat.GHA26$pop, dat.GHA26$id, sep = '_')
dat<-rbind(dat.RMNP, dat.GHA26, fill=T)

dat[ttd1>=0 & ttd2>=0, unique(wolfID)]

dat<-dat[ttd1>=0 & ttd2>=0]

#dat[,'wtd1'] <- as.integer(dat$ttd1/7)

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'))
dat.meta[,'wolfpop'] <- paste(dat.meta$pop, dat.meta$WolfID, sep = '_')

params <- readRDS('data/derived-data/moveParams_all.Rds')

dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')

dat<- merge(dat, dat.meta, by.x = c('id', 'pop'), by.y = c('WolfID', 'pop'))

dat$packDist_end <- ifelse(dat$packDistadj_end >=0, dat$packDistadj_end, 0)


dat[,'propforest_end_adj'] <- dat$propconif_end+dat$propmixed_end +dat$propdecid_end +dat$propshrub_end # remove shrub?
dat[,'propopen_end_adj'] <- dat$propopen_end + dat$propurban_end 

dat[,'wolf_step_id'] <- paste(dat$wolfID, dat$step_id_, sep = '_')

dat$COD <- factor(dat$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
dat$COD[is.na(dat$COD)] <- "control"

dat$ToD_start <- as.factor(dat$ToD_start)


#### only wolves with packmates ####
dat[,uniqueN(step_id_), by=.(wolfID)]
dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo <- dat[ttd1<=31 & !is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo.cod <- merge(dat.wnn.lastmo, dat.meta[,.(wolfpop, pop, COD)], by.x = 'wolfID', by.y = 'wolfpop', all.x = T)

dat.wnn.lastmo.cod[,.(N=uniqueN(wolfID)), by=.(pop, COD)]
dat.meta[,.(N=uniqueN(wolfpop)), by=.(pop, COD)]

dat[,'nnDist_end'] <- ifelse(dat$distance2<=30000, dat$distance2, NA)
dat[,'packDist_end_5'] <- ifelse(dat$packDist_end<=50000, dat$packDist_end, NA)

#### everyone ####

everyone <- glmmTMB(case_ ~# pop + 
                      log_sl:ToD_start +
                      # sl_ + 
                      #log_sl:cos_ta +
                      # log_sl:land_end_adj +
                      log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                      log_sl:COD + cos_ta:COD + 
                      I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                      (1|wolf_step_id) +
                      (0 + (log_sl)|wolfID) +
                      # (0 + (sl_)|wolfID) +
                      #   (0 + (log_sl:cos_ta)|wolfID) +
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
                    map = list(theta=factor(c(NA,1:16))), start = list(theta=c(log(1000),seq(0,0, length.out = 16))))

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
everyone.move <- glmmTMB(case_ ~ #log_sl:ToD_start +
                           log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                           #log_sl:land_end_adj +
                           log_sl:COD + cos_ta:COD +
                           # log_sl:cos_ta +
                           I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                           (1|wolf_step_id) +
                           (0 + (log_sl)|wolfID) +
                           (0 + (cos_ta)|wolfID) +
                           #   (0 + (log_sl:cos_ta)|wolfID) +
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



everyone.habitat <- glmmTMB(case_ ~# log_sl:ToD_start +
                              log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                              # log_sl:COD + cos_ta:COD + 
                              # I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                              #log_sl:cos_ta +
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

everyone.social <- glmmTMB(case_ ~# log_sl:ToD_start +
                             log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                             # log_sl:cos_ta +
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


everyone.pack <- glmmTMB(case_ ~ #log_sl:ToD_start +
                           log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                           # log_sl:COD + cos_ta:COD + 
                           # I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                           #log_sl:cos_ta +
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





#### RSS ####
#### h2 RSS ####
### CDV ####
# forest
CDV.forest75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                        log_sl = mean(log_sl),
                                                        cos_ta = mean(cos_ta),
                                                        # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                        propforest_end_adj = 0.75,
                                                        propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                        propwet_end= mean(propwet_end, na.rm = T),      
                                                        roadDist_end = median(roadDist_end, na.rm = T),
                                                        distance2 = median(distance2, na.rm = T),
                                                        # nnDist_end = median(nnDist_end, na.rm = T),
                                                        packDist_end = median(packDist_end, na.rm = T),
                                                        COD = factor('CDV', levels = levels(COD)),
                                                        ttd1 = 60,
                                                        wolf_step_id = NA, wolfID= NA)]
p.CDV.forest75.2 <- predict(everyone, newdata = CDV.forest75.2, type='link', re.form = NA)
p.CDV.forest75.2.habitat <- predict(everyone.habitat, newdata = CDV.forest75.2, type='link', re.form = NA)
p.CDV.forest75.2.social <- predict(everyone.social, newdata = CDV.forest75.2, type='link', re.form = NA)
p.CDV.forest75.2.pack <- predict(everyone.pack, newdata = CDV.forest75.2, type='link', re.form = NA)



CDV.forest25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = 0.25,
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= mean(propwet_end, na.rm = T),      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('CDV', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.CDV.forest25.2 <- predict(everyone, newdata = CDV.forest25.2, type='link', re.form = NA)
p.CDV.forest25.2.habitat <- predict(everyone.habitat, newdata = CDV.forest25.2, type='link', re.form = NA)
p.CDV.forest25.2.social <- predict(everyone.social, newdata = CDV.forest25.2, type='link', re.form = NA)
p.CDV.forest25.2.pack <- predict(everyone.pack, newdata = CDV.forest25.2, type='link', re.form = NA)


# open
CDV.open75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                            propopen_end_adj = 0.75,
                                                            propwet_end= mean(propwet_end, na.rm = T),      
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            # nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('CDV', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID= NA)]
p.CDV.open75.2 <- predict(everyone, newdata = CDV.open75.2, type='link', re.form = NA)
p.CDV.open75.2.habitat <- predict(everyone.habitat, newdata = CDV.open75.2, type='link', re.form = NA)
p.CDV.open75.2.social <- predict(everyone.social, newdata = CDV.open75.2, type='link', re.form = NA)
p.CDV.open75.2.pack <- predict(everyone.pack, newdata = CDV.open75.2, type='link', re.form = NA)



CDV.open25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                            propopen_end_adj = 0.25,
                                                            propwet_end= mean(propwet_end, na.rm = T),      
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('CDV', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID= NA)]
p.CDV.open25.2 <- predict(everyone, newdata = CDV.open25.2, type='link', re.form = NA)
p.CDV.open25.2.habitat <- predict(everyone.habitat, newdata = CDV.open25.2, type='link', re.form = NA)
p.CDV.open25.2.social <- predict(everyone.social, newdata = CDV.open25.2, type='link', re.form = NA)
p.CDV.open25.2.pack <- predict(everyone.pack, newdata = CDV.open25.2, type='link', re.form = NA)



# wet
CDV.wet75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                          propwet_end= 0.75,      
                                                          roadDist_end = median(roadDist_end, na.rm = T),
                                                          distance2 = median(distance2, na.rm = T),
                                                          # nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('CDV', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.CDV.wet75.2 <- predict(everyone, newdata = CDV.wet75.2, type='link', re.form = NA)
p.CDV.wet75.2.habitat <- predict(everyone.habitat, newdata = CDV.wet75.2, type='link', re.form = NA)
p.CDV.wet75.2.social <- predict(everyone.social, newdata = CDV.wet75.2, type='link', re.form = NA)
p.CDV.wet75.2.pack <- predict(everyone.pack, newdata = CDV.wet75.2, type='link', re.form = NA)



CDV.wet25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj =mean(propopen_end_adj, na.rm = T),
                                                          propwet_end= 0.25,      
                                                          roadDist_end = median(roadDist_end, na.rm = T),
                                                          distance2 = median(distance2, na.rm = T),
                                                          nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('CDV', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.CDV.wet25.2 <- predict(everyone, newdata = CDV.wet25.2, type='link', re.form = NA)
p.CDV.wet25.2.habitat <- predict(everyone.habitat, newdata = CDV.wet25.2, type='link', re.form = NA)
p.CDV.wet25.2.social <- predict(everyone.social, newdata = CDV.wet25.2, type='link', re.form = NA)
p.CDV.wet25.2.pack <- predict(everyone.pack, newdata = CDV.wet25.2, type='link', re.form = NA)


# road
CDV.roadclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                          propwet_end= mean(propwet_end, na.rm = T),      
                                                          roadDist_end = 250,
                                                          distance2 = median(distance2, na.rm = T),
                                                          # nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('CDV', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.CDV.roadclose.2 <- predict(everyone, newdata = CDV.roadclose.2, type='link', re.form = NA)
p.CDV.roadclose.2.habitat <- predict(everyone.habitat, newdata = CDV.roadclose.2, type='link', re.form = NA)
p.CDV.roadclose.2.social <- predict(everyone.social, newdata = CDV.roadclose.2, type='link', re.form = NA)
p.CDV.roadclose.2.pack <- predict(everyone.pack, newdata = CDV.roadclose.2, type='link', re.form = NA)



CDV.roadfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                          propwet_end= mean(propwet_end, na.rm = T),      
                                                          roadDist_end = 1000,
                                                          distance2 = median(distance2, na.rm = T),
                                                          nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('CDV', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.CDV.roadfar.2 <- predict(everyone, newdata = CDV.roadfar.2, type='link', re.form = NA)
p.CDV.roadfar.2.habitat <- predict(everyone.habitat, newdata = CDV.roadfar.2, type='link', re.form = NA)
p.CDV.roadfar.2.social <- predict(everyone.social, newdata = CDV.roadfar.2, type='link', re.form = NA)
p.CDV.roadfar.2.pack <- predict(everyone.pack, newdata = CDV.roadfar.2, type='link', re.form = NA)


# NN
CDV.nnclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),      
                                                             roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = 250,
                                                             # nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID= NA)]
p.CDV.nnclose.2 <- predict(everyone, newdata = CDV.nnclose.2, type='link', re.form = NA)
p.CDV.nnclose.2.habitat <- predict(everyone.habitat, newdata = CDV.nnclose.2, type='link', re.form = NA)
p.CDV.nnclose.2.social <- predict(everyone.social, newdata = CDV.nnclose.2, type='link', re.form = NA)
p.CDV.nnclose.2.pack <- predict(everyone.pack, newdata = CDV.nnclose.2, type='link', re.form = NA)



CDV.nnfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = 1000,
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('CDV', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.CDV.nnfar.2 <- predict(everyone, newdata = CDV.nnfar.2, type='link', re.form = NA)
p.CDV.nnfar.2.habitat <- predict(everyone.habitat, newdata = CDV.nnfar.2, type='link', re.form = NA)
p.CDV.nnfar.2.social <- predict(everyone.social, newdata = CDV.nnfar.2, type='link', re.form = NA)
p.CDV.nnfar.2.pack <- predict(everyone.pack, newdata = CDV.nnfar.2, type='link', re.form = NA)



# pack
CDV.packclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = median(distance2, na.rm = T),
                                                           # nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = 250,
                                                           COD = factor('CDV', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.CDV.packclose.2 <- predict(everyone, newdata = CDV.packclose.2, type='link', re.form = NA)
p.CDV.packclose.2.habitat <- predict(everyone.habitat, newdata = CDV.packclose.2, type='link', re.form = NA)
p.CDV.packclose.2.social <- predict(everyone.social, newdata = CDV.packclose.2, type='link', re.form = NA)
p.CDV.packclose.2.pack <- predict(everyone.pack, newdata = CDV.packclose.2, type='link', re.form = NA)



CDV.packfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= mean(propwet_end, na.rm = T),      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = 1000,
                                                         COD = factor('CDV', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.CDV.packfar.2 <- predict(everyone, newdata = CDV.packfar.2, type='link', re.form = NA)
p.CDV.packfar.2.habitat <- predict(everyone.habitat, newdata = CDV.packfar.2, type='link', re.form = NA)
p.CDV.packfar.2.social <- predict(everyone.social, newdata = CDV.packfar.2, type='link', re.form = NA)
p.CDV.packfar.2.pack <- predict(everyone.pack, newdata = CDV.packfar.2, type='link', re.form = NA)



#### human ####
# forest
human.forest75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = 0.75,
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T),      
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            # nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('human', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID= NA)]
p.human.forest75.2 <- predict(everyone, newdata = human.forest75.2, type='link', re.form = NA)
p.human.forest75.2.habitat <- predict(everyone.habitat, newdata = human.forest75.2, type='link', re.form = NA)
p.human.forest75.2.social <- predict(everyone.social, newdata = human.forest75.2, type='link', re.form = NA)
p.human.forest75.2.pack <- predict(everyone.pack, newdata = human.forest75.2, type='link', re.form = NA)



human.forest25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = 0.25,
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T),      
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('human', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID= NA)]
p.human.forest25.2 <- predict(everyone, newdata = human.forest25.2, type='link', re.form = NA)
p.human.forest25.2.habitat <- predict(everyone.habitat, newdata = human.forest25.2, type='link', re.form = NA)
p.human.forest25.2.social <- predict(everyone.social, newdata = human.forest25.2, type='link', re.form = NA)
p.human.forest25.2.pack <- predict(everyone.pack, newdata = human.forest25.2, type='link', re.form = NA)


# open
human.open75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = 0.75,
                                                          propwet_end= mean(propwet_end, na.rm = T),      
                                                          roadDist_end = median(roadDist_end, na.rm = T),
                                                          distance2 = median(distance2, na.rm = T),
                                                          # nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('human', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.human.open75.2 <- predict(everyone, newdata = human.open75.2, type='link', re.form = NA)
p.human.open75.2.habitat <- predict(everyone.habitat, newdata = human.open75.2, type='link', re.form = NA)
p.human.open75.2.social <- predict(everyone.social, newdata = human.open75.2, type='link', re.form = NA)
p.human.open75.2.pack <- predict(everyone.pack, newdata = human.open75.2, type='link', re.form = NA)



human.open25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = 0.25,
                                                          propwet_end= mean(propwet_end, na.rm = T),      
                                                          roadDist_end = median(roadDist_end, na.rm = T),
                                                          distance2 = median(distance2, na.rm = T),
                                                          nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('human', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.human.open25.2 <- predict(everyone, newdata = human.open25.2, type='link', re.form = NA)
p.human.open25.2.habitat <- predict(everyone.habitat, newdata = human.open25.2, type='link', re.form = NA)
p.human.open25.2.social <- predict(everyone.social, newdata = human.open25.2, type='link', re.form = NA)
p.human.open25.2.pack <- predict(everyone.pack, newdata = human.open25.2, type='link', re.form = NA)



# wet
human.wet75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= 0.75,      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         # nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('human', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.human.wet75.2 <- predict(everyone, newdata = human.wet75.2, type='link', re.form = NA)
p.human.wet75.2.habitat <- predict(everyone.habitat, newdata = human.wet75.2, type='link', re.form = NA)
p.human.wet75.2.social <- predict(everyone.social, newdata = human.wet75.2, type='link', re.form = NA)
p.human.wet75.2.pack <- predict(everyone.pack, newdata = human.wet75.2, type='link', re.form = NA)



human.wet25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj =mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= 0.25,      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('human', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.human.wet25.2 <- predict(everyone, newdata = human.wet25.2, type='link', re.form = NA)
p.human.wet25.2.habitat <- predict(everyone.habitat, newdata = human.wet25.2, type='link', re.form = NA)
p.human.wet25.2.social <- predict(everyone.social, newdata = human.wet25.2, type='link', re.form = NA)
p.human.wet25.2.pack <- predict(everyone.pack, newdata = human.wet25.2, type='link', re.form = NA)


# road
human.roadclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),      
                                                             roadDist_end = 250,
                                                             distance2 = median(distance2, na.rm = T),
                                                             # nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('human', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID= NA)]
p.human.roadclose.2 <- predict(everyone, newdata = human.roadclose.2, type='link', re.form = NA)
p.human.roadclose.2.habitat <- predict(everyone.habitat, newdata = human.roadclose.2, type='link', re.form = NA)
p.human.roadclose.2.social <- predict(everyone.social, newdata = human.roadclose.2, type='link', re.form = NA)
p.human.roadclose.2.pack <- predict(everyone.pack, newdata = human.roadclose.2, type='link', re.form = NA)



human.roadfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = 1000,
                                                           distance2 = median(distance2, na.rm = T),
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('human', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.human.roadfar.2 <- predict(everyone, newdata = human.roadfar.2, type='link', re.form = NA)
p.human.roadfar.2.habitat <- predict(everyone.habitat, newdata = human.roadfar.2, type='link', re.form = NA)
p.human.roadfar.2.social <- predict(everyone.social, newdata = human.roadfar.2, type='link', re.form = NA)
p.human.roadfar.2.pack <- predict(everyone.pack, newdata = human.roadfar.2, type='link', re.form = NA)


# NN
human.nnclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = 250,
                                                           # nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('human', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.human.nnclose.2 <- predict(everyone, newdata = human.nnclose.2, type='link', re.form = NA)
p.human.nnclose.2.habitat <- predict(everyone.habitat, newdata = human.nnclose.2, type='link', re.form = NA)
p.human.nnclose.2.social <- predict(everyone.social, newdata = human.nnclose.2, type='link', re.form = NA)
p.human.nnclose.2.pack <- predict(everyone.pack, newdata = human.nnclose.2, type='link', re.form = NA)



human.nnfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= mean(propwet_end, na.rm = T),      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = 1000,
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('human', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.human.nnfar.2 <- predict(everyone, newdata = human.nnfar.2, type='link', re.form = NA)
p.human.nnfar.2.habitat <- predict(everyone.habitat, newdata = human.nnfar.2, type='link', re.form = NA)
p.human.nnfar.2.social <- predict(everyone.social, newdata = human.nnfar.2, type='link', re.form = NA)
p.human.nnfar.2.pack <- predict(everyone.pack, newdata = human.nnfar.2, type='link', re.form = NA)



# pack
human.packclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),      
                                                             roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = median(distance2, na.rm = T),
                                                             # nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = 250,
                                                             COD = factor('human', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID= NA)]
p.human.packclose.2 <- predict(everyone, newdata = human.packclose.2, type='link', re.form = NA)
p.human.packclose.2.habitat <- predict(everyone.habitat, newdata = human.packclose.2, type='link', re.form = NA)
p.human.packclose.2.social <- predict(everyone.social, newdata = human.packclose.2, type='link', re.form = NA)
p.human.packclose.2.pack <- predict(everyone.pack, newdata = human.packclose.2, type='link', re.form = NA)



human.packfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = median(distance2, na.rm = T),
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = 1000,
                                                           COD = factor('human', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.human.packfar.2 <- predict(everyone, newdata = human.packfar.2, type='link', re.form = NA)
p.human.packfar.2.habitat <- predict(everyone.habitat, newdata = human.packfar.2, type='link', re.form = NA)
p.human.packfar.2.social <- predict(everyone.social, newdata = human.packfar.2, type='link', re.form = NA)
p.human.packfar.2.pack <- predict(everyone.pack, newdata = human.packfar.2, type='link', re.form = NA)


### control ####
# forest
control.forest75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = 0.75,
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T),      
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            # nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('control', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID= NA)]
p.control.forest75.2 <- predict(everyone, newdata = control.forest75.2, type='link', re.form = NA)
p.control.forest75.2.habitat <- predict(everyone.habitat, newdata = control.forest75.2, type='link', re.form = NA)
p.control.forest75.2.social <- predict(everyone.social, newdata = control.forest75.2, type='link', re.form = NA)
p.control.forest75.2.pack <- predict(everyone.pack, newdata = control.forest75.2, type='link', re.form = NA)



control.forest25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                            log_sl = mean(log_sl),
                                                            cos_ta = mean(cos_ta),
                                                            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                            propforest_end_adj = 0.25,
                                                            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                            propwet_end= mean(propwet_end, na.rm = T),      
                                                            roadDist_end = median(roadDist_end, na.rm = T),
                                                            distance2 = median(distance2, na.rm = T),
                                                            nnDist_end = median(nnDist_end, na.rm = T),
                                                            packDist_end = median(packDist_end, na.rm = T),
                                                            COD = factor('control', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID= NA)]
p.control.forest25.2 <- predict(everyone, newdata = control.forest25.2, type='link', re.form = NA)
p.control.forest25.2.habitat <- predict(everyone.habitat, newdata = control.forest25.2, type='link', re.form = NA)
p.control.forest25.2.social <- predict(everyone.social, newdata = control.forest25.2, type='link', re.form = NA)
p.control.forest25.2.pack <- predict(everyone.pack, newdata = control.forest25.2, type='link', re.form = NA)


# open
control.open75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = 0.75,
                                                          propwet_end= mean(propwet_end, na.rm = T),      
                                                          roadDist_end = median(roadDist_end, na.rm = T),
                                                          distance2 = median(distance2, na.rm = T),
                                                          # nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('control', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.control.open75.2 <- predict(everyone, newdata = control.open75.2, type='link', re.form = NA)
p.control.open75.2.habitat <- predict(everyone.habitat, newdata = control.open75.2, type='link', re.form = NA)
p.control.open75.2.social <- predict(everyone.social, newdata = control.open75.2, type='link', re.form = NA)
p.control.open75.2.pack <- predict(everyone.pack, newdata = control.open75.2, type='link', re.form = NA)



control.open25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                          log_sl = mean(log_sl),
                                                          cos_ta = mean(cos_ta),
                                                          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                          propopen_end_adj = 0.25,
                                                          propwet_end= mean(propwet_end, na.rm = T),      
                                                          roadDist_end = median(roadDist_end, na.rm = T),
                                                          distance2 = median(distance2, na.rm = T),
                                                          nnDist_end = median(nnDist_end, na.rm = T),
                                                          packDist_end = median(packDist_end, na.rm = T),
                                                          COD = factor('control', levels = levels(COD)),
                                                          ttd1 = 60,
                                                          wolf_step_id = NA, wolfID= NA)]
p.control.open25.2 <- predict(everyone, newdata = control.open25.2, type='link', re.form = NA)
p.control.open25.2.habitat <- predict(everyone.habitat, newdata = control.open25.2, type='link', re.form = NA)
p.control.open25.2.social <- predict(everyone.social, newdata = control.open25.2, type='link', re.form = NA)
p.control.open25.2.pack <- predict(everyone.pack, newdata = control.open25.2, type='link', re.form = NA)



# wet
control.wet75.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= 0.75,      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         # nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('control', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.control.wet75.2 <- predict(everyone, newdata = control.wet75.2, type='link', re.form = NA)
p.control.wet75.2.habitat <- predict(everyone.habitat, newdata = control.wet75.2, type='link', re.form = NA)
p.control.wet75.2.social <- predict(everyone.social, newdata = control.wet75.2, type='link', re.form = NA)
p.control.wet75.2.pack <- predict(everyone.pack, newdata = control.wet75.2, type='link', re.form = NA)



control.wet25.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj =mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= 0.25,      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = median(distance2, na.rm = T),
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('control', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.control.wet25.2 <- predict(everyone, newdata = control.wet25.2, type='link', re.form = NA)
p.control.wet25.2.habitat <- predict(everyone.habitat, newdata = control.wet25.2, type='link', re.form = NA)
p.control.wet25.2.social <- predict(everyone.social, newdata = control.wet25.2, type='link', re.form = NA)
p.control.wet25.2.pack <- predict(everyone.pack, newdata = control.wet25.2, type='link', re.form = NA)


# road
control.roadclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),      
                                                             roadDist_end = 250,
                                                             distance2 = median(distance2, na.rm = T),
                                                             # nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = median(packDist_end, na.rm = T),
                                                             COD = factor('control', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID= NA)]
p.control.roadclose.2 <- predict(everyone, newdata = control.roadclose.2, type='link', re.form = NA)
p.control.roadclose.2.habitat <- predict(everyone.habitat, newdata = control.roadclose.2, type='link', re.form = NA)
p.control.roadclose.2.social <- predict(everyone.social, newdata = control.roadclose.2, type='link', re.form = NA)
p.control.roadclose.2.pack <- predict(everyone.pack, newdata = control.roadclose.2, type='link', re.form = NA)



control.roadfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = 1000,
                                                           distance2 = median(distance2, na.rm = T),
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('control', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.control.roadfar.2 <- predict(everyone, newdata = control.roadfar.2, type='link', re.form = NA)
p.control.roadfar.2.habitat <- predict(everyone.habitat, newdata = control.roadfar.2, type='link', re.form = NA)
p.control.roadfar.2.social <- predict(everyone.social, newdata = control.roadfar.2, type='link', re.form = NA)
p.control.roadfar.2.pack <- predict(everyone.pack, newdata = control.roadfar.2, type='link', re.form = NA)


# NN
control.nnclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = 250,
                                                           # nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = median(packDist_end, na.rm = T),
                                                           COD = factor('control', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.control.nnclose.2 <- predict(everyone, newdata = control.nnclose.2, type='link', re.form = NA)
p.control.nnclose.2.habitat <- predict(everyone.habitat, newdata = control.nnclose.2, type='link', re.form = NA)
p.control.nnclose.2.social <- predict(everyone.social, newdata = control.nnclose.2, type='link', re.form = NA)
p.control.nnclose.2.pack <- predict(everyone.pack, newdata = control.nnclose.2, type='link', re.form = NA)



control.nnfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                         log_sl = mean(log_sl),
                                                         cos_ta = mean(cos_ta),
                                                         # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                         propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                         propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                         propwet_end= mean(propwet_end, na.rm = T),      
                                                         roadDist_end = median(roadDist_end, na.rm = T),
                                                         distance2 = 1000,
                                                         nnDist_end = median(nnDist_end, na.rm = T),
                                                         packDist_end = median(packDist_end, na.rm = T),
                                                         COD = factor('control', levels = levels(COD)),
                                                         ttd1 = 60,
                                                         wolf_step_id = NA, wolfID= NA)]
p.control.nnfar.2 <- predict(everyone, newdata = control.nnfar.2, type='link', re.form = NA)
p.control.nnfar.2.habitat <- predict(everyone.habitat, newdata = control.nnfar.2, type='link', re.form = NA)
p.control.nnfar.2.social <- predict(everyone.social, newdata = control.nnfar.2, type='link', re.form = NA)
p.control.nnfar.2.pack <- predict(everyone.pack, newdata = control.nnfar.2, type='link', re.form = NA)



# pack
control.packclose.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             propforest_end_adj =  mean(propforest_end_adj, na.rm = T),
                                                             propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                             propwet_end= mean(propwet_end, na.rm = T),      
                                                             roadDist_end = median(roadDist_end, na.rm = T),
                                                             distance2 = median(distance2, na.rm = T),
                                                             # nnDist_end = median(nnDist_end, na.rm = T),
                                                             packDist_end = 250,
                                                             COD = factor('control', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID= NA)]
p.control.packclose.2 <- predict(everyone, newdata = control.packclose.2, type='link', re.form = NA)
p.control.packclose.2.habitat <- predict(everyone.habitat, newdata = control.packclose.2, type='link', re.form = NA)
p.control.packclose.2.social <- predict(everyone.social, newdata = control.packclose.2, type='link', re.form = NA)
p.control.packclose.2.pack <- predict(everyone.pack, newdata = control.packclose.2, type='link', re.form = NA)



control.packfar.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                           log_sl = mean(log_sl),
                                                           cos_ta = mean(cos_ta),
                                                           # land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                           propforest_end_adj = mean(propforest_end_adj, na.rm = T),
                                                           propopen_end_adj = mean(propopen_end_adj, na.rm = T),
                                                           propwet_end= mean(propwet_end, na.rm = T),      
                                                           roadDist_end = median(roadDist_end, na.rm = T),
                                                           distance2 = median(distance2, na.rm = T),
                                                           nnDist_end = median(nnDist_end, na.rm = T),
                                                           packDist_end = 1000,
                                                           COD = factor('control', levels = levels(COD)),
                                                           ttd1 = 60,
                                                           wolf_step_id = NA, wolfID= NA)]
p.control.packfar.2 <- predict(everyone, newdata = control.packfar.2, type='link', re.form = NA)
p.control.packfar.2.habitat <- predict(everyone.habitat, newdata = control.packfar.2, type='link', re.form = NA)
p.control.packfar.2.social <- predict(everyone.social, newdata = control.packfar.2, type='link', re.form = NA)
p.control.packfar.2.pack <- predict(everyone.pack, newdata = control.packfar.2, type='link', re.form = NA)

### INDIVs ###
p.h2.indiv <- function(ids, DT, mod, death, var, value){
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
            propforest_end_adj = ifelse(var == 'forest', value, mean(propforest_end_adj, na.rm = T)),
            propopen_end_adj = ifelse(var == 'open', value, mean(propopen_end_adj, na.rm = T)),
            propwet_end= ifelse(var == 'wet', value, mean(propwet_end, na.rm = T)), 
            roadDist_end = ifelse(var == 'road', value, median(roadDist_end, na.rm = T)),
            distance2 = ifelse(var == 'nn', value, median(distance2, na.rm = T)),
            nnDist_end = median(nnDist_end, na.rm = T),
            packDist_end = ifelse(var == 'pack', value, median(packDist_end, na.rm = T)),
            COD = factor(death, levels = levels(COD)),
            ttd1 = 60,
            wolf_step_id = NA,
            wolfID = i
          )],
          type = "link",
          re.form = NULL
        ), wolfID = i)]
    )
  })}

### CDV ####
CDV.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='CDV', wolfID])

# forest 
p.CDV.forest75.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'forest', value = 0.75)
p.CDV.forest75.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'forest', value = 0.75)

p.CDV.forest25.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'forest', value = 0.25)
p.CDV.forest25.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'forest', value = 0.25)

# open
p.CDV.open75.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'open', value = 0.75)
p.CDV.open75.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'open', value = 0.75)

p.CDV.open25.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'open', value = 0.25)
p.CDV.open25.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'open', value = 0.25)

# wet 
p.CDV.wet75.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'wet', value = 0.75)
p.CDV.wet75.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'wet', value = 0.75)

p.CDV.wet25.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'wet', value = 0.25)
p.CDV.wet25.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'wet', value = 0.25)


# road 
p.CDV.roadclose.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'road', value = 250)
p.CDV.roadclose.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'road', value = 250)

p.CDV.roadfar.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'road', value = 1000)
p.CDV.roadfar.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'road', value = 1000)

# nn 
p.CDV.nnclose.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'nn', value = 250)
p.CDV.nnclose.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'nn', value = 250)

p.CDV.nnfar.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'nn', value = 1000)
p.CDV.nnfar.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'nn', value = 1000)


# pack 
p.CDV.packclose.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'pack', value = 250)
p.CDV.packclose.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'pack', value = 250)
p.CDV.packclose.2.indiv.pack <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', var = 'pack', value = 250)

p.CDV.packfar.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'pack', value = 1000)
p.CDV.packfar.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'pack', value = 1000)
p.CDV.packfar.2.indiv.pack <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', var = 'pack', value = 1000)


### human ####
human.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='human', wolfID])

# forest 
p.human.forest75.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'forest', value = 0.75)
p.human.forest75.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'forest', value = 0.75)

p.human.forest25.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'forest', value = 0.25)
p.human.forest25.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'forest', value = 0.25)

# open
p.human.open75.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'open', value = 0.75)
p.human.open75.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'open', value = 0.75)

p.human.open25.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'open', value = 0.25)
p.human.open25.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'open', value = 0.25)

# wet 
p.human.wet75.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'wet', value = 0.75)
p.human.wet75.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'wet', value = 0.75)

p.human.wet25.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'wet', value = 0.25)
p.human.wet25.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'wet', value = 0.25)


# road 
p.human.roadclose.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'road', value = 250)
p.human.roadclose.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'road', value = 250)

p.human.roadfar.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'road', value = 1000)
p.human.roadfar.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'road', value = 1000)

# nn 
p.human.nnclose.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'nn', value = 250)
p.human.nnclose.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'nn', value = 250)

p.human.nnfar.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'nn', value = 1000)
p.human.nnfar.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'nn', value = 1000)


# pack 
p.human.packclose.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'pack', value = 250)
p.human.packclose.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'pack', value = 250)
p.human.packclose.2.indiv.pack <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', var = 'pack', value = 250)

p.human.packfar.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'pack', value = 1000)
p.human.packfar.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'pack', value = 1000)
p.human.packfar.2.indiv.pack <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', var = 'pack', value = 1000)


### control ####
control.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='control', wolfID])

# forest 
p.control.forest75.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'forest', value = 0.75)
p.control.forest75.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'forest', value = 0.75)

p.control.forest25.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'forest', value = 0.25)
p.control.forest25.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'forest', value = 0.25)

# open
p.control.open75.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'open', value = 0.75)
p.control.open75.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'open', value = 0.75)

p.control.open25.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'open', value = 0.25)
p.control.open25.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'open', value = 0.25)

# wet 
p.control.wet75.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'wet', value = 0.75)
p.control.wet75.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'wet', value = 0.75)

p.control.wet25.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'wet', value = 0.25)
p.control.wet25.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'wet', value = 0.25)


# road 
p.control.roadclose.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'road', value = 250)
p.control.roadclose.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'road', value = 250)

p.control.roadfar.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'road', value = 1000)
p.control.roadfar.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'road', value = 1000)

# nn 
p.control.nnclose.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'nn', value = 250)
p.control.nnclose.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'nn', value = 250)

p.control.nnfar.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'nn', value = 1000)
p.control.nnfar.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'nn', value = 1000)


# pack 
p.control.packclose.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'pack', value = 250)
p.control.packclose.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'pack', value = 250)
p.control.packclose.2.indiv.pack <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', var = 'pack', value = 250)

p.control.packfar.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'pack', value = 1000)
p.control.packfar.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'pack', value = 1000)
p.control.packfar.2.indiv.pack <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', var = 'pack', value = 1000)


#### h1 RSS equations ####
p.h1 <- function(DT, mod, death, var, value){
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          propforest_end_adj = ifelse(var == 'forest', value, mean(propforest_end_adj, na.rm = T)),
          propopen_end_adj = ifelse(var == 'open', value, mean(propopen_end_adj, na.rm = T)),
          propwet_end= ifelse(var == 'wet', value, mean(propwet_end, na.rm = T)), 
          roadDist_end = ifelse(var == 'road', value, median(roadDist_end, na.rm = T)),
          distance2 = ifelse(var == 'nn', value, median(distance2, na.rm = T)),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = ifelse(var == 'pack', value, median(packDist_end, na.rm = T)),
          COD = factor(death, levels = levels(COD)),
          ttd1 = 0:61,
          wolf_step_id = NA,
          wolfID = NA
        )],
        type = "link",
        re.form = NA
      ), 
      COD = death, var = var, value = value, ttd = 0:61)]
    # )
}


p.h1.indiv <- function(ids, DT, mod, death, var, value){
  lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          propforest_end_adj = ifelse(var == 'forest', value, mean(propforest_end_adj, na.rm = T)),
          propopen_end_adj = ifelse(var == 'open', value, mean(propopen_end_adj, na.rm = T)),
          propwet_end= ifelse(var == 'wet', value, mean(propwet_end, na.rm = T)), 
          roadDist_end = ifelse(var == 'road', value, median(roadDist_end, na.rm = T)),
          distance2 = ifelse(var == 'nn', value, median(distance2, na.rm = T)),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = ifelse(var == 'pack', value, median(packDist_end, na.rm = T)),
          COD = factor(death, levels = levels(COD)),
          ttd1 = 0:61,
          wolf_step_id = NA,
          wolfID = i
        )],
        type = "link",
        re.form = NULL
      ), 
      wolfID = i, COD = death, var = var, value = value, ttd = 0:61)]
    # )
  })
}

#### h1 RSS ####
#### POP ####
#### CDV ####
# forest 
p.CDV.forest75.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'forest', value = 0.75)
p.CDV.forest75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'forest', value = 0.75)

p.CDV.forest25.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'forest', value = 0.25)
p.CDV.forest25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'forest', value = 0.25)

# open
p.CDV.open75.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'open', value = 0.75)
p.CDV.open75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'open', value = 0.75)

p.CDV.open25.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'open', value = 0.25)
p.CDV.open25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'open', value = 0.25)

# wet 
p.CDV.wet75.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'wet', value = 0.75)
p.CDV.wet75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'wet', value = 0.75)

p.CDV.wet25.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'wet', value = 0.25)
p.CDV.wet25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'wet', value = 0.25)


# road 
p.CDV.roadclose.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'road', value = 250)
p.CDV.roadclose.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'road', value = 250)

p.CDV.roadfar.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'road', value = 1000)
p.CDV.roadfar.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'CDV', var = 'road', value = 1000)

# nn 
p.CDV.nnclose.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'nn', value = 250)
p.CDV.nnclose.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'CDV', var = 'nn', value = 250)

p.CDV.nnfar.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'nn', value = 1000)
p.CDV.nnfar.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'CDV', var = 'nn', value = 1000)


# pack 
p.CDV.packclose.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'pack', value = 250)
p.CDV.packclose.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'CDV', var = 'pack', value = 250)
p.CDV.packclose.1.pack <-p.h1(DT = dat, mod = everyone.pack, death = 'CDV', var = 'pack', value = 250)

p.CDV.packfar.1 <-p.h1(DT = dat, mod = everyone, death = 'CDV', var = 'pack', value = 1000)
p.CDV.packfar.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'CDV', var = 'pack', value = 1000)
p.CDV.packfar.1.pack <-p.h1(DT = dat, mod = everyone.pack, death = 'CDV', var = 'pack', value = 1000)


### human ####
# forest 
p.human.forest75.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'forest', value = 0.75)
p.human.forest75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'forest', value = 0.75)

p.human.forest25.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'forest', value = 0.25)
p.human.forest25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'forest', value = 0.25)

# open
p.human.open75.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'open', value = 0.75)
p.human.open75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'open', value = 0.75)

p.human.open25.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'open', value = 0.25)
p.human.open25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'open', value = 0.25)

# wet 
p.human.wet75.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'wet', value = 0.75)
p.human.wet75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'wet', value = 0.75)

p.human.wet25.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'wet', value = 0.25)
p.human.wet25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'wet', value = 0.25)


# road 
p.human.roadclose.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'road', value = 250)
p.human.roadclose.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'road', value = 250)

p.human.roadfar.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'road', value = 1000)
p.human.roadfar.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'human', var = 'road', value = 1000)

# nn 
p.human.nnclose.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'nn', value = 250)
p.human.nnclose.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'human', var = 'nn', value = 250)

p.human.nnfar.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'nn', value = 1000)
p.human.nnfar.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'human', var = 'nn', value = 1000)


# pack 
p.human.packclose.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'pack', value = 250)
p.human.packclose.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'human', var = 'pack', value = 250)
p.human.packclose.1.pack <-p.h1(DT = dat, mod = everyone.pack, death = 'human', var = 'pack', value = 250)

p.human.packfar.1 <-p.h1(DT = dat, mod = everyone, death = 'human', var = 'pack', value = 1000)
p.human.packfar.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'human', var = 'pack', value = 1000)
p.human.packfar.1.pack <-p.h1(DT = dat, mod = everyone.pack, death = 'human', var = 'pack', value = 1000)


### control ####
# forest 
p.control.forest75.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'forest', value = 0.75)
p.control.forest75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'forest', value = 0.75)

p.control.forest25.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'forest', value = 0.25)
p.control.forest25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'forest', value = 0.25)

# open
p.control.open75.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'open', value = 0.75)
p.control.open75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'open', value = 0.75)

p.control.open25.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'open', value = 0.25)
p.control.open25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'open', value = 0.25)

# wet 
p.control.wet75.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'wet', value = 0.75)
p.control.wet75.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'wet', value = 0.75)

p.control.wet25.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'wet', value = 0.25)
p.control.wet25.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'wet', value = 0.25)


# road 
p.control.roadclose.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'road', value = 250)
p.control.roadclose.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'road', value = 250)

p.control.roadfar.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'road', value = 1000)
p.control.roadfar.1.habitat <-p.h1(DT = dat, mod = everyone.habitat, death = 'control', var = 'road', value = 1000)

# nn 
p.control.nnclose.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'nn', value = 250)
p.control.nnclose.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'control', var = 'nn', value = 250)

p.control.nnfar.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'nn', value = 1000)
p.control.nnfar.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'control', var = 'nn', value = 1000)


# pack 
p.control.packclose.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'pack', value = 250)
p.control.packclose.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'control', var = 'pack', value = 250)
p.control.packclose.1.pack <-p.h1(DT = dat, mod = everyone.pack, death = 'control', var = 'pack', value = 250)

p.control.packfar.1 <-p.h1(DT = dat, mod = everyone, death = 'control', var = 'pack', value = 1000)
p.control.packfar.1.social <-p.h1(DT = dat, mod = everyone.social, death = 'control', var = 'pack', value = 1000)
p.control.packfar.1.pack <-p.h1(DT = dat, mod = everyone.pack, death = 'control', var = 'pack', value = 1000)



#### INDIVS ####
#### CDV ####
# forest 
p.CDV.forest75.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'forest', value = 0.75)
p.CDV.forest75.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'forest', value = 0.75)

p.CDV.forest25.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'forest', value = 0.25)
p.CDV.forest25.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'forest', value = 0.25)

# open
p.CDV.open75.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'open', value = 0.75)
p.CDV.open75.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'open', value = 0.75)

p.CDV.open25.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'open', value = 0.25)
p.CDV.open25.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'open', value = 0.25)

# wet 
p.CDV.wet75.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'wet', value = 0.75)
p.CDV.wet75.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'wet', value = 0.75)

p.CDV.wet25.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'wet', value = 0.25)
p.CDV.wet25.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'wet', value = 0.25)


# road 
p.CDV.roadclose.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'road', value = 250)
p.CDV.roadclose.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'road', value = 250)

p.CDV.roadfar.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'road', value = 1000)
p.CDV.roadfar.1.indiv.habitat <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', var = 'road', value = 1000)

# nn 
p.CDV.nnclose.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'nn', value = 250)
p.CDV.nnclose.1.indiv.social <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'nn', value = 250)

p.CDV.nnfar.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'nn', value = 1000)
p.CDV.nnfar.1.indiv.social <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'nn', value = 1000)


# pack 
p.CDV.packclose.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'pack', value = 250)
p.CDV.packclose.1.indiv.social <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'pack', value = 250)
p.CDV.packclose.1.indiv.pack <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', var = 'pack', value = 250)

p.CDV.packfar.1.indiv <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', var = 'pack', value = 1000)
p.CDV.packfar.1.indiv.social <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', var = 'pack', value = 1000)
p.CDV.packfar.1.indiv.pack <-p.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.pack, death = 'CDV', var = 'pack', value = 1000)


### human ####
human.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='human', wolfID])

# forest 
p.human.forest75.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'forest', value = 0.75)
p.human.forest75.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'forest', value = 0.75)

p.human.forest25.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'forest', value = 0.25)
p.human.forest25.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'forest', value = 0.25)

# open
p.human.open75.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'open', value = 0.75)
p.human.open75.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'open', value = 0.75)

p.human.open25.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'open', value = 0.25)
p.human.open25.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'open', value = 0.25)

# wet 
p.human.wet75.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'wet', value = 0.75)
p.human.wet75.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'wet', value = 0.75)

p.human.wet25.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'wet', value = 0.25)
p.human.wet25.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'wet', value = 0.25)


# road 
p.human.roadclose.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'road', value = 250)
p.human.roadclose.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'road', value = 250)

p.human.roadfar.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'road', value = 1000)
p.human.roadfar.1.indiv.habitat <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', var = 'road', value = 1000)

# nn 
p.human.nnclose.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'nn', value = 250)
p.human.nnclose.1.indiv.social <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'nn', value = 250)

p.human.nnfar.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'nn', value = 1000)
p.human.nnfar.1.indiv.social <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'nn', value = 1000)


# pack 
p.human.packclose.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'pack', value = 250)
p.human.packclose.1.indiv.social <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'pack', value = 250)
p.human.packclose.1.indiv.pack <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', var = 'pack', value = 250)

p.human.packfar.1.indiv <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', var = 'pack', value = 1000)
p.human.packfar.1.indiv.social <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', var = 'pack', value = 1000)
p.human.packfar.1.indiv.pack <-p.h1.indiv(ids = human.wolfID, DT = dat, mod = everyone.pack, death = 'human', var = 'pack', value = 1000)


### control ####
control.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='control', wolfID])

# forest 
p.control.forest75.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'forest', value = 0.75)
p.control.forest75.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'forest', value = 0.75)

p.control.forest25.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'forest', value = 0.25)
p.control.forest25.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'forest', value = 0.25)

# open
p.control.open75.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'open', value = 0.75)
p.control.open75.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'open', value = 0.75)

p.control.open25.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'open', value = 0.25)
p.control.open25.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'open', value = 0.25)

# wet 
p.control.wet75.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'wet', value = 0.75)
p.control.wet75.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'wet', value = 0.75)

p.control.wet25.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'wet', value = 0.25)
p.control.wet25.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'wet', value = 0.25)


# road 
p.control.roadclose.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'road', value = 250)
p.control.roadclose.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'road', value = 250)

p.control.roadfar.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'road', value = 1000)
p.control.roadfar.1.indiv.habitat <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', var = 'road', value = 1000)

# nn 
p.control.nnclose.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'nn', value = 250)
p.control.nnclose.1.indiv.social <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'nn', value = 250)

p.control.nnfar.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'nn', value = 1000)
p.control.nnfar.1.indiv.social <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'nn', value = 1000)


# pack 
p.control.packclose.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'pack', value = 250)
p.control.packclose.1.indiv.social <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'pack', value = 250)
p.control.packclose.1.indiv.pack <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', var = 'pack', value = 250)

p.control.packfar.1.indiv <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', var = 'pack', value = 1000)
p.control.packfar.1.indiv.social <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', var = 'pack', value = 1000)
p.control.packfar.1.indiv.pack <-p.h1.indiv(ids = control.wolfID, DT = dat, mod = everyone.pack, death = 'control', var = 'pack', value = 1000)



#### CALC RSS ####
### POP ####
### CDV ####
# forest
logRSS.CDV.forest75 <- p.CDV.forest75.1[, .(COD, var, value, h1, h2 = p.CDV.forest75.2, rss = h1 - p.CDV.forest75.2, ttd)]
logRSS.CDV.forest25 <- p.CDV.forest25.1[, .(COD, var, value, h1, h2 = p.CDV.forest25.2, rss = h1 - p.CDV.forest25.2, ttd)]

logRSS.CDV.forest75.habitat <- p.CDV.forest75.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.forest75.2.habitat, rss = h1 - p.CDV.forest75.2.habitat, ttd)]
logRSS.CDV.forest25.habitat <- p.CDV.forest25.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.forest25.2.habitat, rss = h1 - p.CDV.forest25.2.habitat, ttd)]


# open
logRSS.CDV.open75 <- p.CDV.open75.1[, .(COD, var, value, h1, h2 = p.CDV.open75.2, rss = h1 - p.CDV.open75.2, ttd)]
logRSS.CDV.open25 <- p.CDV.open25.1[, .(COD, var, value, h1, h2 = p.CDV.open25.2, rss = h1 - p.CDV.open25.2, ttd)]

logRSS.CDV.open75.habitat <- p.CDV.open75.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.open75.2.habitat, rss = h1 - p.CDV.open75.2.habitat, ttd)]
logRSS.CDV.open25.habitat <- p.CDV.open25.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.open25.2.habitat, rss = h1 - p.CDV.open25.2.habitat, ttd)]


# wet
logRSS.CDV.wet75 <- p.CDV.wet75.1[, .(COD, var, value, h1, h2 = p.CDV.wet75.2, rss = h1 - p.CDV.wet75.2, ttd)]
logRSS.CDV.wet25 <- p.CDV.wet25.1[, .(COD, var, value, h1, h2 = p.CDV.wet25.2, rss = h1 - p.CDV.wet25.2, ttd)]

logRSS.CDV.wet75.habitat <- p.CDV.wet75.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.wet75.2.habitat, rss = h1 - p.CDV.wet75.2.habitat, ttd)]
logRSS.CDV.wet25.habitat <- p.CDV.wet25.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.wet25.2.habitat, rss = h1 - p.CDV.wet25.2.habitat, ttd)]


# road
logRSS.CDV.roadclose <- p.CDV.roadclose.1[, .(COD, var, value, h1, h2 = p.CDV.roadclose.2, rss = h1 - p.CDV.roadclose.2, ttd)]
logRSS.CDV.roadfar <- p.CDV.roadfar.1[, .(COD, var, value, h1, h2 = p.CDV.roadfar.2, rss = h1 - p.CDV.roadfar.2, ttd)]

logRSS.CDV.roadclose.habitat <- p.CDV.roadclose.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.roadclose.2.habitat, rss = h1 - p.CDV.roadclose.2.habitat, ttd)]
logRSS.CDV.roadfar.habitat <- p.CDV.roadfar.1.habitat[, .(COD, var, value, h1, h2 = p.CDV.roadfar.2.habitat, rss = h1 - p.CDV.roadfar.2.habitat, ttd)]


# nn
logRSS.CDV.nnclose <- p.CDV.nnclose.1[, .(COD, var, value, h1, h2 = p.CDV.nnclose.2, rss = h1 - p.CDV.nnclose.2, ttd)]
logRSS.CDV.nnfar <- p.CDV.nnfar.1[, .(COD, var, value, h1, h2 = p.CDV.nnfar.2, rss = h1 - p.CDV.nnfar.2, ttd)]

logRSS.CDV.nnclose.social <- p.CDV.nnclose.1.social[, .(COD, var, value, h1, h2 = p.CDV.nnclose.2.social, rss = h1 - p.CDV.nnclose.2.social, ttd)]
logRSS.CDV.nnfar.social <- p.CDV.nnfar.1.social[, .(COD, var, value, h1, h2 = p.CDV.nnfar.2.social, rss = h1 - p.CDV.nnfar.2.social, ttd)]


# pack
logRSS.CDV.packclose <- p.CDV.packclose.1[, .(COD, var, value, h1, h2 = p.CDV.packclose.2, rss = h1 - p.CDV.packclose.2, ttd)]
logRSS.CDV.packfar <- p.CDV.packfar.1[, .(COD, var, value, h1, h2 = p.CDV.packfar.2, rss = h1 - p.CDV.packfar.2, ttd)]

logRSS.CDV.packclose.social <- p.CDV.packclose.1.social[, .(COD, var, value, h1, h2 = p.CDV.packclose.2.social, rss = h1 - p.CDV.packclose.2.social, ttd)]
logRSS.CDV.packfar.social <- p.CDV.packfar.1.social[, .(COD, var, value, h1, h2 = p.CDV.packfar.2.social, rss = h1 - p.CDV.packfar.2.social, ttd)]

logRSS.CDV.packclose.pack <- p.CDV.packclose.1.pack[, .(COD, var, value, h1, h2 = p.CDV.packclose.2.pack, rss = h1 - p.CDV.packclose.2.pack, ttd)]
logRSS.CDV.packfar.pack <- p.CDV.packfar.1.pack[, .(COD, var, value, h1, h2 = p.CDV.packfar.2.pack, rss = h1 - p.CDV.packfar.2.pack, ttd)]




### human ####
# forest
logRSS.human.forest75 <- p.human.forest75.1[, .(COD, var, value, h1, h2 = p.human.forest75.2, rss = h1 - p.human.forest75.2, ttd)]
logRSS.human.forest25 <- p.human.forest25.1[, .(COD, var, value, h1, h2 = p.human.forest25.2, rss = h1 - p.human.forest25.2, ttd)]

logRSS.human.forest75.habitat <- p.human.forest75.1.habitat[, .(COD, var, value, h1, h2 = p.human.forest75.2.habitat, rss = h1 - p.human.forest75.2.habitat, ttd)]
logRSS.human.forest25.habitat <- p.human.forest25.1.habitat[, .(COD, var, value, h1, h2 = p.human.forest25.2.habitat, rss = h1 - p.human.forest25.2.habitat, ttd)]


# open
logRSS.human.open75 <- p.human.open75.1[, .(COD, var, value, h1, h2 = p.human.open75.2, rss = h1 - p.human.open75.2, ttd)]
logRSS.human.open25 <- p.human.open25.1[, .(COD, var, value, h1, h2 = p.human.open25.2, rss = h1 - p.human.open25.2, ttd)]

logRSS.human.open75.habitat <- p.human.open75.1.habitat[, .(COD, var, value, h1, h2 = p.human.open75.2.habitat, rss = h1 - p.human.open75.2.habitat, ttd)]
logRSS.human.open25.habitat <- p.human.open25.1.habitat[, .(COD, var, value, h1, h2 = p.human.open25.2.habitat, rss = h1 - p.human.open25.2.habitat, ttd)]


# wet
logRSS.human.wet75 <- p.human.wet75.1[, .(COD, var, value, h1, h2 = p.human.wet75.2, rss = h1 - p.human.wet75.2, ttd)]
logRSS.human.wet25 <- p.human.wet25.1[, .(COD, var, value, h1, h2 = p.human.wet25.2, rss = h1 - p.human.wet25.2, ttd)]

logRSS.human.wet75.habitat <- p.human.wet75.1.habitat[, .(COD, var, value, h1, h2 = p.human.wet75.2.habitat, rss = h1 - p.human.wet75.2.habitat, ttd)]
logRSS.human.wet25.habitat <- p.human.wet25.1.habitat[, .(COD, var, value, h1, h2 = p.human.wet25.2.habitat, rss = h1 - p.human.wet25.2.habitat, ttd)]


# road
logRSS.human.roadclose <- p.human.roadclose.1[, .(COD, var, value, h1, h2 = p.human.roadclose.2, rss = h1 - p.human.roadclose.2, ttd)]
logRSS.human.roadfar <- p.human.roadfar.1[, .(COD, var, value, h1, h2 = p.human.roadfar.2, rss = h1 - p.human.roadfar.2, ttd)]

logRSS.human.roadclose.habitat <- p.human.roadclose.1.habitat[, .(COD, var, value, h1, h2 = p.human.roadclose.2.habitat, rss = h1 - p.human.roadclose.2.habitat, ttd)]
logRSS.human.roadfar.habitat <- p.human.roadfar.1.habitat[, .(COD, var, value, h1, h2 = p.human.roadfar.2.habitat, rss = h1 - p.human.roadfar.2.habitat, ttd)]


# nn
logRSS.human.nnclose <- p.human.nnclose.1[, .(COD, var, value, h1, h2 = p.human.nnclose.2, rss = h1 - p.human.nnclose.2, ttd)]
logRSS.human.nnfar <- p.human.nnfar.1[, .(COD, var, value, h1, h2 = p.human.nnfar.2, rss = h1 - p.human.nnfar.2, ttd)]

logRSS.human.nnclose.social <- p.human.nnclose.1.social[, .(COD, var, value, h1, h2 = p.human.nnclose.2.social, rss = h1 - p.human.nnclose.2.social, ttd)]
logRSS.human.nnfar.social <- p.human.nnfar.1.social[, .(COD, var, value, h1, h2 = p.human.nnfar.2.social, rss = h1 - p.human.nnfar.2.social, ttd)]


# pack
logRSS.human.packclose <- p.human.packclose.1[, .(COD, var, value, h1, h2 = p.human.packclose.2, rss = h1 - p.human.packclose.2, ttd)]
logRSS.human.packfar <- p.human.packfar.1[, .(COD, var, value, h1, h2 = p.human.packfar.2, rss = h1 - p.human.packfar.2, ttd)]

logRSS.human.packclose.social <- p.human.packclose.1.social[, .(COD, var, value, h1, h2 = p.human.packclose.2.social, rss = h1 - p.human.packclose.2.social, ttd)]
logRSS.human.packfar.social <- p.human.packfar.1.social[, .(COD, var, value, h1, h2 = p.human.packfar.2.social, rss = h1 - p.human.packfar.2.social, ttd)]

logRSS.human.packclose.pack <- p.human.packclose.1.pack[, .(COD, var, value, h1, h2 = p.human.packclose.2.pack, rss = h1 - p.human.packclose.2.pack, ttd)]
logRSS.human.packfar.pack <- p.human.packfar.1.pack[, .(COD, var, value, h1, h2 = p.human.packfar.2.pack, rss = h1 - p.human.packfar.2.pack, ttd)]



### control ####
# forest
logRSS.control.forest75 <- p.control.forest75.1[, .(COD, var, value, h1, h2 = p.control.forest75.2, rss = h1 - p.control.forest75.2, ttd)]
logRSS.control.forest25 <- p.control.forest25.1[, .(COD, var, value, h1, h2 = p.control.forest25.2, rss = h1 - p.control.forest25.2, ttd)]

logRSS.control.forest75.habitat <- p.control.forest75.1.habitat[, .(COD, var, value, h1, h2 = p.control.forest75.2.habitat, rss = h1 - p.control.forest75.2.habitat, ttd)]
logRSS.control.forest25.habitat <- p.control.forest25.1.habitat[, .(COD, var, value, h1, h2 = p.control.forest25.2.habitat, rss = h1 - p.control.forest25.2.habitat, ttd)]


# open
logRSS.control.open75 <- p.control.open75.1[, .(COD, var, value, h1, h2 = p.control.open75.2, rss = h1 - p.control.open75.2, ttd)]
logRSS.control.open25 <- p.control.open25.1[, .(COD, var, value, h1, h2 = p.control.open25.2, rss = h1 - p.control.open25.2, ttd)]

logRSS.control.open75.habitat <- p.control.open75.1.habitat[, .(COD, var, value, h1, h2 = p.control.open75.2.habitat, rss = h1 - p.control.open75.2.habitat, ttd)]
logRSS.control.open25.habitat <- p.control.open25.1.habitat[, .(COD, var, value, h1, h2 = p.control.open25.2.habitat, rss = h1 - p.control.open25.2.habitat, ttd)]


# wet
logRSS.control.wet75 <- p.control.wet75.1[, .(COD, var, value, h1, h2 = p.control.wet75.2, rss = h1 - p.control.wet75.2, ttd)]
logRSS.control.wet25 <- p.control.wet25.1[, .(COD, var, value, h1, h2 = p.control.wet25.2, rss = h1 - p.control.wet25.2, ttd)]

logRSS.control.wet75.habitat <- p.control.wet75.1.habitat[, .(COD, var, value, h1, h2 = p.control.wet75.2.habitat, rss = h1 - p.control.wet75.2.habitat, ttd)]
logRSS.control.wet25.habitat <- p.control.wet25.1.habitat[, .(COD, var, value, h1, h2 = p.control.wet25.2.habitat, rss = h1 - p.control.wet25.2.habitat, ttd)]


# road
logRSS.control.roadclose <- p.control.roadclose.1[, .(COD, var, value, h1, h2 = p.control.roadclose.2, rss = h1 - p.control.roadclose.2, ttd)]
logRSS.control.roadfar <- p.control.roadfar.1[, .(COD, var, value, h1, h2 = p.control.roadfar.2, rss = h1 - p.control.roadfar.2, ttd)]

logRSS.control.roadclose.habitat <- p.control.roadclose.1.habitat[, .(COD, var, value, h1, h2 = p.control.roadclose.2.habitat, rss = h1 - p.control.roadclose.2.habitat, ttd)]
logRSS.control.roadfar.habitat <- p.control.roadfar.1.habitat[, .(COD, var, value, h1, h2 = p.control.roadfar.2.habitat, rss = h1 - p.control.roadfar.2.habitat, ttd)]


# nn
logRSS.control.nnclose <- p.control.nnclose.1[, .(COD, var, value, h1, h2 = p.control.nnclose.2, rss = h1 - p.control.nnclose.2, ttd)]
logRSS.control.nnfar <- p.control.nnfar.1[, .(COD, var, value, h1, h2 = p.control.nnfar.2, rss = h1 - p.control.nnfar.2, ttd)]

logRSS.control.nnclose.social <- p.control.nnclose.1.social[, .(COD, var, value, h1, h2 = p.control.nnclose.2.social, rss = h1 - p.control.nnclose.2.social, ttd)]
logRSS.control.nnfar.social <- p.control.nnfar.1.social[, .(COD, var, value, h1, h2 = p.control.nnfar.2.social, rss = h1 - p.control.nnfar.2.social, ttd)]


# pack
logRSS.control.packclose <- p.control.packclose.1[, .(COD, var, value, h1, h2 = p.control.packclose.2, rss = h1 - p.control.packclose.2, ttd)]
logRSS.control.packfar <- p.control.packfar.1[, .(COD, var, value, h1, h2 = p.control.packfar.2, rss = h1 - p.control.packfar.2, ttd)]

logRSS.control.packclose.social <- p.control.packclose.1.social[, .(COD, var, value, h1, h2 = p.control.packclose.2.social, rss = h1 - p.control.packclose.2.social, ttd)]
logRSS.control.packfar.social <- p.control.packfar.1.social[, .(COD, var, value, h1, h2 = p.control.packfar.2.social, rss = h1 - p.control.packfar.2.social, ttd)]

logRSS.control.packclose.pack <- p.control.packclose.1.pack[, .(COD, var, value, h1, h2 = p.control.packclose.2.pack, rss = h1 - p.control.packclose.2.pack, ttd)]
logRSS.control.packfar.pack <- p.control.packfar.1.pack[, .(COD, var, value, h1, h2 = p.control.packfar.2.pack, rss = h1 - p.control.packfar.2.pack, ttd)]



### INDIVS ####
### CDV ####
# forest
h1.CDV.forest75.indiv <- rbindlist(p.CDV.forest75.1.indiv)
h2.CDV.forest75.indiv <- rbindlist(p.CDV.forest75.2.indiv)

h1.CDV.forest25.indiv <- rbindlist(p.CDV.forest25.1.indiv)
h2.CDV.forest25.indiv <- rbindlist(p.CDV.forest25.2.indiv)

logRSS.CDV.forest75.indiv <- merge(h1.CDV.forest75.indiv, h2.CDV.forest75.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.forest75.indiv[,'rss'] <- logRSS.CDV.forest75.indiv$h1 - logRSS.CDV.forest75.indiv$h2

logRSS.CDV.forest25.indiv <- merge(h1.CDV.forest25.indiv, h2.CDV.forest25.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.forest25.indiv[,'rss'] <- logRSS.CDV.forest25.indiv$h1 - logRSS.CDV.forest25.indiv$h2


h1.CDV.forest75.indiv.habitat <- rbindlist(p.CDV.forest75.1.indiv.habitat)
h2.CDV.forest75.indiv.habitat <- rbindlist(p.CDV.forest75.2.indiv.habitat)

h1.CDV.forest25.indiv.habitat <- rbindlist(p.CDV.forest25.1.indiv.habitat)
h2.CDV.forest25.indiv.habitat <- rbindlist(p.CDV.forest25.2.indiv.habitat)

logRSS.CDV.forest75.indiv.habitat <- merge(h1.CDV.forest75.indiv.habitat, h2.CDV.forest75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.forest75.indiv.habitat[,'rss'] <- logRSS.CDV.forest75.indiv.habitat$h1 - logRSS.CDV.forest75.indiv.habitat$h2

logRSS.CDV.forest25.indiv.habitat <- merge(h1.CDV.forest25.indiv.habitat, h2.CDV.forest25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.forest25.indiv.habitat[,'rss'] <- logRSS.CDV.forest25.indiv.habitat$h1 - logRSS.CDV.forest25.indiv.habitat$h2


# open
h1.CDV.open75.indiv <- rbindlist(p.CDV.open75.1.indiv)
h2.CDV.open75.indiv <- rbindlist(p.CDV.open75.2.indiv)

h1.CDV.open25.indiv <- rbindlist(p.CDV.open25.1.indiv)
h2.CDV.open25.indiv <- rbindlist(p.CDV.open25.2.indiv)

logRSS.CDV.open75.indiv <- merge(h1.CDV.open75.indiv, h2.CDV.open75.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.open75.indiv[,'rss'] <- logRSS.CDV.open75.indiv$h1 - logRSS.CDV.open75.indiv$h2

logRSS.CDV.open25.indiv <- merge(h1.CDV.open25.indiv, h2.CDV.open25.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.open25.indiv[,'rss'] <- logRSS.CDV.open25.indiv$h1 - logRSS.CDV.open25.indiv$h2


h1.CDV.open75.indiv.habitat <- rbindlist(p.CDV.open75.1.indiv.habitat)
h2.CDV.open75.indiv.habitat <- rbindlist(p.CDV.open75.2.indiv.habitat)

h1.CDV.open25.indiv.habitat <- rbindlist(p.CDV.open25.1.indiv.habitat)
h2.CDV.open25.indiv.habitat <- rbindlist(p.CDV.open25.2.indiv.habitat)

logRSS.CDV.open75.indiv.habitat <- merge(h1.CDV.open75.indiv.habitat, h2.CDV.open75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.open75.indiv.habitat[,'rss'] <- logRSS.CDV.open75.indiv.habitat$h1 - logRSS.CDV.open75.indiv.habitat$h2

logRSS.CDV.open25.indiv.habitat <- merge(h1.CDV.open25.indiv.habitat, h2.CDV.open25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.open25.indiv.habitat[,'rss'] <- logRSS.CDV.open25.indiv.habitat$h1 - logRSS.CDV.open25.indiv.habitat$h2


# wet
h1.CDV.wet75.indiv <- rbindlist(p.CDV.wet75.1.indiv)
h2.CDV.wet75.indiv <- rbindlist(p.CDV.wet75.2.indiv)

h1.CDV.wet25.indiv <- rbindlist(p.CDV.wet25.1.indiv)
h2.CDV.wet25.indiv <- rbindlist(p.CDV.wet25.2.indiv)

logRSS.CDV.wet75.indiv <- merge(h1.CDV.wet75.indiv, h2.CDV.wet75.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.wet75.indiv[,'rss'] <- logRSS.CDV.wet75.indiv$h1 - logRSS.CDV.wet75.indiv$h2

logRSS.CDV.wet25.indiv <- merge(h1.CDV.wet25.indiv, h2.CDV.wet25.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.wet25.indiv[,'rss'] <- logRSS.CDV.wet25.indiv$h1 - logRSS.CDV.wet25.indiv$h2


h1.CDV.wet75.indiv.habitat <- rbindlist(p.CDV.wet75.1.indiv.habitat)
h2.CDV.wet75.indiv.habitat <- rbindlist(p.CDV.wet75.2.indiv.habitat)

h1.CDV.wet25.indiv.habitat <- rbindlist(p.CDV.wet25.1.indiv.habitat)
h2.CDV.wet25.indiv.habitat <- rbindlist(p.CDV.wet25.2.indiv.habitat)

logRSS.CDV.wet75.indiv.habitat <- merge(h1.CDV.wet75.indiv.habitat, h2.CDV.wet75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.wet75.indiv.habitat[,'rss'] <- logRSS.CDV.wet75.indiv.habitat$h1 - logRSS.CDV.wet75.indiv.habitat$h2

logRSS.CDV.wet25.indiv.habitat <- merge(h1.CDV.wet25.indiv.habitat, h2.CDV.wet25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.wet25.indiv.habitat[,'rss'] <- logRSS.CDV.wet25.indiv.habitat$h1 - logRSS.CDV.wet25.indiv.habitat$h2


# road
h1.CDV.roadclose.indiv <- rbindlist(p.CDV.roadclose.1.indiv)
h2.CDV.roadclose.indiv <- rbindlist(p.CDV.roadclose.2.indiv)

h1.CDV.roadfar.indiv <- rbindlist(p.CDV.roadfar.1.indiv)
h2.CDV.roadfar.indiv <- rbindlist(p.CDV.roadfar.2.indiv)

logRSS.CDV.roadclose.indiv <- merge(h1.CDV.roadclose.indiv, h2.CDV.roadclose.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.roadclose.indiv[,'rss'] <- logRSS.CDV.roadclose.indiv$h1 - logRSS.CDV.roadclose.indiv$h2

logRSS.CDV.roadfar.indiv <- merge(h1.CDV.roadfar.indiv, h2.CDV.roadfar.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.roadfar.indiv[,'rss'] <- logRSS.CDV.roadfar.indiv$h1 - logRSS.CDV.roadfar.indiv$h2


h1.CDV.roadclose.indiv.habitat <- rbindlist(p.CDV.roadclose.1.indiv.habitat)
h2.CDV.roadclose.indiv.habitat <- rbindlist(p.CDV.roadclose.2.indiv.habitat)

h1.CDV.roadfar.indiv.habitat <- rbindlist(p.CDV.roadfar.1.indiv.habitat)
h2.CDV.roadfar.indiv.habitat <- rbindlist(p.CDV.roadfar.2.indiv.habitat)

logRSS.CDV.roadclose.indiv.habitat <- merge(h1.CDV.roadclose.indiv.habitat, h2.CDV.roadclose.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.roadclose.indiv.habitat[,'rss'] <- logRSS.CDV.roadclose.indiv.habitat$h1 - logRSS.CDV.roadclose.indiv.habitat$h2

logRSS.CDV.roadfar.indiv.habitat <- merge(h1.CDV.roadfar.indiv.habitat, h2.CDV.roadfar.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.CDV.roadfar.indiv.habitat[,'rss'] <- logRSS.CDV.roadfar.indiv.habitat$h1 - logRSS.CDV.roadfar.indiv.habitat$h2


# nn
h1.CDV.nnclose.indiv <- rbindlist(p.CDV.nnclose.1.indiv)
h2.CDV.nnclose.indiv <- rbindlist(p.CDV.nnclose.2.indiv)

h1.CDV.nnfar.indiv <- rbindlist(p.CDV.nnfar.1.indiv)
h2.CDV.nnfar.indiv <- rbindlist(p.CDV.nnfar.2.indiv)

logRSS.CDV.nnclose.indiv <- merge(h1.CDV.nnclose.indiv, h2.CDV.nnclose.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.nnclose.indiv[,'rss'] <- logRSS.CDV.nnclose.indiv$h1 - logRSS.CDV.nnclose.indiv$h2

logRSS.CDV.nnfar.indiv <- merge(h1.CDV.nnfar.indiv, h2.CDV.nnfar.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.nnfar.indiv[,'rss'] <- logRSS.CDV.nnfar.indiv$h1 - logRSS.CDV.nnfar.indiv$h2


h1.CDV.nnclose.indiv.social <- rbindlist(p.CDV.nnclose.1.indiv.social)
h2.CDV.nnclose.indiv.social <- rbindlist(p.CDV.nnclose.2.indiv.social)

h1.CDV.nnfar.indiv.social <- rbindlist(p.CDV.nnfar.1.indiv.social)
h2.CDV.nnfar.indiv.social <- rbindlist(p.CDV.nnfar.2.indiv.social)

logRSS.CDV.nnclose.indiv.social <- merge(h1.CDV.nnclose.indiv.social, h2.CDV.nnclose.indiv.social, by = c('wolfID'), all.x = T)
logRSS.CDV.nnclose.indiv.social[,'rss'] <- logRSS.CDV.nnclose.indiv.social$h1 - logRSS.CDV.nnclose.indiv.social$h2

logRSS.CDV.nnfar.indiv.social <- merge(h1.CDV.nnfar.indiv.social, h2.CDV.nnfar.indiv.social, by = c('wolfID'), all.x = T)
logRSS.CDV.nnfar.indiv.social[,'rss'] <- logRSS.CDV.nnfar.indiv.social$h1 - logRSS.CDV.nnfar.indiv.social$h2


# pack
h1.CDV.packclose.indiv <- rbindlist(p.CDV.packclose.1.indiv)
h2.CDV.packclose.indiv <- rbindlist(p.CDV.packclose.2.indiv)

h1.CDV.packfar.indiv <- rbindlist(p.CDV.packfar.1.indiv)
h2.CDV.packfar.indiv <- rbindlist(p.CDV.packfar.2.indiv)

logRSS.CDV.packclose.indiv <- merge(h1.CDV.packclose.indiv, h2.CDV.packclose.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.packclose.indiv[,'rss'] <- logRSS.CDV.packclose.indiv$h1 - logRSS.CDV.packclose.indiv$h2

logRSS.CDV.packfar.indiv <- merge(h1.CDV.packfar.indiv, h2.CDV.packfar.indiv, by = c('wolfID'), all.x = T)
logRSS.CDV.packfar.indiv[,'rss'] <- logRSS.CDV.packfar.indiv$h1 - logRSS.CDV.packfar.indiv$h2


h1.CDV.packclose.indiv.social <- rbindlist(p.CDV.packclose.1.indiv.social)
h2.CDV.packclose.indiv.social <- rbindlist(p.CDV.packclose.2.indiv.social)

h1.CDV.packfar.indiv.social <- rbindlist(p.CDV.packfar.1.indiv.social)
h2.CDV.packfar.indiv.social <- rbindlist(p.CDV.packfar.2.indiv.social)

logRSS.CDV.packclose.indiv.social <- merge(h1.CDV.packclose.indiv.social, h2.CDV.packclose.indiv.social, by = c('wolfID'), all.x = T)
logRSS.CDV.packclose.indiv.social[,'rss'] <- logRSS.CDV.packclose.indiv.social$h1 - logRSS.CDV.packclose.indiv.social$h2

logRSS.CDV.packfar.indiv.social <- merge(h1.CDV.packfar.indiv.social, h2.CDV.packfar.indiv.social, by = c('wolfID'), all.x = T)
logRSS.CDV.packfar.indiv.social[,'rss'] <- logRSS.CDV.packfar.indiv.social$h1 - logRSS.CDV.packfar.indiv.social$h2


h1.CDV.packclose.indiv.pack <- rbindlist(p.CDV.packclose.1.indiv.pack)
h2.CDV.packclose.indiv.pack <- rbindlist(p.CDV.packclose.2.indiv.pack)

h1.CDV.packfar.indiv.pack <- rbindlist(p.CDV.packfar.1.indiv.pack)
h2.CDV.packfar.indiv.pack <- rbindlist(p.CDV.packfar.2.indiv.pack)

logRSS.CDV.packclose.indiv.pack <- merge(h1.CDV.packclose.indiv.pack, h2.CDV.packclose.indiv.pack, by = c('wolfID'), all.x = T)
logRSS.CDV.packclose.indiv.pack[,'rss'] <- logRSS.CDV.packclose.indiv.pack$h1 - logRSS.CDV.packclose.indiv.pack$h2

logRSS.CDV.packfar.indiv.pack <- merge(h1.CDV.packfar.indiv.pack, h2.CDV.packfar.indiv.pack, by = c('wolfID'), all.x = T)
logRSS.CDV.packfar.indiv.pack[,'rss'] <- logRSS.CDV.packfar.indiv.pack$h1 - logRSS.CDV.packfar.indiv.pack$h2



### human ####
# forest
h1.human.forest75.indiv <- rbindlist(p.human.forest75.1.indiv)
h2.human.forest75.indiv <- rbindlist(p.human.forest75.2.indiv)

h1.human.forest25.indiv <- rbindlist(p.human.forest25.1.indiv)
h2.human.forest25.indiv <- rbindlist(p.human.forest25.2.indiv)

logRSS.human.forest75.indiv <- merge(h1.human.forest75.indiv, h2.human.forest75.indiv, by = c('wolfID'), all.x = T)
logRSS.human.forest75.indiv[,'rss'] <- logRSS.human.forest75.indiv$h1 - logRSS.human.forest75.indiv$h2

logRSS.human.forest25.indiv <- merge(h1.human.forest25.indiv, h2.human.forest25.indiv, by = c('wolfID'), all.x = T)
logRSS.human.forest25.indiv[,'rss'] <- logRSS.human.forest25.indiv$h1 - logRSS.human.forest25.indiv$h2


h1.human.forest75.indiv.habitat <- rbindlist(p.human.forest75.1.indiv.habitat)
h2.human.forest75.indiv.habitat <- rbindlist(p.human.forest75.2.indiv.habitat)

h1.human.forest25.indiv.habitat <- rbindlist(p.human.forest25.1.indiv.habitat)
h2.human.forest25.indiv.habitat <- rbindlist(p.human.forest25.2.indiv.habitat)

logRSS.human.forest75.indiv.habitat <- merge(h1.human.forest75.indiv.habitat, h2.human.forest75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.forest75.indiv.habitat[,'rss'] <- logRSS.human.forest75.indiv.habitat$h1 - logRSS.human.forest75.indiv.habitat$h2

logRSS.human.forest25.indiv.habitat <- merge(h1.human.forest25.indiv.habitat, h2.human.forest25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.forest25.indiv.habitat[,'rss'] <- logRSS.human.forest25.indiv.habitat$h1 - logRSS.human.forest25.indiv.habitat$h2


# open
h1.human.open75.indiv <- rbindlist(p.human.open75.1.indiv)
h2.human.open75.indiv <- rbindlist(p.human.open75.2.indiv)

h1.human.open25.indiv <- rbindlist(p.human.open25.1.indiv)
h2.human.open25.indiv <- rbindlist(p.human.open25.2.indiv)

logRSS.human.open75.indiv <- merge(h1.human.open75.indiv, h2.human.open75.indiv, by = c('wolfID'), all.x = T)
logRSS.human.open75.indiv[,'rss'] <- logRSS.human.open75.indiv$h1 - logRSS.human.open75.indiv$h2

logRSS.human.open25.indiv <- merge(h1.human.open25.indiv, h2.human.open25.indiv, by = c('wolfID'), all.x = T)
logRSS.human.open25.indiv[,'rss'] <- logRSS.human.open25.indiv$h1 - logRSS.human.open25.indiv$h2


h1.human.open75.indiv.habitat <- rbindlist(p.human.open75.1.indiv.habitat)
h2.human.open75.indiv.habitat <- rbindlist(p.human.open75.2.indiv.habitat)

h1.human.open25.indiv.habitat <- rbindlist(p.human.open25.1.indiv.habitat)
h2.human.open25.indiv.habitat <- rbindlist(p.human.open25.2.indiv.habitat)

logRSS.human.open75.indiv.habitat <- merge(h1.human.open75.indiv.habitat, h2.human.open75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.open75.indiv.habitat[,'rss'] <- logRSS.human.open75.indiv.habitat$h1 - logRSS.human.open75.indiv.habitat$h2

logRSS.human.open25.indiv.habitat <- merge(h1.human.open25.indiv.habitat, h2.human.open25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.open25.indiv.habitat[,'rss'] <- logRSS.human.open25.indiv.habitat$h1 - logRSS.human.open25.indiv.habitat$h2


# wet
h1.human.wet75.indiv <- rbindlist(p.human.wet75.1.indiv)
h2.human.wet75.indiv <- rbindlist(p.human.wet75.2.indiv)

h1.human.wet25.indiv <- rbindlist(p.human.wet25.1.indiv)
h2.human.wet25.indiv <- rbindlist(p.human.wet25.2.indiv)

logRSS.human.wet75.indiv <- merge(h1.human.wet75.indiv, h2.human.wet75.indiv, by = c('wolfID'), all.x = T)
logRSS.human.wet75.indiv[,'rss'] <- logRSS.human.wet75.indiv$h1 - logRSS.human.wet75.indiv$h2

logRSS.human.wet25.indiv <- merge(h1.human.wet25.indiv, h2.human.wet25.indiv, by = c('wolfID'), all.x = T)
logRSS.human.wet25.indiv[,'rss'] <- logRSS.human.wet25.indiv$h1 - logRSS.human.wet25.indiv$h2


h1.human.wet75.indiv.habitat <- rbindlist(p.human.wet75.1.indiv.habitat)
h2.human.wet75.indiv.habitat <- rbindlist(p.human.wet75.2.indiv.habitat)

h1.human.wet25.indiv.habitat <- rbindlist(p.human.wet25.1.indiv.habitat)
h2.human.wet25.indiv.habitat <- rbindlist(p.human.wet25.2.indiv.habitat)

logRSS.human.wet75.indiv.habitat <- merge(h1.human.wet75.indiv.habitat, h2.human.wet75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.wet75.indiv.habitat[,'rss'] <- logRSS.human.wet75.indiv.habitat$h1 - logRSS.human.wet75.indiv.habitat$h2

logRSS.human.wet25.indiv.habitat <- merge(h1.human.wet25.indiv.habitat, h2.human.wet25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.wet25.indiv.habitat[,'rss'] <- logRSS.human.wet25.indiv.habitat$h1 - logRSS.human.wet25.indiv.habitat$h2


# road
h1.human.roadclose.indiv <- rbindlist(p.human.roadclose.1.indiv)
h2.human.roadclose.indiv <- rbindlist(p.human.roadclose.2.indiv)

h1.human.roadfar.indiv <- rbindlist(p.human.roadfar.1.indiv)
h2.human.roadfar.indiv <- rbindlist(p.human.roadfar.2.indiv)

logRSS.human.roadclose.indiv <- merge(h1.human.roadclose.indiv, h2.human.roadclose.indiv, by = c('wolfID'), all.x = T)
logRSS.human.roadclose.indiv[,'rss'] <- logRSS.human.roadclose.indiv$h1 - logRSS.human.roadclose.indiv$h2

logRSS.human.roadfar.indiv <- merge(h1.human.roadfar.indiv, h2.human.roadfar.indiv, by = c('wolfID'), all.x = T)
logRSS.human.roadfar.indiv[,'rss'] <- logRSS.human.roadfar.indiv$h1 - logRSS.human.roadfar.indiv$h2


h1.human.roadclose.indiv.habitat <- rbindlist(p.human.roadclose.1.indiv.habitat)
h2.human.roadclose.indiv.habitat <- rbindlist(p.human.roadclose.2.indiv.habitat)

h1.human.roadfar.indiv.habitat <- rbindlist(p.human.roadfar.1.indiv.habitat)
h2.human.roadfar.indiv.habitat <- rbindlist(p.human.roadfar.2.indiv.habitat)

logRSS.human.roadclose.indiv.habitat <- merge(h1.human.roadclose.indiv.habitat, h2.human.roadclose.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.roadclose.indiv.habitat[,'rss'] <- logRSS.human.roadclose.indiv.habitat$h1 - logRSS.human.roadclose.indiv.habitat$h2

logRSS.human.roadfar.indiv.habitat <- merge(h1.human.roadfar.indiv.habitat, h2.human.roadfar.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.human.roadfar.indiv.habitat[,'rss'] <- logRSS.human.roadfar.indiv.habitat$h1 - logRSS.human.roadfar.indiv.habitat$h2


# nn
h1.human.nnclose.indiv <- rbindlist(p.human.nnclose.1.indiv)
h2.human.nnclose.indiv <- rbindlist(p.human.nnclose.2.indiv)

h1.human.nnfar.indiv <- rbindlist(p.human.nnfar.1.indiv)
h2.human.nnfar.indiv <- rbindlist(p.human.nnfar.2.indiv)

logRSS.human.nnclose.indiv <- merge(h1.human.nnclose.indiv, h2.human.nnclose.indiv, by = c('wolfID'), all.x = T)
logRSS.human.nnclose.indiv[,'rss'] <- logRSS.human.nnclose.indiv$h1 - logRSS.human.nnclose.indiv$h2

logRSS.human.nnfar.indiv <- merge(h1.human.nnfar.indiv, h2.human.nnfar.indiv, by = c('wolfID'), all.x = T)
logRSS.human.nnfar.indiv[,'rss'] <- logRSS.human.nnfar.indiv$h1 - logRSS.human.nnfar.indiv$h2


h1.human.nnclose.indiv.social <- rbindlist(p.human.nnclose.1.indiv.social)
h2.human.nnclose.indiv.social <- rbindlist(p.human.nnclose.2.indiv.social)

h1.human.nnfar.indiv.social <- rbindlist(p.human.nnfar.1.indiv.social)
h2.human.nnfar.indiv.social <- rbindlist(p.human.nnfar.2.indiv.social)

logRSS.human.nnclose.indiv.social <- merge(h1.human.nnclose.indiv.social, h2.human.nnclose.indiv.social, by = c('wolfID'), all.x = T)
logRSS.human.nnclose.indiv.social[,'rss'] <- logRSS.human.nnclose.indiv.social$h1 - logRSS.human.nnclose.indiv.social$h2

logRSS.human.nnfar.indiv.social <- merge(h1.human.nnfar.indiv.social, h2.human.nnfar.indiv.social, by = c('wolfID'), all.x = T)
logRSS.human.nnfar.indiv.social[,'rss'] <- logRSS.human.nnfar.indiv.social$h1 - logRSS.human.nnfar.indiv.social$h2


# pack
h1.human.packclose.indiv <- rbindlist(p.human.packclose.1.indiv)
h2.human.packclose.indiv <- rbindlist(p.human.packclose.2.indiv)

h1.human.packfar.indiv <- rbindlist(p.human.packfar.1.indiv)
h2.human.packfar.indiv <- rbindlist(p.human.packfar.2.indiv)

logRSS.human.packclose.indiv <- merge(h1.human.packclose.indiv, h2.human.packclose.indiv, by = c('wolfID'), all.x = T)
logRSS.human.packclose.indiv[,'rss'] <- logRSS.human.packclose.indiv$h1 - logRSS.human.packclose.indiv$h2

logRSS.human.packfar.indiv <- merge(h1.human.packfar.indiv, h2.human.packfar.indiv, by = c('wolfID'), all.x = T)
logRSS.human.packfar.indiv[,'rss'] <- logRSS.human.packfar.indiv$h1 - logRSS.human.packfar.indiv$h2


h1.human.packclose.indiv.social <- rbindlist(p.human.packclose.1.indiv.social)
h2.human.packclose.indiv.social <- rbindlist(p.human.packclose.2.indiv.social)

h1.human.packfar.indiv.social <- rbindlist(p.human.packfar.1.indiv.social)
h2.human.packfar.indiv.social <- rbindlist(p.human.packfar.2.indiv.social)

logRSS.human.packclose.indiv.social <- merge(h1.human.packclose.indiv.social, h2.human.packclose.indiv.social, by = c('wolfID'), all.x = T)
logRSS.human.packclose.indiv.social[,'rss'] <- logRSS.human.packclose.indiv.social$h1 - logRSS.human.packclose.indiv.social$h2

logRSS.human.packfar.indiv.social <- merge(h1.human.packfar.indiv.social, h2.human.packfar.indiv.social, by = c('wolfID'), all.x = T)
logRSS.human.packfar.indiv.social[,'rss'] <- logRSS.human.packfar.indiv.social$h1 - logRSS.human.packfar.indiv.social$h2


h1.human.packclose.indiv.pack <- rbindlist(p.human.packclose.1.indiv.pack)
h2.human.packclose.indiv.pack <- rbindlist(p.human.packclose.2.indiv.pack)

h1.human.packfar.indiv.pack <- rbindlist(p.human.packfar.1.indiv.pack)
h2.human.packfar.indiv.pack <- rbindlist(p.human.packfar.2.indiv.pack)

logRSS.human.packclose.indiv.pack <- merge(h1.human.packclose.indiv.pack, h2.human.packclose.indiv.pack, by = c('wolfID'), all.x = T)
logRSS.human.packclose.indiv.pack[,'rss'] <- logRSS.human.packclose.indiv.pack$h1 - logRSS.human.packclose.indiv.pack$h2

logRSS.human.packfar.indiv.pack <- merge(h1.human.packfar.indiv.pack, h2.human.packfar.indiv.pack, by = c('wolfID'), all.x = T)
logRSS.human.packfar.indiv.pack[,'rss'] <- logRSS.human.packfar.indiv.pack$h1 - logRSS.human.packfar.indiv.pack$h2



### control ####
# forest
h1.control.forest75.indiv <- rbindlist(p.control.forest75.1.indiv)
h2.control.forest75.indiv <- rbindlist(p.control.forest75.2.indiv)

h1.control.forest25.indiv <- rbindlist(p.control.forest25.1.indiv)
h2.control.forest25.indiv <- rbindlist(p.control.forest25.2.indiv)

logRSS.control.forest75.indiv <- merge(h1.control.forest75.indiv, h2.control.forest75.indiv, by = c('wolfID'), all.x = T)
logRSS.control.forest75.indiv[,'rss'] <- logRSS.control.forest75.indiv$h1 - logRSS.control.forest75.indiv$h2

logRSS.control.forest25.indiv <- merge(h1.control.forest25.indiv, h2.control.forest25.indiv, by = c('wolfID'), all.x = T)
logRSS.control.forest25.indiv[,'rss'] <- logRSS.control.forest25.indiv$h1 - logRSS.control.forest25.indiv$h2


h1.control.forest75.indiv.habitat <- rbindlist(p.control.forest75.1.indiv.habitat)
h2.control.forest75.indiv.habitat <- rbindlist(p.control.forest75.2.indiv.habitat)

h1.control.forest25.indiv.habitat <- rbindlist(p.control.forest25.1.indiv.habitat)
h2.control.forest25.indiv.habitat <- rbindlist(p.control.forest25.2.indiv.habitat)

logRSS.control.forest75.indiv.habitat <- merge(h1.control.forest75.indiv.habitat, h2.control.forest75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.forest75.indiv.habitat[,'rss'] <- logRSS.control.forest75.indiv.habitat$h1 - logRSS.control.forest75.indiv.habitat$h2

logRSS.control.forest25.indiv.habitat <- merge(h1.control.forest25.indiv.habitat, h2.control.forest25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.forest25.indiv.habitat[,'rss'] <- logRSS.control.forest25.indiv.habitat$h1 - logRSS.control.forest25.indiv.habitat$h2


# open
h1.control.open75.indiv <- rbindlist(p.control.open75.1.indiv)
h2.control.open75.indiv <- rbindlist(p.control.open75.2.indiv)

h1.control.open25.indiv <- rbindlist(p.control.open25.1.indiv)
h2.control.open25.indiv <- rbindlist(p.control.open25.2.indiv)

logRSS.control.open75.indiv <- merge(h1.control.open75.indiv, h2.control.open75.indiv, by = c('wolfID'), all.x = T)
logRSS.control.open75.indiv[,'rss'] <- logRSS.control.open75.indiv$h1 - logRSS.control.open75.indiv$h2

logRSS.control.open25.indiv <- merge(h1.control.open25.indiv, h2.control.open25.indiv, by = c('wolfID'), all.x = T)
logRSS.control.open25.indiv[,'rss'] <- logRSS.control.open25.indiv$h1 - logRSS.control.open25.indiv$h2


h1.control.open75.indiv.habitat <- rbindlist(p.control.open75.1.indiv.habitat)
h2.control.open75.indiv.habitat <- rbindlist(p.control.open75.2.indiv.habitat)

h1.control.open25.indiv.habitat <- rbindlist(p.control.open25.1.indiv.habitat)
h2.control.open25.indiv.habitat <- rbindlist(p.control.open25.2.indiv.habitat)

logRSS.control.open75.indiv.habitat <- merge(h1.control.open75.indiv.habitat, h2.control.open75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.open75.indiv.habitat[,'rss'] <- logRSS.control.open75.indiv.habitat$h1 - logRSS.control.open75.indiv.habitat$h2

logRSS.control.open25.indiv.habitat <- merge(h1.control.open25.indiv.habitat, h2.control.open25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.open25.indiv.habitat[,'rss'] <- logRSS.control.open25.indiv.habitat$h1 - logRSS.control.open25.indiv.habitat$h2


# wet
h1.control.wet75.indiv <- rbindlist(p.control.wet75.1.indiv)
h2.control.wet75.indiv <- rbindlist(p.control.wet75.2.indiv)

h1.control.wet25.indiv <- rbindlist(p.control.wet25.1.indiv)
h2.control.wet25.indiv <- rbindlist(p.control.wet25.2.indiv)

logRSS.control.wet75.indiv <- merge(h1.control.wet75.indiv, h2.control.wet75.indiv, by = c('wolfID'), all.x = T)
logRSS.control.wet75.indiv[,'rss'] <- logRSS.control.wet75.indiv$h1 - logRSS.control.wet75.indiv$h2

logRSS.control.wet25.indiv <- merge(h1.control.wet25.indiv, h2.control.wet25.indiv, by = c('wolfID'), all.x = T)
logRSS.control.wet25.indiv[,'rss'] <- logRSS.control.wet25.indiv$h1 - logRSS.control.wet25.indiv$h2


h1.control.wet75.indiv.habitat <- rbindlist(p.control.wet75.1.indiv.habitat)
h2.control.wet75.indiv.habitat <- rbindlist(p.control.wet75.2.indiv.habitat)

h1.control.wet25.indiv.habitat <- rbindlist(p.control.wet25.1.indiv.habitat)
h2.control.wet25.indiv.habitat <- rbindlist(p.control.wet25.2.indiv.habitat)

logRSS.control.wet75.indiv.habitat <- merge(h1.control.wet75.indiv.habitat, h2.control.wet75.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.wet75.indiv.habitat[,'rss'] <- logRSS.control.wet75.indiv.habitat$h1 - logRSS.control.wet75.indiv.habitat$h2

logRSS.control.wet25.indiv.habitat <- merge(h1.control.wet25.indiv.habitat, h2.control.wet25.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.wet25.indiv.habitat[,'rss'] <- logRSS.control.wet25.indiv.habitat$h1 - logRSS.control.wet25.indiv.habitat$h2


# road
h1.control.roadclose.indiv <- rbindlist(p.control.roadclose.1.indiv)
h2.control.roadclose.indiv <- rbindlist(p.control.roadclose.2.indiv)

h1.control.roadfar.indiv <- rbindlist(p.control.roadfar.1.indiv)
h2.control.roadfar.indiv <- rbindlist(p.control.roadfar.2.indiv)

logRSS.control.roadclose.indiv <- merge(h1.control.roadclose.indiv, h2.control.roadclose.indiv, by = c('wolfID'), all.x = T)
logRSS.control.roadclose.indiv[,'rss'] <- logRSS.control.roadclose.indiv$h1 - logRSS.control.roadclose.indiv$h2

logRSS.control.roadfar.indiv <- merge(h1.control.roadfar.indiv, h2.control.roadfar.indiv, by = c('wolfID'), all.x = T)
logRSS.control.roadfar.indiv[,'rss'] <- logRSS.control.roadfar.indiv$h1 - logRSS.control.roadfar.indiv$h2


h1.control.roadclose.indiv.habitat <- rbindlist(p.control.roadclose.1.indiv.habitat)
h2.control.roadclose.indiv.habitat <- rbindlist(p.control.roadclose.2.indiv.habitat)

h1.control.roadfar.indiv.habitat <- rbindlist(p.control.roadfar.1.indiv.habitat)
h2.control.roadfar.indiv.habitat <- rbindlist(p.control.roadfar.2.indiv.habitat)

logRSS.control.roadclose.indiv.habitat <- merge(h1.control.roadclose.indiv.habitat, h2.control.roadclose.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.roadclose.indiv.habitat[,'rss'] <- logRSS.control.roadclose.indiv.habitat$h1 - logRSS.control.roadclose.indiv.habitat$h2

logRSS.control.roadfar.indiv.habitat <- merge(h1.control.roadfar.indiv.habitat, h2.control.roadfar.indiv.habitat, by = c('wolfID'), all.x = T)
logRSS.control.roadfar.indiv.habitat[,'rss'] <- logRSS.control.roadfar.indiv.habitat$h1 - logRSS.control.roadfar.indiv.habitat$h2


# nn
h1.control.nnclose.indiv <- rbindlist(p.control.nnclose.1.indiv)
h2.control.nnclose.indiv <- rbindlist(p.control.nnclose.2.indiv)

h1.control.nnfar.indiv <- rbindlist(p.control.nnfar.1.indiv)
h2.control.nnfar.indiv <- rbindlist(p.control.nnfar.2.indiv)

logRSS.control.nnclose.indiv <- merge(h1.control.nnclose.indiv, h2.control.nnclose.indiv, by = c('wolfID'), all.x = T)
logRSS.control.nnclose.indiv[,'rss'] <- logRSS.control.nnclose.indiv$h1 - logRSS.control.nnclose.indiv$h2

logRSS.control.nnfar.indiv <- merge(h1.control.nnfar.indiv, h2.control.nnfar.indiv, by = c('wolfID'), all.x = T)
logRSS.control.nnfar.indiv[,'rss'] <- logRSS.control.nnfar.indiv$h1 - logRSS.control.nnfar.indiv$h2


h1.control.nnclose.indiv.social <- rbindlist(p.control.nnclose.1.indiv.social)
h2.control.nnclose.indiv.social <- rbindlist(p.control.nnclose.2.indiv.social)

h1.control.nnfar.indiv.social <- rbindlist(p.control.nnfar.1.indiv.social)
h2.control.nnfar.indiv.social <- rbindlist(p.control.nnfar.2.indiv.social)

logRSS.control.nnclose.indiv.social <- merge(h1.control.nnclose.indiv.social, h2.control.nnclose.indiv.social, by = c('wolfID'), all.x = T)
logRSS.control.nnclose.indiv.social[,'rss'] <- logRSS.control.nnclose.indiv.social$h1 - logRSS.control.nnclose.indiv.social$h2

logRSS.control.nnfar.indiv.social <- merge(h1.control.nnfar.indiv.social, h2.control.nnfar.indiv.social, by = c('wolfID'), all.x = T)
logRSS.control.nnfar.indiv.social[,'rss'] <- logRSS.control.nnfar.indiv.social$h1 - logRSS.control.nnfar.indiv.social$h2


# pack
h1.control.packclose.indiv <- rbindlist(p.control.packclose.1.indiv)
h2.control.packclose.indiv <- rbindlist(p.control.packclose.2.indiv)

h1.control.packfar.indiv <- rbindlist(p.control.packfar.1.indiv)
h2.control.packfar.indiv <- rbindlist(p.control.packfar.2.indiv)

logRSS.control.packclose.indiv <- merge(h1.control.packclose.indiv, h2.control.packclose.indiv, by = c('wolfID'), all.x = T)
logRSS.control.packclose.indiv[,'rss'] <- logRSS.control.packclose.indiv$h1 - logRSS.control.packclose.indiv$h2

logRSS.control.packfar.indiv <- merge(h1.control.packfar.indiv, h2.control.packfar.indiv, by = c('wolfID'), all.x = T)
logRSS.control.packfar.indiv[,'rss'] <- logRSS.control.packfar.indiv$h1 - logRSS.control.packfar.indiv$h2


h1.control.packclose.indiv.social <- rbindlist(p.control.packclose.1.indiv.social)
h2.control.packclose.indiv.social <- rbindlist(p.control.packclose.2.indiv.social)

h1.control.packfar.indiv.social <- rbindlist(p.control.packfar.1.indiv.social)
h2.control.packfar.indiv.social <- rbindlist(p.control.packfar.2.indiv.social)

logRSS.control.packclose.indiv.social <- merge(h1.control.packclose.indiv.social, h2.control.packclose.indiv.social, by = c('wolfID'), all.x = T)
logRSS.control.packclose.indiv.social[,'rss'] <- logRSS.control.packclose.indiv.social$h1 - logRSS.control.packclose.indiv.social$h2

logRSS.control.packfar.indiv.social <- merge(h1.control.packfar.indiv.social, h2.control.packfar.indiv.social, by = c('wolfID'), all.x = T)
logRSS.control.packfar.indiv.social[,'rss'] <- logRSS.control.packfar.indiv.social$h1 - logRSS.control.packfar.indiv.social$h2


h1.control.packclose.indiv.pack <- rbindlist(p.control.packclose.1.indiv.pack)
h2.control.packclose.indiv.pack <- rbindlist(p.control.packclose.2.indiv.pack)

h1.control.packfar.indiv.pack <- rbindlist(p.control.packfar.1.indiv.pack)
h2.control.packfar.indiv.pack <- rbindlist(p.control.packfar.2.indiv.pack)

logRSS.control.packclose.indiv.pack <- merge(h1.control.packclose.indiv.pack, h2.control.packclose.indiv.pack, by = c('wolfID'), all.x = T)
logRSS.control.packclose.indiv.pack[,'rss'] <- logRSS.control.packclose.indiv.pack$h1 - logRSS.control.packclose.indiv.pack$h2

logRSS.control.packfar.indiv.pack <- merge(h1.control.packfar.indiv.pack, h2.control.packfar.indiv.pack, by = c('wolfID'), all.x = T)
logRSS.control.packfar.indiv.pack[,'rss'] <- logRSS.control.packfar.indiv.pack$h1 - logRSS.control.packfar.indiv.pack$h2

#### gather RSSs ####
logRSS <- rbind(logRSS.CDV.forest25, logRSS.CDV.forest75, logRSS.CDV.open25, logRSS.CDV.open75, logRSS.CDV.wet25, logRSS.CDV.wet75, 
                logRSS.CDV.roadclose, logRSS.CDV.roadfar, logRSS.CDV.nnclose, logRSS.CDV.nnfar, logRSS.CDV.packclose, logRSS.CDV.packfar,
                logRSS.human.forest25, logRSS.human.forest75, logRSS.human.open25, logRSS.human.open75, logRSS.human.wet25, logRSS.human.wet75, 
                logRSS.human.roadclose, logRSS.human.roadfar, logRSS.human.nnclose, logRSS.human.nnfar, logRSS.human.packclose, logRSS.human.packfar,
                logRSS.control.forest25, logRSS.control.forest75, logRSS.control.open25, logRSS.control.open75, logRSS.control.wet25, logRSS.control.wet75, 
                logRSS.control.roadclose, logRSS.control.roadfar, logRSS.control.nnclose, logRSS.control.nnfar, logRSS.control.packclose, logRSS.control.packfar)

logRSS.indiv <- rbind(logRSS.CDV.forest25.indiv, logRSS.CDV.forest75.indiv, logRSS.CDV.open25.indiv, logRSS.CDV.open75.indiv, logRSS.CDV.wet25.indiv, logRSS.CDV.wet75.indiv, 
                logRSS.CDV.roadclose.indiv, logRSS.CDV.roadfar.indiv, logRSS.CDV.nnclose.indiv, logRSS.CDV.nnfar.indiv, logRSS.CDV.packclose.indiv, logRSS.CDV.packfar.indiv,
                logRSS.human.forest25.indiv, logRSS.human.forest75.indiv, logRSS.human.open25.indiv, logRSS.human.open75.indiv, logRSS.human.wet25.indiv, logRSS.human.wet75.indiv, 
                logRSS.human.roadclose.indiv, logRSS.human.roadfar.indiv, logRSS.human.nnclose.indiv, logRSS.human.nnfar.indiv, logRSS.human.packclose.indiv, logRSS.human.packfar.indiv,
                logRSS.control.forest25.indiv, logRSS.control.forest75.indiv, logRSS.control.open25.indiv, logRSS.control.open75.indiv, logRSS.control.wet25.indiv, logRSS.control.wet75.indiv, 
                logRSS.control.roadclose.indiv, logRSS.control.roadfar.indiv, logRSS.control.nnclose.indiv, logRSS.control.nnfar.indiv, logRSS.control.packclose.indiv, logRSS.control.packfar.indiv)

# logRSS.indiv.model <- rbind(logRSS.CDV.forest25.indiv.habitat, logRSS.CDV.forest75.indiv.habitat, 
#                       logRSS.CDV.open25.indiv.habitat, logRSS.CDV.open75.indiv.habitat, 
#                       logRSS.CDV.wet25.indiv.habitat, logRSS.CDV.wet75.indiv.habitat, 
#                       logRSS.CDV.roadclose.indiv.habitat, logRSS.CDV.roadfar.indiv.habitat, 
#                       logRSS.CDV.nnclose.indiv.social, logRSS.CDV.nnfar.indiv.social, 
#                       logRSS.CDV.packclose.indiv, logRSS.CDV.packfar.indiv) ## CHECK WORKING MODEL

#### graph colors ####

cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")
gcolors <- c("deepskyblue", "purple", "dark green")

logRSS.indiv$COD <- factor(logRSS.indiv$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))
#### GRAPHS ####
forest.75 <- ggplot(data=setDT(logRSS.indiv)[var == 'forest'& value ==0.25], aes(-ttd, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("High proportion of forest") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-10,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
forest.75


open.75 <- ggplot(data=setDT(logRSS.indiv)[var == 'open'& value ==0.25], aes(-ttd, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("High proportion of open") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-10,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
open.75

wet.75 <- ggplot(data=setDT(logRSS.indiv)[var == 'wet'& value ==0.25], aes(-ttd, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("High proportion of wet") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-10,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75

road.close <- ggplot(data=setDT(logRSS.indiv)[var == 'road'& value ==250], aes(-ttd, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("Close to road") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-10,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
road.close

nn.close <- ggplot(data=setDT(logRSS.indiv)[var == 'nn'& value ==250], aes(-ttd, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("Close to nearest neighbor") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-10,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.close

pack.close <- ggplot(data=setDT(logRSS.indiv)[var == 'pack'& value ==250], aes(-ttd, rss, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'pack'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("Close to pack boundary") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #ylim(-10,8) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
pack.close

