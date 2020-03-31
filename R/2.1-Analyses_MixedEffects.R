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

summary(dat$land_end)

# dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet', 
#                                ifelse(dat$land_end == 'mixed'|dat$land_end == 'deciduous', 'forest',
#                                       ifelse(dat$land_end == 'coniferous', 'coniferous', 'open')))
dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet',
                               ifelse(dat$land_end == 'mixed'|dat$land_end == 'deciduous'|dat$land_end == 'coniferous', 'forest','open'))

#dat[,'forest']
dat[,'propforest_end_adj'] <- dat$propconif_end+dat$propmixed_end +dat$propdecid_end
dat[,'propopen_end_adj'] <- dat$propopen_end + dat$propurban_end #+dat$propshrub_end 

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


######
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


dat[COD=='cdv',.(last(packDistadj_end), min(packDistadj_end), max(packDistadj_end)), by=.(wolfID)]


#### full model #### 
#### all time ####
dat[,uniqueN(step_id_), by=.(wolfID)]
dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo <- dat[ttd1<=31 & !is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo.cod <- merge(dat.wnn.lastmo, dat.meta[,.(wolfpop, pop, COD)], by.x = 'wolfID', by.y = 'wolfpop', all.x = T)

dat.wnn.lastmo.cod[,.(N=uniqueN(wolfID)), by=.(pop, COD)]
dat.meta[,.(N=uniqueN(wolfpop)), by=.(pop, COD)]

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
                      log_sl:ToD_start +
                       # log_sl:land_end_adj +
                      log_sl:propforest_end_adj + log_sl:propopen_end_adj + log_sl:propwet_end +
                        log_sl:COD + cos_ta:COD + 
                        I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                        (1|wolf_step_id) +
                        (0 + (log_sl)|wolfID) +
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
# everyone.ran_pars <-broom.mixed::tidy(everyone, effect= 'ran_pars')
everyone.se <-setDT(everyone.ran_vals)[group=='wolfID']


# everyone.all.indiv <- coef(everyone)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
#   pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
#   mutate(method = "ME")



#### by model ####
everyone.move <- glmmTMB(case_ ~ log_sl:ToD_start +
                              log_sl:land_end_adj +
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

AICtab(everyone.move, everyone.habitat, everyone, everyone.social)


#### GATHERING RESULTS ####

cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")


#### everyone ####
everyone.indiv<- merge(everyone.se[,.(wolfID= level, term, estimate, se=std.error)], dat.meta[,.(wolfpop, COD)], by.x ='wolfID', by.y= 'wolfpop', all.x=T)
everyone.indiv$COD <- factor(everyone.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))
everyone.indiv$COD[is.na(everyone.indiv$COD)] <- "control"

unique(everyone.indiv$term)
everyone.all.betas.names <- c("log_sl", "cos_ta", 
                             "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "I(log(1 + roadDist_end))",
                             "I(log(1 + distance2))", "I(log(1 + packDistadj_end))",
                             "I(log(ttd1 + 1)):log_sl", "I(log(ttd1 + 1)):cos_ta", 
                             "I(log(ttd1 + 1)):land_end_adjforest", "I(log(ttd1 + 1)):land_end_adjopen", "I(log(ttd1 + 1)):land_end_adjwet", "I(log(ttd1 + 1)):I(log(1 + roadDist_end))",
                             "I(log(ttd1 + 1)):I(log(1 + distance2))", "I(log(ttd1 + 1)):I(log(1 + packDistadj_end))")
everyone.indiv$term <- factor(everyone.indiv$term, levels = everyone.all.betas.names, labels = c("log_sl", "cos_ta",'forest', "open", "wet", "roadDist",
                                                                                                                     "nnDist", "boundaryDist",
                                                                                                                     "log_sl-ttd", "cos_ta-ttd", "forest-ttd", "open-ttd", "wet-ttd", "roadDist-ttd",
                                                                                                                     "nnDist-ttd", "boundaryDist-ttd"))


#saveRDS(everyone.indiv, 'data/derived-data/everyone_betas_COD.Rds')
#everyone.indiv<-readRDS('data/derived-data/everyone_betas_lastmoNN.Rds')


everyone.ttd <- everyone.indiv[term %like% "ttd", ]

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


everyone.main <- everyone.indiv[!(term %like% "ttd"), ]

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
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-.35,.3)

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
  scale_color_manual(values = cbPalette) + ylim(-8,12)

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
  scale_color_manual(values = cbPalette) + ylim(-1,2.5)




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
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + ylim(-.35,.3)

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
  scale_color_manual(values = cbPalette) + ylim(-8,12)

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
  scale_color_manual(values = cbPalette) + ylim(-1,2.5)


main.move/ttd.move
main.hab/ttd.hab
main.dist/ttd.dist

##### last mo NEED TO UPDATE NAMES ####
everyone.lastmo.all.indiv<- merge(everyone.lastmo.all.indiv, dat.meta[,.(wolfpop, COD)], by.x ='wolfID', by.y= 'wolfpop', all.x=T)
everyone.lastmo.all.indiv$COD <- factor(everyone.lastmo.all.indiv$COD, levels = c('none','human','cdv'), labels = c('control','human','CDV'))

minmax <- setDT(everyone.lastmo.all.indiv)[,.(min= min(estimate), max=max(estimate)), by = .(term, COD)]
everyone.lastmo.all.betas <- minmax[min!= max]
everyone.lastmo.all.betas.names <- unique(everyone.lastmo.all.betas$term)

everyone.lastmo.all.betas.names <- c("log_sl", "cos_ta", 
                              "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
                              "log(1 + distance2)", "log(1 + packDistadj_end)",
                              "log(ttd1 + 1):log_sl", "log(ttd1 + 1):cos_ta", 
                              "log(ttd1 + 1):land_end_adjforest", "log(ttd1 + 1):land_end_adjopen", "log(ttd1 + 1):land_end_adjwet", "log(ttd1 + 1):log(1 + roadDist_end)", 
                              "log(ttd1 + 1):log(1 + distance2)", "log(ttd1 + 1):log(1 + packDistadj_end)")

everyone.lastmo.all.indiv.betas <- everyone.lastmo.all.indiv[term %chin% everyone.lastmo.all.betas.names]
# everyone.lastmo.betas <- c("log_sl", "cos_ta", "land_end_adjforest", "land_end_adjopen", "land_end_adjwet", "log(1 + roadDist_end)",
#                 "log(1 + distance2)", "log(1 + packDistadj_end)")

everyone.lastmo.all.indiv.betas$term <- factor(everyone.lastmo.all.indiv.betas$term, levels = everyone.lastmo.all.betas.names, labels = c("log_sl", "cos_ta",'forest', "open", "wet", "roadDist",
                                                                                                                     "nnDist", "boundaryDist",
                                                                                                                     "log_sl-ttd", "cos_ta-ttd", "forest-ttd", "open-ttd", "wet-ttd", "roadDist-ttd",
                                                                                                                     "nnDist-ttd", "boundaryDist-ttd"))

everyone.lastmo.all.indiv.betas$COD[is.na(everyone.lastmo.all.indiv.betas$COD)] <- "control"

#saveRDS(everyone.lastmo.all.indiv.betas, 'data/derived-data/everyone_lastmo_betas.Rds')

everyone.lastmo.ttd <- everyone.lastmo.all.indiv.betas[term %like% "ttd", ]

ttd.vars.lastmo <- ggplot(everyone.lastmo.ttd, aes(term, (estimate), fill = COD)) +
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


ggplot(everyone.lastmo.ttd[term== 'wet-ttd' | term== 'boundaryDist-ttd'], aes(term, (estimate), fill = COD)) +
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


everyone.lastmo.main <- everyone.lastmo.all.indiv.betas[!(term %like% "ttd"), ]

main.vars.lastmo <- ggplot(everyone.lastmo.main, aes(term, (estimate), fill = COD)) +
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

ttd.vars.lastmo
main.vars.lastmo



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
                                                           wolf_step_id = NA, wolfID= 'RMNP_W02')]
p.CDV.1day.2 <- predict(everyone, newdata = CDV.1day.2, type='link', re.form = NA)
p.CDV.1day.2.habitat <- predict(everyone.habitat, newdata = CDV.1day.2, type='link', re.form = NA)
p.CDV.1day.2.social <- predict(everyone.social, newdata = CDV.1day.2, type='link', re.form = NA)



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
                                                        wolf_step_id = NA, wolfID= 'RMNP_W02')]
p.CDV.60day.2 <- predict(everyone, newdata = CDV.60day.2, type='link', re.form = NA)
p.CDV.60day.2.habitat <- predict(everyone.habitat, newdata = CDV.60day.2, type='link', re.form = NA)
p.CDV.60day.2.social <- predict(everyone.social, newdata = CDV.60day.2, type='link', re.form = NA)


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
                                                        wolf_step_id = NA, wolfID= 'RMNP_W02')]
p.human.1day.2 <- predict(everyone, newdata = human.1day.2, type='link', re.form = NA)
p.human.1day.2.habitat <- predict(everyone.habitat, newdata = human.1day.2, type='link', re.form = NA)
p.human.1day.2.social <- predict(everyone.social, newdata = human.1day.2, type='link', re.form = NA)



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
                                                         wolf_step_id = NA, wolfID= 'RMNP_W02')]
p.human.60day.2 <- predict(everyone, newdata = human.60day.2, type='link', re.form = NA)
p.human.60day.2.habitat <- predict(everyone.habitat, newdata = human.60day.2, type='link', re.form = NA)
p.human.60day.2.social <- predict(everyone.social, newdata = human.60day.2, type='link', re.form = NA)


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
                                                        wolf_step_id = NA, wolfID= 'RMNP_W02')]
p.control.1day.2 <- predict(everyone, newdata = control.1day.2, type='link', re.form = NA)
p.control.1day.2.habitat <- predict(everyone.habitat, newdata = control.1day.2, type='link', re.form = NA)
p.control.1day.2.social <- predict(everyone.social, newdata = control.1day.2, type='link', re.form = NA)


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
                                                         wolf_step_id = NA, wolfID= 'RMNP_W02')]
p.control.60day.2 <- predict(everyone, newdata = control.60day.2, type='link', re.form = NA)
p.control.60day.2.habitat <- predict(everyone.habitat, newdata = control.60day.2, type='link', re.form = NA)
p.control.60day.2.social <- predict(everyone.social, newdata = control.60day.2, type='link', re.form = NA)


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
  
p.CDV.60day.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 60)
p.CDV.60day.2.indiv.habitat <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.habitat, death = 'CDV', t2death = 60)
p.CDV.60day.2.indiv.social <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone.social, death = 'CDV', t2death = 60)


### human ###
human.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='human', wolfID])

p.human.1day.2.indiv <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 1)
p.human.1day.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 1)
p.human.1day.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 1)

p.human.60day.2.indiv <- p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone, death = 'human', t2death = 60)
p.human.60day.2.indiv.habitat <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.habitat, death = 'human', t2death = 60)
p.human.60day.2.indiv.social <-p.h2.indiv(ids = human.wolfID, DT = dat, mod = everyone.social, death = 'human', t2death = 60)


### control ###
control.wolfID <- unique(dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD=='control', wolfID])

p.control.1day.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 1)
p.control.1day.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 1)
p.control.1day.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 1)

p.control.60day.2.indiv <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone, death = 'control', t2death = 60)
p.control.60day.2.indiv.habitat <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.habitat, death = 'control', t2death = 60)
p.control.60day.2.indiv.social <-p.h2.indiv(ids = control.wolfID, DT = dat, mod = everyone.social, death = 'control', t2death = 60)




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
                                                             wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                  wolf_step_id = NA, wolfID = 'RMNP_W02')]

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





##### graphs ####
logRSS.forest.indiv <- setDT(logRSS.forest.indiv)
logRSS.forest.indiv[,'COD'] <- as.factor(logRSS.forest.indiv$COD)
logRSS.forest.indiv[,'ttd'] <- as.factor(logRSS.forest.indiv$ttd)

logRSS.forest.indiv.se <- unique(logRSS.forest.indiv[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.forest.pop <- merge(logRSS.forest, logRSS.forest.indiv.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.forest.pop)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Forest") +
  ggtitle("1 days") 
# 
# ggplot(data=setDT(logRSS.forest.indiv)[ttd=='1 day'], aes(x, rss, colour=COD)) +
#   #geom_line() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
#   geom_smooth()

ggplot(data=setDT(logRSS.forest.pop)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2)) +
  ylab("logRSS") + xlab("Forest") +
  ggtitle("60 days") 



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
                                                               wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                  wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                    wolf_step_id = NA, wolfID = 'RMNP_W02')]

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





##### graphs ####
logRSS.open.indiv <- setDT(logRSS.open.indiv)
logRSS.open.indiv[,'COD'] <- as.factor(logRSS.open.indiv$COD)
logRSS.open.indiv[,'ttd'] <- as.factor(logRSS.open.indiv$ttd)

logRSS.open.indiv.se <- unique(logRSS.open.indiv[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.open.pop <- merge(logRSS.open, logRSS.open.indiv.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.open.pop)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Open") +
  ggtitle("1 days") 
# 
# ggplot(data=setDT(logRSS.open.indiv)[ttd=='1 day'], aes(x, rss, colour=COD)) +
#   #geom_line() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
#   geom_smooth()

ggplot(data=setDT(logRSS.open.pop)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2)) +
  ylab("logRSS") + xlab("open") +
  ggtitle("60 days") 




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
                                                             wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                  wolf_step_id = NA, wolfID = 'RMNP_W02')]

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





##### graphs ####
logRSS.wet.indiv <- setDT(logRSS.wet.indiv)
logRSS.wet.indiv[,'COD'] <- as.factor(logRSS.wet.indiv$COD)
logRSS.wet.indiv[,'ttd'] <- as.factor(logRSS.wet.indiv$ttd)

logRSS.wet.indiv.se <- unique(logRSS.wet.indiv[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.wet.pop <- merge(logRSS.wet, logRSS.wet.indiv.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.wet.pop)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("wet") +
  ggtitle("1 days") 
# 
# ggplot(data=setDT(logRSS.wet.indiv)[ttd=='1 day'], aes(x, rss, colour=COD)) +
#   #geom_line() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
#   geom_smooth()

ggplot(data=setDT(logRSS.wet.pop)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2)) +
  ylab("logRSS") + xlab("wet") +
  ggtitle("60 days") 

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
                          wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                              wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
    wolfID = 'RMNP_W02'
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
                                                                wolf_step_id = NA, wolfID = 'RMNP_W02')]

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





##### graphs ####
logRSS.road.indiv <- setDT(logRSS.road.indiv)
logRSS.road.indiv[,'COD'] <- as.factor(logRSS.road.indiv$COD)
logRSS.road.indiv[,'ttd'] <- as.factor(logRSS.road.indiv$ttd)

logRSS.road.indiv.se <- unique(logRSS.road.indiv[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.road.pop <- merge(logRSS.road, logRSS.road.indiv.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.road.pop)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("1 days") 
# 
# ggplot(data=setDT(logRSS.road.indiv)[ttd=='1 day'], aes(x, rss, colour=COD)) +
#   #geom_line() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
#   geom_smooth()

ggplot(data=setDT(logRSS.road.pop)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2)) +
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("60 days") 

  



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
                                                             wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
                                                              wolf_step_id = NA, wolfID = 'RMNP_W02')]
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
                                                           wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
                                                            wolf_step_id = NA, wolfID = 'RMNP_W02')]
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
                                                           wolf_step_id = NA, wolfID = 'RMNP_W02')]

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
                                                            wolf_step_id = NA, wolfID = 'RMNP_W02')]
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


##### graphs ####
logRSS.nn.indiv <- setDT(logRSS.nn.indiv)
logRSS.nn.indiv[,'COD'] <- as.factor(logRSS.nn.indiv$COD)
logRSS.nn.indiv[,'ttd'] <- as.factor(logRSS.nn.indiv$ttd)

logRSS.nn.indiv.se <- unique(logRSS.nn.indiv[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.nn.pop <- merge(logRSS.nn, logRSS.nn.indiv.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.nn.pop)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))


ggplot(data=setDT(logRSS.nn.pop)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))



logRSS.nn.indiv.social <- setDT(logRSS.nn.indiv.social)
logRSS.nn.indiv.social[,'COD'] <- as.factor(logRSS.nn.indiv.social$COD)
logRSS.nn.indiv.social[,'ttd'] <- as.factor(logRSS.nn.indiv.social$ttd)

logRSS.nn.indiv.social.se <- unique(logRSS.nn.indiv.social[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.nn.pop.social <- merge(logRSS.social, logRSS.nn.indiv.social.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.nn.pop.social)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to NN (m)") +
  ggtitle("1 day") 

# ggplot(data=setDT(logRSS.nn.indiv.social)[ttd=='1 day'], aes(x, rss, colour=COD)) +
#   geom_point() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)# +
#   #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))


ggplot(data=setDT(logRSS.nn.pop.social)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to NN (m)") +
  ggtitle("60 days") 

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
                                                           packDist_end = seq(0, max(packDist_end), length.out = 150),
                                                           COD = factor('CDV', levels = levels(COD)),
                                                           ttd1 = 1,
                                                           wolf_step_id = NA, wolfID = 'RMNP_W02')]

p.pack.CDV.1day.1 <- predict(everyone, newdata = pack.CDV.1day.1, type='link', re.form = NA)
p.pack.CDV.1day.1.social <- predict(everyone.social, newdata = pack.CDV.1day.1, type='link', re.form = NA)


logRSS.pack.CDV.1day<- p.pack.CDV.1day.1 - p.CDV.1day.2
logRSS.pack.CDV.1day.social<- p.pack.CDV.1day.1.social - p.CDV.1day.2.social


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
                                                            packDist_end = seq(0, max(packDist_end), length.out = 150),
                                                            COD = factor('CDV', levels = levels(COD)),
                                                            ttd1 = 60,
                                                            wolf_step_id = NA, wolfID = 'RMNP_W02')]
p.pack.CDV.60day.1 <- predict(everyone, newdata = pack.CDV.60day.1, type='link', re.form=NA)
p.pack.CDV.60day.1.social <- predict(everyone.social, newdata = pack.CDV.60day.1, type='link', re.form=NA)


logRSS.pack.CDV.60day<- p.pack.CDV.60day.1 - p.CDV.60day.2
logRSS.pack.CDV.60day.social<- p.pack.CDV.60day.1.social - p.CDV.60day.2.social


logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'CDV', ttd = '1 day', var = 'boundary', rss = logRSS.pack.CDV.1day, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'CDV', ttd = '60 days', var = 'boundary', rss = logRSS.pack.CDV.60day, x = seq(from = 0, to = 15000, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '1 day', var = 'boundary', rss = logRSS.pack.CDV.1day.social, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'CDV', ttd = '60 days', var = 'boundary', rss = logRSS.pack.CDV.60day.social, x = seq(from = 0, to = 15000, length.out = 150)))

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
          packDist_end = seq(0, max(packDist_end), length.out = 150),
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
                                                             packDist_end = seq(0, max(packDist_end), length.out = 150),
                                                             COD = factor('human', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA, wolfID = 'RMNP_W02')]

p.pack.human.1day.1 <- predict(everyone, newdata = pack.human.1day.1, type='link', re.form = NA)
p.pack.human.1day.1.social <- predict(everyone.social, newdata = pack.human.1day.1, type='link', re.form = NA)


logRSS.pack.human.1day<- p.pack.human.1day.1 - p.human.1day.2
logRSS.pack.human.1day.social<- p.pack.human.1day.1.social - p.human.1day.2.social


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
                                                              packDist_end = seq(0, max(packDist_end), length.out = 150),
                                                              COD = factor('human', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = 'RMNP_W02')]
p.pack.human.60day.1 <- predict(everyone, newdata = pack.human.60day.1, type='link', re.form=NA)
p.pack.human.60day.1.social <- predict(everyone.social, newdata = pack.human.60day.1, type='link', re.form=NA)



logRSS.pack.human.60day<- p.pack.human.60day.1 - p.human.60day.2
logRSS.pack.human.60day.social<- p.pack.human.60day.1.social - p.human.60day.2.social

logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'human', ttd = '1 day', var = 'boundary', rss = logRSS.pack.human.1day, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'human', ttd = '60 days', var = 'boundary', rss = logRSS.pack.human.60day, x = seq(from = 0, to = 15000, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '1 day', var = 'boundary', rss = logRSS.pack.human.1day.social, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'human', ttd = '60 days', var = 'boundary', rss = logRSS.pack.human.60day.social, x = seq(from = 0, to = 15000, length.out = 150)))

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
                                                               packDist_end = seq(0, max(packDist_end), length.out = 150),
                                                               COD = factor('control', levels = levels(COD)),
                                                               ttd1 = 1,
                                                               wolf_step_id = NA, wolfID = 'RMNP_W02')]

p.pack.control.1day.1 <- predict(everyone, newdata = pack.control.1day.1, type='link', re.form = NA)
p.pack.control.1day.1.social <- predict(everyone.social, newdata = pack.control.1day.1, type='link', re.form = NA)


logRSS.pack.control.1day<- p.pack.control.1day.1 - p.control.1day.2
logRSS.pack.control.1day.social<- p.pack.control.1day.1.social - p.control.1day.2.social


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
                                                                packDist_end = seq(0, max(packDist_end), length.out = 150),
                                                                COD = factor('control', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = 'RMNP_W02')]
p.pack.control.60day.1 <- predict(everyone, newdata = pack.control.60day.1, type='link', re.form=NA)
p.pack.control.60day.1.social <- predict(everyone.social, newdata = pack.control.60day.1, type='link', re.form=NA)



logRSS.pack.control.60day<- p.pack.control.60day.1 - p.control.60day.2
logRSS.pack.control.60day.social<- p.pack.control.60day.1.social - p.control.60day.2.social

logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'control', ttd = '1 day', var = 'boundary', rss = logRSS.pack.control.1day, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.nn <- rbind(logRSS.nn, data.frame(COD= 'control', ttd = '60 days', var = 'boundary', rss = logRSS.pack.control.60day, x = seq(from = 0, to = 15000, length.out = 150)))

logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '1 day', var = 'boundary', rss = logRSS.pack.control.1day.social, x = seq(from = 0, to = 15000, length.out = 150)))
logRSS.social <- rbind(logRSS.social, data.frame(COD= 'control', ttd = '60 days', var = 'boundary', rss = logRSS.pack.control.60day.social, x = seq(from = 0, to = 15000, length.out = 150)))


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




##### gathering RSS
logRSS.pack.indiv <- rbind(logRSS.pack.control.1day.indiv, logRSS.pack.control.60day.indiv, 
                         logRSS.pack.human.1day.indiv, logRSS.pack.human.60day.indiv, 
                         logRSS.pack.CDV.1day.indiv, logRSS.pack.CDV.60day.indiv)

logRSS.pack.indiv.social <- rbind(logRSS.pack.control.1day.indiv.social, logRSS.pack.control.60day.indiv.social, 
                                logRSS.pack.human.1day.indiv.social, logRSS.pack.human.60day.indiv.social, 
                                logRSS.pack.CDV.1day.indiv.social, logRSS.pack.CDV.60day.indiv.social)


#### all ###

logRSS.indiv <- rbind(logRSS.forest.indiv, logRSS.open.indiv, logRSS.wet.indiv, logRSS.road.indiv, 
                      logRSS.nn.indiv, logRSS.pack.indiv)
logRSS.indiv.model <- rbind(logRSS.forest.indiv.habitat, logRSS.open.indiv.habitat, logRSS.wet.indiv.habitat, logRSS.road.indiv.habitat, 
                      logRSS.nn.indiv.social, logRSS.pack.indiv.social)

logRSS <- rbind(logRSS.forest, logRSS.open, logRSS.wet, logRSS.road, logRSS.nn)
logRSS.model <- rbind(logRSS.forest.habitat, logRSS.open.habitat, logRSS.wet.habitat, logRSS.road.habitat, logRSS.social)

#saveRDS(logRSS.indiv, 'data/derived-data/logRSS_indiv2.Rds')
#saveRDS(logRSS.indiv.model, 'data/derived-data/logRSS_indiv_model2.Rds')


#saveRDS(logRSS, 'data/derived-data/logRSS2.Rds')
#saveRDS(logRSS.model, 'data/derived-data/logRSS_model2.Rds')


##### graphs ####
logRSS.pack.indiv <- setDT(logRSS.pack.indiv)
logRSS.pack.indiv[,'COD'] <- as.factor(logRSS.pack.indiv$COD)
logRSS.pack.indiv[,'ttd'] <- as.factor(logRSS.pack.indiv$ttd)

logRSS.pack.indiv.se <- unique(logRSS.pack.indiv[,.(se=se(rss)), by = .(x, COD, ttd, var)])

logRSS.pack.pop <- merge(setDT(logRSS.nn)[var=='boundary'], logRSS.pack.indiv.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.pack.pop)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))


ggplot(data=setDT(logRSS.pack.pop)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))



logRSS.pack.indiv.social <- setDT(logRSS.pack.indiv.social)
logRSS.pack.indiv.social[,'COD'] <- as.factor(logRSS.pack.indiv.social$COD)
logRSS.pack.indiv.social[,'ttd'] <- as.factor(logRSS.pack.indiv.social$ttd)

logRSS.pack.indiv.social.se <- unique(logRSS.pack.indiv.social[,.(se=se(rss), var), by = .(x, COD, ttd)])

logRSS.pack.pop.social <- merge(setDT(logRSS.social)[var=='boundary'], logRSS.pack.indiv.social.se, by = c('x', 'COD', 'ttd','var'))


ggplot(data=setDT(logRSS.pack.pop.social)[ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to NN (m)") +
  ggtitle("1 day") 

# ggplot(data=setDT(logRSS.pack.indiv.social)[ttd=='1 day'], aes(x, rss, colour=COD)) +
#   geom_point() +
#   geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)# +
#   #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))


ggplot(data=setDT(logRSS.pack.pop.social)[ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_line() +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to NN (m)") +
  ggtitle("60 days") 


#### GRAPHS ####

logRSS.indiv <- setDT(logRSS.indiv)
logRSS.indiv[,'COD'] <- factor(logRSS.indiv$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))
logRSS.indiv[,'ttd'] <- as.factor(logRSS.indiv$ttd)

logRSS.indiv.se <- unique(logRSS.indiv[,.(se=se(rss)), by = .(x, var, COD, ttd)])

logRSS.pop <- merge(setDT(logRSS), logRSS.indiv.se, by = c('x', 'COD', 'ttd','var'))


logRSS.indiv.model <- setDT(logRSS.indiv.model)
logRSS.indiv.model[,'COD'] <- factor(logRSS.indiv.model$COD, levels = c('control','human','CDV'), labels = c('control','human','CDV'))
logRSS.indiv.model[,'ttd'] <- as.factor(logRSS.indiv.model$ttd)

logRSS.indiv.model.se <- unique(logRSS.indiv.model[,.(se=se(rss)), by = .(x, var, COD, ttd)])

logRSS.pop.model <- merge(setDT(logRSS.model), logRSS.indiv.model.se, by = c('x', 'COD', 'ttd','var'))


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
  ylim(-0.7,1.3) +
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
  ylim(-0.7,1.3) +
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
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion forest") +
  ggtitle("a) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-1.2,1) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

forest.60.b <- ggplot(data=setDT(logRSS.indiv)[var == 'forest'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion forest") +
  ggtitle("b) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-1.2,1) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

forest.1.b|forest.60.b

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
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion open") +
  ggtitle("a) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-2,2) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

open.60.b <- ggplot(data=setDT(logRSS.indiv)[var == 'open'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion open") +
  ggtitle("b) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-2,2) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

open.1.b|open.60.b


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
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion wet") +
  ggtitle("a) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-1,2.5) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

wet.60.b <- ggplot(data=setDT(logRSS.indiv)[var == 'wet'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'glm') +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Proportion wet") +
  ggtitle("b) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-1,2.5) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

wet.1.b|wet.60.b


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


road.1.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("a) 1 day") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,2.5) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

road.60.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'road'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'road'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to road (m)") +
  ggtitle("b) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,2.5) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

road.1.model.b|road.60.model.b

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

nn.1.model|nn.60.model


nn.1.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
 # geom_line(data=logRSS.pop.model[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to nearest neighbor (m)") +
  ggtitle("a) 1 day") +
  theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-.5,6) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

nn.60.model.b <- ggplot(data=setDT(logRSS.indiv.model)[var == 'nn'& ttd=='60 days'], aes(x, rss, colour=COD)) +
  geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), method = 'loess') +
  # geom_line(data=logRSS.pop.model[var == 'nn'& ttd=='1 day'], aes(x, rss, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss - 1.96*se), ymax = (rss + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Distance to nearest neighbor (m)") +
  ggtitle("b) 60 days") +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
 ylim(-.5,6) +
  scale_colour_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  scale_fill_manual("", values = c("deepskyblue", "purple", "dark green"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))

nn.1.model.b|nn.60.model.b







