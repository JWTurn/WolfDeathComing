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

dat[,'land_end_adj'] <- ifelse(dat$land_end == 'wet', 'wet', 
                               ifelse(dat$land_end == 'mixed'|dat$land_end == 'deciduous'|dat$land_end == 'coniferous', 'forest','open'))

#dat[,'forest']
# dat[,'propforest_end_adj'] <- dat$propconif_end+dat$propmixed_end +dat$propdecid_end 
# dat[,'propopen_end_adj'] <- dat$propopen_end +dat$propshrub_end + dat$propurban_end
# 
# 
# dat[,mean(propwet_end), by=.(wolfID,ua)]
# dat[,mean(propforest_end_adj), by=.(wolfID,ua)]
# dat[,mean(propopen_end_adj), by=.(wolfID,ua)]
# dat<-dat[COD == 'cdv'|COD == 'human'|COD == 'none']

dat[,'wolf_step_id'] <- paste(dat$wolfID, dat$step_id_, sep = '_')

dat[ua=='used',.(road=mean(roadDist_end, na.rm = T), nn= mean(distance2, na.rm = T), pack= mean(packDistadj_end, na.rm = T)), by= COD]

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


#### everyone ####

everyone <- glmmTMB(case_ ~ log_sl:ToD_start +
                       # log_sl:land_end_adj +
                        log_sl:COD + cos_ta:COD + 
                        I(log(ttd1 + 1)):log_sl:COD + I(log(ttd1 + 1)):cos_ta:COD +
                        (1|wolf_step_id) +
                        (0 + (log_sl)|wolfID) +
                        (0 + (cos_ta)|wolfID) +
                        (0 + (I(log(ttd1 + 1)):log_sl)|wolfID) +
                        (0 + (I(log(ttd1 + 1)):cos_ta)|wolfID) +
                        COD:land_end_adj + I(log(1+roadDist_end)):COD +
                        COD:I(log(ttd1 + 1)):land_end_adj +  I(log(ttd1 + 1)):I(log(1+roadDist_end)):COD +
                        
                        (0 + land_end_adj|wolfID) + (0 + (I(log(ttd1 + 1)):land_end_adj)|wolfID) +
                        (0 + I(log(1+roadDist_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+roadDist_end)))|wolfID) +
                        
                        I(log(1+distance2)):COD + I(log(1+packDistadj_end)):COD +
                        I(log(ttd1 + 1)):I(log(1+distance2)):COD + I(log(ttd1 + 1)):I(log(1+packDistadj_end)):COD +
                        
                        (0 + I(log(1+distance2))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+distance2)))|wolfID) +
                        (0 + I(log(1+packDistadj_end))|wolfID) + (0 + (I(log(ttd1 + 1)):I(log(1+packDistadj_end)))|wolfID)
                      , family=poisson(),
                      data = dat[wolfID %chin% dat.wnn.lastmo$wolfID], 
                    map = list(theta=factor(c(NA,1:22))), start = list(theta=c(log(1000),seq(0,0, length.out = 22))))

# everyone$parameters$theta[1] <- log(1e3)
# nvar_parm <- length(everyone$parameters$theta)
# everyone$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
# everyone <- glmmTMB:::fitTMB(everyone)
summary(everyone)

summary(everyone)$coef$cond[-1, "Estimate"]
popeveryone<- summary(everyone)$coef$cond[-1, 1:2]
#saveRDS(popeveryone, 'data/derived-data/popeveryone_COD.Rds')
sum.everyone<- summary(everyone)$coef$cond
#saveRDS(sum.everyone, 'data/derived-data/summarypopeveryone_COD.Rds')

everyone.ran_vals <-broom.mixed::tidy(everyone, effect= 'ran_vals')
# everyone.ran_pars <-broom.mixed::tidy(everyone, effect= 'ran_pars')
everyone.se <-setDT(everyone.ran_vals)[group=='wolfID']


# everyone.all.indiv <- coef(everyone)$cond$wolfID %>% rownames_to_column("wolfID") %>% 
#   pivot_longer(-wolfID, names_to = "term", values_to = "estimate") %>% 
#   mutate(method = "ME")


everyone.lastmo <- glmmTMB(case_ ~ log_sl:ToD_start +
                      log_sl:land_end_adj +
                      log(ttd1+1):log_sl + cos_ta + log(ttd1+1):cos_ta +
                      (1|wolf_step_id) +
                      (0 + (log_sl)|wolfID) +
                      (0 + (cos_ta)|wolfID) +
                      (0 + (log(ttd1+1):log_sl)|wolfID) +
                      (0 + (log(ttd1+1):cos_ta)|wolfID) +
                      land_end_adj + log(1+roadDist_end) +
                      log(ttd1+1):land_end_adj +  log(ttd1+1):log(1+roadDist_end) +
                      
                      (0 + land_end_adj|wolfID) + (0 + (log(ttd1+1):land_end_adj)|wolfID) +
                      (0 + (log(1+roadDist_end))|wolfID) + (0 + (log(ttd1+1):log(1+roadDist_end))|wolfID) +
                      
                      log(1+distance2) + log(1+packDistadj_end) +
                      log(ttd1+1):log(1+distance2) + log(ttd1+1):log(1+packDistadj_end) +
                      
                      (0 + (log(1+distance2))|wolfID) + (0 + (log(ttd1+1):log(1+distance2))|wolfID) +
                      (0 + (log(1+packDistadj_end))|wolfID) + (0 + (log(ttd1+1):log(1+packDistadj_end))|wolfID)
                    , family=poisson(),
                    data = dat[ttd1<=31 & wolfID %chin% dat.wnn.lastmo$wolfID], doFit=FALSE)

everyone.lastmo$parameters$theta[1] <- log(1e3)
nvar_parm <- length(everyone.lastmo$parameters$theta)
everyone.lastmo$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))
everyone.lastmo <- glmmTMB:::fitTMB(everyone.lastmo)
summary(everyone.lastmo)

summary(everyone.lastmo)$coef$cond[-1, "Estimate"]
popeveryone.lastmo<- summary(everyone.lastmo)$coef$cond[-1, 1:2]
saveRDS(popeveryone.lastmo, 'data/derived-data/popeveryone_lastmo.Rds')
summary(everyone.lastmo)$varcor

everyone.lastmo.ran_vals <-broom.mixed::tidy(everyone.lastmo, effect= 'ran_vals')
everyone.lastmo.se <-setDT(everyone.lastmo.ran_vals)[group=='wolfID']



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
#### road ####
#### CDV ####
road.CDV.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                          log_sl = mean(log_sl),
                          cos_ta = mean(cos_ta),
                          land_end_adj = factor('forest', levels = levels(land_end_adj)),
                          roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                          distance2 = mean(distance2, na.rm = T),
                          packDistadj_end = mean(packDistadj_end, na.rm = T),
                          COD = factor('CDV', levels = levels(COD)),
                          ttd1 = 1,
                          wolf_step_id = NA,
                          wolfID = NA)]

p.road.CDV.1day.1 <- predict(everyone, newdata = road.CDV.1day.1, type='link')


road.CDV.1day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             roadDist_end = mean(roadDist_end, na.rm = T),
                                                             distance2 = mean(distance2, na.rm = T),
                                                             packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA, wolfID = NA)]
p.road.CDV.1day.2 <- predict(everyone, newdata = road.CDV.1day.2, type='link')

logRSS.road.CDV.1day<- p.road.CDV.1day.1 - p.road.CDV.1day.2


road.CDV.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                             distance2 = mean(distance2, na.rm = T),
                                                             packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID = NA)]
p.road.CDV.60day.1 <- predict(everyone, newdata = road.CDV.60day.1, type='link')


road.CDV.60day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             roadDist_end = mean(roadDist_end, na.rm = T),
                                                             distance2 = mean(distance2, na.rm = T),
                                                             packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 60,
                                                             wolf_step_id = NA, wolfID = NA)]
p.road.CDV.60day.2 <- predict(everyone, newdata = road.CDV.60day.2, type='link')

logRSS.road.CDV.60day<- p.road.CDV.60day.1 - p.road.CDV.60day.2

logRSS <- data.frame(COD= 'CDV', ttd = '1 day', var = 'road', rss = logRSS.road.CDV.1day, x = seq(from = 0, to = 3000, length.out = 100))
logRSS <- rbind(logRSS, data.frame(COD= 'CDV', ttd = '60 days', var = 'road', rss = logRSS.road.CDV.60day, x = seq(from = 0, to = 3000, length.out = 100)))


#### indivs ####
road.CDV.1day.1.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'CDV',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                             distance2 = mean(distance2, na.rm = T),
                                                             packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                             COD = factor('CDV', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA),
                                                             by=.(wolfID) ]

road.CDV.1day.1.indiv[,'wolfID2'] <- road.CDV.1day.1.indiv$wolfID
road.CDV.1day.1.indiv[,'wolfID'] <- NA

p.road.CDV.1day.1.indiv <- road.CDV.1day.1.indiv[,.(h1 = predict(everyone,
                              newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = "link"),
                              x = seq(from = 0, to = 3000, length.out = 100)), 
                              by=.(wolfID2)]


road.CDV.1day.2.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'CDV',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                   log_sl = mean(log_sl),
                                                                   cos_ta = mean(cos_ta),
                                                                   land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                   roadDist_end = mean(roadDist_end, na.rm = T),
                                                                   distance2 = mean(distance2, na.rm = T),
                                                                   packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                   COD = factor('CDV', levels = levels(COD)),
                                                                   ttd1 = 1,
                                                                   wolf_step_id = NA),
                             by=.(wolfID) ]

road.CDV.1day.2.indiv[,'wolfID2'] <- road.CDV.1day.2.indiv$wolfID
road.CDV.1day.2.indiv[,'wolfID'] <- NA

p.road.CDV.1day.2.indiv <- road.CDV.1day.2.indiv[,.(h2 = predict(everyone,
                                                          newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link')), 
                                                 by=.(wolfID2)]
p.road.CDV.1day.1.indiv<- p.road.CDV.1day.1.indiv[,.(wolfID = wolfID2, h1, x, ttd = '1 day', COD = 'CDV', var = 'road')]
p.road.CDV.1day.2.indiv<- p.road.CDV.1day.2.indiv[,.(wolfID = wolfID2, h2, ttd = '1 day', COD = 'CDV', var = 'road')]

logRSS.road.CDV.indiv <- merge(p.road.CDV.1day.1.indiv, p.road.CDV.1day.2.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.CDV.indiv[,'rss'] <- logRSS.road.CDV.indiv$h1 - logRSS.road.CDV.indiv$h2



road.CDV.60day.1.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'CDV',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                                           log_sl = mean(log_sl),
                                                                                           cos_ta = mean(cos_ta),
                                                                                           land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                                           roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                                           distance2 = mean(distance2, na.rm = T),
                                                                                           packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                                           COD = factor('CDV', levels = levels(COD)),
                                                                                           ttd1 = 60,
                                                                                           wolf_step_id = NA),
                                  by=.(wolfID) ]

road.CDV.60day.1.indiv[,'wolfID2'] <- road.CDV.60day.1.indiv$wolfID
road.CDV.60day.1.indiv[,'wolfID'] <- NA

p.road.CDV.60day.1.indiv <- road.CDV.60day.1.indiv[,.(h1 = predict(everyone,
                                                                           newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link'),
                                                              x = seq(from = 0, to = 3000, length.out = 100)), 
                                                           by=.(wolfID2)]


road.CDV.60day.2.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'CDV',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                                           log_sl = mean(log_sl),
                                                                                           cos_ta = mean(cos_ta),
                                                                                           land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                                           roadDist_end = mean(roadDist_end, na.rm = T),
                                                                                           distance2 = mean(distance2, na.rm = T),
                                                                                           packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                                           COD = factor('CDV', levels = levels(COD)),
                                                                                           ttd1 = 60,
                                                                                           wolf_step_id = NA),
                                  by=.(wolfID) ]

road.CDV.60day.2.indiv[,'wolfID2'] <- road.CDV.60day.2.indiv$wolfID
road.CDV.60day.2.indiv[,'wolfID'] <- NA

p.road.CDV.60day.2.indiv <- road.CDV.60day.2.indiv[,.(h2 = predict(everyone,
                                                                           newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link')), 
                                                           by=.(wolfID2)]
p.road.CDV.60day.1.indiv<- p.road.CDV.60day.1.indiv[,.(wolfID = wolfID2, h1, x, ttd = '60 days', COD = 'CDV', var = 'road')]
p.road.CDV.60day.2.indiv<- p.road.CDV.60day.2.indiv[,.(wolfID = wolfID2, h2, ttd = '60 days', COD = 'CDV', var = 'road')]

logRSS.road.CDV.60day.indiv <- merge(p.road.CDV.60day.1.indiv, p.road.CDV.60day.2.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.CDV.60day.indiv[,'rss'] <- logRSS.road.CDV.60day.indiv$h1 - logRSS.road.CDV.60day.indiv$h2



#### human ####

road.human.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                             distance2 = mean(distance2, na.rm = T),
                                                             packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                             COD = factor('human', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA,
                                                             wolfID = NA)]

p.road.human.1day.1 <- predict(everyone, newdata = road.human.1day.1, type='link')


road.human.1day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                             log_sl = mean(log_sl),
                                                             cos_ta = mean(cos_ta),
                                                             land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                             roadDist_end = mean(roadDist_end, na.rm = T),
                                                             distance2 = mean(distance2, na.rm = T),
                                                             packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                             COD = factor('human', levels = levels(COD)),
                                                             ttd1 = 1,
                                                             wolf_step_id = NA, wolfID = NA)]
p.road.human.1day.2 <- predict(everyone, newdata = road.human.1day.2, type='link')

logRSS.road.human.1day<- p.road.human.1day.1 - p.road.human.1day.2


road.human.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                              log_sl = mean(log_sl),
                                                              cos_ta = mean(cos_ta),
                                                              land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                              roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                              distance2 = mean(distance2, na.rm = T),
                                                              packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                              COD = factor('human', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]
p.road.human.60day.1 <- predict(everyone, newdata = road.human.60day.1, type='link')


road.human.60day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                              log_sl = mean(log_sl),
                                                              cos_ta = mean(cos_ta),
                                                              land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                              roadDist_end = mean(roadDist_end, na.rm = T),
                                                              distance2 = mean(distance2, na.rm = T),
                                                              packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                              COD = factor('human', levels = levels(COD)),
                                                              ttd1 = 60,
                                                              wolf_step_id = NA, wolfID = NA)]
p.road.human.60day.2 <- predict(everyone, newdata = road.human.60day.2, type='link')

logRSS.road.human.60day<- p.road.human.60day.1 - p.road.human.60day.2

logRSS <- rbind(logRSS, data.frame(COD= 'human', ttd = '1 day', var = 'road', rss = logRSS.road.human.1day, x = seq(from = 0, to = 3000, length.out = 100)))
logRSS <- rbind(logRSS, data.frame(COD= 'human', ttd = '60 days', var = 'road', rss = logRSS.road.human.60day, x = seq(from = 0, to = 3000, length.out = 100)))


#### indivs ####
road.human.1day.1.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'human',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                   log_sl = mean(log_sl),
                                                                   cos_ta = mean(cos_ta),
                                                                   land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                   roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                   distance2 = mean(distance2, na.rm = T),
                                                                   packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                   COD = factor('human', levels = levels(COD)),
                                                                   ttd1 = 1,
                                                                   wolf_step_id = NA),
                             by=.(wolfID) ]

road.human.1day.1.indiv[,'wolfID2'] <- road.human.1day.1.indiv$wolfID
road.human.1day.1.indiv[,'wolfID'] <- NA

p.road.human.1day.1.indiv <- road.human.1day.1.indiv[,.(h1 = predict(everyone,
                                                                 newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = "link"),
                                                    x = seq(from = 0, to = 3000, length.out = 100)), 
                                                 by=.(wolfID2)]


road.human.1day.2.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'human',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                   log_sl = mean(log_sl),
                                                                   cos_ta = mean(cos_ta),
                                                                   land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                   roadDist_end = mean(roadDist_end, na.rm = T),
                                                                   distance2 = mean(distance2, na.rm = T),
                                                                   packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                   COD = factor('human', levels = levels(COD)),
                                                                   ttd1 = 1,
                                                                   wolf_step_id = NA),
                             by=.(wolfID) ]

road.human.1day.2.indiv[,'wolfID2'] <- road.human.1day.2.indiv$wolfID
road.human.1day.2.indiv[,'wolfID'] <- NA

p.road.human.1day.2.indiv <- road.human.1day.2.indiv[,.(h2 = predict(everyone,
                                                                 newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link')), 
                                                 by=.(wolfID2)]

p.road.human.1day.1.indiv<- p.road.human.1day.1.indiv[,.(wolfID = wolfID2, h1, x, ttd = '1 day', COD = 'human', var = 'road')]
p.road.human.1day.2.indiv<- p.road.human.1day.2.indiv[,.(wolfID = wolfID2, h2, ttd = '1 day', COD = 'human', var = 'road')]

logRSS.road.human.indiv <- merge(p.road.human.1day.1.indiv, p.road.human.1day.2.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.human.indiv[,'rss'] <- logRSS.road.human.indiv$h1 - logRSS.road.human.indiv$h2



road.human.60day.1.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'human',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                                           log_sl = mean(log_sl),
                                                                                           cos_ta = mean(cos_ta),
                                                                                           land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                                           roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                                           distance2 = mean(distance2, na.rm = T),
                                                                                           packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                                           COD = factor('human', levels = levels(COD)),
                                                                                           ttd1 = 60,
                                                                                           wolf_step_id = NA),
                                  by=.(wolfID) ]

road.human.60day.1.indiv[,'wolfID2'] <- road.human.60day.1.indiv$wolfID
road.human.60day.1.indiv[,'wolfID'] <- NA

p.road.human.60day.1.indiv <- road.human.60day.1.indiv[,.(h1 = predict(everyone,
                                                                           newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link'),
                                                              x = seq(from = 0, to = 3000, length.out = 100)), 
                                                           by=.(wolfID2)]


road.human.60day.2.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'human',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                                           log_sl = mean(log_sl),
                                                                                           cos_ta = mean(cos_ta),
                                                                                           land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                                           roadDist_end = mean(roadDist_end, na.rm = T),
                                                                                           distance2 = mean(distance2, na.rm = T),
                                                                                           packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                                           COD = factor('human', levels = levels(COD)),
                                                                                           ttd1 = 60,
                                                                                           wolf_step_id = NA),
                                  by=.(wolfID) ]

road.human.60day.2.indiv[,'wolfID2'] <- road.human.60day.2.indiv$wolfID
road.human.60day.2.indiv[,'wolfID'] <- NA

p.road.human.60day.2.indiv <- road.human.60day.2.indiv[,.(h2 = predict(everyone,
                                                                           newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link')), 
                                                           by=.(wolfID2)]
p.road.human.60day.1.indiv<- p.road.human.60day.1.indiv[,.(wolfID = wolfID2, h1, x, ttd = '60 days', COD = 'human', var = 'road')]
p.road.human.60day.2.indiv<- p.road.human.60day.2.indiv[,.(wolfID = wolfID2, h2, ttd = '60 days', COD = 'human', var = 'road')]

logRSS.road.human.60day.indiv <- merge(p.road.human.60day.1.indiv, p.road.human.60day.2.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.human.60day.indiv[,'rss'] <- logRSS.road.human.60day.indiv$h1 - logRSS.road.human.60day.indiv$h2



#### control ####


road.control.1day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                               log_sl = mean(log_sl),
                                                               cos_ta = mean(cos_ta),
                                                               land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                               roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                               distance2 = mean(distance2, na.rm = T),
                                                               packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                               COD = factor('control', levels = levels(COD)),
                                                               ttd1 = 1,
                                                               wolf_step_id = NA,
                                                               wolfID = NA)]

p.road.control.1day.1 <- predict(everyone, newdata = road.control.1day.1, type='link')


road.control.1day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                               log_sl = mean(log_sl),
                                                               cos_ta = mean(cos_ta),
                                                               land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                               roadDist_end = mean(roadDist_end, na.rm = T),
                                                               distance2 = mean(distance2, na.rm = T),
                                                               packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                               COD = factor('control', levels = levels(COD)),
                                                               ttd1 = 1,
                                                               wolf_step_id = NA, wolfID = NA)]
p.road.control.1day.2 <- predict(everyone, newdata = road.control.1day.2, type='link')

logRSS.road.control.1day<- p.road.control.1day.1 - p.road.control.1day.2


road.control.60day.1 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                distance2 = mean(distance2, na.rm = T),
                                                                packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                COD = factor('control', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]
p.road.control.60day.1 <- predict(everyone, newdata = road.control.60day.1, type='link')


road.control.60day.2 <- dat[wolfID %chin% dat.wnn.lastmo$wolfID,.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                log_sl = mean(log_sl),
                                                                cos_ta = mean(cos_ta),
                                                                land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                roadDist_end = mean(roadDist_end, na.rm = T),
                                                                distance2 = mean(distance2, na.rm = T),
                                                                packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                COD = factor('control', levels = levels(COD)),
                                                                ttd1 = 60,
                                                                wolf_step_id = NA, wolfID = NA)]
p.road.control.60day.2 <- predict(everyone, newdata = road.control.60day.2, type='link')

logRSS.road.control.60day<- p.road.control.60day.1 - p.road.control.60day.2

logRSS <- rbind(logRSS, data.frame(COD= 'control', ttd = '1 day', var = 'road', rss = logRSS.road.control.1day, x = seq(from = 0, to = 3000, length.out = 100)))
logRSS <- rbind(logRSS, data.frame(COD= 'control', ttd = '60 days', var = 'road', rss = logRSS.road.control.60day, x = seq(from = 0, to = 3000, length.out = 100)))

#### indivs ####
road.control.1day.1.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'control',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                   log_sl = mean(log_sl),
                                                                   cos_ta = mean(cos_ta),
                                                                   land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                   roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                   distance2 = mean(distance2, na.rm = T),
                                                                   packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                   COD = factor('control', levels = levels(COD)),
                                                                   ttd1 = 1,
                                                                   wolf_step_id = NA),
                             by=.(wolfID) ]

road.control.1day.1.indiv[,'wolfID2'] <- road.control.1day.1.indiv$wolfID
road.control.1day.1.indiv[,'wolfID'] <- NA

p.road.control.1day.1.indiv <- road.control.1day.1.indiv[,.(h1 = predict(everyone,
                                                                 newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link'),
                                                    x = seq(from = 0, to = 3000, length.out = 100)), 
                                                 by=.(wolfID2)]


road.control.1day.2.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'control',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                   log_sl = mean(log_sl),
                                                                   cos_ta = mean(cos_ta),
                                                                   land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                   roadDist_end = mean(roadDist_end, na.rm = T),
                                                                   distance2 = mean(distance2, na.rm = T),
                                                                   packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                   COD = factor('control', levels = levels(COD)),
                                                                   ttd1 = 1,
                                                                   wolf_step_id = NA),
                             by=.(wolfID) ]

road.control.1day.2.indiv[,'wolfID2'] <- road.control.1day.2.indiv$wolfID
road.control.1day.2.indiv[,'wolfID'] <- NA

p.road.control.1day.2.indiv <- road.control.1day.2.indiv[,.(h2 = predict(everyone,
                                                                 newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link')), 
                                                 by=.(wolfID2)]
p.road.control.1day.1.indiv<- p.road.control.1day.1.indiv[,.(wolfID = wolfID2, h1, x, ttd = '1 day', COD = 'control', var = 'road')]
p.road.control.1day.2.indiv<- p.road.control.1day.2.indiv[,.(wolfID = wolfID2, h2, ttd = '1 day', COD = 'control', var = 'road')]

logRSS.road.control.1day.indiv <- merge(p.road.control.1day.1.indiv, p.road.control.1day.2.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.control.1day.indiv[,'rss'] <- logRSS.road.control.1day.indiv$h1 - logRSS.road.control.1day.indiv$h2



road.control.60day.1.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'control',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                                          log_sl = mean(log_sl),
                                                                                          cos_ta = mean(cos_ta),
                                                                                          land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                                          roadDist_end = seq(from = 0, to = 3000, length.out = 100),
                                                                                          distance2 = mean(distance2, na.rm = T),
                                                                                          packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                                          COD = factor('control', levels = levels(COD)),
                                                                                          ttd1 = 60,
                                                                                          wolf_step_id = NA),
                                 by=.(wolfID) ]

road.control.60day.1.indiv[,'wolfID2'] <- road.control.60day.1.indiv$wolfID
road.control.60day.1.indiv[,'wolfID'] <- NA

p.road.control.60day.1.indiv <- road.control.60day.1.indiv[,.(h1 = predict(everyone,
                                                                         newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link'),
                                                            x = seq(from = 0, to = 3000, length.out = 100)), 
                                                         by=.(wolfID2)]


road.control.60day.2.indiv <- dat[wolfID %chin% dat.wnn.lastmo$wolfID & COD == 'control',.(ToD_start = factor('day', levels = levels(ToD_start)),
                                                                                          log_sl = mean(log_sl),
                                                                                          cos_ta = mean(cos_ta),
                                                                                          land_end_adj = factor('forest', levels = levels(land_end_adj)),
                                                                                          roadDist_end = mean(roadDist_end, na.rm = T),
                                                                                          distance2 = mean(distance2, na.rm = T),
                                                                                          packDistadj_end = mean(packDistadj_end, na.rm = T),
                                                                                          COD = factor('control', levels = levels(COD)),
                                                                                          ttd1 = 60,
                                                                                          wolf_step_id = NA),
                                 by=.(wolfID) ]

road.control.60day.2.indiv[,'wolfID2'] <- road.control.60day.2.indiv$wolfID
road.control.60day.2.indiv[,'wolfID'] <- NA

p.road.control.60day.2.indiv <- road.control.60day.2.indiv[,.(h2 = predict(everyone,
                                                                         newdata = .SD[,.(ToD_start, log_sl, cos_ta, land_end_adj, roadDist_end, distance2, packDistadj_end, COD, ttd1, wolf_step_id, wolfID)], type = 'link')), 
                                                         by=.(wolfID2)]
p.road.control.60day.1.indiv<- p.road.control.60day.1.indiv[,.(wolfID = wolfID2, h1, x, ttd = '60 days', COD = 'control', var = 'road')]
p.road.control.60day.2.indiv<- p.road.control.60day.2.indiv[,.(wolfID = wolfID2, h2, ttd = '60 days', COD = 'control', var = 'road')]

logRSS.road.control.60day.indiv <- merge(p.road.control.60day.1.indiv, p.road.control.60day.2.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)
logRSS.road.control.60day.indiv[,'rss'] <- logRSS.road.control.60day.indiv$h1 - logRSS.road.control.60day.indiv$h2










#####
logRSS.road.indiv <- rbind()

ggplot(setDT(logRSS)[ttd=='1 day'], aes(x, rss)) +
  geom_line(aes(colour = COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) 

ggplot(setDT(logRSS)[ttd=='60 days'], aes(x, rss)) +
  geom_line(aes(colour = COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) 



