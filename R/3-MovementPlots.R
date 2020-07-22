### Movement graphs ====
# Julie Turner
# Started: June 1 2020


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

dat <- readRDS('data/derived-data/everyone_betas_COD.Rds')
dat.pop <- readRDS('data/derived-data/summarypopeveryone_COD.Rds')
dat.pop <- setDT(dat.pop)[effect =='fixed' & (term %like% 'log_sl:COD' |term %like%  'cos_ta')]
params <- readRDS('data/derived-data/moveParams_all.Rds')

#dat <- everyone.indiv

dat.wide <- dcast(data =dat, wolfID + COD ~ term, value.var = 'estimate')

dat.wide <- setDT(merge(dat.wide, params[,.(wolfID,shape, scale, kappa)], by = 'wolfID', all.x = T))

cod.params <- dat.wide[,.(mean.shp = mean(shape, na.rm = T), se.shp = se(shape),
                          mean.scl = mean(scale, na.rm = T), se.scl = se(scale), 
                          mean.kap = mean(kappa, na.rm = T), se.kap = se(kappa)), by =.(COD)]

#dat.pop <- sum.everyone
dat.pop$term <- gsub('[[:punct:]]', '', dat.pop$term)
dat.pop$term <- gsub(' ', '', dat.pop$term)
dat.pop$term <- gsub('1', '', dat.pop$term)
dat.pop[,'COD'] <- ifelse(dat.pop$term %like% 'control', 'control',
                          ifelse(dat.pop$term %like% 'human', 'human', 'CDV'))
dat.pop$term <- gsub('CODcontrol', '', dat.pop$term)
dat.pop$term <- gsub('CODhuman', '', dat.pop$term)
dat.pop$term <- gsub('CODCDV', '', dat.pop$term)
dat.pop$term <- gsub('Ilog', '_', dat.pop$term)

dat.pop.wide <- dcast(data =setDT(dat.pop)[is.na(group)], COD ~ term, value.var = c('estimate', 'std.error'))
dat.pop.wide <- setDT(merge(dat.pop.wide, cod.params, by = 'COD', all.x = T))



t2d <- seq(0, 61, length.out = 100)

intercept <- 6.576
forest <- 0.059931*0.7858649
open <- 0.245925*0.05676634
wet <- 0.125641*0.1573687

dat.wide[, spd:= list(list((shape+log_sl + intercept + forest + open + wet +(`log_sl-ttd`*t2d))*(scale))), by=.(wolfID)]
dat.wide[, dir:= list(list(kappa + cos_ta + (`cos_ta-ttd`*t2d))), by=.(wolfID)]
dat.wide[, ttd:= list(list(seq(0,61,length.out = 100))), by=.(wolfID)]

move <- dat.wide[, .(spd = unlist(spd), dir=unlist(dir), ttd= unlist(ttd)), by=.(wolfID,COD)]
move[,spd_hr :=spd/5000]

move[spd_hr<0, unique(wolfID)]
move[,pop:=gsub('_.*$','',wolfID)]

gcolors <- c("deepskyblue", "purple", "dark green")
speed <- ggplot(data=move[spd_hr >=0], aes(x=-ttd, y=(spd_hr), color = COD)) + 
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = FALSE, show.legend = F, method = 'lm')+
  theme_classic() +
  theme(text = element_text(size=15)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x =  element_text(size = 15)) + 
  #  theme(legend.position = "none") +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) +
  ggtitle("a)") +
  xlab("Time to death (days)") + ylab("Speed (km/hour)")
speed 

speed.disease <- ggplot(data=move[spd_hr >=0 & COD != 'human'], aes(x=-ttd, y=(spd_hr), color = COD)) + 
  geom_line(aes(group = wolfID),alpha = .5, linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = FALSE, show.legend = T, method = 'lm')+
  theme_classic() +
  theme(text = element_text(size=15)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x =  element_text(size = 15)) + 
  #  theme(legend.position = "none") +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) +
 # ggtitle("a)") +
  xlab("Time to death (days)") + ylab("Speed (km/hour)")
speed.disease 


direction <- ggplot(data=move, aes(x=-ttd, y=dir, color = COD)) + 
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = FALSE)+
  theme_classic() +
  theme(text = element_text(size=15)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x =  element_text(size = 15)) + 
  #  theme(legend.position = "none") +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) +
  ggtitle("b)") +
  xlab("Time to death (days)") + ylab("Concentration of turn angle")
direction

speed|direction


speed2 <- ggplot(data=move[spd_hr >=0 ], aes(x=-ttd, y=(spd_hr), color = COD)) + 
  #geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = T, show.legend = F, method = 'lm')+
  theme_classic() +
  theme(text = element_text(size=15)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #theme(axis.text.x =  element_text(size = 15)) +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  #  theme(legend.position = "none") +
 # ylim(-1,7) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) + theme(legend.text = element_text(size = 10)) +
  ggtitle("a) Speed") +
  xlab("Time to death (days)") + ylab("Speed (km/hour)")
speed2 


direction2 <- ggplot(data=move, aes(x=-ttd, y=dir, color = COD)) + 
  #geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = T)+
  theme_classic() +
  theme(text = element_text(size=15)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x =  element_text(size = 15)) + 
  #  theme(legend.position = "none") +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) +
  ggtitle("b)") +
  xlab("Time to death (days)") + ylab("Concentration of turn angle")
direction2

speed2|direction2


ggplot(data=move[spd_hr >=0], aes(x=-ttd, y=(spd_hr), color = pop)) + 
  #geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = pop), se = T, show.legend = T, method = 'lm')+
  theme_classic() +
  theme(text = element_text(size=15)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x =  element_text(size = 15)) + 
  #  theme(legend.position = "none") +
  # ylim(-1,7) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) +
  ggtitle("a)") +
  xlab("Time to death (days)") + ylab("Speed (km/hour)")

move.pop <- move[,.(mean.spd = mean(spd, na.rm = T), se.spd = se(spd), 
                    mean.dir = mean(dir, na.rm = T), se.dir = se(dir)), by=.(ttd,COD)]

ggplot(move.pop, aes( x = -ttd, y= mean.spd, color = COD))+
  geom_line() +
  geom_ribbon(aes(-ttd, ymin = (mean.spd - 1.96*se.spd), ymax = (mean.spd + 1.96*se.spd), fill=COD, alpha = .2), show.legend = F)




####pop####

dat.pop.wide[, spd:= list(list((mean.shp+estimate_logsl+(estimate_logsl_ttd*t2d))*(mean.scl))), by=.(COD)]
dat.pop.wide[, spd.upr:= list(list((mean.shp+(estimate_logsl + (1.96*std.error_logsl))+((estimate_logsl_ttd + (1.96*std.error_logsl_ttd))*t2d))*(mean.scl))), by=.(COD)]
dat.pop.wide[, spd.lwr:= list(list((mean.shp+(estimate_logsl - (1.96*std.error_logsl))+((estimate_logsl_ttd - (1.96*std.error_logsl_ttd))*t2d))*(mean.scl))), by=.(COD)]

dat.pop.wide[, dir:= list(list(kappa + cos_ta + (`cos_ta-ttd`*t2d))), by=.(wolfID)]
dat.pop.wide[, ttd:= list(list(seq(1,61,length.out = 100))), by=.(wolfID)]

dat.pop.wide[, ttd:= list(list(seq(1,61,length.out = 100))), by=.(COD)]

move.pop <- dat.pop.wide[, .(spd = unlist(spd), spd.upr=unlist(spd.upr), spd.lwr=unlist(spd.lwr), ttd= unlist(ttd)), by=.(COD)]
move.pop[,spd_hr :=spd/500]
move.pop[,spd.upr_hr :=spd.upr/500]
move.pop[,spd.lwr_hr :=spd.lwr/500]

ggplot(move.pop, aes( x = -ttd, y= spd, color = COD))+
  geom_line() +
  geom_ribbon(aes(-ttd, ymin = (spd.lwr), ymax = (spd.upr), fill=COD, alpha = .2), show.legend = F)



