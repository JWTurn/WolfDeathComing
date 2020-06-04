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


dat.wide <- dcast(data =dat, wolfID + COD ~ term, value.var = 'estimate')

dat.wide <- setDT(merge(dat.wide, params[,.(wolfID,shape, scale, kappa)], by = 'wolfID', all.x = T))

t2d <- seq(0, 61, length.out = 100)

dat.wide[, spd:= list(list((shape+log_sl+(`log_sl-ttd`*t2d))*(scale))), by=.(wolfID)]
dat.wide[, dir:= list(list(kappa + cos_ta + (`cos_ta-ttd`*t2d))), by=.(wolfID)]
dat.wide[, ttd:= list(list(seq(1,61,length.out = 100))), by=.(wolfID)]

move <- dat.wide[, .(spd = unlist(spd), dir=unlist(dir), ttd= unlist(ttd)), by=.(wolfID,COD)]
move[,spd_hr :=spd/500]

gcolors <- c("deepskyblue", "purple", "dark green")
speed <- ggplot(data=move[wolfID != 'GHA26_W35'], aes(x=-ttd, y=spd_hr, color = COD)) + 
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = FALSE, show.legend = F)+
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


speed2 <- ggplot(data=move, aes(x=-ttd, y=spd_hr, color = COD)) + 
  #geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_hline(yintercept=790.9842, linetype='dashed', size = 1) +
  geom_smooth(size = 1.5, aes(fill = COD), se = T, show.legend = F)+
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


move.pop <- move[,.(mean.spd = mean(spd, na.rm = T), se.spd = se(spd), 
                    mean.dir = mean(dir, na.rm = T), se.dir = se(dir)), by=.(ttd,COD)]

ggplot(move.pop, aes( x = -ttd, y= mean.spd, color = COD))+
  geom_line() +
  geom_ribbon(aes(-ttd, ymin = (mean.spd - 1.96*se.spd), ymax = (mean.spd + 1.96*se.spd), fill=COD, alpha = .2), show.legend = F)


