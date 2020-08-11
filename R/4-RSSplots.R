### RSS ====
# Julie Turner
# Revised July 6 2020


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('dplyr', 'tidyr', 'ggplot2', 'patchwork','data.table')
lapply(libs, require, character.only = TRUE)

### Function ----
se <- function(x){
  sd(x, na.rm = T)/ sqrt(length(na.omit(x)))
}


### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

logRSS.prop <- readRDS('data/derived-data/logRSS_indiv.Rds')
logRSS.model.prop <- readRDS('data/derived-data/logRSS_indiv_model.Rds')

logRSS.hab <- logRSS.prop[ttd == '60 days' & (var == 'forest'|var=='open'|var=='wet'|var=='road')]
logRSS.hab <- logRSS.hab[CJ(wolfID = wolfID, var, x = c(0.25, 0.75, 250, 3000), unique = TRUE)
   , on = .(wolfID, var, x)
   , .(wolfID, var, value = x.x, rss)
   , roll = "nearest"]
logRSS.hab<-logRSS.hab[value!= 0 & value !=1]
logRSS.hab$value[logRSS.hab$value%like%0.25] <- 0.25
logRSS.hab$value[logRSS.hab$value%like%0.7] <- 0.75
logRSS.hab$value[logRSS.hab$value%like%242] <- 250
logRSS.hab$value[logRSS.hab$value%like%3000] <- 3000

logRSS.dist <- logRSS.model.prop[ttd == '60 days' & (var == 'nn'|var=='boundary')]
logRSS.dist <- logRSS.dist[CJ(wolfID = wolfID, var, x = c(250, 3000), unique = TRUE)
                         , on = .(wolfID, var, x)
                         , .(wolfID, var, value = x.x, rss)
                         , roll = "nearest"]
logRSS.dist$value[logRSS.dist$value<500] <- 250
logRSS.dist$value[logRSS.dist$value>500] <- 3000
logRSS.dist$var[logRSS.dist$var %like% 'boundary'] <- 'pack'
unique(logRSS.dist$var)

logRSS.60 <- rbind(logRSS.hab, logRSS.dist)
names(logRSS.60)[names(logRSS.60)=="rss"] <- "rss.60"

logRSS.indiv <- readRDS('data/derived-data/logRSS_indiv_ttd.Rds')
logRSS.indiv.model <- readRDS('data/derived-data/logRSS_indiv_model_ttd.Rds')

logRSS <- rbind(logRSS.indiv[var != 'nn' & var != 'pack'], logRSS.indiv.model[var == 'nn' | var == 'pack'])

logRSS.60$value <- as.numeric(logRSS.60$value)

logRSS <- merge(logRSS, logRSS.60, by = c('wolfID','var', 'value'))
logRSS[,rss.ttd := rss + rss.60]



#### graph colors ####

cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")
gcolors <- c("deepskyblue", "purple", "dark green")


#### GRAPHS ####

forest.25 <- ggplot(data=setDT(logRSS)[var == 'forest'& value ==0.25], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("Low proportion of forest") + 
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
 # ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
forest.25

forest.75 <- ggplot(data=setDT(logRSS)[var == 'forest'& value ==0.75], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'forest'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
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
 # ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
forest.75


open.25 <- ggplot(data=setDT(logRSS)[var == 'open'& value ==0.25], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("Low proportion of open") + 
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
 # ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
open.25

open.75 <- ggplot(data=setDT(logRSS)[var == 'open'& value ==0.75], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("a) High proportion of open areas") + 
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
 # ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
open.75


wet.25 <- ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.25], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("Low proportion of wet") + 
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
  #ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.25

wet.75 <- ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.75], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("b) High proportion of wet areas") + 
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
  #ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75



road.close <- ggplot(data=setDT(logRSS)[var == 'road'& value ==250], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'road'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("c) Close to road") + 
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
  #ylim(-5,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
road.close

road.far <- ggplot(data=setDT(logRSS)[var == 'road'& value ==3000], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'road'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("far to road") + 
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
  #ylim(-5,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
road.far


nn.close <- ggplot(data=setDT(logRSS)[var == 'nn'& value ==250], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("b) Close to nearest neighbor") + 
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
 # ylim(-3,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.close

nn.far <- ggplot(data=setDT(logRSS)[var == 'nn'& value ==3000], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("far to nearest neighbor") + 
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
 # ylim(-3,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.far


pack.close <- ggplot(data=setDT(logRSS)[var == 'pack'& value ==250], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'pack'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
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
  #ylim(-4,4) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
pack.close

pack.far <- ggplot(data=setDT(logRSS)[var == 'pack'& value ==3000], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'pack'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("far to pack boundary") + 
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
  #ylim(-4,4) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
pack.far

# to show
open.75|wet.75|road.close
speed|nn.close

### others 
open.25|open.75
road.close|road.far
wet.25|wet.75
nn.close|nn.far

forest.25|forest.75
pack.close|pack.far

##### Disease graphs ###########
wet.disease <- ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.75 & COD != 'human'], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID),alpha = 0.5, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), show.legend = T) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("High proportion of wet areas") + 
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
  #ylim(-6,15) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.disease

nn.disease <- ggplot(data=setDT(logRSS)[var == 'nn'& value ==250 & COD != 'human'], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID),alpha = .5, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to death (days)") +
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
  # ylim(-3,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.disease

wet.disease
nn.disease
speed.disease

#########

ggplot(data=setDT(logRSS)[var == 'open'& value ==0.25], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop)) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
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

ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.25], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop)) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
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

ggplot(data=setDT(logRSS)[var == 'road'& value ==250], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop)) +
  # geom_line(data=logRSS.pop[var == 'road'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("High proportion of road") + 
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

ggplot(data=setDT(logRSS)[var == 'nn'& value ==250], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop)) +
  # geom_line(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  ylab("logRSS") + xlab("Time to Death") +
  ggtitle("High proportion of nn") + 
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

