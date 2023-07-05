### RSS ====
# Julie Turner
# Revised July 6 2020


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('dplyr', 'tidyr', 'ggplot2', 'patchwork','data.table', 'viridis')
lapply(libs, require, character.only = TRUE)

### Function ----
se <- function(x){
  sd(x, na.rm = T)/ sqrt(length(na.omit(x)))
}


### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

gen.h2.indiv <- readRDS('data/derived-data/genh2_indiv_ttd_scaled.Rds')


logRSS <- readRDS('data/derived-data/logRSS_indiv_ttd.Rds')
logRSS <- merge(logRSS, gen.h2.indiv, by = 'wolfID')
logRSS[, rss.60 := h2 - h2.60]


logRSS[,rss.ttd := rss + rss.60]

dat.meta <- fread(paste0(raw, 'wolf_metadata_all.csv'))
dat.meta[,'wolfpop'] <- paste(dat.meta$pop, dat.meta$WolfID, sep = '_')

logRSS <- merge(logRSS, dat.meta[,.(wolfpop, PackID)], by.x = 'wolfID', by.y = 'wolfpop')


logRSS[, COD2 := ifelse(COD == 'human', 'human-associated', as.character(COD))]
logRSS$COD2 <- factor(logRSS$COD2, 
                      levels = c('control', 'human-associated', 'CDV'))
summary(logRSS$COD2)

#### graph colors ####

cbPalette = c("#A95AA1", "#85C0F9", "#0F2080")
gcolors <- c("deepskyblue", "purple", "dark green")


#### GRAPHS ####
open.75 <- ggplot(data=setDT(logRSS)[var == 'open'& value ==0.75], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("a) High proportion of open areas") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  #+  
  #theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
open.75

wet.75 <- ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.75], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop), show.legend = T) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("b) High proportion of wet areas") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) + theme(legend.text = element_text(size = 10)) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  #+  
  #theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75

road.close <- ggplot(data=setDT(logRSS)[var == 'road'& value ==250], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("c) Close to road") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) + theme(legend.text = element_text(size = 10)) +
  ylim(-10,10) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
road.close


nn.close <- ggplot(data=setDT(logRSS)[var == 'nn'& value ==250], aes(-ttd, rss.ttd, colour=pop)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = pop), se = F, show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'nn'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("b) Proximity to nearest neighbor") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  # ylim(-3,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors) # +  
  #theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.close

# to show
open.75+wet.75+road.close +plot_layout(ncol = 1)
speed + nn.close +plot_layout(ncol = 1)


#### by pack ####
open.75.pack <- ggplot(data=setDT(logRSS)[var == 'open'& value ==0.75], aes(-ttd, rss.ttd, colour=PackID)) +
  geom_line(aes(group = wolfID, linetype = COD, colour = PackID),alpha = .5, show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = PackID), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  #facet_wrap(~COD) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("a) High proportion of open areas") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_viridis(discrete = T)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
open.75.pack

wet.75.pack <- ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.75], aes(-ttd, rss.ttd, colour=PackID)) +
  geom_line(aes(group = wolfID, linetype = COD, colour = PackID),alpha = .5, show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = PackID), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  #facet_wrap(~COD) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("b) High proportion of wet areas") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_viridis(discrete = T)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75.pack

road.close.pack <- ggplot(data=setDT(logRSS)[var == 'road'& value ==250], aes(-ttd, rss.ttd, colour=PackID)) +
  geom_line(aes(group = wolfID, linetype = COD, colour = PackID),alpha = .5, show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = PackID), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  #facet_wrap(~COD) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("c) Close to road") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_viridis(discrete = T)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
road.close.pack

nn.close.pack <- ggplot(data=setDT(logRSS)[var == 'nn'& value ==250], aes(-ttd, rss.ttd, colour=PackID)) +
  geom_line(aes(group = wolfID, linetype = COD, colour = PackID),alpha = .5, show.legend = T) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = PackID), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'open'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  #facet_wrap(~COD) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("d) Proximity to nearest neighbor") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylim(-10,10) +
  scale_colour_viridis(discrete = T)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.close.pack

open.75.pack + wet.75.pack + road.close.pack + nn.close.pack


#### RMNP only ----
#logRSS.rmnp <- setDT(logRSS)[pop == 'RMNP' & COD != 'human']
logRSS.rmnp <- setDT(logRSS)[COD != 'human']
logRSS.rmnp[,COD3:=ifelse(COD == 'control' & pop == 'GHA26', 'control-noCDV', as.character(COD))]
logRSS.rmnp$COD3 <- factor(logRSS.rmnp$COD3, 
                    levels = c('control-noCDV', 'control', 'CDV'))

wet.75.rmnp <- ggplot(data=logRSS.rmnp[var == 'wet'& value ==0.75], aes(-ttd, rss.ttd, color = COD2)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD2), show.legend = T) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, linewidth = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
 # facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("High proportion of wet areas") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) + theme(legend.text = element_text(size = 10)) +
  ylim(-10,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75.rmnp

nn.close <- ggplot(data=logRSS.rmnp[var == 'nn'& value ==250], aes(-ttd, rss.ttd, colour=COD2)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD2), show.legend = T) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, linewidth = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  #facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("Proximity to nearest neighbor") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  # ylim(-3,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors) # +  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.close

nn.close2 <- ggplot(data=logRSS.rmnp[var == 'nn'& value ==250], aes(-ttd, rss.ttd, colour=COD3)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD3), show.legend = T) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, linewidth = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  #facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("Proximity to nearest neighbor") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  # ylim(-3,6) +
  scale_colour_manual("", values = gcolors2)  +  
  scale_fill_manual("", values = gcolors2) # +  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
nn.close2

wet.75.rmnp2 <- ggplot(data=logRSS.rmnp[var == 'wet'& value ==0.75], aes(-ttd, rss.ttd, color = COD3)) +
  geom_line(aes(group = wolfID),alpha = .95, linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD3), show.legend = T) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, linewidth = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  # facet_wrap(~COD2) +
  ylab("logRSS") + xlab("Time to death (days)") +
  ggtitle("High proportion of wet areas") + 
  theme_classic() +
  theme(plot.title=element_text(size=12,hjust = 0.05),axis.text.x = element_text(size=12), axis.title = element_text(size=15),axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) + theme(legend.text = element_text(size = 10)) +
  ylim(-10,6) +
  scale_colour_manual("", values = gcolors2)  +  
  scale_fill_manual("", values = gcolors2)  #+  
#theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75.rmnp2

####
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

wet.75b <- ggplot(data=setDT(logRSS)[var == 'wet'& value ==0.75], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD), show.legend = F) +
  # geom_line(data=logRSS.pop[var == 'wet'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  facet_wrap(~pop) +
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
  ylim(-10,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
wet.75b



road.closeb <- ggplot(data=setDT(logRSS)[var == 'road'& value ==250], aes(-ttd, rss.ttd, colour=COD)) +
  geom_line(aes(group = wolfID,alpha = .0001), linetype ='twodash', show.legend = F) +
  #geom_point(shape = 1, aes(alpha = .001), show.legend = F) +
  geom_smooth(size = 1.5, aes(fill = COD)) +
  # geom_line(data=logRSS.pop[var == 'road'& ttd=='1 day'], aes(x, rss.ttd, colour=COD)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) +
  #geom_ribbon(aes(x, ymin = (rss.ttd - 1.96*se), ymax = (rss.ttd + 1.96*se), fill=COD, alpha = .2))+
  facet_wrap(~pop) +
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
  ylim(-10,6) +
  scale_colour_manual("", values = gcolors)  +  
  scale_fill_manual("", values = gcolors)  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 10))
road.closeb

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

