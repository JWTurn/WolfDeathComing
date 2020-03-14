# RSS mixed effects output
# 17 February 2020
# created by Julie Turner

### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival','forcats', 'patchwork')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

#indivOUT <- readRDS('data/derived-data/full_betas.Rds')
fullOUT <- readRDS('data/derived-data/everyone_betas_lastmoNN.Rds')


### Function ----
se <- function(x){
  sd(x, na.rm = T)/ sqrt(length(na.omit(x)))
}


### simpler RSS for intx terms -- not log-transformed betas Δhi ⋅(βi + βij ⋅hj (x1 )] 
rss.intx <- function(x, xintx, delta, h){
  delta*(x + (xintx*h))
}



### create CIs ----
fullOUT[,'var'] <- ifelse(fullOUT$term %like% 'ttd', 'intx','var')
fullOUT[,'lwr'] <- fullOUT$estimate - (fullOUT$se*1.96)
fullOUT[,'upr'] <- fullOUT$estimate + (fullOUT$se*1.96)

beta.se <- fullOUT[,.(mean = mean(estimate), se= se(estimate)),.(term, COD)]
beta.se[,'lwr'] <- beta.se$mean - (beta.se$se*1.96)
beta.se[,'upr'] <- beta.se$mean + (beta.se$se*1.96)
beta.se[,'var'] <- ifelse(beta.se$term %like% 'ttd', 'intx','var')
beta.se[,'term2'] <- gsub("-ttd", "", beta.se$term)

pd <- position_dodge(0.5) # move them .05 to the left and right
main <- ggplot(beta.se[var=='var'], aes(x=term, y=mean, colour=COD)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
    geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)

ttd<-ggplot(beta.se[var=='intx'], aes(x=term, y=mean, colour=COD)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)

main/ttd

fullOUT$wolfID<-as.factor(fullOUT$wolfID)
fullOUT$var<-as.factor(fullOUT$var)


none.nn.main <- ggplot(fullOUT[var=='var'& COD == 'control'], aes(x=term, y=estimate, colour=wolfID)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)

none.nn.ttd <- ggplot(fullOUT[var=='intx'& COD == 'control'], aes(x=term, y=estimate, colour=wolfID)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)
none.nn.main/none.nn.ttd

cdv.nn.main <- ggplot(fullOUT[var=='var'& COD == 'CDV'], aes(x=term, y=estimate, colour=wolfID)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)

cdv.nn.ttd <- ggplot(fullOUT[var=='intx'& COD == 'CDV'], aes(x=term, y=estimate, colour=wolfID)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)
cdv.nn.main/cdv.nn.ttd
######
### RSS for two log-transformed interaction betas
# transform data
beta.se.wide <- dcast(beta.se, COD + term2 ~ var, value.var = c('mean'))
beta.se.wide <- plyr::rename(beta.se.wide, c('var'='beta', 'intx'='betaintx'))


# 1day
hj.1 <- 1
# 2wks
hj.2 <- 14
# 1mo
hj.3 <- 30
# 2mo
hj.4 <- 60


# 1day
log.hj.1 <- log(1 + 1)
# 2wks
log.hj.2 <- log(1 + 14)
# 1mo
log.hj.3 <- log(1 + 30)
# 2mo 
log.hj.4 <- log(1 + 60)

# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100
hi.nn <- 5001
delta.hi.nn <- 1:5000

#### pop RSS nn ----

humanrss.nn <- beta.se.wide[COD=='human' & term2=='nnDist']

humanrss.nn.1 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(humanrss.nn$beta + (humanrss.nn$betaintx*hj.1))
humanrss.nn.2 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(humanrss.nn$beta + (humanrss.nn$betaintx*hj.2))
humanrss.nn.3 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(humanrss.nn$beta + (humanrss.nn$betaintx*hj.3))
humanrss.nn.4 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(humanrss.nn$beta + (humanrss.nn$betaintx*hj.4))

lnrss.nn.pop.human =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.nn),y=(humanrss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.nn),y=(humanrss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.nn),y=(humanrss.nn.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.nn),y=(humanrss.nn.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to NN (m)") +
  ggtitle("b) human") +
  ylim(-0.01,20) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.nn.pop.human



CDVrss.nn <- beta.se.wide[COD=='CDV' & term2=='nnDist']

CDVrss.nn.1 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(CDVrss.nn$beta + (CDVrss.nn$betaintx*hj.1))
CDVrss.nn.2 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(CDVrss.nn$beta + (CDVrss.nn$betaintx*hj.2))
CDVrss.nn.3 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(CDVrss.nn$beta + (CDVrss.nn$betaintx*hj.3))
CDVrss.nn.4 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(CDVrss.nn$beta + (CDVrss.nn$betaintx*hj.4))

lnrss.nn.pop.CDV =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.nn),y=(CDVrss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.nn),y=(CDVrss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.nn),y=(CDVrss.nn.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.nn),y=(CDVrss.nn.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to NN (m)") +
  ggtitle("c) CDV") +
  ylim(-0.01,20) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.nn.pop.CDV



nonerss.nn <- beta.se.wide[COD=='control' & term2=='nnDist']

nonerss.nn.1 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(nonerss.nn$beta + (nonerss.nn$betaintx*hj.1))
nonerss.nn.2 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(nonerss.nn$beta + (nonerss.nn$betaintx*hj.2))
nonerss.nn.3 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(nonerss.nn$beta + (nonerss.nn$betaintx*hj.3))
nonerss.nn.4 <- (log(hi.nn/(hi.nn-delta.hi.nn)))^(nonerss.nn$beta + (nonerss.nn$betaintx*hj.4))

lnrss.nn.pop.none =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.nn),y=(nonerss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.nn),y=(nonerss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.nn),y=(nonerss.nn.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.nn),y=(nonerss.nn.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to NN (m)") +
  ggtitle("a) control") +
  ylim(-0.01,20) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.nn.pop.none

lnrss.nn.pop.none/lnrss.nn.pop.human/lnrss.nn.pop.CDV

#### pop RSS road dist ----

humanrss.road <- beta.se.wide[COD=='human' & term2=='roadDist']

# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100
hi.road <- 3001
delta.hi.road <- 1:3000

humanrss.road.1 <- (log(hi.road/(hi.road-delta.hi.road)))^(humanrss.road$beta + (humanrss.road$betaintx*hj.1))
humanrss.road.2 <- (log(hi.road/(hi.road-delta.hi.road)))^(humanrss.road$beta + (humanrss.road$betaintx*hj.2))
humanrss.road.3 <- (log(hi.road/(hi.road-delta.hi.road)))^(humanrss.road$beta + (humanrss.road$betaintx*hj.3))
humanrss.road.4 <- (log(hi.road/(hi.road-delta.hi.road)))^(humanrss.road$beta + (humanrss.road$betaintx*hj.4))

lnrss.road.pop.human =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.road),y=(humanrss.road.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.road),y=(humanrss.road.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.road),y=(humanrss.road.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.road),y=(humanrss.road.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to road (m)") +
  ggtitle("b) human") +
  ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.road.pop.human



CDVrss.road <- beta.se.wide[COD=='CDV' & term2=='roadDist']

CDVrss.road.1 <- (log(hi.road/(hi.road-delta.hi.road)))^(CDVrss.road$beta + (CDVrss.road$betaintx*hj.1))
CDVrss.road.2 <- (log(hi.road/(hi.road-delta.hi.road)))^(CDVrss.road$beta + (CDVrss.road$betaintx*hj.2))
CDVrss.road.3 <- (log(hi.road/(hi.road-delta.hi.road)))^(CDVrss.road$beta + (CDVrss.road$betaintx*hj.3))
CDVrss.road.4 <- (log(hi.road/(hi.road-delta.hi.road)))^(CDVrss.road$beta + (CDVrss.road$betaintx*hj.4))

lnrss.road.pop.CDV =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.road),y=(CDVrss.road.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.road),y=(CDVrss.road.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.road),y=(CDVrss.road.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.road),y=(CDVrss.road.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to road (m)") +
  ggtitle("c) CDV") +
  ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.road.pop.CDV



nonerss.road <- beta.se.wide[COD=='control' & term2=='roadDist']

nonerss.road.1 <- (log(hi.road/(hi.road-delta.hi.road)))^(nonerss.road$beta + (nonerss.road$betaintx*hj.1))
nonerss.road.2 <- (log(hi.road/(hi.road-delta.hi.road)))^(nonerss.road$beta + (nonerss.road$betaintx*hj.2))
nonerss.road.3 <- (log(hi.road/(hi.road-delta.hi.road)))^(nonerss.road$beta + (nonerss.road$betaintx*hj.3))
nonerss.road.4 <- (log(hi.road/(hi.road-delta.hi.road)))^(nonerss.road$beta + (nonerss.road$betaintx*hj.4))

lnrss.road.pop.none =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.road),y=(nonerss.road.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.road),y=(nonerss.road.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.road),y=(nonerss.road.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.road),y=(nonerss.road.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to road (m)") +
  ggtitle("a) control") +
  ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.road.pop.none

lnrss.road.pop.none/lnrss.road.pop.human/lnrss.road.pop.CDV


#### pop RSS pack dist ----

humanrss.pack <- beta.se.wide[COD=='human' & term2=='boundaryDist']

# delta hi should be based off avg/median nndist
# hi.pack <- 2:5000
# delta.hi.pack <- 1
hi.pack <- 3001
delta.hi.pack <- 0:3000

humanrss.pack.1 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(humanrss.pack$beta + (humanrss.pack$betaintx*hj.1))
humanrss.pack.2 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(humanrss.pack$beta + (humanrss.pack$betaintx*hj.2))
humanrss.pack.3 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(humanrss.pack$beta + (humanrss.pack$betaintx*hj.3))
humanrss.pack.4 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(humanrss.pack$beta + (humanrss.pack$betaintx*hj.4))

lnrss.pack.pop.human =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.pack),y=(humanrss.pack.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.pack),y=(humanrss.pack.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.pack),y=(humanrss.pack.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.pack),y=(humanrss.pack.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to pack (m)") +
  ggtitle("b) human") +
  ylim(-0.01,3) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.pack.pop.human



CDVrss.pack <- beta.se.wide[COD=='CDV' & term2=='boundaryDist']

CDVrss.pack.1 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(CDVrss.pack$beta + (CDVrss.pack$betaintx*hj.1))
CDVrss.pack.2 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(CDVrss.pack$beta + (CDVrss.pack$betaintx*hj.2))
CDVrss.pack.3 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(CDVrss.pack$beta + (CDVrss.pack$betaintx*hj.3))
CDVrss.pack.4 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(CDVrss.pack$beta + (CDVrss.pack$betaintx*hj.4))

lnrss.pack.pop.CDV =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.pack),y=(CDVrss.pack.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.pack),y=(CDVrss.pack.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.pack),y=(CDVrss.pack.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.pack),y=(CDVrss.pack.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to pack (m)") +
  ggtitle("c) CDV") +
  ylim(-0.01,3) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.pack.pop.CDV



nonerss.pack <- beta.se.wide[COD=='control' & term2=='boundaryDist']

nonerss.pack.1 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(nonerss.pack$beta + (nonerss.pack$betaintx*hj.1))
nonerss.pack.2 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(nonerss.pack$beta + (nonerss.pack$betaintx*hj.2))
nonerss.pack.3 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(nonerss.pack$beta + (nonerss.pack$betaintx*hj.3))
nonerss.pack.4 <- (log(hi.pack/(hi.pack-delta.hi.pack)))^(nonerss.pack$beta + (nonerss.pack$betaintx*hj.4))

lnrss.pack.pop.none =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.pack),y=(nonerss.pack.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.pack),y=(nonerss.pack.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.pack),y=(nonerss.pack.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.pack),y=(nonerss.pack.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to pack (m)") +
  ggtitle("a) control") +
  ylim(-0.01,3) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.pack.pop.none

(lnrss.road.pop.none|lnrss.road.pop.human|lnrss.road.pop.CDV)/(lnrss.nn.pop.none|lnrss.nn.pop.human|lnrss.nn.pop.CDV)/(lnrss.pack.pop.none|lnrss.pack.pop.human|lnrss.pack.pop.CDV)



#### habitat RSS ####

#### pop RSS road dist ----

humanrss.forest <- beta.se.wide[COD=='human' & term2=='forest']

# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100
hi.forest <- 1
delta.hi.forest <- seq(0, 1, .0001)

humanrss.forest.1 <- delta.hi.forest*(humanrss.forest$beta + (humanrss.forest$betaintx*log.hj.1))
humanrss.forest.2 <- delta.hi.forest*(humanrss.forest$beta + (humanrss.forest$betaintx*log.hj.2))
humanrss.forest.3 <- delta.hi.forest*(humanrss.forest$beta + (humanrss.forest$betaintx*log.hj.3))
humanrss.forest.4 <- delta.hi.forest*(humanrss.forest$beta + (humanrss.forest$betaintx*log.hj.4))


# humanrss.forest.1 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(humanrss.forest$beta + (humanrss.forest$betaintx*hj.1))
# humanrss.forest.2 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(humanrss.forest$beta + (humanrss.forest$betaintx*hj.2))
# humanrss.forest.3 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(humanrss.forest$beta + (humanrss.forest$betaintx*hj.3))
# humanrss.forest.4 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(humanrss.forest$beta + (humanrss.forest$betaintx*hj.4))

rss.forest.pop.human =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.forest),y=(humanrss.forest.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.forest),y=(humanrss.forest.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.forest),y=(humanrss.forest.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.forest),y=(humanrss.forest.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Forest") +
  ggtitle("b) human") +
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'none') + theme(legend.text = element_text(size = 20))

rss.forest.pop.human



CDVrss.forest <- beta.se.wide[COD=='CDV' & term2=='forest']

CDVrss.forest.1 <- delta.hi.forest*(CDVrss.forest$beta + (CDVrss.forest$betaintx*log.hj.1))
CDVrss.forest.2 <- delta.hi.forest*(CDVrss.forest$beta + (CDVrss.forest$betaintx*log.hj.2))
CDVrss.forest.3 <- delta.hi.forest*(CDVrss.forest$beta + (CDVrss.forest$betaintx*log.hj.3))
CDVrss.forest.4 <- delta.hi.forest*(CDVrss.forest$beta + (CDVrss.forest$betaintx*log.hj.4))

# 
# CDVrss.forest.1 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(CDVrss.forest$beta + (CDVrss.forest$betaintx*hj.1))
# CDVrss.forest.2 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(CDVrss.forest$beta + (CDVrss.forest$betaintx*hj.2))
# CDVrss.forest.3 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(CDVrss.forest$beta + (CDVrss.forest$betaintx*hj.3))
# CDVrss.forest.4 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(CDVrss.forest$beta + (CDVrss.forest$betaintx*hj.4))

rss.forest.pop.CDV =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.forest),y=(CDVrss.forest.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.forest),y=(CDVrss.forest.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.forest),y=(CDVrss.forest.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.forest),y=(CDVrss.forest.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Forest") +
  ggtitle("c) CDV") +
  
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = "right") + theme(legend.text = element_text(size = 20))

rss.forest.pop.CDV



nonerss.forest <- beta.se.wide[COD=='control' & term2=='forest']

nonerss.forest.1 <- delta.hi.forest*(nonerss.forest$beta + (nonerss.forest$betaintx*log.hj.1))
nonerss.forest.2 <- delta.hi.forest*(nonerss.forest$beta + (nonerss.forest$betaintx*log.hj.2))
nonerss.forest.3 <- delta.hi.forest*(nonerss.forest$beta + (nonerss.forest$betaintx*log.hj.3))
nonerss.forest.4 <- delta.hi.forest*(nonerss.forest$beta + (nonerss.forest$betaintx*log.hj.4))

# nonerss.forest.1 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(nonerss.forest$beta + (nonerss.forest$betaintx*hj.1))
# nonerss.forest.2 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(nonerss.forest$beta + (nonerss.forest$betaintx*hj.2))
# nonerss.forest.3 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(nonerss.forest$beta + (nonerss.forest$betaintx*hj.3))
# nonerss.forest.4 <- (log(hi.forest/(hi.forest-delta.hi.forest)))^(nonerss.forest$beta + (nonerss.forest$betaintx*hj.4))

rss.forest.pop.none =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.forest),y=(nonerss.forest.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.forest),y=(nonerss.forest.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.forest),y=(nonerss.forest.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.forest),y=(nonerss.forest.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Forest") +
  ggtitle("a) control") +
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'none') + theme(legend.text = element_text(size = 20))

rss.forest.pop.none





humanrss.wet <- beta.se.wide[COD=='human' & term2=='wet']

# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100
hi.wet <- 1
delta.hi.wet <- seq(0, 1, .0001)

humanrss.wet.1 <- delta.hi.wet*(humanrss.wet$beta + (humanrss.wet$betaintx*log.hj.1))
humanrss.wet.2 <- delta.hi.wet*(humanrss.wet$beta + (humanrss.wet$betaintx*log.hj.2))
humanrss.wet.3 <- delta.hi.wet*(humanrss.wet$beta + (humanrss.wet$betaintx*log.hj.3))
humanrss.wet.4 <- delta.hi.wet*(humanrss.wet$beta + (humanrss.wet$betaintx*log.hj.4))


# humanrss.wet.1 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(humanrss.wet$beta + (humanrss.wet$betaintx*hj.1))
# humanrss.wet.2 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(humanrss.wet$beta + (humanrss.wet$betaintx*hj.2))
# humanrss.wet.3 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(humanrss.wet$beta + (humanrss.wet$betaintx*hj.3))
# humanrss.wet.4 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(humanrss.wet$beta + (humanrss.wet$betaintx*hj.4))

rss.wet.pop.human =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.wet),y=(humanrss.wet.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.wet),y=(humanrss.wet.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.wet),y=(humanrss.wet.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.wet),y=(humanrss.wet.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Wet") +
  ggtitle("b) human") +
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'none') + theme(legend.text = element_text(size = 20))

rss.wet.pop.human



CDVrss.wet <- beta.se.wide[COD=='CDV' & term2=='wet']

CDVrss.wet.1 <- delta.hi.wet*(CDVrss.wet$beta + (CDVrss.wet$betaintx*log.hj.1))
CDVrss.wet.2 <- delta.hi.wet*(CDVrss.wet$beta + (CDVrss.wet$betaintx*log.hj.2))
CDVrss.wet.3 <- delta.hi.wet*(CDVrss.wet$beta + (CDVrss.wet$betaintx*log.hj.3))
CDVrss.wet.4 <- delta.hi.wet*(CDVrss.wet$beta + (CDVrss.wet$betaintx*log.hj.4))

# 
# CDVrss.wet.1 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(CDVrss.wet$beta + (CDVrss.wet$betaintx*hj.1))
# CDVrss.wet.2 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(CDVrss.wet$beta + (CDVrss.wet$betaintx*hj.2))
# CDVrss.wet.3 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(CDVrss.wet$beta + (CDVrss.wet$betaintx*hj.3))
# CDVrss.wet.4 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(CDVrss.wet$beta + (CDVrss.wet$betaintx*hj.4))

rss.wet.pop.CDV =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.wet),y=(CDVrss.wet.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.wet),y=(CDVrss.wet.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.wet),y=(CDVrss.wet.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.wet),y=(CDVrss.wet.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Wet") +
  ggtitle("c) CDV") +
  
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 20))

rss.wet.pop.CDV



nonerss.wet <- beta.se.wide[COD=='control' & term2=='wet']

nonerss.wet.1 <- delta.hi.wet*(nonerss.wet$beta + (nonerss.wet$betaintx*log.hj.1))
nonerss.wet.2 <- delta.hi.wet*(nonerss.wet$beta + (nonerss.wet$betaintx*log.hj.2))
nonerss.wet.3 <- delta.hi.wet*(nonerss.wet$beta + (nonerss.wet$betaintx*log.hj.3))
nonerss.wet.4 <- delta.hi.wet*(nonerss.wet$beta + (nonerss.wet$betaintx*log.hj.4))

# nonerss.wet.1 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(nonerss.wet$beta + (nonerss.wet$betaintx*hj.1))
# nonerss.wet.2 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(nonerss.wet$beta + (nonerss.wet$betaintx*hj.2))
# nonerss.wet.3 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(nonerss.wet$beta + (nonerss.wet$betaintx*hj.3))
# nonerss.wet.4 <- (log(hi.wet/(hi.wet-delta.hi.wet)))^(nonerss.wet$beta + (nonerss.wet$betaintx*hj.4))

rss.wet.pop.none =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.wet),y=(nonerss.wet.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.wet),y=(nonerss.wet.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.wet),y=(nonerss.wet.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.wet),y=(nonerss.wet.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Wet") +
  ggtitle("a) control") +
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'none') + theme(legend.text = element_text(size = 20))

rss.wet.pop.none





humanrss.open <- beta.se.wide[COD=='human' & term2=='open']

# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100
hi.open <- 1
delta.hi.open <- seq(0, 1, .0001)

humanrss.open.1 <- delta.hi.open*(humanrss.open$beta + (humanrss.open$betaintx*log.hj.1))
humanrss.open.2 <- delta.hi.open*(humanrss.open$beta + (humanrss.open$betaintx*log.hj.2))
humanrss.open.3 <- delta.hi.open*(humanrss.open$beta + (humanrss.open$betaintx*log.hj.3))
humanrss.open.4 <- delta.hi.open*(humanrss.open$beta + (humanrss.open$betaintx*log.hj.4))


# humanrss.open.1 <- (log(hi.open/(hi.open-delta.hi.open)))^(humanrss.open$beta + (humanrss.open$betaintx*hj.1))
# humanrss.open.2 <- (log(hi.open/(hi.open-delta.hi.open)))^(humanrss.open$beta + (humanrss.open$betaintx*hj.2))
# humanrss.open.3 <- (log(hi.open/(hi.open-delta.hi.open)))^(humanrss.open$beta + (humanrss.open$betaintx*hj.3))
# humanrss.open.4 <- (log(hi.open/(hi.open-delta.hi.open)))^(humanrss.open$beta + (humanrss.open$betaintx*hj.4))

rss.open.pop.human =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.open),y=(humanrss.open.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.open),y=(humanrss.open.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.open),y=(humanrss.open.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.open),y=(humanrss.open.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Open") +
  ggtitle("b) human") +
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'none') + theme(legend.text = element_text(size = 20))

rss.open.pop.human



CDVrss.open <- beta.se.wide[COD=='CDV' & term2=='open']

CDVrss.open.1 <- delta.hi.open*(CDVrss.open$beta + (CDVrss.open$betaintx*log.hj.1))
CDVrss.open.2 <- delta.hi.open*(CDVrss.open$beta + (CDVrss.open$betaintx*log.hj.2))
CDVrss.open.3 <- delta.hi.open*(CDVrss.open$beta + (CDVrss.open$betaintx*log.hj.3))
CDVrss.open.4 <- delta.hi.open*(CDVrss.open$beta + (CDVrss.open$betaintx*log.hj.4))

# 
# CDVrss.open.1 <- (log(hi.open/(hi.open-delta.hi.open)))^(CDVrss.open$beta + (CDVrss.open$betaintx*hj.1))
# CDVrss.open.2 <- (log(hi.open/(hi.open-delta.hi.open)))^(CDVrss.open$beta + (CDVrss.open$betaintx*hj.2))
# CDVrss.open.3 <- (log(hi.open/(hi.open-delta.hi.open)))^(CDVrss.open$beta + (CDVrss.open$betaintx*hj.3))
# CDVrss.open.4 <- (log(hi.open/(hi.open-delta.hi.open)))^(CDVrss.open$beta + (CDVrss.open$betaintx*hj.4))

rss.open.pop.CDV =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.open),y=(CDVrss.open.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.open),y=(CDVrss.open.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.open),y=(CDVrss.open.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.open),y=(CDVrss.open.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Open") +
  ggtitle("c) CDV") +
  
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'right') + theme(legend.text = element_text(size = 20))

rss.open.pop.CDV



nonerss.open <- beta.se.wide[COD=='control' & term2=='open']

nonerss.open.1 <- delta.hi.open*(nonerss.open$beta + (nonerss.open$betaintx*log.hj.1))
nonerss.open.2 <- delta.hi.open*(nonerss.open$beta + (nonerss.open$betaintx*log.hj.2))
nonerss.open.3 <- delta.hi.open*(nonerss.open$beta + (nonerss.open$betaintx*log.hj.3))
nonerss.open.4 <- delta.hi.open*(nonerss.open$beta + (nonerss.open$betaintx*log.hj.4))

# nonerss.open.1 <- (log(hi.open/(hi.open-delta.hi.open)))^(nonerss.open$beta + (nonerss.open$betaintx*hj.1))
# nonerss.open.2 <- (log(hi.open/(hi.open-delta.hi.open)))^(nonerss.open$beta + (nonerss.open$betaintx*hj.2))
# nonerss.open.3 <- (log(hi.open/(hi.open-delta.hi.open)))^(nonerss.open$beta + (nonerss.open$betaintx*hj.3))
# nonerss.open.4 <- (log(hi.open/(hi.open-delta.hi.open)))^(nonerss.open$beta + (nonerss.open$betaintx*hj.4))

rss.open.pop.none =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi.open),y=(nonerss.open.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi.open),y=(nonerss.open.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi.open),y=(nonerss.open.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi.open),y=(nonerss.open.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Open") +
  ggtitle("a) control") +
  #ylim(-0.01,2.5) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = 'none') + theme(legend.text = element_text(size = 20))

rss.open.pop.none


(rss.forest.pop.none|rss.forest.pop.human|rss.forest.pop.CDV)/(rss.wet.pop.none|rss.wet.pop.human|rss.wet.pop.CDV)/(rss.open.pop.none|rss.open.pop.human|rss.open.pop.CDV)

#########
beta <- merge(fullOUT, beta.se[,.(term, COD, se)], by = c('term', 'COD'))
beta[,'lwr'] <- beta$estimate - (beta$se*1.96)
beta[,'upr'] <- beta$estimate + (beta$se*1.96)
beta[,'var'] <- ifelse(beta$term %like% 'ttd', 'intx','var')
beta[,'term2'] <- gsub("-ttd", "", beta$term)

beta.wide <- dcast(beta, wolfID + term2 ~ var, value.var = c('estimate'))
beta.wide <- plyr::rename(beta.wide, c('var'='beta', 'intx'='betaintx'))

lwr.wide <- dcast(beta, wolfID + term2 ~ var, value.var = c('lwr'))
lwr.wide <- plyr::rename(lwr.wide, c('var'='lwr', 'intx'='lwrintx'))

upr.wide <- dcast(beta, wolfID + term2 ~ var, value.var = c('upr'))
upr.wide <- plyr::rename(upr.wide, c('var'='upr', 'intx'='uprintx'))

### RSS nndist ----


delta.hi <- log(1:5000)

#ttd 1mo, 2 weeks, 1 day
# 1day
hj.1 <- log(1 + 1)
# 2wks
hj.2 <- log(1 + 14)
# 1mo
hj.3 <- log(1 + 30)
# 2mo 
hj.4 <- log(1 + 60)


### simpler RSS for not log-transformed betas Δhi ⋅(βi + βij ⋅hj (x1 )] 
rss.nn <- beta.wide[term2=='nnDist']
# rss.nn<- rss.nn[, .(oneD = rss.intx(beta, betaintx, delta.hi, hj.1), twoW = rss.intx(beta, betaintx, delta.hi, hj.2),
#                     oneM = rss.intx(beta, betaintx, delta.hi, hj.3), twoM = rss.intx(beta, betaintx, delta.hi, hj.4)), by = .(wolfID)]

W02<- rss.nn[wolfID=='RMNP_W02']
W26<- rss.nn[wolfID=='RMNP_W26']
rss.nn.1 <- delta.hi*(W26$beta + (W26$betaintx*hj.1))
rss.nn.2 <- delta.hi*(W26$beta + (W26$betaintx*hj.2))
rss.nn.3 <- delta.hi*(W26$beta + (W26$betaintx*hj.3))
rss.nn.4 <- delta.hi*(W26$beta + (W26$betaintx*hj.4))

# I've negated RSS value here so it's more intuitive to understand distance to selecting to be farther
RMNPW26.nn =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(rss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(rss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(rss.nn.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi),y=(rss.nn.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("log Distance to NN (m)") +
  #r = r + ylim(-0.01,3)
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

RMNPW26.nn



### RSS for two log-transformed interaction betas


hj.1 <- 1
# 2wks
hj.2 <- 14
# 1mo
hj.3 <- 30
# 2mo
hj.4 <- 60

# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100
hi <- 5001
delta.hi <- 1:5000


lnrss.nn.1 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.1))
lnrss.nn.2 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.2))
lnrss.nn.3 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.3))
lnrss.nn.4 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.4))

lnrss.nn.W02 =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(lnrss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(lnrss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(lnrss.nn.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi),y=(lnrss.nn.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to NN (m)") +
  ylim(-0.01,10) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.nn.W02


#### pop RSS ----
beta.se[,'term2'] <- gsub("-ttd", "", beta.se$term)

beta.se.wide <- dcast(beta.se, COD + term2 ~ var, value.var = c('mean'))
beta.se.wide <- plyr::rename(beta.se.wide, c('var'='beta', 'intx'='betaintx'))

beta.se.cdv <- beta.se.wide[COD=='CDV']
poprss.nn <- beta.se.cdv[term2=='nnDist']

pop.lnrss.nn.1 <- (log(hi/(hi-delta.hi)))^(poprss.nn$beta + (poprss.nn$betaintx*hj.1))
pop.lnrss.nn.2 <- (log(hi/(hi-delta.hi)))^(poprss.nn$beta + (poprss.nn$betaintx*hj.2))
pop.lnrss.nn.3 <- (log(hi/(hi-delta.hi)))^(poprss.nn$beta + (poprss.nn$betaintx*hj.3))
pop.lnrss.nn.4 <- (log(hi/(hi-delta.hi)))^(poprss.nn$beta + (poprss.nn$betaintx*hj.4))

lnrss.nn.pop =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(pop.lnrss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(pop.lnrss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(pop.lnrss.nn.3), colour = "30 days"),size = 1) +
  geom_line(aes(x=(delta.hi),y=(pop.lnrss.nn.4), colour = "60 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to NN (m)") +
  ylim(-0.01,20) +
  scale_colour_manual("", values = c("gray", "black", "gray33", 'blue'))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

lnrss.nn.pop

######

popcdv <- readRDS('data/derived-data/popcdv.Rds')
pop.cdv <- as.data.frame(popcdv)
popcdv.nn <- data.frame(term=rownames(pop.cdv), beta=pop.cdv[,1], se=pop.cdv[,2], row.names=NULL)
popcdv.nn[,'var'] <- ifelse(popcdv.nn$term%like% 'ttd', 'intx','var')
popcdv.nn <-setDT(popcdv.nn)[term %like% 'distance2']
popcdv.nn[,'sel'] <- 'nnDist'

popcdv.wide <- dcast(popcdv.nn, sel ~ var, value.var = c('beta'))
popcdv.wide <- popcdv.wide[,.(sel, beta= var, betaintx=intx)]

rss.nn.1 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.1))
rss.nn.2 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.2))
rss.nn.3 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.3))
rss.nn.4 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.4))

rss.nn =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(-rss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(-rss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(-rss.nn.3), colour = "30 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("log Distance to NN (m)") +
  #r = r + ylim(-0.01,3)
  scale_colour_manual("", values = c("gray", "black", "gray33"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

rss.nn



## pop cdv last mo
pop.cdv.lastmo <- as.data.frame(popcdv.lastmo)
pop.cdv.lastmo <- data.frame(term=rownames(pop.cdv.lastmo), beta=pop.cdv.lastmo[,1], se=pop.cdv.lastmo[,2], row.names=NULL)
pop.cdv.lastmo[,'var'] <- ifelse(pop.cdv.lastmo$term%like% 'ttd', 'intx','var')
popcdv.lastmo.nn <-setDT(pop.cdv.lastmo)[term %like% 'distance2']
popcdv.lastmo.nn[,'sel'] <- 'nnDist'

popcdv.lastmo.wide <- dcast(popcdv.lastmo.nn, sel ~ var, value.var = c('beta'))
popcdv.lastmo.wide <- popcdv.lastmo.wide[,.(sel, beta= var, betaintx=intx)]

rss.nn.1 <- delta.hi*(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.1))
rss.nn.2 <- delta.hi*(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.2))
rss.nn.3 <- delta.hi*(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.3))


rss.nn.lastmo =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(rss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(rss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(rss.nn.3), colour = "30 days"),size = 1) +
  theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
              axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("log Distance to NN (m)") +
#r = r + ylim(-0.01,3)
 scale_colour_manual("", values = c("gray", "black", "gray33"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

rss.nn.lastmo





popcdv.lastmo.pack <-setDT(pop.cdv.lastmo)[term %like% 'pack']
popcdv.lastmo.pack[,'sel'] <- 'packDist'

popcdv.lastmo.pack.wide <- dcast(popcdv.lastmo.pack, sel ~ var, value.var = c('beta'))
popcdv.lastmo.pack.wide <- popcdv.lastmo.pack.wide[,.(sel, beta= var, betaintx=intx)]

rss.pack.1 <- delta.hi*(popcdv.lastmo.pack.wide$beta + (popcdv.lastmo.pack.wide$betaintx*hj.1))
rss.pack.2 <- delta.hi*(popcdv.lastmo.pack.wide$beta + (popcdv.lastmo.pack.wide$betaintx*hj.2))
rss.pack.3 <- delta.hi*(popcdv.lastmo.pack.wide$beta + (popcdv.lastmo.pack.wide$betaintx*hj.3))


rss.pack.lastmo =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(rss.pack.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(rss.pack.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(rss.pack.3), colour = "30 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("log Distance to pack boundary (m)") +
  #r = r + ylim(-0.01,3)
  scale_colour_manual("", values = c("gray", "black", "gray33"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

rss.pack.lastmo

### RSS for two log-transformed interaction betas

# delta hi should be based off avg/median nndist
hi <- 5001
delta.hi <- 1:5000

hj.1 <- 1
# 2wks
hj.2 <- 14
# 1mo
hj.3 <- 30

rss.nn.1 <- (log(hi/(hi-delta.hi)))^(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.1))
rss.nn.2 <- (log(hi/(hi-delta.hi)))^(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.2))
rss.nn.3 <- (log(hi/(hi-delta.hi)))^(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.3))


rss.log.nn.lastmo =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7) + 
  geom_line(aes(x=(delta.hi),y=(rss.nn.1), colour = "1 day"), size = 1) + 
  geom_line(aes(x=(delta.hi),y=(rss.nn.2), colour = "14 days"), size = 1) +
  geom_line(aes(x=(delta.hi),y=(rss.nn.3), colour = "30 days"),size = 1) +
  theme_bw()  + theme(
    #panel.background =element_rect(colour = "black", fill=NA, size=1),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20)) +
  theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) +
  ylab("RSS") + xlab("Distance to NN (m)") +
  ylim(-0.01,3) +
  scale_colour_manual("", values = c("gray", "black", "gray33"))  +  
  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

print(rss.log.nn.lastmo)
