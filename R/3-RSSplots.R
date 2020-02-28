# RSS mixed effects output
# 17 February 2020
# created by Julie Turner

### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival','forcats')
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

beta.se <- fullOUT[,.(mean = mean(estimate), se= se(estimate)),.(term, COD)]
beta.se[,'lwr'] <- beta.se$mean - (beta.se$se*1.96)
beta.se[,'upr'] <- beta.se$mean + (beta.se$se*1.96)
beta.se[,'var'] <- ifelse(beta.se$term %like% 'ttd', 'intx','var')

pd <- position_dodge(0.5) # move them .05 to the left and right
ggplot(beta.se[var=='var'], aes(x=term, y=mean, colour=COD)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
    geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)

ggplot(beta.se[var=='intx'], aes(x=term, y=mean, colour=COD)) + 
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.1, position=pd) +
  geom_point(position=pd) +
  theme_bw()  + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = .7)) +
  geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)

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


#### NOT WORKING ####
### RSS for two log-transformed interaction betas

# delta hi should be based off avg/median nndist
hi <- 101:1000
delta.hi <- 100

hj.1 <- 1
# 2wks
hj.2 <- 14
# 1mo
hj.3 <- 30

rss.nn.1 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.1))
rss.nn.2 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.2))
rss.nn.3 <- (log(hi/(hi-delta.hi)))^(W02$beta + (W02$betaintx*hj.3))

r =  ggplot() + geom_hline(yintercept = 0,colour = "black",lty = 2, size = .7)
r = r + geom_line(aes(x=(delta.hi),y=(rss.nn.1), colour = "1 day"), size = 1) 
#r = r + geom_line(aes(x=(hj-250),y=log(rssroadd_lo), colour = "Day"), size = 1, lty = 3) 
#r = r + geom_line(aes(x=(hj-250),y=log(rssroadd_hi), colour = "Day"), size = 1, lty = 3) 
r = r + geom_line(aes(x=(delta.hi),y=(rss.nn.2), colour = "14 days"), size = 1) 
# r = r + geom_line(aes(x=(hj-250),y=log(rssroadt_lo), colour = "Twilight"), size = 1, lty = 3) 
# r = r + geom_line(aes(x=(hj-250),y=log(rssroadt_hi), colour = "Twilight"), size = 1, lty = 3) 
r = r + geom_line(aes(x=(delta.hi),y=(rss.nn.3), colour = "30 days"),size = 1) 
# r = r + geom_line(aes(x=(hj-250),y=log(rssroadn_lo), colour = "Night"),size = 1, lty = 3) 
# r = r + geom_line(aes(x=(hj-250),y=log(rssroadn_hi), colour = "Night"),size = 1, lty = 3) 
r = r + theme_bw()  + theme(
  #panel.background =element_rect(colour = "black", fill=NA, size=1),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = .7))
r = r + theme(plot.title=element_text(size=20,hjust = 0.05),axis.text.x = element_text(size=20), axis.title = element_text(size=25),axis.text.y = element_text(size=20))
r = r + theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
              axis.text.y = element_text(margin=margin(10,10,10,10,"pt")))+ theme(axis.ticks.length = unit(-0.25, "cm")) 
r = r + ylab("RSS") + xlab("log Distance to NN (m)")
#r = r + ylim(-0.01,3)
r = r +scale_colour_manual("", 
                           values = c("gray", "black", "gray33"))  
r = r +  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

print(r)


#### pop RSS ----
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
