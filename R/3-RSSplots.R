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

fullOUT <- readRDS('data/derived-data/full_betas.Rds')


### Function ----
se <- function(x){
  sd(x, na.rm = T)/ sqrt(length(na.omit(x)))
}

### simpler RSS for intx terms -- not log-transformed betas Δhi ⋅(βi + βij ⋅hj (x1 )] 
rss.intx <- function(x, xintx, delta, h){
  delta*(x + (xintx*h))
}



### create CIs ----

beta.se <- fullOUT[,.(se= se(estimate)),.(term)]

beta <- merge(fullOUT, beta.se, by = 'term')
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
r = r + geom_line(aes(x=(delta.hi),y=(rss.nn.4), colour = "60 days"),size = 1) 
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
                           values = c("gray", "black", "gray33", 'purple'))  
r = r +  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

print(r)


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
popcdv.nn <- data.frame(term=names(popcdv), beta=popcdv, row.names=NULL)
popcdv.nn[,'var'] <- ifelse(popcdv.nn$term%like% 'ttd', 'intx','var')
popcdv.nn <-setDT(popcdv.nn)[term %like% 'distance2']
popcdv.nn[,'sel'] <- 'nnDist'

popcdv.wide <- dcast(popcdv.nn, sel ~ var, value.var = c('beta'))
popcdv.wide <- popcdv.wide[,.(sel, beta= var, betaintx=intx)]

rss.nn.1 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.1))
rss.nn.2 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.2))
rss.nn.3 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.3))
rss.nn.4 <- delta.hi*(popcdv.wide$beta + (popcdv.wide$betaintx*hj.4))

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
r = r + geom_line(aes(x=(delta.hi),y=(rss.nn.4), colour = "60 days"),size = 1) 
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
                           values = c("gray", "black", "gray33", 'purple'))  
r = r +  theme(legend.key = element_blank()) + theme(legend.position = c(.75,.9)) + theme(legend.text = element_text(size = 20))

print(r)



## pop cdv last mo
popcdv.lastmo.nn <- data.frame(term=names(popcdv.lastmo), beta=popcdv.lastmo, row.names=NULL)
popcdv.lastmo.nn[,'var'] <- ifelse(popcdv.lastmo.nn$term%like% 'ttd', 'intx','var')
popcdv.lastmo.nn <-setDT(popcdv.lastmo.nn)[term %like% 'distance2']
popcdv.lastmo.nn[,'sel'] <- 'nnDist'

popcdv.lastmo.wide <- dcast(popcdv.lastmo.nn, sel ~ var, value.var = c('beta'))
popcdv.lastmo.wide <- popcdv.lastmo.wide[,.(sel, beta= var, betaintx=intx)]

rss.nn.1 <- delta.hi*(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.1))
rss.nn.2 <- delta.hi*(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.2))
rss.nn.3 <- delta.hi*(popcdv.lastmo.wide$beta + (popcdv.lastmo.wide$betaintx*hj.3))


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


