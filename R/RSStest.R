# RSS test
# 15 Oct. 2019

### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'amt', 'lubridate', 'tidyr', 'ggplot2','survival','forcats')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'


#  ln[hi(x1)/(hi(x1)-Δhi)][βi +βij*hj(x1)]
# varying x1

# nndist <- socpropland1moOUT[,.(wolfID, term, betann = nndist, betaintx=`nndist:ttd`)]
# nndist.w <- dcast(nndist,wolfID~ term, value.var = c( 'betann', 'betaintx'))
# nndist.w[,'betann_lwr'] <- nndist.w$betann_coef - (nndist.w$betann_se*1.96)
# nndist.w[,'betann_upr'] <- nndist.w$betann_coef + (nndist.w$betann_se*1.96)
# nndist.w[,'betaintx_lwr'] <- nndist.w$betaintx_coef - (nndist.w$betaintx_se*1.96)
# nndist.w[,'betaintx_upr'] <- nndist.w$betaintx_coef + (nndist.w$betaintx_se*1.96)
# 
# saveRDS(nndist.w, 'data/derived-data/socprop1moCoef.Rds')
nndist.w <- readRDS('data/derived-data/socprop1moCoef.Rds')

W02 <- nndist.w[wolfID=='RMNP_W02',.(betann_coef, betann_lwr, betann_upr, betaintx_coef, betaintx_lwr, betaintx_upr)]

W26 <- nndist.w[wolfID=='RMNP_W26',.(betann_coef, betann_lwr, betann_upr, betaintx_coef, betaintx_lwr, betaintx_upr)]

W06 <- nndist.w[wolfID=='RMNP_W06',.(betann_coef, betann_lwr, betann_upr, betaintx_coef, betaintx_lwr, betaintx_upr)]

W09 <- nndist.w[wolfID=='GHA26_W09',.(betann_coef, betann_lwr, betann_upr, betaintx_coef, betaintx_lwr, betaintx_upr)]
W27 <- nndist.w[wolfID=='RMNP_W27',.(betann_coef, betann_lwr, betann_upr, betaintx_coef, betaintx_lwr, betaintx_upr)]

# Population RSS at three levels of road availability
# delta hi should be based off avg/median nndist
# hi <- 101:1000
# delta.hi <- 100

#### go back and change to simpler equation and log everything in it to make scales same

# hi <- 1001:5000 # territories max ~5km probably closer 2-3km (find avg diameter of HR)

delta.hi <- log(1:5000)

#ttd 1mo, 2 weeks, 1 day
# 1day
hj.1 <- log(1 + 1)
# 2wks
hj.2 <- log(1 + 14)
# 1mo
hj.3 <- log(1 + 30)

# # 1day
# hj.1 <- 1
# # 2wks
# hj.2 <- 14
# # 1mo
# hj.3 <- 30

# ### RSS for two log-transformed interaction betas
# rss.nn.1 <- (log(hi/(hi-delta.hi)))^(W26$betann_coef + (W26$betaintx_coef*hj.1))
# rss.nn.2 <- (log(hi/(hi-delta.hi)))^(W26$betann_coef + (W26$betaintx_coef*hj.2))
# rss.nn.3 <- (log(hi/(hi-delta.hi)))^(W26$betann_coef + (W26$betaintx_coef*hj.3))

### simpler RSS for not log-transformed betas Δhi ⋅(βi + βij ⋅hj (x1 )] 
rss.nn.1 <- delta.hi*(W27$betann_coef + (W27$betaintx_coef*hj.1))
rss.nn.2 <- delta.hi*(W27$betann_coef + (W27$betaintx_coef*hj.2))
rss.nn.3 <- delta.hi*(W27$betann_coef + (W27$betaintx_coef*hj.3))


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





