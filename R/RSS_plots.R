# The input files here are RSF model outputs from broom::tidy:
# 1) Three different population level models: 'gFR1', 'HRgFR1', and 'pop'
# 2) A data table with individual level models, including 47 individuals: 'indivdataset'

# These RSS plots are relative selection for forest at two locations differing 
# in their availability of roads in elk home ranges

# Each model was run 120 times (something you may be doing if bootstrapping) so 
# the first part is just getting the means and CIs of the selection coefficients
# I wanted in order to plot single lines on the RSS plots.

######################################################
### 1 Mean selection coefficients for RSS calculations
######################################################

# Subset out GFR1 deciduous selection coefficient
hx <- gFR1[term=="deciduous"]
meanhx <- mean(hx$estimate)
upperhx <- meanhx+1.96*stderror(hx$estimate)
lowerhx <- meanhx-1.96*stderror(hx$estimate)
# Subset out GFR1 deciduous x road selection coefficient
hxij <- gFR1[term=="deciduous:avgDistToRoadLn"]
meanhxij <- mean(hxij$estimate)
upperhxij <- meanhxij+1.96*stderror(hxij$estimate)
lowerhxij <- meanhxij-1.96*stderror(hxij$estimate)
# Subset out HRGFR1 deciduous selection coefficient
HRhx <- HRgFR1[term=="deciduous"]
meanHRhx <- mean(HRhx$estimate)
upperHRhx <- meanHRhx+1.96*stderror(HRhx$estimate)
lowerHRhx <- meanHRhx-1.96*stderror(HRhx$estimate)
# Subset out HRGFR1 deciduous x road selection coefficient
HRhxij <- HRgFR1[term=="deciduous:avgDistToRoadLn"]
meanHRhxij <- mean(HRhxij$estimate)
upperHRhxij <- meanHRhxij+1.96*stderror(HRhxij$estimate)
lowerHRhxij <- meanHRhxij-1.96*stderror(HRhxij$estimate)

# Subset out population deciduous selection coefficient
hxpop <- pop[term=="deciduous"]
meanhxpop <- mean(hxpop$estimate)

##########################################################################
### 2 Subset road covariates and deciduous area from the individual models
##########################################################################

# Loop through all IDs to form data.table summary of mean betas and HR cover by ID
meanByID <- data.table()
for (i in unique(indivdataset$idNumb)){
  h <- subset(indivdataset, idNumb==i)
  j <- h[,!c(2:5)]
  k.mean <- j[, lapply(.SD, mean), by = covariate]
  l.mean <- cbind(i,k.mean)
  meanByID <- rbind(meanByID, l.mean)
}
# Loop through all IDs to form data.table summary of standard errors by ID
seByID <- data.table()
for (i in unique(indivdataset$idNumb)){
  h <- subset(indivdataset, idNumb==i)
  j <- h[,!c(2:5)]
  k.se <- j[, lapply(.SD, stderror), by = covariate]
  k.se$beta <- k.se$beta*1.96
  l.se <- cbind(i,k.se)
  seByID <- rbind(seByID, l.se)
}

# merge mean and SE
seByID <- seByID[,!c(4,5)]
colnames(seByID) <- c('i', 'covariate', 'confint')
meanByID <- merge(meanByID, seByID, by=c('i','covariate'))

# Selection by mean
# Combine average road distance with deciduous selection covariate by iteration
avgRoad <- meanByID[covariate=="distToRoadLn"]
avgRoad <- log(avgRoad[,4])
avgRoad$iteration <- rep(1:43)
betaDeciduous <- meanByID[covariate=="deciduous"]
betaDeciduous <- betaDeciduous[,3]
betaDeciduous$iteration <- rep(1:43)
seDeciduous <- meanByID[covariate=="deciduous"]
seDeciduous <- seDeciduous[,6]
seDeciduous$iteration <- rep(1:43)
# Combine into dataset
indivData <- merge(avgRoad, betaDeciduous, by="iteration")
indivData <- merge(indivData, seDeciduous, by="iteration")
# create upper and lower CI bounds
indivData$lwr <- indivData$beta-indivData$confint
indivData$upr <- indivData$beta+indivData$confint

###########################################
### 3 Calculate relative selection strength
###########################################

# Population RSS at three levels of road availability
deciduous.range <- seq(1,0,length.out=100)
delta.hi <- 1-deciduous.range

# road availability in home ranges 
# (these are just the max, min, and median distance
# to roads in my home ranges)
# high
hj.1 <- 8.7
# intermediate
hj.2 <- 7.3
# low
hj.3 <- 6.5

# This is the RSS equation I'm using: Δhi⋅(βi +βij⋅hj(x1)
# because my coefficients are interactions between roads and forest

### A) Population level RSS:
# GFR rss at high, medium, low deciduous forest
rss.high <- delta.hi*(meanhx+meanhxij*hj.1)
rss.mid <- delta.hi*(meanhx+meanhxij*hj.2)
rss.low <- delta.hi*(meanhx+meanhxij*hj.3)

# HRGFR rss at high, medium, low deciduous forest
rss.high.HR <- delta.hi*(meanHRhx+meanHRhxij*hj.1)
rss.mid.HR <- delta.hi*(meanHRhx+meanHRhxij*hj.2)
rss.low.HR <- delta.hi*(meanHRhx+meanHRhxij*hj.3)

# Upper/lower bounds of GFR rss at high, medium, low deciduous forest
rss.high.upr <- delta.hi*(upperhx+upperhxij*hj.1)
rss.mid.upr <- delta.hi*(upperhx+upperhxij*hj.2)
rss.low.upr <- delta.hi*(upperhx+upperhxij*hj.3)
rss.high.lwr <- delta.hi*(lowerhx+lowerhxij*hj.1)
rss.mid.lwr <- delta.hi*(lowerhx+lowerhxij*hj.2)
rss.low.lwr <- delta.hi*(lowerhx+lowerhxij*hj.3)

# Upper/lower bounds of HRGFR rss at high, medium, low deciduous forest
rss.high.HR.upr <- delta.hi*(upperHRhx+upperHRhxij*hj.1)
rss.mid.HR.upr <- delta.hi*(upperHRhx+upperHRhxij*hj.2)
rss.low.HR.upr <- delta.hi*(upperHRhx+upperHRhxij*hj.3)
rss.high.HR.lwr <- delta.hi*(lowerHRhx+lowerHRhxij*hj.1)
rss.mid.HR.lwr <- delta.hi*(lowerHRhx+lowerHRhxij*hj.2)
rss.low.HR.lwr <- delta.hi*(lowerHRhx+lowerHRhxij*hj.3)

# predicted rss from population model
rss.pop <- delta.hi*meanhxpop

### B) Individual level RSS:
# Relative selection strength for each individual as they move from low to high NDVI
rss.indivs <- data.table()
for (i in row(indivData)){
  rss <- delta.hi*indivData$beta[i]
  rss.lwr <- delta.hi*indivData$lwr[i]
  rss.upr <- delta.hi*indivData$upr[i]
  cov <- rep(indivData$avgHR[i], length(rss))
  elk <- as.data.table(cbind(i, delta.hi, rss, rss.lwr, rss.upr, cov))
  rss.indivs <- rbind(elk, rss.indivs)
}

######################################
### 4 Plot relative selection strength
######################################

# Preparing data frames for plots
#rename columns
colnames(rss.indivs) <- c('idNumb', 'deltaDeciduous', 'rss', 'rss.lwr', 'rss.upr', 'roadsavgHR')

#create data table for population rss
#GFR
rss.avgs<- as.data.table(cbind(rss.high, rss.high.upr, rss.high.lwr, rss.mid, rss.mid.upr,
                               rss.mid.lwr, rss.low, rss.low.upr, rss.low.lwr, rss.pop, delta.hi))
#HRGFR
rss.avgs.HR<- as.data.table(cbind(rss.high.HR, rss.high.HR.upr, rss.high.HR.lwr, rss.mid.HR, rss.mid.HR.upr,
                                  rss.mid.HR.lwr, rss.low.HR, rss.low.HR.upr, rss.low.HR.lwr, rss.pop, delta.hi))
#optionally melt data.tables
#rss.avgs <- melt(rss.avgs, id.vars='delta.hi')
#colnames(rss.avgs) <- c('delta.hi','cover.scale', 'value')

# Plot OFR (one of the population level models)
gfr1aplot <- ggplot(data=rss.indivs, aes(x=deltaDeciduous, y=rss))
# plot individual lines and CIs
gfr1aplot <- gfr1aplot + geom_line(aes(group=idNumb, colour=roadsavgHR), size=0.2, alpha=0.5, col="#3D4849")
gfr1aplot <- gfr1aplot + geom_ribbon(aes(group=idNumb, ymin=rss.lwr,ymax=rss.upr), alpha=0.2, fill="#3D4849")
# plot pop nmodel (another population level model)
gfr1aplot <- gfr1aplot + geom_line(data=rss.avgs.HR, aes(x=delta.hi, y=rss.pop), size=2, linetype="solid", col="#56B4E9")
# plot model lines and CIs
gfr1aplot <- gfr1aplot + geom_line(data=rss.avgs, aes(x=delta.hi, y=rss.high), size=1, linetype="dashed", col="#CC79A7")
gfr1aplot <- gfr1aplot + geom_line(data=rss.avgs, aes(x=delta.hi, y=rss.mid), size=1, linetype="dotdash", col="#CC79A7")
gfr1aplot <- gfr1aplot + geom_line(data=rss.avgs, aes(x=delta.hi, y=rss.low), size=1, linetype="dotted", col="#CC79A7")
gfr1aplot <- gfr1aplot + geom_ribbon(data=rss.avgs.HR, aes(x=delta.hi, ymin=rss.high.lwr,ymax=rss.high.upr),alpha=0.2, fill="#CC79A7")
gfr1aplot <- gfr1aplot + geom_ribbon(data=rss.avgs.HR, aes(x=delta.hi, ymin=rss.mid.lwr,ymax=rss.mid.upr),alpha=0.2, fill="#CC79A7")
gfr1aplot <- gfr1aplot + geom_ribbon(data=rss.avgs.HR, aes(x=delta.hi, ymin=rss.low.lwr,ymax=rss.low.upr),alpha=0.2, fill="#CC79A7")
# zero line
gfr1aplot <- gfr1aplot + geom_abline(intercept=0, slope=0, col='black')
# plot aesthetics
gfr1aplot <- gfr1aplot + theme(panel.grid.major = element_blank())  + theme(panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill="white")) + theme(axis.line= element_line(colour="black"))
gfr1aplot <- gfr1aplot + theme(axis.text = element_text(size=15, colour="black")) +
  theme(axis.title.x = element_blank()) + theme(axis.title.y=element_text(size=18, colour='black'))
theme(axis.text = element_text(size=18, colour="black"))
# axis labels
gfr1aplot <- gfr1aplot + ylab("RSS for forest at x1")
# annotate letter
gfr1aplot <- gfr1aplot + annotate("text", x=0.1, y=6.7, label="a", size=10)
# axis limits
gfr1aplot <- gfr1aplot + ylim(-4,7)
# add legend
gfr1aplot <- gfr1aplot + annotate("text", x=0.25, y=5.3, label="Avg. dist. to road:", size=5) 
gfr1aplot <- gfr1aplot + geom_segment(aes(x = 0, y = 4.5, xend = .1, yend = 4.5), linetype='dashed')
gfr1aplot <- gfr1aplot + annotate("text", x=0.25, y=4.5, label="6 km", size=5)
gfr1aplot <- gfr1aplot + geom_segment(aes(x = 0, y = 4, xend = .1, yend = 4), linetype='dotdash')
gfr1aplot <- gfr1aplot + annotate("text", x=0.25, y=4, label="1.5 km", size=5)
gfr1aplot <- gfr1aplot + geom_segment(aes(x = 0, y = 3.5, xend = .1, yend = 3.5), linetype='dotted')
gfr1aplot <- gfr1aplot + annotate("text", x=0.25, y=3.5, label="0.5 km", size=5)
