#### adding prop landcover, but move to data prep ####
# Generate spatial points object from UTM coordinates
locs <- sp::SpatialPoints(data.frame(dat$x1_, dat$y1_))
locs_end <- sp::SpatialPoints(data.frame(dat$x2_, dat$y2_))

# Read in the landcover map with different habitat types
lc <-  raster(paste0(raw, 'RMNPlandcover.tif'), )
cland <- fread(paste0(raw, 'rcl_cowu.csv'))
lc <- raster::reclassify(lc, cland)
# To get rid of NAs, I make them a new category (10), and then add it to the legend
lc[is.na(lc)] <- 5

closed <- lc == 1
names(closed) <- "closed"
open <- lc == 2
names(open) <- "open"
wet <- lc == 3
names(wet) <- "wet"

## This creates an object which can be used to make a layer of specified diameter
# The d value is what determines the buffer size if you want to change it.
## If you're doing multiple landcover classes, you only need to run this line once, as long as each of the habitat variables has the same resolution
# (e.g., the "Wetland" here could be any cover type)

Buff100 <- focalWeight(closed, d=100, type = 'circle')

## This generates a new raster where each cell corresponds to the mean wetland within a 100m buffer.
# Since it's all 1s and 0s, this is the same as the proportion of wetland surrounding the focal variable
wetBuff100 <- focal(wet, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
closedBuff100 <- focal(closed, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)
openBuff100 <- focal(open, Buff100, na.rm = TRUE, pad = TRUE, padValue = 0)

## You can then use the "extract" function as you would on any habitat layer and add it to your data:
dat[,'propwet_start'] <- extract(wetBuff100, locs)
dat[,'propclosed_start'] <- extract(closedBuff100, locs)
dat[,'propopen_start'] <- extract(openBuff100, locs)

dat[,'propwet_end'] <- extract(wetBuff100, locs_end)
dat[,'propclosed_end'] <- extract(closedBuff100, locs_end)
dat[,'propopen_end'] <- extract(openBuff100, locs_end)

#saveRDS(dat, 'data/derived-data/ssfAllCov.Rds')


dat.avail <- dat[ua=='avail', .(ttd1, wetavail= sum(wet=='wet')/9, closedavail= sum(land_end=='closed')/9, openavail= sum(land_end=='open')/9, 
                                propwetavail= sum(propwet_end)/9, propclosedavail= sum(propclosed_end)/9, propopenavail= sum(propopen_end)/9), by=.(id, step_id_)]

dat.avail[,cor(wetavail,closedavail), by = .(id)]
dat.avail[,cor(propwetavail,propclosedavail), by = .(id)]


dat.wet.probs <- dat[id =='W06' | id =='W09' | id =='W10' | id =='W11']
data.probs.avail <- dat.avail[id =='W06' | id =='W09' | id =='W10' | id =='W11']
par(mfrow = c(1,3))
plot(wetavail ~ ttd1, data = dat.avail[id=='W06'], main = 'W06') #low R2
abline(lm(wetavail ~ ttd1, dat.avail[id=='W06']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W06']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W06'], main = 'W06')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W06']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W06']))

plot(openavail ~ ttd1, data = dat.avail[id=='W06'], main = 'W06')
abline(lm(openavail ~ ttd1, dat.avail[id=='W06']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W06']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W09'], main = 'W09') #low R2
abline(lm(wetavail ~ ttd1, dat.avail[id=='W09']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W09']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W09'], main = 'W09')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W09']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W09']))

plot(openavail ~ ttd1, data = dat.avail[id=='W09'], main = 'W09')
abline(lm(openavail ~ ttd1, dat.avail[id=='W09']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W09']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W10'], main = 'W10') # R2 below .15
abline(lm(wetavail ~ ttd1, dat.avail[id=='W10']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W10']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W10'], main = 'W10')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W10']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W10']))

plot(openavail ~ ttd1, data = dat.avail[id=='W10'], main = 'W10')
abline(lm(openavail ~ ttd1, dat.avail[id=='W10']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W10']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W11'], main = 'W11') # low R2
abline(lm(wetavail ~ ttd1, dat.avail[id=='W11']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W11']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W11'], main = 'W11')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W11']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W11']))

plot(openavail ~ ttd1, data = dat.avail[id=='W11'], main = 'W11')
abline(lm(openavail ~ ttd1, dat.avail[id=='W11']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W11']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W05'], main = 'W05') # low R2
abline(lm(wetavail ~ ttd1, dat.avail[id=='W05']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W05']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W05'], main = 'W05')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W05']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W05']))

plot(openavail ~ ttd1, data = dat.avail[id=='W05'], main = 'W05')
abline(lm(openavail ~ ttd1, dat.avail[id=='W05']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W05']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W22'], main = 'W22') #R2 below .15
abline(lm(wetavail ~ ttd1, dat.avail[id=='W22']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W22']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W22'], main = 'W22')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W22']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W22']))

plot(openavail ~ ttd1, data = dat.avail[id=='W22'], main = 'W22')
abline(lm(openavail ~ ttd1, dat.avail[id=='W22']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W22']))

plot(wetavail ~ ttd1, data = dat.avail[id=='W04'], main = 'W04') ## R2 below .1
abline(lm(wetavail ~ ttd1, dat.avail[id=='W04']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W04']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W04'], main = 'W04')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W04']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W04']))

plot(openavail ~ ttd1, data = dat.avail[id=='W04'], main = 'W04')
abline(lm(openavail ~ ttd1, dat.avail[id=='W04']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W04']))

plot(wetavail ~ ttd1, data = dat.avail[id=='W02'], main = 'W02') ## R2 below .15
abline(lm(wetavail ~ ttd1, dat.avail[id=='W02']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W02']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W02'], main = 'W02')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W02']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W02']))

plot(openavail ~ ttd1, data = dat.avail[id=='W02'], main = 'W02')
abline(lm(openavail ~ ttd1, dat.avail[id=='W02']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W02']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W13'], main = 'W13') ## low R2 
abline(lm(wetavail ~ ttd1, dat.avail[id=='W13']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W13']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W13'], main = 'W13')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W13']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W13']))

plot(openavail ~ ttd1, data = dat.avail[id=='W13'], main = 'W13')
abline(lm(openavail ~ ttd1, dat.avail[id=='W13']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W13']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W12'], main = 'W12') ## R2 below .1
abline(lm(wetavail ~ ttd1, dat.avail[id=='W12']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W12']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W12'], main = 'W12')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W12']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W12']))

plot(openavail ~ ttd1, data = dat.avail[id=='W12'], main = 'W12')
abline(lm(openavail ~ ttd1, dat.avail[id=='W12']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W12']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W26'], main = 'W26') ## R2 below .1
abline(lm(wetavail ~ ttd1, dat.avail[id=='W26']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W26']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W26'], main = 'W26')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W26']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W26']))

plot(openavail ~ ttd1, data = dat.avail[id=='W26'], main = 'W26')
abline(lm(openavail ~ ttd1, dat.avail[id=='W26']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W26']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W27'], main = 'W27') ## R2 below .1, has more time in open
abline(lm(wetavail ~ ttd1, dat.avail[id=='W27']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W27']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W27'], main = 'W27')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W27']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W27']))

plot(openavail ~ ttd1, data = dat.avail[id=='W27'], main = 'W27')
abline(lm(openavail ~ ttd1, dat.avail[id=='W27']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W27']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W19'], main = 'W19') ## R2 below .1, slightly more open
abline(lm(wetavail ~ ttd1, dat.avail[id=='W19']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W19']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W19'], main = 'W19')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W19']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W19']))

plot(openavail ~ ttd1, data = dat.avail[id=='W19'], main = 'W19')
abline(lm(openavail ~ ttd1, dat.avail[id=='W19']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W19']))


plot(wetavail ~ ttd1, data = dat.avail[id=='W15'], main = 'W15') ## R2 below .1, more open
abline(lm(wetavail ~ ttd1, dat.avail[id=='W15']))
summary(lm(wetavail ~ ttd1, dat.avail[id=='W15']))

plot(closedavail ~ ttd1, data = dat.avail[id=='W15'], main = 'W15')
abline(lm(closedavail ~ ttd1, dat.avail[id=='W15']))
summary(lm(closedavail ~ ttd1, dat.avail[id=='W15']))

plot(openavail ~ ttd1, data = dat.avail[id=='W15'], main = 'W15')
abline(lm(openavail ~ ttd1, dat.avail[id=='W15']))
summary(lm(openavail ~ ttd1, dat.avail[id=='W15']))



dat.wet.avail <- dat[ua=='avail', .(ttd1, wetavail= sum(wet=='wet')/9, propwetavail= sum(propwet_end)/9), by=.(id, step_id_)]
dat.wet.probs.avail <- dat.wet.probs[ua=='avail', .(ttd1, wetavail= sum(wet=='wet')/9, propwetavail= sum(propwet_end)/9), by=.(id, step_id_)]
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail, groups = id)
abline(lm(wetavail ~ ttd1, dat.wet.probs.avail))
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W06'])

lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W09'],  type ='l')
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W10'],  type ='l')
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W11'],  type ='l')

lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W06'],  type ='l')
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W09'],  type ='l')
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W10'],  type ='l')
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W11'],  type ='l')

lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W04'],  type ='l') #ok
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W05'],  type ='l') #ok
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W02'],  type ='l') # days missing 30-40
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W13'],  type ='l') #missing periods
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W12'],  type ='l') # ok with a peek ~20 days
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W26'],  type ='l') # ok but generally low prop
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W22'],  type ='l') # ok
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W27'],  type ='l') # much less closer to death
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W19'],  type ='l') # ok
lattice::xyplot(propwetavail ~ ttd1, data = dat.wet.probs.avail[id=='W15'],  type ='l') # missing, low props

lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W04'],  type ='l') #ok
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W05'],  type ='l') #ok
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W02'],  type ='l') # days missing 30-40
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W13'],  type ='l') #missing periods
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W12'],  type ='l') # ok with a peak ~20 days, higer ttd ~lower avail prop
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W26'],  type ='l') # ok but generally low prop
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W22'],  type ='l') # ok
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W27'],  type ='l') # much less closer to death
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W19'],  type ='l') # ok, kind of lower props
lattice::xyplot(wetavail ~ ttd1, data = dat.wet.probs.avail[id=='W15'],  type ='l') # missing, low props


dat.wet.probs[ua == 'used',cor(propwet_end, ttd1), by=.(id)]