CreateTrack <- function(x.col, y.col, date.col, crs, ID, NumbRandSteps, sl_distr, ta_distr) {
  #print(ID)
  #create track from dataset
  trk <- track(x.col, y.col, date.col, ID, crs) %>%
    #function turns locs into steps
    steps()
  #remove any steps that span more than 2hr
  trk$dt_ <- difftime(trk$t2_, trk$t1_, unit='hours')
  trk <- subset(trk, trk$dt_ > 1.9 & trk$dt_ < 2.1, drop = T)
  #generate random steps
  trk %>%
    random_steps(n = NumbRandSteps, sl_distr, ta_distr) %>%
    sl_params()
}

#run function by wolf ID
WParams <- W2018[, CreateTrack(10, x.col = EASTING, y.col = NORTHING, date.col = ts, ID = Collar_ID, crs = utm14N,
                               sl_distr = "gamma", ta_distr = "vonmises"),
                 by = Collar_ID]

W_SLparams <- WParams[!is.na(V1)]