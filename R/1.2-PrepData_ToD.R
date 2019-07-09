### Data - Time of Day ====
# Julie Turner
# Started: July 3 2019




### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'lubridate', 'tidyr')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'



# Read in Sunset/Sunrise times
day <- readDS('data/derived-data/sunsetsunriseRMNP_2016-2017.Rds')

daymst <- day[,.(Date, TwiStartDate, SunriseDate, SunsetDate, TwiEndDate)]

# Read in SSF
ssf.soc <- readRDS("data/derived-data/ssfAll.Rds")

ssf.soc[, 't1mst'] <- ssf.soc$t1_ - hours(6)
ssf.soc[, 't2mst'] <- ssf.soc$t2_ - hours(6)
ssf.soc[,'date1mst'] <- as.POSIXct(format(ssf.soc$t1mst, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")



# Merge based on Date
fullmst <- merge(ssf.soc, daymst, by.x = 'date1mst', by.y = 'Date', all.x = T)



# Specify phase

fullmst[,'ToD_start'] <- ifelse((fullmst$SunriseDate + hours(1)) < fullmst$t1mst & fullmst$t1mst < (fullmst$SunsetDate - hours(1)), 'day',
                                ifelse(fullmst$TwiStart > fullmst$t1mst | fullmst$t1mst > fullmst$TwiEnd, 'night', 'twilight'))

# get rid of now unneeded columns
ssf.wolf <- fullmst[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta, ToD_start, tod_amt_start = tod_start_, 
                           parkYN_start, parkYN_end, roadDist_start, roadDist_end, lnparkdist_start, lnparkdist_end, land_start, land_end, 
                           id, nn1, nn2, distance1, distance2, timegroup1, timegroup2, packYN_start, packYN_end, packDist_start, packDist_end, 
                           ttd1, ttd2)]

saveRDS(ssf.wolf, 'data/derived-data/ssfAllCov.Rds')
