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
day.RMNP <- readRDS('data/derived-data/sunsetsunriseRMNP_2016-2017.Rds')
day.GHA26 <- readRDS('data/derived-data/sunsetsunriseGHA26_2014-2018.Rds')


daymst.RMNP <- day.RMNP[,.(Date, TwiStartDate, SunriseDate, SunsetDate, TwiEndDate)]
daymst.GHA26 <- day.GHA26[,.(Date, TwiStartDate, SunriseDate, SunsetDate, TwiEndDate)]

# Read in SSF
ssf.RMNP <- readRDS("data/derived-data/ssfAll.Rds")
ssf.GHA26 <- readRDS("data/derived-data/ssfAll_GHA26.Rds")

# not needed for GHA26
ssf.RMNP[, 't1mst'] <- ssf.RMNP$t1_ - hours(6)
ssf.RMNP[, 't2mst'] <- ssf.RMNP$t2_ - hours(6)
ssf.RMNP[,'date1mst'] <- as.POSIXct(format(ssf.RMNP$t1mst, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")


ssf.GHA26[,'date1mst'] <- as.POSIXct(format(ssf.GHA26$t1_, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")



# Merge based on Date
fullmst.RMNP <- merge(ssf.RMNP, daymst.RMNP, by.x = 'date1mst', by.y = 'Date', all.x = T)

fullmst.GHA26 <- merge(ssf.GHA26, daymst.GHA26, by.x = 'date1mst', by.y = 'Date', all.x = T)


# Specify phase

fullmst.RMNP[,'ToD_start'] <- ifelse((fullmst.RMNP$SunriseDate + hours(1)) < fullmst.RMNP$t1mst & fullmst.RMNP$t1mst < (fullmst.RMNP$SunsetDate - hours(1)), 'day',
                                ifelse(fullmst.RMNP$TwiStart > fullmst.RMNP$t1mst | fullmst.RMNP$t1mst > fullmst.RMNP$TwiEnd, 'night', 'twilight'))

fullmst.GHA26[,'ToD_start'] <- ifelse((fullmst.GHA26$SunriseDate + hours(1)) < fullmst.GHA26$t1_ & fullmst.GHA26$t1_ < (fullmst.GHA26$SunsetDate - hours(1)), 'day',
                                     ifelse(fullmst.GHA26$TwiStart > fullmst.GHA26$t1_ | fullmst.GHA26$t1_ > fullmst.GHA26$TwiEnd, 'night', 'twilight'))


# get rid of now unneeded columns
ssf.wolf.RMNP <- fullmst.RMNP[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta, ToD_start,
                           parkYN_start, parkYN_end, roadDist_start, roadDist_end, parkDist_end, lnparkdist_start, lnparkdist_end, parkDistadj_end,
                       land_start, land_end, propwet_end, propopen_end, propconif_end, propmixed_end, wet_end, open_end, conif_end, mixed_end,
                       id, nn1, nn2, distance1, distance2, timegroup1, timegroup2, packYN_start, packYN_end, packDist_start, packDist_end, packDistadj_end,
                           ttd1, ttd2)]

ssf.wolf.GHA26 <- fullmst.GHA26[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta, ToD_start,
                       roadDist_start, roadDist_end, 
                       land_start, land_end, propwet_end, propopen_end, propconif_end, propmixed_end, wet_end, open_end, conif_end, mixed_end,
                       id, nn1, nn2, distance1, distance2, timegroup1, timegroup2, packYN_start, packYN_end, packDist_start, packDist_end, packDistadj_end,
                       ttd1, ttd2)]

saveRDS(ssf.wolf.RMNP, 'data/derived-data/ssfAllCov_RMNP.Rds')
saveRDS(ssf.wolf.GHA26, 'data/derived-data/ssfAllCov_GHA26.Rds')



