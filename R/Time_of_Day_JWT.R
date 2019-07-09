### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'lubridate', 'tidyr')
lapply(libs, require, character.only = TRUE)




### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'



# Read in Sunset/Sunrise times
ss2016 <- fread(paste0(raw,"sunset_sunrise2016.csv"), header=TRUE)
ss2016[,'year'] <- 2016
ss2017 <- fread(paste0(raw,"sunset_sunrise2017.csv"),header=TRUE)
ss2017[,'year'] <- 2017
ss2018 <- fread(paste0(raw,"sunset_sunrise2018.csv"),header=TRUE)
ss2018[,'year'] <- 2018
day <- rbind(ss2016, ss2017, ss2018)
day[,'dateyear'] <- paste(day$year, day$Date, sep = '-')
day <- day[, .(Date = dateyear, TwiStart = `Civil Twilight Start`, Sunrise = Sunrise, Sunset = Sunset, TwiEnd = `Civil Twilight End`)]
day[,'Date'] <- as.POSIXct(day$Date, tz = 'UTC', "%Y-%m-%d")
day[,'TwiStartDate'] <- as.POSIXct(paste(day$Date, day$TwiStart, sep = ' '), tz = 'UTC', "%Y-%m-%d %H:%M")
day[,'SunriseDate'] <- as.POSIXct(paste(day$Date, day$Sunrise, sep = ' '), tz = 'UTC', "%Y-%m-%d %H:%M")
day[,'SunsetDate'] <- as.POSIXct(paste(day$Date, day$Sunset, sep = ' '), tz = 'UTC', "%Y-%m-%d %H:%M")
day[,'TwiEndDate'] <- as.POSIXct(paste(day$Date, day$TwiEnd, sep = ' '), tz = 'UTC', "%Y-%m-%d %H:%M")


day[,'gmtTwiStartDate'] <- day$TwiStartDate + hours(6)
day[,'gmtSunriseDate'] <- day$SunriseDate + hours(6)
day[,'gmtSunsetDate'] <- day$SunsetDate + hours(6)
day[,'gmtTwiEndDate'] <- day$TwiEndDate + hours(6)
day[,'gmtDate'] <-  as.POSIXct(format(day$gmtTwiStartDate, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")

# saveRDS(day, 'data/derived-data/sunsetsunriseRMNP_2016-2017.Rds')



daygmt <- day[,.(gmtTwiStartDate, gmtSunriseDate, gmtSunsetDate, gmtTwiEndDate)]
daygmt[,'date'] <-  as.POSIXct(format(daygmt$gmtTwiStartDate, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")

daymst <- day[,.(Date, TwiStartDate, SunriseDate, SunsetDate, TwiEndDate)]

# Read in SSF
ssf.soc <- readRDS("data/derived-data/ssfAll.Rds")
ssf.soc[,'date1'] <- as.POSIXct(format(ssf.soc$t1_, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")
# ssf.soc[,'time1'] <- as.ITime(ssf.soc$t1_)
ssf.soc[,'date2'] <- as.POSIXct(format(ssf.soc$t2_, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")
# ssf.soc[,'time2'] <- as.ITime(ssf.soc$t2_)

ssf.soc[, 't1mst'] <- ssf.soc$t1_ - hours(6)
ssf.soc[, 't2mst'] <- ssf.soc$t2_ - hours(6)
ssf.soc[,'date1mst'] <- as.POSIXct(format(ssf.soc$t1mst, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")


locs <- ssf.soc[,.(id, t1_, t2_, date1, date2, time1, time2)]
locs <- unique(locs)

locsmst <- ssf.soc[,.(id, t1mst,t2mst)]
locsmst <- unique(locsmst)
locsmst[,'date'] <- as.POSIXct(format(locsmst$t1mst, "%Y-%m-%d"), tz = 'UTC', "%Y-%m-%d")

# Merge based on Date
full <- merge(locs, daygmt, by.x = 'date1', by.y = 'gmtDate')
fullgmt <- merge(ssf.soc, daygmt, by.x = 'date1', by.y = 'gmtDate', all.x = T)


fullmst <- merge(ssf.soc, daymst, by.x = 'date1mst', by.y = 'Date', all.x = T)



# Specify phase
# fullmst[, c('Day', 'Night', 'Twilight') := .(
#   ifelse((Sunrise + 3600) < time1 & time1 < (Sunset - 3600), 1, 0),
#   ifelse(TwiStart > time1 | time1 > TwiEnd, 1, 0),
#   ifelse(TwiStart < time1 & time1 < (Sunrise + 3600) | time1 > (Sunset - 3600) & time1 < TwiEnd, 1, 0))]
# head(full)
# 
# ssf.soc[,'ToD_start'] <- ifelse((daymst$SunriseDate + hours(1)) < ssf.soc$t1mst & ssf.soc$t1mst < (daymst$SunsetDate - hours(1)), 'day',
#                                 ifelse(daymst$TwiStart > ssf.soc$t1mst | ssf.soc$t1mst > daymst$TwiEnd, 'night', 'twilight'))
#   

#### This one works
fullmst[,'ToD_start'] <- ifelse((fullmst$SunriseDate + hours(1)) < fullmst$t1mst & fullmst$t1mst < (fullmst$SunsetDate - hours(1)), 'day',
                                ifelse(fullmst$TwiStart > fullmst$t1mst | fullmst$t1mst > fullmst$TwiEnd, 'night', 'twilight'))

ssf.wolf <- fullmst[,.(burst_, step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, dt_, sl_, log_sl, ta_, cos_ta, ToD_start, tod_amt_start = tod_start_, 
                           parkYN_start, parkYN_end, roadDist_start, roadDist_end, lnparkdist_start, lnparkdist_end, land_start, land_end, 
                           id, nn1, nn2, distance1, distance2, timegroup1, timegroup2, packYN_start, packYN_end, packDist_start, packDist_end, 
                           ttd1, ttd2)]

saveRDS(ssf.wolf, 'data/derived-data/ssfAllCov.Rds')
# full.ToD <- fullmst[,.(id, t1mst, ToD_start)]
# 
# ssf.soc.ToD <- merge(ssf.soc, full.ToD, by = c('id', 't1mst'), all.x = T, allow.cartesian = T)

# full[,'ToD_start'] <- ifelse(between(full$time1, (daygmt$gmtTwiEndDate), (daygmt$gmtTwiStartDate)), 'night', 
#                                      ifelse(between(full$time1, (daygmt$gmtSunriseDate + hours(1)), (daygmt$gmtSunsetDate - hours(1))), 'day', 'twilight'))
# 
# 
# ssf.soc[,'ToD_start'] <- ifelse(ssf.soc$t1_ %within% interval(daygmt$gmtTwiEndDate, daygmt$gmtTwiStartDate), 'night', 
#                              ifelse(ssf.soc$t1_ %within% 
#                               interval((daygmt$gmtSunriseDate), (daygmt$gmtSunsetDate)), 'day', 'twilight'))
# 
# ssf.soc[,'ToD_start'] <- ifelse(ssf.soc$t1mst %within% interval(daymst$TwiEndDate, daymst$TwiStartDate), 'night', 
#                                 ifelse(ssf.soc$t1mst %within% 
#                                          interval((daymst$SunriseDate), (daymst$SunsetDate)), 'day', 'twilight'))
# 


# ssf.soc[,'ToD_start'] <- ifelse(ssf.soc$t1_ %between% interval(daygmt$gmtTwiEndDate, daygmt$gmtTwiStartDate), 'night', 
#                                 ifelse(ssf.soc$t1_ %between% 
#                                          interval((daygmt$gmtSunriseDate), (daygmt$gmtSunsetDate)), 'day', 'twilight'))
# 
# 
# locs[,'ToD_start'] <- ifelse(locs$t1_ %within% (daygmt$gmtTwiEndDate %--% daygmt$gmtTwiStartDate), 'night', 
#                                 ifelse(locs$t1_ %within% 
#                                          ((daygmt$gmtSunriseDate) %--% (daygmt$gmtSunsetDate)), 'day', 'twilight'))
# 



# fwrite(full, '5138_ToD.csv')

