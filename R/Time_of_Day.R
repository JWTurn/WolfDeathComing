library(data.table)
setwd("~/Documents/Masters/Data/Time of Day")

# Read in Sunset/Sunrise times
day <- fread("Sunset-Sunrise_2016.csv", header=TRUE, skip=3)
day <- day[, .(Date = V2, TwiStart = V4, Sunrise = V5, Sunset = V7, TwiEnd = V8)]

# Read in Locs
locs <- fread("5138_Locs.csv")

# Convert to date
timeFields <- c('TwiStart', 'Sunrise',
                'Sunset', 'TwiEnd')
day[, (timeFields) := lapply(.SD, as.ITime), .SDcols = timeFields]

locs <- locs[, Time := as.ITime(Time)]

# Merge based on Date
full <- merge(day, locs)

# Specify phase
full[, c('Day', 'Night', 'Twilight') := .(
  ifelse((Sunrise + 3600) < Time & Time < (Sunset - 3600), 1, 0),
  ifelse(TwiStart > Time | Time > TwiEnd, 1, 0),
  ifelse(TwiStart < Time & Time < (Sunrise + 3600) | Time > (Sunset - 3600) & Time < TwiEnd, 1, 0))]
head(full)

fwrite(full, '5138_ToD.csv')

