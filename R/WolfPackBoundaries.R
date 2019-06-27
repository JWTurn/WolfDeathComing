### Making pack boundaries ====
# Julie Turner
# Started: June 20 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
libs <- c('data.table', 'adehabitatHR', 'sp', 'rgdal', 'raster')
lapply(libs, require, character.only = TRUE)

### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

utm14N <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolfpts_boundary.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")

#dat$PackID <- factor(dat$PackID)

#### 2016 ####
dat.2016 <- dat[Year == 2016]
dat.2016$PackID <- as.factor(dat.2016$PackID)



dat.2016[, c('EASTING', 'NORTHING') := as.data.table(project(cbind(Longitude, Latitude), utm14N))] # not sure I needed to do this

pack.2016 <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016$X, dat.2016$Y)), data=dat.2016, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
dat.2016[,.(.N), PackID]

# loc.2016 <- data.frame("x"=dat.2016$X,"y"=dat.2016$Y)
# dat.pack.2016 <- dat.2016[,.( PackID)]
# proj4string <- sp::CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# packs.2016.other<- SpatialPointsDataFrame(loc.2016, dat.pack.2016, proj4string = proj4string)
# plot(packs.2016.other)

plot(pack.2016)

pack.kde.2016<- kernelUD(pack.2016[,'PackID'], h="href", grid=30)
image(pack.kde.2016)
packrange.2016 <- getverticeshr(pack.kde.2016, percent = 95)
plot(packrange.2016, col = 1:4)

pack.mcp.2016 <- mcp(pack.2016[, "PackID"], percent = 95, unin = "m", unout = "km2")
plot(pack.mcp.2016, col = c(1:4))

#### by pack ####

dat.2016.GL<-dat.2016[PackID == 'GL']
dat.2016.GL$PackID <- factor(dat.2016.GL$PackID)
dat.2016.WW<-dat.2016[PackID == 'WW']
dat.2016.WW$PackID <- factor(dat.2016.WW$PackID)
dat.2016.BD<-dat.2016[PackID == 'BD']
dat.2016.BD$PackID <- factor(dat.2016.BD$PackID)
dat.2016.RC<-dat.2016[PackID == 'RC']
dat.2016.RC$PackID <- factor(dat.2016.RC$PackID)

pack.2016.GL <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.GL$X, dat.2016.GL$Y)), data=dat.2016.GL, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.GL<- kernelUD(pack.2016.GL[,'PackID'], h="href", grid=30)
packrange.2016.GL <- getverticeshr(pack.kde.2016.GL, percent = 95)


pack.2016.WW <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.WW$X, dat.2016.WW$Y)), data=dat.2016.WW, 
                                       proj4string = 
                                         CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.WW<- kernelUD(pack.2016.WW[,'PackID'], h="href", grid=30)
packrange.2016.WW <- getverticeshr(pack.kde.2016.WW, percent = 95)


pack.2016.BD <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.BD$X, dat.2016.BD$Y)), data=dat.2016.BD, 
                                       proj4string = 
                                         CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.BD<- kernelUD(pack.2016.BD[,'PackID'], h="href", grid=30)
packrange.2016.BD <- getverticeshr(pack.kde.2016.BD, percent = 95)


pack.2016.RC <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.RC$X, dat.2016.RC$Y)), data=dat.2016.RC, 
                                       proj4string = 
                                         CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.RC<- kernelUD(pack.2016.RC[,'PackID'], h="href", grid=30)
packrange.2016.RC <- getverticeshr(pack.kde.2016.RC, percent = 95)


#### 2017 ####
dat.2017 <- dat[Year == 2017]
dat.2017$PackID <- factor(dat.2017$PackID)
dat.2017[,uniqueN(fixID), PackID]
dat.2017[, c('EASTING', 'NORTHING') := as.data.table(project(cbind(Longitude, Latitude), utm14N))] # not sure I needed to do this

pack.2017 <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2017$EASTING, dat.2017$NORTHING)), data=dat.2017, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

plot(pack.2017)

pack.kde.2017<- kernelUD(pack.2017[,'PackID'], h="href", grid=30)
image(pack.kde.2017)
packrange.2017 <- getverticeshr(pack.kde.2017, percent = 95)
plot(packrange.2017, col = 1:6)

pack.mcp.2017 <- mcp(pack.2017[, "PackID"], percent = 95, unin = "m", unout = "km2")
plot(pack.mcp.2017, col = c(1:6))


#### by pack ####

dat.2017.BD<-dat.2017[PackID == 'BD']
dat.2017.RC<-dat.2017[PackID == 'RC']
dat.2017.AD<-dat.2017[PackID == 'AD']
dat.2017.BL<-dat.2017[PackID == 'BL']
dat.2017.BT<-dat.2017[PackID == 'BT']
dat.2017.SL<-dat.2017[PackID == 'SL']

