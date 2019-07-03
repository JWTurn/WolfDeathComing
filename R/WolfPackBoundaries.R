### Making pack boundaries ====
# Julie Turner
# Started: June 20 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
libs <- c('data.table', 'adehabitatHR', 'sp', 'rgdal', 'raster', 'dplyr')
lapply(libs, require, character.only = TRUE)

### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'
maps_out <- 'data/map-output/'

utm14N <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolfpts_boundary.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")


### clusters data ----
dat.clusters <- fread(paste0(raw, 'InvestigatedPoints0131_cleaned20180614.csv'))
dat.clusters[,'packID'] <- factor(dat.clusters$Pack_ID, c('Baldy Lake', 'No Collar', 'Whitewater Lake', 'Gunn Lake',
                           'Deep Lake', 'Ranck Creek', 'Lake Audy', 'Spruce Lake', 'Birdtail Valley', 'Block'),
                           labels = c('BD', 'NA', 'WW', 'GL', 'DL', 'RC', 'AD', 'SL', 'BT', 'BL'))
dat.clusters$Clstr_Strt <- as.POSIXct(dat.clusters$Clstr_Strt, tz = 'UTC', "%Y-%m-%d")
dat.clusters$Clstr_End <- as.POSIXct(dat.clusters$Clstr_End, tz = 'UTC', "%Y-%m-%d")
dat.clusters[,'year'] <- lubridate::year(dat.clusters$Clstr_Strt)
dat.clusters <- dat.clusters[packID != 'NA']
dat.clusters$packID <- factor(dat.clusters$packID)

#dat$PackID <- factor(dat$PackID)

#### 2016-2017 ####
pack <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.clusters$X_proj, dat.clusters$Y_proj)), data=dat.clusters, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

dat.clusters[,.(.N), packID]


plot(pack, col = c(1:8))



pack.kde<- kernelUD(pack[,'packID'], h="href", grid=30)
image(pack.kde)
packrange <- getverticeshr(pack.kde, percent = 95)
plot(packrange, col = 1:8)

writeOGR(packrange, dsn = 'data/map-output/RMNPkde', layer = 'RMNPkde', driver = "ESRI Shapefile")


pack.mcp <- mcp(pack[, "packID"], percent = 95, unin = "m", unout = "km2")
plot(pack.mcp, col = c(1:8))
paste0(raw, 'RMNPwolfpts_boundary.csv')
writeOGR(pack.mcp, dsn = paste0(maps_out, 'RMNPmcp'), layer = 'RMNPmcp', driver = "ESRI Shapefile")


#### by pack ####

# pack = 'BD'

makePackRange <- function(dat, pack){
  dat.sub<-dat[packID == pack]
  dat.sub$packID <- factor(dat.sub$packID)
  
  spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.sub$X_proj, dat.sub$Y_proj)), data=dat.sub, 
                                 proj4string = 
                                   CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
  kde <- kernelUD(spdf[,'packID'], h="href", grid=30)
  image(kde)
  packrange <- getverticeshr(kde, percent = 95)
  plot(packrange)
  
  mcp <- mcp(spdf[, "packID"], percent = 95, unin = "m", unout = "km2")
  plot(mcp)
  
  writeOGR(packrange, dsn = paste0(maps_out, paste(pack, 'kde', sep = '_')), layer = paste(pack, 'kde', sep = '_'), driver = "ESRI Shapefile")
  writeOGR(mcp, dsn = paste0(maps_out, paste(pack, 'mcp', sep = '_')), layer = paste(pack, 'mcp', sep = '_'), driver = "ESRI Shapefile")
}

makePackRange(dat.clusters, 'AD')
makePackRange(dat.clusters, 'BD')
makePackRange(dat.clusters, 'BL')
makePackRange(dat.clusters, 'BT')
makePackRange(dat.clusters, 'DL')
makePackRange(dat.clusters, 'GL')
makePackRange(dat.clusters, 'SL')
makePackRange(dat.clusters, 'WW')





#### 2016 ####
dat.2016 <- dat.clusters[year == 2016]
dat.2016$packID <- as.factor(dat.2016$packID)
dat.2016 <- dat.2016[packID != 'NA']



dat.2016[, c('EASTING', 'NORTHING') := as.data.table(project(cbind(X_proj, Y_proj), utm14N))] # not sure I needed to do this
# didn't like this

pack.2016 <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016$X_proj, dat.2016$Y_proj)), data=dat.2016, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
dat.2016[,.(.N), packID]

# loc.2016 <- data.frame("x"=dat.2016$X,"y"=dat.2016$Y)
# dat.pack.2016 <- dat.2016[,.( PackID)]
# proj4string <- sp::CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# packs.2016.other<- SpatialPointsDataFrame(loc.2016, dat.pack.2016, proj4string = proj4string)
# plot(packs.2016.other)

plot(pack.2016, col = c(1:4))

pack.kde.2016<- kernelUD(pack.2016[,'packID'], h="href", grid=30)
image(pack.kde.2016)
packrange.2016 <- getverticeshr(pack.kde.2016, percent = 95)
plot(packrange.2016, col = 1:4)

pack.mcp.2016 <- mcp(pack.2016[, "packID"], percent = 95, unin = "m", unout = "km2")
plot(pack.mcp.2016, col = c(1:4))

#### by pack ####

dat.2016.GL<-dat.2016[packID == 'GL']
dat.2016.GL$packID <- factor(dat.2016.GL$packID)
dat.2016.WW<-dat.2016[packID == 'WW']
dat.2016.WW$packID <- factor(dat.2016.WW$packID)
dat.2016.BD<-dat.2016[packID == 'BD']
dat.2016.BD$packID <- factor(dat.2016.BD$packID)
dat.2016.DL<-dat.2016[packID == 'DL']
dat.2016.DL$packID <- factor(dat.2016.DL$packID)

pack.2016.GL <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.GL$X, dat.2016.GL$Y)), data=dat.2016.GL, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.GL<- kernelUD(pack.2016.GL[,'packID'], h="href", grid=30)
packrange.2016.GL <- getverticeshr(pack.kde.2016.GL, percent = 95)
image(pack.kde.2016.GL)
plot(packrange.2016.GL)

pack.2016.WW <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.WW$X, dat.2016.WW$Y)), data=dat.2016.WW, 
                                       proj4string = 
                                         CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.WW<- kernelUD(pack.2016.WW[,'packID'], h="href", grid=30)
packrange.2016.WW <- getverticeshr(pack.kde.2016.WW, percent = 95)
image(pack.kde.2016.WW)
plot(packrange.2016.WW)


pack.2016.BD <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.BD$X, dat.2016.BD$Y)), data=dat.2016.BD, 
                                       proj4string = 
                                         CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.BD<- kernelUD(pack.2016.BD[,'packID'], h="href", grid=30)
packrange.2016.BD <- getverticeshr(pack.kde.2016.BD, percent = 95)
image(pack.kde.2016.BD)
plot(packrange.2016.BD)


pack.2016.DL <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2016.DL$X, dat.2016.DL$Y)), data=dat.2016.DL, 
                                       proj4string = 
                                         CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pack.kde.2016.DL<- kernelUD(pack.2016.DL[,'packID'], h="href", grid=30)
packrange.2016.DL <- getverticeshr(pack.kde.2016.DL, percent = 95)
image(pack.kde.2016.DL)
plot(packrange.2016.DL)


#### 2017 ####
dat.2017 <- dat.clusters[year == 2017]
dat.2017$packID <- factor(dat.2017$packID)
dat.2017 <- dat.2017[packID != 'NA']
dat.2017[,.(.N), packID]
#dat.2017[, c('EASTING', 'NORTHING') := as.data.table(project(cbind(Longitude, Latitude), utm14N))] # not sure I needed to do this

pack.2017 <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.2017$X_proj, dat.2017$Y_proj)), data=dat.2017, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

plot(pack.2017)

pack.kde.2017<- kernelUD(pack.2017[,'packID'], h="href", grid=30)
image(pack.kde.2017)
packrange.2017 <- getverticeshr(pack.kde.2017, percent = 95)
plot(packrange.2017, col = 1:6)

pack.mcp.2017 <- mcp(pack.2017[, "PackID"], percent = 95, unin = "m", unout = "km2")
plot(pack.mcp.2017, col = c(1:6))


#### by pack ####

dat.2017.BD<-dat.2017[PackID == 'BD']
#dat.2017.RC<-dat.2017[PackID == 'RC']
dat.2017.AD<-dat.2017[PackID == 'AD']
dat.2017.BL<-dat.2017[PackID == 'BL']
dat.2017.BT<-dat.2017[PackID == 'BT']
dat.2017.SL<-dat.2017[PackID == 'SL']

