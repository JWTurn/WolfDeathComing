# remotes::install_github('bbolker/broom.mixed')
libs <- c('data.table', 'adehabitatHR', 'sp', 'rgdal', 'raster', 'dplyr')
lapply(libs, require, character.only = TRUE)

### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'
maps_out <- 'data/map-output/'

utm14N <- "+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
utm15N <- "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

### rarified data ----
dat.RMNP <- fread(paste0(raw, 'RMNPwolfpts_boundary.csv'))
dat.RMNP$datetime <- paste(dat.RMNP$gmtDate, dat.RMNP$gmtTime)
dat.RMNP$datetime <- as.POSIXct(dat.RMNP$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")


### clusters data ----
dat.clusters.RMNP <- fread(paste0(raw, 'InvestigatedPoints0131_cleaned20180614.csv'))
dat.clusters.RMNP[,'packID'] <- factor(dat.clusters.RMNP$Pack_ID, c('Baldy Lake', 'No Collar', 'Whitewater Lake', 'Gunn Lake',
                                                                    'Deep Lake', 'Ranck Creek', 'Lake Audy', 'Spruce Lake', 'Birdtail Valley', 'Block'),
                                       labels = c('BD', 'NA', 'WW', 'GL', 'DL', 'RC', 'AD', 'SL', 'BT', 'BL'))
dat.clusters.RMNP$Clstr_Strt <- as.POSIXct(dat.clusters.RMNP$Clstr_Strt, tz = 'UTC', "%Y-%m-%d")
dat.clusters.RMNP$Clstr_End <- as.POSIXct(dat.clusters.RMNP$Clstr_End, tz = 'UTC', "%Y-%m-%d")
dat.clusters.RMNP[,'year'] <- lubridate::year(dat.clusters.RMNP$Clstr_Strt)
dat.clusters.RMNP <- dat.clusters.RMNP[packID != 'NA']
dat.clusters.RMNP$packID <- factor(dat.clusters.RMNP$packID)

#### 2016-2017 ####
pack.RMNP <- SpatialPointsDataFrame(coords=as.data.frame(cbind(dat.clusters.RMNP$X_proj, dat.clusters.RMNP$Y_proj)), data=dat.clusters.RMNP, 
                                    proj4string = 
                                      CRS("+proj=utm +zone=14 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


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

## RMNP
makePackRange(dat.clusters.RMNP, 'AD')
makePackRange(dat.clusters.RMNP, 'BD')
makePackRange(dat.clusters.RMNP, 'BL')
makePackRange(dat.clusters.RMNP, 'BT')
makePackRange(dat.clusters.RMNP, 'DL')
makePackRange(dat.clusters.RMNP, 'GL')
makePackRange(dat.clusters.RMNP, 'SL')
makePackRange(dat.clusters.RMNP, 'WW')