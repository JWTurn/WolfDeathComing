############################
##### 1 - Load packages ----

libs <- c('lubridate', 'amt', 'data.table', 'tidyr', 'raster')
lapply(libs, require, character.only=TRUE)

###################################
### 2 - Generate Grid function ----

GenerateGrid <- function(spacing, crs, boundary, DT = NULL, xCol = NULL, yCol = NULL){
  
  mcpExtent <- as(extent(boundary), 'SpatialPolygons')
  extentDT <- data.table(mcpExtent@bbox, keep.rownames = TRUE)
  
  extentDT[, dif := abs(min - max)]
  
  
  if(extentDT[order(-dif)][1, rn] == "x"){
    xSeq <- extentDT[rn == "x", seq(min, max, by = spacing)]
    ySeq <- extentDT[rn == "y", seq(min, min + (spacing * (length(xSeq) - 1)),
                                    by = spacing)]
  } else {
    ySeq <- extentDT[rn == "y", seq(min, max, by = spacing)]
    xSeq <- extentDT[rn == "x", seq(min, min + (spacing * (length(ySeq) - 1)),
                                    by = spacing)]
  }
  
  r <- raster::extent(head(xSeq, n = 1), tail(xSeq, n = 1),
                      head(ySeq, n = 1), tail(ySeq, n = 1))
  ra <- raster(r)
  res(ra) <- c(spacing, spacing)
  
  projection(ra) <- CRS(crs)
  
  rSP <- rasterToPoints(ra, spatial = TRUE)
  
  data.table(rSP@coords)
  
}



###############################
#### 3 - Load spatial data ----

# A) land cover rasters
ag_cov <- raster('spatial/agriculture_raster.tif')
names(ag_cov) <- 'agricultural'
grass_cov <- raster('spatial/grassland_raster.tif')
names(grass_cov) <- 'grassland'
conif_cov <- raster('spatial/coniferous_raster.tif')
names(conif_cov) <- 'coniferous'
mw_cov <- raster('spatial/mixedwood_raster.tif')
names(mw_cov) <- 'mixedwood'
water_cov <- raster('spatial/water_raster.tif')
names(water_cov) <- 'water'
# combine landcover into raster brick
habitat_cov <- brick(ag_cov, grass_cov, conif_cov, mw_cov, water_cov)
# C) density raster
yrs <- c('2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2015', '2016')
for(i in yrs){
  j <- raster(paste0('spatial/', paste('elk_density_', '.tif', sep=i)))
  assign(paste0('elk_density_', i, sep=''), j)
}

################################
#### 2 - Load location data ----

dat <- readRDS('input/elkdatbyEYs.rds') 

#########################################
#### 3 - Clean NAs and add ID column ----

dat <- dat[, elkyear := .GRP, 
           by = c('EarTag', 'intyear')] %>% filter(!is.na('X')) %>% 
  amt::select(x='X', y='Y', t = 'DateTime', id = 'elkyear')

# make datetime column POSIX
dat$t <- as.POSIXct(dat$t, tz = 'America/Winnipeg')

#### SKIP 4, THIS IS FOR LN'S DENSITY STUFF ####

#####################################################
#####################################################
#####                                           #####
##### 4 - Extract mean/sd density values        #####
##### from density rasters according to habitat #####
#####                                           #####
#####################################################
#####################################################

# Match projection of habitat raster to density
habitat_cov <- raster::projectRaster(habitat_cov, crs = crs(elk_density_2003))

# Create extract points grid to sample at resolution of habitat raster
extract_pts <- GenerateGrid(spacing = 30, boundary = elk_density_2003, crs = '+proj=utm +zone=14 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

# Extract values from habitat rasters
hab_vals <- as.data.table(extract(habitat_cov, extract_pts))

# Extract values from density rasters
densitylist <- list(elk_density_2003, elk_density_2004, 
                    elk_density_2005, elk_density_2006, elk_density_2007, elk_density_2008, 
                    elk_density_2009, elk_density_2010, elk_density_2011, elk_density_2012, elk_density_2015, elk_density_2016)
densitytable <- data.table()
for(drast in densitylist){
  dtab <- as.data.table(extract(drast, extract_pts))
  colnames(dtab) <- drast@data@names
  densitytable <- cbind(densitytable, dtab)
}

# Combine density and habitat data
densitytable <- cbind(hab_vals, densitytable)

# Pull rows where agriculture = 1
density_ag <- subset(densitytable, agriculture==1)
# Pull rows where conifer = 1
density_conif <- subset(densitytable, coniferous==1)
# Pull rows where mixedwood = 1
density_mw <- subset(densitytable, mixedwood==1)

# Mean and SD density of agriculture cells
mean_ag_density <- mean(as.numeric(density_ag[, lapply(.SD, mean, na.rm=TRUE), .SDcols=c(3:14)]))
sd_ag_density <- mean(as.numeric(density_ag[, lapply(.SD, sd, na.rm=TRUE), .SDcols=c(3:14)]))

# Mean and SD of density conifer cells
mean_conif_density <- mean(as.numeric(density_conif[, lapply(.SD, mean, na.rm=TRUE), .SDcols=c(3:14)]))
sd_conif_density <- mean(as.numeric(density_conif[, lapply(.SD, sd, na.rm=TRUE), .SDcols=c(3:14)]))

# Mean and SD density of mixedwood cells
mean_mw_density <- mean(as.numeric(density_mw[, lapply(.SD, mean, na.rm=TRUE), .SDcols=c(3:14)]))
sd_mw_density <- mean(as.numeric(density_mw[, lapply(.SD, sd, na.rm=TRUE), .SDcols=c(3:14)]))

# Melt data tables to plot distribution of density values in each habitat
agplot <- melt(density_ag, id.vars=c('coniferous', 'agriculture', 'mixedwood'))
mwplot <- melt(density_mw, id.vars=c('coniferous', 'agriculture', 'mixedwood'))
conifplot <- melt(density_conif, id.vars=c('coniferous', 'agriculture', 'mixedwood'))

# Plot distributions
ggplot(agplot, aes(x=value, group=variable)) + geom_density() + geom_vline(xintercept = (mean_ag_density+(2*sd_ag_density)), linetype='dashed')
ggplot(mwplot, aes(x=value, group=variable)) + geom_density() + geom_vline(xintercept = (mean_mw_density+(2*sd_mw_density)), linetype='dashed')
ggplot(conifplot, aes(x=value, group=variable)) + geom_density() + geom_vline(xintercept = (mean_conif_density+(2*sd_conif_density)), linetype='dashed')

############################################################
############################################################
#####                                                  #####
##### 5 - Loop through to create table of proportion   #####
##### from density rasters according to habitat        #####
#####                                                  #####
############################################################
############################################################

##########################################################
#### 4 - Loop through to create table of proportion 
# of each habitat in available steps for each elk ----

mean_by_habs <- data.table()

for(i in unique(dat$id)){
  tryCatch({
    print(i)
    dat_1 <- dat %>% filter(id==i)
    
    # A) Make track with x, y, time
    dat_1 <- mk_track(dat_1, .x=x, .y=y, .t=t, crs = sp::CRS("+init=epsg:26914"))
 ###SKIP THIS ###   
    # # Calculate time interval between fix rates
    # start <- ymd_hms(dat_1$t_[1])
    # end <- ymd_hms(dat_1$t_[2])
    # time.interval <- start %--% end
    # time.duration <- as.duration(time.interval)
    # dur <- as.numeric(time.duration, "hours")
    # 
    # # If fix is on the half hour (i.e. not divisible by 1),
    # # optionally round every second date time (this will make the fixes every 2 
    # # hours for bursts of three steps:
    # 
    # # for(j in dat_1[seq(1,nrow(dat_1),2) ,]){
    # # if(dur > 2.4 & dur < 2.6){
    # #   dat_1$t_ <- ceiling_date(dat_1$t_, unit = "hour", change_on_boundary = NULL,
    # #                week_start = getOption("lubridate.week.start", 7))
    # # }
    # # }
    # 
    # if(dur > 1.7 & dur < 2.6){
    #   dur <- 2
    # }
    # 
    # if(!dur %% 1 == 0){
    #   dur <- ceiling(dur)}
    # 
    # if(dur > 2.6){
    #   dur <- dur + 0.82
    # }
    # 
    # # B) Specify which density raster to use
    # # Each year has a different raster including elk counts 
    # # from a 3-year window surrounding the current year
    # 
    # # choose raster:
    # if(all(dat_1$t_ %within% interval(ymd("2003-10-31"), ymd("2004-04-02")))){
    #   elk_density <- elk_density_2004
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2004-10-31"), ymd("2004-05-02")))){
    #   elk_density <- elk_density_2005
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2005-10-31"), ymd("2006-04-02")))){
    #   elk_density <- elk_density_2006
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2006-10-31"), ymd("2007-04-02")))){
    #   elk_density <- elk_density_2007
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2007-10-31"), ymd("2008-04-02")))){
    #   elk_density <- elk_density_2008
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2008-10-31"), ymd("2009-04-02")))){
    #   elk_density <- elk_density_2009
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2009-10-31"), ymd("2010-04-02")))){
    #   elk_density <- elk_density_2010
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2010-10-31"), ymd("2011-04-02")))){
    #   elk_density <- elk_density_2011
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2011-10-31"), ymd("2012-04-02")))){
    #   elk_density <- elk_density_2012
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2014-10-31"), ymd("2015-04-02")))){
    #   elk_density <- elk_density_2015
    # }
    # if(all(dat_1$t_ %within% interval(ymd("2015-10-31"), ymd("2016-04-02")))){
    #   elk_density <- elk_density_2016
    # }
    # 
    # # rename raster values to common value
    # names(elk_density) <- 'elk_density'
    # 
    # C) Resample track to constant fix rate with a tolerance around the fix rate
    # to accommodate differences in timing of fix rate:
    # fix rate of 2 hours with a tolerance of 30 MINUTES
    # to accommodate half hour fix rates
    stps <- track_resample(dat_1, rate=hours(2), tolerance=minutes(30)) %>%
      
      # D) Filter to bursts with at least 3 consecutive points (required to calculate TA)
      filter_min_n_burst(min_n=3) %>% steps_by_burst() %>%
      
      # E) calculates whether the location was taken at day or night
      time_of_day(include.crepuscule=FALSE, where='both')
    
    hab_by_step <- stps %>% random_steps(n=9, sl_distr = "gamma", ta_distr = "vonmises") %>%
      extract_covariates(habitat_cov, where='end') %>%
      extract_covariates(elk_density, where='both') 
    
    # subset out available steps only
    hab_by_step <- hab_by_step[hab_by_step$case_ == FALSE ,]
    
    # divide number of available steps in each habitat by total available steps
    avg_covs <- cbind(sum(hab_by_step$agricultural, na.rm=TRUE)/nrow(hab_by_step), sum(hab_by_step$grassland, na.rm=TRUE)/nrow(hab_by_step),
                      sum(hab_by_step$coniferous, na.rm=TRUE)/nrow(hab_by_step), sum(hab_by_step$mixedwood, na.rm=TRUE)/nrow(hab_by_step), 
                      sum(hab_by_step$water, na.rm=TRUE)/nrow(hab_by_step), sum(hab_by_step$elk_density_end, na.rm=TRUE), 
                      sum(hab_by_step$elk_density_start, na.rm=TRUE))
  
    # bind elk name
    avg_covs <- cbind(i, avg_covs)
    
    # bind to data table
    mean_by_habs <- rbind(mean_by_habs, avg_covs)

}, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})}

colnames(mean_by_habs) <- c('elkID', colnames(hab_by_step[,15:21]))

# Subset elk years with at least 2% agricultural end steps and with end steps in density
agany <- mean_by_habs[agricultural>0.02 & elk_density_start>0]
agany <- as.character(unique(agany$elkID))

# Subset elk years with at least 2% conifer end steps and with end steps in density
conifany <- mean_by_habs[coniferous>0.02 & elk_density_start>0]
conifany <- as.character(unique(conifany$elkID))

# Subset elk years with end steps in density
densany <- mean_by_habs[elk_density_start>0]
densany <- as.character(unique(densany$elkID))

#####################
#### 5 - Save outputs

# Save vectors of elk IDs to go in each model
# Elk with agriculture
saveRDS(agany, 'input/elkIDs_ag_any.rds')
# Elk with coniferous
saveRDS(conifany, 'input/elkIDs_conif_any.rds')    
# Elk with density
saveRDS(densany, 'input/elkIDs_density_any.rds') 

# Save means and SDs
# Agriculture
saveRDS(mean_ag_density, 'input/mean_density_in_ag.rds')
saveRDS(sd_ag_density, 'input/sd_density_in_ag.rds')
# Coniferous
saveRDS(mean_conif_density, 'input/mean_density_in_conif.rds')
saveRDS(sd_conif_density, 'input/sd_density_in_conif.rds')
# Mixedwood
saveRDS(mean_mw_density, 'input/mean_density_in_mw.rds')
saveRDS(sd_mw_density, 'input/sd_density_in_mw.rds')


    