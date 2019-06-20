### Data ====
# Julie Turner
# Started: June 20 2019


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
libs <- c('data.table', 'spatsoc')
lapply(libs, require, character.only = TRUE)

### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### rarified data ----
dat <- fread(paste0(raw, 'RMNPwolf_rarified.csv'))
dat$datetime <- paste(dat$gmtDate, dat$gmtTime)
dat$datetime <- as.POSIXct(dat$datetime, tz = 'UTC', "%Y-%m-%d %H:%M:%S")
group_times(dat, datetime = 'datetime', threshold = '20 minutes')

### create steps first, then do nn, but only compare to original nn for each random step, 
### not whoever may be the new closest  

dat_nn <- edge_nn(dat, id = 'WolfID', coords = c('X', 'Y'),
                 timegroup = 'timegroup', splitBy = 'PackID')
