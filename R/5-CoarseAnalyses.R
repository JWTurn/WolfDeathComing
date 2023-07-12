### GLM ====
# Julie Turner
# started July 27 2020


### Packages ----
# remotes::install_github('bbolker/broom.mixed')
# remotes::install_github('ropensci/spatsoc')
libs <- c('data.table', 'dplyr', 'tidyr', 'ggplot2','survival', 'glmmTMB', 'tibble', 'bbmle', 'patchwork', 
          'broom.mixed', 'performance', 'lme4', 'ggthemes')
lapply(libs, require, character.only = TRUE)

### Function ----
se <- function(x){
  sd(x, na.rm = T)/ sqrt(length(na.omit(x)))
}


### Input data ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### SSF data with all covariates ----
set.seed(53)
dat <- readRDS('data/derived-data/ssfAllCov_pub.Rds')

dat[,'ua'] <- ifelse(dat$case_ == T, 'used', 'avail')


## limit to wolves with other collared individuals in the pack at the same time 
dat[,uniqueN(step_id_), by=.(wolfID)]
dat.wnn <- dat[!is.na(distance2),uniqueN(step_id_), by=.(wolfID)]
dat.wnn.lastmo <- dat[ttd1<=31 & !is.na(distance2),uniqueN(step_id_), by=.(wolfID)]

# binary of if they were in their pack boundary or not
dat[,inout := ifelse(packYN_end== 'pack', 0, 1)]

DT <- dat[wolfID %in% dat.wnn.lastmo$wolfID & case_==TRUE]

# proportion of time they were outside their pack boundary
DT[,propout:=inout/.N, by = .(wolfID)]

DT.sub <- unique(DT[wolfID!='GHA26_W38',.( propout=sum(inout)/.N), by =.(wolfID,pop,COD)])
DT.sub[, mean(propout), COD]

#### ANALYSES ----
ssa <- glmmTMB(inout~COD + COD:log(ttd1+1) + (1|wolf_step_id), data = dat[wolfID %in% dat.wnn.lastmo$wolfID], family = 'poisson',
             map = list(theta=factor(c(NA))), start = list(theta=c(log(1000))))
tidy(ssa)


binom <- glmer(inout~COD + COD:log(ttd1+1) + (1|wolfID), data = DT, family = 'binomial')
tidy(binom)
performance(binom)
check_model(binom)

prop <- glm(propout~COD , data = DT.sub, family = 'binomial')
tidy(prop)
model_performance(prop)
check_model(prop)

#### GRAPHS ----
ggplot(DT.sub, aes(COD, propout, color= COD))+
  geom_boxplot(aes( fill = COD), alpha = 0.1, outlier.shape = NA)+
  geom_jitter() +
  theme_classic() +
  theme(text = element_text(size=15)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x =  element_text(size = 15)) + 
  theme(plot.margin = margin(0.1, 1, .1, .1, "cm")) +
  xlab('') + ylab('Proportion of steps outside of the pack boundary') +
  #theme(legend.position = "none") +
  scale_color_colorblind() + scale_fill_colorblind()



