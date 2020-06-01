# RSS code for going through individuals
# JWT 

# function for predicting the road h1 part of the RSS for individuals
# this specifies the specific type of death and individuals, t2death is so that I can categorize it as 60 days and 1 day, but you may not need to do this type of thing

p.road.h1.indiv <- function(ids, DT, mod, death, t2death){
  lapply(ids, function(i) {
    #unique(
    DT[#wolfID == i,
      ,.(h1 = predict(
        mod,
        newdata = .SD[, .(
          ToD_start = factor('day', levels = levels(ToD_start)),
          log_sl = mean(log_sl),
          cos_ta = mean(cos_ta),
          # land_end_adj = factor('forest', levels = levels(land_end_adj)),
          propforest_end_adj = mean(propforest_end_adj, na.rm = T),
          propopen_end_adj = mean(propopen_end_adj, na.rm = T),
          propwet_end= mean(propwet_end, na.rm = T), 
          roadDist_end = seq(from = 0, to = 3000, length.out = 100),
          distance2 = median(distance2, na.rm = T),
          nnDist_end = median(nnDist_end, na.rm = T),
          packDist_end = median(packDist_end, na.rm = T),
          COD = factor(death, levels = levels(COD)),
          ttd1 = t2death,
          wolf_step_id = NA, # this is for the mixed effects model
          wolfID = i
        )],
        type = "link",
        re.form = NULL # this is for the mixed effects model
      ), wolfID = i)]
    # )
  })
}

# run funcition on the model of interest
p.road.CDV.1day.1.indiv <- p.road.h1.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)

# saving everything together and labeled so I can merge all RSS values together later
h1.road.CDV.1day.indiv <- data.table(rbindlist(p.road.CDV.1day.1.indiv),
                                     ttd = '1 day', COD = 'CDV', var = 'road', x = seq(from = 0, to = 3000, length.out = 100))


# function for h2 (h2 is the same for all similar models basically)
p.h2.indiv <- function(ids, DT, mod, death, t2death){
  lapply(ids, function(i) {
    unique(
      DT[#wolfID == i,
        ,.(h2 = predict(
          mod,
          newdata = .SD[, .(
            ToD_start = factor('day', levels = levels(ToD_start)),
            log_sl = mean(log_sl),
            cos_ta = mean(cos_ta),
            # land_end_adj = factor('forest', levels = levels(land_end_adj)),
            propforest_end_adj = mean(propforest_end_adj, na.rm = T),
            propopen_end_adj = mean(propopen_end_adj, na.rm = T),
            propwet_end= mean(propwet_end, na.rm = T), 
            roadDist_end = median(roadDist_end, na.rm = T),
            distance2 = median(distance2, na.rm = T),
            nnDist_end = median(nnDist_end, na.rm = T),
            packDist_end = median(packDist_end, na.rm = T),
            COD = factor(death, levels = levels(COD)),
            ttd1 = t2death,
            wolf_step_id = NA, #mixed effects
            wolfID = i
          )],
          type = "link",
          re.form = NULL # mixed effects
        ), wolfID = i)]
    )
  })}


# run predict model
p.CDV.1day.2.indiv <-p.h2.indiv(ids = CDV.wolfID, DT = dat, mod = everyone, death = 'CDV', t2death = 1)

# saving everything together and labeled so I can merge all RSS values together later
h2.road.CDV.1day.indiv <- data.table(rbindlist(p.CDV.1day.2.indiv),
                                     ttd = '1 day', COD = 'CDV', var = 'road')


# merge h1 and h2
logRSS.road.CDV.1day.indiv <- merge(h1.road.CDV.1day.indiv, h2.road.CDV.1day.indiv, by = c('wolfID', 'ttd', 'COD', 'var'), all.x = T)

# calculate RSS
logRSS.road.CDV.1day.indiv[,'rss'] <- logRSS.road.CDV.1day.indiv$h1 - logRSS.road.CDV.1day.indiv$h2


