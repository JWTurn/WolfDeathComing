W06 <- fread('W06_troubleshoot_dat.csv')
W06.model <- W06[ttd1>31, clogit(case_ ~ log_sl:ToD_start + land_end +log_sl:land_end +
                                   log(ttd1+1):log_sl + log(ttd1+1):cos_ta + strata(stepjum))]