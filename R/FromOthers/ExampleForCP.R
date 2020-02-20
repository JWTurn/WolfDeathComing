require(survival)
ModelList<-list()
AllResultsNoIntNoSL<-data.frame()
options(warn=1)
AnimalList<-unique(Basic15$IdSeason)
for (i in AnimalList){
  subdat <- subset(Basic15, IdSeason == i)
  PropWidei <- subset(PropWide, IdSeason == i)#prop wide is a table of the proportion of availability for each habitat type per indiviual per season
  print("Identify Use")
  lis<-PropWidei[,'PropLIS']
  cli<-PropWidei[,'PropCLI']
  pipetrans<-PropWidei[, 'PropPipeTrans']
  transport<-PropWidei[,'PropTransp']
  poly<-PropWidei[,'PropPoly']
  water<-PropWidei[,'PropWater']
  lis[is.na(lis)] <- 0
  cli[is.na(cli)] <- 0
  pipetrans[is.na(pipetrans)] <- 0
  transport[is.na(transport)] <- 0
  poly[is.na(poly)] <- 0
  water[is.na(water)] <- 0
  print("Begin fitting models")
  ##Try fitting the full model
  print("Try Full Model")
  if((lis >= 0.01) &
     (cli >= 0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly >=0.01) &
     (water >=0.01)){
    print("Full model")
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + Transport:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } else {model<-NA}
  
  ##If that didn't work, then fit the model without transport
  print("Try No Transport")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    print("No Transport")
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS
  print("Try no LIS")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    print("No LISl")
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + Transport:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without Poly
  print("Try No Poly")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without PipeTrans
  print("Try no PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >=0.01) &
     (poly >= 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + Transport:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without CLI
  print("Try no CLI")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans >= 0.01) &
     (transport >=0.01) &
     (poly >= 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + LIS:LogSL + Transport:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try No Transport or LIS")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or Poly
  print("Try No Transport or Poly")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or Pipetrans
  print("Try No Transport or PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS or Poly
  print("Try No LIS or Poly")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS or PipeTrans
  print("Try No LIS or PT")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly >= 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Transport:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA} 
  
  ##If that didn't work, then fit the model without Poly or PipeTrans
  print("Try No Poly or PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA} 
  
  ##If that didn't work, then fit the model without transport or CLI
  print("Try no C:I or Transport")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans >= 0.01) &
     (transport < 0.01) &
     (poly >= 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + Poly:LogSL + PipeTrans:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, LIS, or poly
  print("Try No Transport LIS or Poly")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, LIS or Pipetrans
  print("Try No Transport LIS or PT")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, poly or Pipetrans
  print("Try No Transport Poly or PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS, Poly or Pipetrans
  print("Try No LIS Poly or PT")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try No Transport CLI or LIS")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli < 0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only PT and Transport")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli < 0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only CLI and Transport")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >= 0.01) &
     (pipetrans < 0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only LIS and Transport")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only LIS and PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + PipeTrans:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only LIS and Poly")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >= 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, LIS, Poly or Pipetrans
  print("Try Only CLI")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only poly
  print("Try only Poly")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + Poly:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only PT
  print("Try only PT")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli < 0.01) &
     (pipetrans >= 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only LIS
  print("Try only PT")
  if((is.na(model)[[0.01]]) &
     (lis >=  0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only Transport
  print("Try only PT")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + Transport:LogSL + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without disturbance at all 
  print("Try No Disturbance")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water >=0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + StepWater:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##GOING INTO NO WATER
  print("Going into no water")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly >=0.01) &
     (water < 0.01)){
    print("Full model")
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + Transport:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport
  print("Try No Transport")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    print("No Transport")
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS
  print("Try no LIS")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    print("No LISl")
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + Transport:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without Poly
  print("Try No Poly")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without PipeTrans
  print("Try no PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >=0.01) &
     (poly >= 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + Transport:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without CLI
  print("Try no CLI")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans >= 0.01) &
     (transport >=0.01) &
     (poly >= 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + LIS:LogSL + Transport:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try No Transport or LIS")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or Poly
  print("Try No Transport or Poly")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + PipeTrans:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or Pipetrans
  print("Try No Transport or PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS or Poly
  print("Try No LIS or Poly")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS or PipeTrans
  print("Try No LIS or PT")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly >= 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Transport:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA} 
  
  ##If that didn't work, then fit the model without Poly or PipeTrans
  print("Try No Poly or PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA} 
  
  ##If that didn't work, then fit the model without transport or CLI
  print("Try no C:I or Transport")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans >= 0.01) &
     (transport < 0.01) &
     (poly >= 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + Poly:LogSL + PipeTrans:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, LIS, or poly
  print("Try No Transport LIS or Poly")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + PipeTrans:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, LIS or Pipetrans
  print("Try No Transport LIS or PT")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, poly or Pipetrans
  print("Try No Transport Poly or PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + LIS:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without LIS, Poly or Pipetrans
  print("Try No LIS Poly or PT")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try No Transport CLI or LIS")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli < 0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only PT and Transport")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli < 0.01) &
     (pipetrans >=0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only CLI and Transport")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >= 0.01) &
     (pipetrans < 0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only LIS and Transport")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport >=0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only LIS and PT")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans >=0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + PipeTrans:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport or LIS
  print("Try Only LIS and Poly")
  if((is.na(model)[[0.01]]) &
     (lis >= 0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >= 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without transport, LIS, Poly or Pipetrans
  print("Try Only CLI")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli >=0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + ConLowIce:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only poly
  print("Try only Poly")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly >=0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + Poly:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only PT
  print("Try only PT")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli < 0.01) &
     (pipetrans >= 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + PipeTrans:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only LIS
  print("Try only PT")
  if((is.na(model)[[0.01]]) &
     (lis >=  0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + LIS:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model with only Transport
  print("Try only PT")
  if((is.na(model)[[0.01]]) &
     (lis <  0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport >= 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + Transport:LogSL + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##If that didn't work, then fit the model without disturbance at all 
  print("Try No Disturbance")
  if((is.na(model)[[0.01]]) &
     (lis < 0.01) &
     (cli < 0.01) &
     (pipetrans < 0.01) &
     (transport < 0.01) &
     (poly < 0.01) &
     (water < 0.01)){ 
    subdat<-as.data.frame(subdat)
    model<-clogit(Case ~ LC3 + LogD2LF + LogD2Poly + LogD2Water + LogSL + LogSL:ToD2 + CosTA + LogSL:CosTA + strata(StartId), method = "approximate", data=subdat)
  } #else {model<-NA}
  
  ##Now extract the model information
  coefficient <- coef(summary(model))[,1]
  SE <- coef(summary(model))[,3]
  P<- coef(summary(model))[,5]
  ModelList[[i]] <- data.frame(ID = i, Coeff = round(coefficient, digits =8), SE = SE, P = P)
  ModelList[[i]][,5]<-row.names(ModelList[[i]])
  names(ModelList[[i]])[5]<-"Variable"
}
AllResultsNoIntNoSL = do.call(rbind, ModelList)

##EXAMPLE AVERAGING:
ModelList<-list()
AvePolyLogSL<-data.frame()
oPolyions(warn=1)
SSList<-unique(ScaledAvail$SS)
SSList<-as.factor(SSList)
for (i in SSList){
  subdat <- subset(ScaledAvail, ScaledAvail$SS == i)
  subdat <- subdat[!is.na(subdat$PolyLogSL),]
  sample <- nrow(subdat)
  if ((sample > 2)){
    wei_WSFglm_PolyLogSL = glm(PolyLogSL ~ ScaleLogSL + ScalePercPoly, data=subdat, weight=PolyLogSLWeight)
    coef(summary(wei_WSFglm_PolyLogSL))
    coefficient <- coef(summary(wei_WSFglm_PolyLogSL))[,1]
    SE <- coef(summary(wei_WSFglm_PolyLogSL))[,2]
    P<- coef(summary(wei_WSFglm_PolyLogSL))[,4]
    ModelList[[i]] <- data.frame(ID = i, Coeff = round(coefficient, digits =8), SE = SE, P = P)
    ModelList[[i]][,5]<-row.names(ModelList[[i]])
    names(ModelList[[i]])[5]<-"Variable"
  } else {
    ModelList[[i]] <- data.frame(ID = i, Coeff = NA, SE = NA, P = NA)
    #ModelList[[i]][,5]<-row.names(ModelList[[i]])
    ModelList[[i]][,5]<-NA
    names(ModelList[[i]])[5]<-"Variable"
  }
}
AvePolyLogSL = do.call(rbind, ModelList)