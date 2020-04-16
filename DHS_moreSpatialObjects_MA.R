source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/hydro.middle_arm.sp.RData')
load(file='data/processed/hydro.middle_arm.df.RData')
load(file='data/processed/waves.middle_arm.sp.RData')
load(file='data/processed/waves_ubot.middle_arm.df.RData')
load(file='data/processed/waves_ubot.middle_arm.sp.RData')
load(file='data/processed/waves.middle_arm.df.RData')
load(file='data/processed/sediment.data_MA.RData')
load(file='data/processed/sites_MA.RData')
load(file='data/processed/DH.zone_MA.sp.RData')
load(file='data/processed/DH.zone_MA.df.RData')

## Dissolve into a single polygon
DH.wh_MA.sp  <- unionSpatialPolygons(DH.zone_MA.sp ,ID=rep(1, length(DH.zone_MA.sp$Name)))
DH.wh_MA.df  <- broom::tidy(DH.wh_MA.sp)

save(DH.wh_MA.sp, file='data/processed/DH.wh_MA.sp.RData')
save(DH.wh_MA.df, file='data/processed/DH.wh_MA.df.RData')

