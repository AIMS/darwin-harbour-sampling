source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/hydro.outer_east.sp.RData')
load(file='data/processed/hydro.outer_east.df.RData')
load(file='data/processed/waves.outer_east.sp.RData')
load(file='data/processed/waves_ubot.outer_east.df.RData')
load(file='data/processed/waves_ubot.outer_east.sp.RData')
load(file='data/processed/waves.outer_east.df.RData')
load(file='data/processed/sediment.data.RData')
load(file='data/processed/outer_sites.RData')
load(file='data/processed/DH.zone.sp.RData')
load(file='data/processed/DH.zone.df.RData')

## Dissolve into a single polygon
DH.wh.sp  <- unionSpatialPolygons(DH.zone.sp ,ID=rep(1, length(DH.zone.sp$Name)))
DH.wh.df  <- broom::tidy(DH.wh.sp)

save(DH.wh.sp, file='data/processed/DH.wh.sp.RData')
save(DH.wh.df, file='data/processed/DH.wh.df.RData')


East.arm.sp = subset(DH.zone.sp, Name=='East Arm')
East.arm.df = broom::tidy(East.arm.sp)
save(East.arm.sp, file='data/processed/East.arm.sp.RData')
save(East.arm.df, file='data/processed/East.arm.df.RData')
wch=point.in.polygon(sediment.data$Longitude, sediment.data$Latitude, East.arm.df$long, East.arm.df$lat)
sediment.data.east = sediment.data[wch==1,]
save(sediment.data.east, file='data/processed/sediment.data.east.RData')

Outer.sp = subset(DH.zone.sp, Name=='Outer Harbour')
Outer.df = broom::tidy(Outer.sp)
save(Outer.sp, file='data/processed/Outer.sp.RData')
save(Outer.df, file='data/processed/Outer.df.RData')
wch=point.in.polygon(sediment.data$Longitude, sediment.data$Latitude, Outer.df$long, Outer.df$lat)
sediment.data.outer = sediment.data[wch==1,]
save(sediment.data.outer, file='data/processed/sediment.data.outer.RData')

wch=point.in.polygon(outer_sites$Longitude, outer_sites$Latitude, Outer.df$long, Outer.df$lat)
outer_sites.outer = outer_sites[wch==1,]
save(outer_sites.outer, file='data/processed/outer_sites.outer.RData')


## ggplot() +
##   geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill='white', color='black') +
##   geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude))

## Now crop the hydro and wave data to East Arm
a1=subset(DH.zone.sp, Name=='East Arm')
hydro_east.sp=crop(mask(hydro.outer_east.sp, a1), a1)
save(hydro_east.sp, file='data/processed/hydro_east.sp.RData')
waves_east.sp=crop(mask(waves.outer_east.sp, a1), a1)
save(waves_east.sp, file='data/processed/waves_east.sp.RData')
waves_ubot_east.sp=crop(mask(waves_ubot.outer_east.sp, a1), a1)
save(waves_ubot_east.sp, file='data/processed/waves_ubot_east.sp.RData')

a1=subset(DH.zone.sp, Name=='Outer Harbour')
hydro_outer.sp=crop(mask(hydro.outer_east.sp, a1), a1)
save(hydro_outer.sp, file='data/processed/hydro_outer.sp.RData')
waves_outer.sp=crop(mask(waves.outer_east.sp, a1), a1)
save(waves_outer.sp, file='data/processed/waves_outer.sp.RData')
waves_ubot_outer.sp=crop(mask(waves_ubot.outer_east.sp, a1), a1)
save(waves_ubot_outer.sp, file='data/processed/waves_ubot_outer.sp.RData')
waves_ubot.outer.df = rasterToPoints(waves_ubot.outer.sp) %>% as.data.frame
save(waves_ubot_outer.sp, file='data/processed/waves_ubot_outer.df.RData')
