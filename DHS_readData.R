source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

## ---- Read in the GIS data (shapefiles)

##,----------------
##| Blackmore River
##`----------------
blackmore.sp<- maptools:::readShapeSpatial("data/GIS/blackmore_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
blackmore.sp[['Name']] <- 'Blackmore River'
blackmore.sp.utm = blackmore.sp
blackmore.sp <- spTransform(blackmore.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
blackmore.fort <- broom::tidy(blackmore.sp,region='Name')
g = ggplot(blackmore.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black') + coord_map()
#g

##,---------
##| East Arm
##`---------
east_arm.sp<- maptools:::readShapeSpatial("data/GIS/east_arm_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
east_arm.sp[['Name']] <- 'East Arm'
east_arm.sp.utm=east_arm.sp
east_arm.sp <- spTransform(east_arm.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
east_arm.fort <- broom::tidy(east_arm.sp)
g=ggplot(east_arm.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black')+coord_map()
##g

##,----------------
##| Elizabeth River
##`----------------
elizabeth.sp<- maptools:::readShapeSpatial("data/GIS/elizabeth_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
elizabeth.sp[['Name']] <- 'Elizabeth River'
elizabeth.sp.utm=elizabeth.sp
elizabeth.sp <- spTransform(elizabeth.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
elizabeth.fort <- broom::tidy(elizabeth.sp)
g=ggplot(elizabeth.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black')+coord_map()
##g

##,---------------
##| Middle Harbour
##`---------------
middleharbour.sp<- maptools:::readShapeSpatial("data/GIS/middleharbour_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
middleharbour.sp[['Name']] <- 'Middle Harbour'
middleharbour.sp.utm=middleharbour.sp
middleharbour.sp <- spTransform(middleharbour.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
middleharbour.fort <- broom::tidy(middleharbour.sp)
g=ggplot(middleharbour.fort, aes(y=lat, x=long, group=group,fill=hole)) + geom_polygon(color='black',show.legend=FALSE)+scale_fill_manual(values=c('grey40','white'))+coord_map()
##g

##,--------------
##| Outer Harbour
##`--------------
outerharbour.sp<- maptools:::readShapeSpatial("data/GIS/outerharbour_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
outerharbour.sp[['Name']] <- 'Outer Harbour'
outerharbour.sp.utm=outerharbour.sp
outerharbour.sp <- spTransform(outerharbour.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
outerharbour.fort <- broom::tidy(outerharbour.sp)
g=ggplot(outerharbour.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black')+coord_map()
##g

##,----------
##| Shoal Bay
##`----------
shoalbay.sp<- maptools:::readShapeSpatial("data/GIS/shoalbay_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
shoalbay.sp[['Name']] <- 'Shoal Bay'
shoalbay.sp.utm=shoalbay.sp
shoalbay.sp <- spTransform(shoalbay.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
shoalbay.fort <- broom::tidy(shoalbay.sp)
g=ggplot(shoalbay.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black')+coord_map()
##g

##,------------
##| Western Arm
##`------------
westarm.sp<- maptools:::readShapeSpatial("data/GIS/westarm_mask.shp",
                                           proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
westarm.sp[['Name']] <- 'West Arm'
westarm.sp.utm=westarm.sp
westarm.sp <- spTransform(westarm.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
westarm.fort <- broom::tidy(westarm.sp)
g=ggplot(westarm.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black')+coord_map()
##g

##################################################################
## Combine these all into a single Zone chunked Spatial Polygon ##
## Note, the Zones must be added alphebetically                 ##
##################################################################
DH.zone.sp <- rbind(blackmore.sp,east_arm.sp,elizabeth.sp,middleharbour.sp,outerharbour.sp,shoalbay.sp,westarm.sp,makeUniqueIDs = TRUE)
save(DH.zone.sp, file='data/processed/DH.zone.sp.RData')
DH.zone.df <- broom::tidy(DH.zone.sp, region='Name')
spatial <- read.csv('config/spatial.csv', strip.white=TRUE)
centroids = data.frame(ZoneName=DH.zone.sp$Name, gCentroid(DH.zone.sp,byid=TRUE))
spatial = spatial %>% left_join(centroids)

sf.poly=st_as_sf(DH.zone.sp)
                                        #sf.point=st_as_sf(spatial, coords=c()
#register_tile_source(osmn = "http://a.tiles.wmflabs.org/osm-no-labels/${z}/${x}/${y}.png")
g=ggplot() +
    #annotation_map_tile(type='osm', zoom=10) +
    #annotation_map_tile(type='thunderforestoutdoors', zoom=10) +
    #annotation_map_tile(type='osmn', zoom=10) +
    #annotation_map_tile(type='stamenwatercolor', zoom=10) +
    geom_sf(data=sf.poly, alpha=0.8, aes(fill=Name), show.legend=FALSE) +
    scale_fill_manual(breaks=sf.poly$Name, values=c('grey90', 'grey40', 'grey90', 'grey90', 'grey40', 'grey90', 'grey90')) +
    geom_point(data=spatial, aes(y=y, x=x),show.legend=FALSE) +
    geom_segment(data=spatial, aes(y=y, x=x, yend=Lab_lat, xend=Lab_long)) +
    geom_text(data=spatial, aes(y=Lab_lat, x=Lab_long, label=ZoneName, hjust=hjust, vjust=vjust)) +
    coord_sf() +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
ggsave('output/Map1.pdf', g, width=7, height=8)
ggsave('output/Map1.png', g, width=7, height=8, units='in', dpi=300)

save(DH.zone.df, file='data/processed/DH.zone.df.RData')

Zones.df = data.frame(ZoneName=unique(DH.zone.df$id),rgeos:::gCentroid(DH.zone.sp, byid=TRUE))
spatial <- read.csv('config/spatial.csv', strip.white=TRUE)

## the following generates a temporary figure of the same dimensions
## as the intended map while it is working out how much space the
## labels will take up on a map of that size (given the font size)
pdf(paste0('junk.pdf'), width=7, height=7)
Zones.df = Zones.df %>% full_join(spatial,by='ZoneName') %>% mutate(w = strwidth(ZoneName, cex=0.7,units='figure'),h = strheight(ZoneName, cex=0.7,units='figure'))
dev.off()
## Remove the temp figure
system('rm junk.pdf')
save(Zones.df, file='data/processed/Zones.df.RData')


####################################################################
## Combine all the Zones into Regions and produce a single Region ##
## chunked Spatial Polygon.                                       ##
####################################################################
Region1 = rgeos:::gUnion(blackmore.sp,east_arm.sp,byid=TRUE,id='Upper')
Region1 = rgeos:::gUnion(Region1,elizabeth.sp, byid=TRUE,id='Upper')
Region1 = rgeos:::gUnion(Region1,westarm.sp, byid=TRUE,id='Upper')
Region2 = rgeos:::gUnion(middleharbour.sp,middleharbour.sp, byid=TRUE,id='Middle') #need to union with itself to allow the rbind below to work
Region3 = rgeos:::gUnion(outerharbour.sp,shoalbay.sp,byid=TRUE,id='Outer')
DH.region.sp <- rbind(Region1,Region2,Region3,makeUniqueIDs = TRUE)
DH.region.fort <- broom::tidy(DH.region.sp, region='ID')
ggplot(DH.region.fort, aes(y=lat, x=long, group=group,fill=id)) + geom_polygon(color='black')+coord_map()
save(DH.region.fort, file='data/processed/DH.region.fort.RData')

Regions.df = data.frame(RegionName=unique(DH.region.fort$id),rgeos:::gCentroid(DH.region.sp, byid=TRUE))
spatial <- read.csv('config/spatial.csv', strip.white=TRUE)
## the following generates a temporary figure of the same dimensions
## as the intended map while it is working out how much space the
## labels will take up on a map of that size (given the font size)
pdf(paste0('junk.pdf'), width=7, height=7)
Regions.df = Regions.df %>% full_join(spatial,by='RegionName') %>% mutate(w = strwidth(RegionName, cex=0.7,units='figure'),h = strheight(RegionName, cex=0.7,units='figure'))
dev.off()
## Remove the temp figure
system('rm junk.pdf')
save(Regions.df, file='data/processed/Regions.df.RData')

## Get a google map
#map <- ggmap:::get_googlemap(center=c(mean(bbox(DH.zone.sp)[1,]),mean(bbox(DH.zone.sp)[2,])), scale=2,zoom=10, maptype="roadmap", style='feature:road|visibility:off&style=element:labels|visibility:off')
#save(map, file='data/processed/map.RData')

## ----end


## ---- Read in the site data
sites = read_csv('data/primary/Munksgaard_sites.csv', trim_ws=TRUE,
                 col_types = cols(Sample=col_factor(NULL),
                                  Zone=col_factor(NULL),
                                  Date=col_date(format='%d/%m/%y'),
                                  Datum=col_factor(NULL)))
head(sites)

g=ggplot(DH.region.fort) +
    geom_polygon(color='black',aes(y=lat, x=long, group=group,fill=id))+
    geom_point(data=sites, aes(y=Latitude, x=Longitude)) +
    coord_map()
##g

chem = read_csv('data/primary/Munksgaard_chem.csv', trim_ws=TRUE,
                 col_types = cols(Sample=col_factor(NULL)))
chem
## some of the fields have < entries and so are considered character vectors.
## the > values indicate that the values are less than the limit of detection.
## I will replace them with half the next smallest value for that chemical.
chem = chem %>% mutate_if(is.character, function(x) {
    xx = x[x!='<'];
    Min=min(as.numeric(xx), na.rm=TRUE);
    x[x=='<'] = 0.5*Min;
    as.numeric(x)
    })

## Join this to the sites data
sediment.data = sites %>% full_join(chem)
save(sediment.data, file='data/processed/sediment.data.RData')

g=ggplot() +
    geom_sf(data=sf.poly, alpha=0.8) +
    coord_sf() +
    geom_point(data=sediment.data, aes(y=Latitude, x=Longitude)) +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())

## East Arm samples
East.arm.sp = subset(DH.zone.sp, Name=='East Arm')
East.arm.df = broom::tidy(East.arm.sp)
wch=point.in.polygon(sediment.data$Longitude, sediment.data$Latitude, East.arm.df$long, East.arm.df$lat)
sediment.data.east = sediment.data[wch==1,]
sf.poly=st_as_sf(East.arm.sp)
g1=ggplot() +
    geom_sf(data=sf.poly, alpha=0.8) +
    coord_sf() +
    geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude)) +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
ggsave('output/Map2a.pdf', g1, width=7, height=6)
ggsave('output/Map2a.png', g1, width=7, height=6, units='in', dpi=300)

## East Arm samples
Outer.sp = subset(DH.zone.sp, Name=='Outer Harbour')
Outer.df = broom::tidy(Outer.sp)
wch=point.in.polygon(sediment.data$Longitude, sediment.data$Latitude, Outer.df$long, Outer.df$lat)
sediment.data.outer = sediment.data[wch==1,]
sf.poly=st_as_sf(Outer.sp)
g2=ggplot() +
    geom_sf(data=sf.poly, alpha=0.8) +
    coord_sf() +
    geom_point(data=sediment.data.outer, aes(y=Latitude, x=Longitude)) +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
ggsave('output/Map2b.pdf', g2, width=7, height=6)
ggsave('output/Map2b.png', g2, width=7, height=6, units='in', dpi=300)

grid.arrange(g1,g2, nrow=1)

sf.poly=st_as_sf(rbind(Outer.sp, East.arm.sp))
s.d = rbind(data.frame(Name='Outer Harbour', sediment.data.outer), data.frame(Name='East Arm', sediment.data.east))
g=ggplot() +
   #annotate(geom='polygon', x=blackmore.fort$long,y=blackmore.fort$lat, group=blackmore.fort$group, fill=NA, color='black') +
    geom_sf(data=st_as_sf(rbind(Outer.sp, East.arm.sp, westarm.sp, middleharbour.sp,blackmore.sp,elizabeth.sp)), fill='grey99') +
    geom_point(data=sediment.data, aes(y=Latitude, x=Longitude), color='grey40', fill=NA, shape=21) +
    geom_sf(data=sf.poly, alpha=0.8) +
    coord_sf(xlim=range(c(sediment.data$Longitude,outerharbour.fort$long)),ylim=range(c(sediment.data$Latitude,outerharbour.fort$lat))) +
    geom_point(data=s.d, aes(y=Latitude, x=Longitude)) +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
ggsave('output/Map2.pdf', g, width=7, height=6)
ggsave('output/Map2.png', g, width=7, height=6, units='in', dpi=300)

## ----end


## ---- Hydrodynamic data
bedShear_50p <- raster('data/Hydrodynamics/bedShear_50p.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
bedShear_75p <- raster('data/Hydrodynamics/bedShear_75p.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
velocity_50p <- raster('data/Hydrodynamics/velocity_50p.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
velocity_75p <- raster('data/Hydrodynamics/velocity_75p.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
hydro = brick(bedShear_50p,bedShear_75p,velocity_50p,velocity_75p)

##change the projection into +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
hydro = projectRaster(hydro,
              crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84') 
## now we need to apply a mask to narrow this to the focal area

outer_east.sp = rbind(east_arm.sp,outerharbour.sp,makeUniqueIDs = TRUE)
hydro.outer_east.sp = b = crop(mask(hydro, outer_east.sp), outer_east.sp)
hydro.outer_east.df = rasterToPoints(hydro.outer_east.sp) %>% as.data.frame
save(hydro.outer_east.sp, file='data/processed/hydro.outer_east.sp.RData')
save(hydro.outer_east.df, file='data/processed/hydro.outer_east.df.RData')

hydro.east.sp = crop(mask(hydro, east_arm.sp), east_arm.sp)
hydro.east.df = rasterToPoints(hydro.east.sp) %>% as.data.frame
save(hydro.east.sp, file='data/processed/hydro.east.sp.RData')
save(hydro.east.df, file='data/processed/hydro.east.df.RData')

plot_func = function (df, name) {
    ggplot() +
        geom_tile(data=df, aes(y=y, x=x, fill=Value)) +
        coord_sf(crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84') +
        scale_fill_gradientn(name,colors=tim.colors(12)) +
        theme_bw()+
        theme(axis.title.y=element_blank(), axis.title.x=element_blank())
}
dat = hydro.east.df %>% gather(key=Variable, value=Value,-x,-y)
dat = dat %>% group_by(Variable) %>% nest() %>%
    mutate(g=map2(data, Variable, plot_func))
grid.arrange(grobs=dat$g)
ggsave(filename='output/hydro_east.pdf', grid.arrange(grobs=dat$g), width=10, height=7)
ggsave(filename='output/hydro_east.png', grid.arrange(grobs=dat$g), width=10, height=7, units='in', dpi=300)

hydro.outer.sp = crop(mask(hydro, outerharbour.sp), outerharbour.sp)
hydro.outer.df = rasterToPoints(hydro.outer.sp) %>% as.data.frame
save(hydro.outer.sp, file='data/processed/hydro.outer.sp.RData')
save(hydro.outer.df, file='data/processed/hydro.outer.df.RData')
dat = hydro.outer.df %>% gather(key=Variable, value=Value,-x,-y)
dat = dat %>% group_by(Variable) %>% nest() %>%
    mutate(g=map2(data, Variable, plot_func))
grid.arrange(grobs=dat$g)
ggsave(filename='output/hydro_outer.pdf', grid.arrange(grobs=dat$g), width=10, height=5)
ggsave(filename='output/hydro_outer.png', grid.arrange(grobs=dat$g), width=10, height=5, units='in', dpi=300)

## ----end

## ---- Wave data
beagle_10_00_Ubot <- raster('data/Waves/beagle_10_00_Ubot.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
beagle_10_90_Ubot <- raster('data/Waves/beagle_10_90_Ubot.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
beagle_10_140_Ubot <- raster('data/Waves/beagle_10_140_Ubot.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
beagle_10_270_Ubot <- raster('data/Waves/beagle_10_270_Ubot.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
beagle_10_315_Ubot <- raster('data/Waves/beagle_10_315_Ubot.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))

waveShear_10_00 <- raster('data/Waves/waveShear_10_00.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
waveShear_10_90 <- raster('data/Waves/waveShear_10_90.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
waveShear_10_140 <- raster('data/Waves/waveShear_10_140.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
waveShear_10_270 <- raster('data/Waves/waveShear_10_270.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))
waveShear_10_315 <- raster('data/Waves/waveShear_10_315.tif', crs=CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'))

waves_ubot = brick(beagle_10_00_Ubot,beagle_10_90_Ubot,beagle_10_140_Ubot,beagle_10_270_Ubot,beagle_10_315_Ubot)
waves = brick(waveShear_10_00, waveShear_10_90, waveShear_10_140, waveShear_10_270, waveShear_10_315)

##change the projection into +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
waves_ubot = projectRaster(waves_ubot,
                      crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84')
waves = projectRaster(waves,
              crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84') 
## now we need to apply a mask to narrow this to the focal area

outer_east.sp = rbind(east_arm.sp,outerharbour.sp,makeUniqueIDs = TRUE)
waves_ubot.outer_east.sp = b = crop(mask(waves_ubot, outer_east.sp), outer_east.sp)
waves_ubot.outer_east.df = rasterToPoints(waves_ubot.outer_east.sp) %>% as.data.frame
save(waves_ubot.outer_east.sp, file='data/processed/waves_ubot.outer_east.sp.RData')
save(waves_ubot.outer_east.df, file='data/processed/waves_ubot.outer_east.df.RData')
waves.outer_east.sp = b = crop(mask(waves, outer_east.sp), outer_east.sp)
waves.outer_east.df = rasterToPoints(waves.outer_east.sp) %>% as.data.frame
save(waves.outer_east.sp, file='data/processed/waves.outer_east.sp.RData')
save(waves.outer_east.df, file='data/processed/waves.outer_east.df.RData')

waves_ubot.east.sp = b = crop(mask(waves_ubot, east_arm.sp), east_arm.sp)
waves_ubot.east.df = rasterToPoints(waves_ubot.east.sp) %>% as.data.frame
save(waves_ubot.east.sp, file='data/processed/waves_ubot.east.sp.RData')
save(waves_ubot.east.df, file='data/processed/waves_ubot.east.df.RData')

dat = waves_ubot.east.df %>% gather(key=Variable, value=Value,-x,-y)
dat = dat %>% group_by(Variable) %>% nest() %>%
    mutate(g=map2(data, Variable, plot_func))
grid.arrange(grobs=dat$g)
ggsave(filename='output/waves_ubot_east.pdf', grid.arrange(grobs=dat$g), width=10, height=10.5)
ggsave(filename='output/waves_ubot_east.png', grid.arrange(grobs=dat$g), width=10, height=10, units='in', dpi=300)


waves_ubot.outer.sp = b = crop(mask(waves_ubot, outerharbour.sp), outerharbour.sp)
waves_ubot.outer.df = rasterToPoints(waves_ubot.outer.sp) %>% as.data.frame
dat = waves_ubot.outer.df %>% gather(key=Variable, value=Value,-x,-y)
dat = dat %>% group_by(Variable) %>% nest() %>%
    mutate(g=map2(data, Variable, plot_func))
grid.arrange(grobs=dat$g)
ggsave(filename='output/waves_ubot_outer.pdf', grid.arrange(grobs=dat$g), width=11, height=7)
ggsave(filename='output/waves_ubot_outer.png', grid.arrange(grobs=dat$g), width=11, height=7, units='in', dpi=300)



## ----end

## ---- Grain size
grainsize = read_csv('data/primary/grainsize.csv', trim_ws=TRUE)
head(grainsize)

nms=names(grainsize)[-1]
grainsize = grainsize %>% gather(key=Size,value=Value, -SampleName) %>%
    mutate(Size=gsub('\\(','\n(',Size),
           Size=gsub('%','',Size),
           Size=factor(Size, levels=unique(Size)))
grainsize %>% group_by(Size) %>% summarise(mean(Value, na.rm=TRUE))

g=ggplot() +
    geom_boxplot(data=grainsize, aes(y=Value, x=Size)) +
    scale_y_continuous('% Abundance in sediments') +
    scale_x_discrete('Grain size category')+
    theme_bw()
ggsave(filename='output/particle_sizes.pdf', g, width=8, height=3)
ggsave(filename='output/particle_sizes.png', g, width=6, height=2, units='in', dpi=300)

    #geom_density(data=grainsize, aes(x=Value, fill=Size))
## ----end


## waves_ubot.outer_east.df %>%
##     ggplot() +
##     geom_tile(aes(y=y, x=x, fill=beagle_10_00_Ubot)) + scale_fill_gradientn(colors=rev(terrain.colors(10))) +
##     geom_point(data=sediment.data, aes(y=Latitude, x=Longitude))


## ggplot(bb, aes(y=y, x=x)) + geom_tile(aes(fill=bedShear_50p)) + scale_fill_gradientn(colors=rev(terrain.colors(10)))
## ggplot(bb, aes(y=y, x=x)) + geom_tile(aes(fill=bedShear_50p)) + scale_fill_gradientn(colors=(heat.colors(10)))

## cc = bb %>% gather(key=Variable, value=Value, -x,-y)
## g1 = cc %>% filter(Variable=='bedShear_50p') %>%
##     ggplot(aes(y=y, x=x)) + geom_tile(aes(fill=Value)) + scale_fill_gradientn(colors=(heat.colors(10)))

## g2 = cc %>% filter(Variable=='bedShear_75p') %>%
##     ggplot(aes(y=y, x=x)) + geom_tile(aes(fill=Value)) + scale_fill_gradientn(colors=(heat.colors(10)))

## g3 = cc %>% filter(Variable=='velocity_50p') %>%
##     ggplot(aes(y=y, x=x)) + geom_tile(aes(fill=Value)) + scale_fill_gradientn(colors=(heat.colors(10)))

## g4 = cc %>% filter(Variable=='velocity_75p') %>%
##     ggplot(aes(y=y, x=x)) + geom_tile(aes(fill=Value)) + scale_fill_gradientn(colors=(heat.colors(10)))

## gridExtra::grid.arrange(g1,g2, g3, g4, ncol=2)
## ----end
