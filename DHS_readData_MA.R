source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')


## ---- Read in the GIS data (shapefiles)
## Lynda Radke has provided a Middle Arm shapefile
middle_arm.sp<- maptools:::readShapeSpatial("data/GIS/MA_WA_StudyArea.shp",
                                            proj4string = CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'),
                                            repair=TRUE,force_ring=T,verbose=TRUE)
temp_data = middle_arm.sp@data %>% summarize(OBJECTID = 1, DH=1, Shape_Leng=mean(Shape_Leng), Shape_Area=mean(Shape_Area),Name='Middle Arm')
middle_arm.sp  <- unionSpatialPolygons(middle_arm.sp ,ID=rep(1, length(middle_arm.sp$Zone_Name)))
middle_arm.sp = SpatialPolygonsDataFrame(middle_arm.sp, data=temp_data)
middle_arm.sp <- spTransform(middle_arm.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
middle_arm.fort <- broom::tidy(middle_arm.sp,region='Name')
g = ggplot(middle_arm.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black') #+ coord_map()
save(middle_arm.sp, file='data/processed/middle_arm.sp.RData')
middle_arm.df = broom::tidy(middle_arm.sp, region='Name')
save(middle_arm.df, file='data/processed/middle_arm.df.RData')

## The previous version had a map of the entire harbour that outlined the focal areas
## So we need to add in the East Arm
load(file='data/processed/east_arm.sp.RData')
DH.zone_MA.sp <- rbind(middle_arm.sp, east_arm.sp, makeUniqueIDs = TRUE)
save(DH.zone_MA.sp, file='data/processed/DH.zone_MA.sp.RData')
DH.zone_MA.df <- broom::tidy(DH.zone_MA.sp, region='Name')
save(DH.zone_MA.df, file='data/processed/DH.zone_MA.df.RData')

spatial <- read.csv('config/spatial.csv', strip.white=TRUE)
centroids = data.frame(ZoneName=DH.zone.sp$Name, gCentroid(DH.zone.sp,byid=TRUE))
spatial = spatial %>% left_join(centroids)

sf.poly_MA=st_as_sf(DH.zone_MA.sp)
spatial = spatial %>% filter(RegionName!='Outer', ZoneName != 'Shoal Bay', ZoneName != 'Elizabeth River') %>%
    mutate(ZoneName=ifelse(ZoneName=='Outer Harbour', 'Outer Harbour\n and Shoal Bay',ZoneName)) %>%
    mutate(Lab_lat = ifelse(ZoneName=='Blackmore River', -12.7,
                     ifelse(ZoneName=='West Arm', -12.65,
                     ifelse(ZoneName=='Middle Harbour', -12.4, Lab_lat))),
           Lab_long = ifelse(ZoneName=='West Arm', 130.75,
                      ifelse(ZoneName=='Middle Harbour', 130.75, Lab_long)),
           hjust = ifelse(ZoneName=='West Arm', 0.5,
                   ifelse(ZoneName=='Middle Harbour', 0.5, hjust)),
           vjust = ifelse(ZoneName=='West Arm', 1,
                      ifelse(ZoneName=='Middle Harbour', 0, vjust))
           )
g = ggplot() +
    geom_sf(data=sf.poly_MA, alpha=0.8, aes(fill=Name), show.legend=FALSE) +
    #scale_fill_manual(breaks=sf.poly$Name, values=c('grey40', 'grey90', 'grey40', 'grey90', 'grey40')) +
    geom_point(data=spatial, aes(y=y, x=x),show.legend=FALSE) +
    geom_segment(data=spatial, aes(y=y, x=x, yend=Lab_lat, xend=Lab_long)) +
    geom_text(data=spatial[,], aes(y=Lab_lat, x=Lab_long, label=ZoneName, hjust=hjust, vjust=vjust)) +
    coord_sf() +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
g
ggsave('output/Map_MA1.pdf', g, width=7, height=8)
ggsave('output/Map_MA1.png', g, width=7, height=8, units='in', dpi=300)

## ----end

## ---- Read in the GIS data (shapefiles)
## There are some areas that have a hard seabed that are not ideal for sampling, these should be masked out
hard_seabed_MA.sp<- maptools:::readShapeSpatial("data/GIS/Hard_Seabed_MA.shp",
                                          proj4string = CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
hard_seabed_MA.sp  <- unionSpatialPolygons(hard_seabed_MA.sp ,ID=rep(1, length(hard_seabed_MA.sp$gridcode)))
hard_seabed_MA.sp <- spTransform(hard_seabed_MA.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
hard_seabed_MA.fort <- broom::tidy(hard_seabed_MA.sp,region='Name')
g = ggplot(hard_seabed_MA.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black') + coord_map()
save(hard_seabed_MA.sp, file='data/processed/hard_seabed_MA.sp.RData')

## Shipping Chanels Merged
shipping_MA.sp<- maptools:::readShapeSpatial("data/GIS/All_MiddleArm_Channels.shp",
                                          proj4string = CRS('+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
shipping_MA.sp  <- unionSpatialPolygons(shipping_MA.sp ,ID=rep(1, length(shipping_MA.sp$Id)))
shipping_MA.sp <- spTransform(shipping_MA.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
shipping_MA.fort <- broom::tidy(shipping_MA.sp,region='Name')
g = ggplot(shipping_MA.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black') + coord_map()
save(shipping_MA.sp, file='data/processed/shipping_MA.sp.RData')

## Exclusion Zones Merged
exclusion_MA.sp<- maptools:::readShapeSpatial("data/GIS/DLNG_Marine_Traffic_Exclusion_Zones.shp",
                                          proj4string = CRS('+proj=utm +zone=52 +south +datum=WGS84 +units=m +no_defs'),repair=TRUE,force_ring=T,verbose=TRUE)
exclusion_MA.sp  <- unionSpatialPolygons(exclusion_MA.sp ,ID=rep(1, length(exclusion_MA.sp$Id)))
exclusion_MA.sp <- spTransform(exclusion_MA.sp, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
exclusion_MA.fort <- broom::tidy(exclusion_MA.sp,region='Name')
g = ggplot(exclusion_MA.fort, aes(y=lat, x=long, group=group)) + geom_polygon(color='black') + coord_map()
save(exclusion_MA.sp, file='data/processed/exclusion_MA.sp.RData')

## ----end

# The above wont be masked out until after the modelling stage.  The reason for this is that the rest of the
## polygons are used to define the barrier models, howver these exclusions are not barriers to physical processes -
## rather they are exclusion zones for sampling.  Therefore, they will only be considered in the determination
## of sampling domains in the DHS_sampling.R routines.


## Munksgaard data ==========================================================
## ---- Read in the site data
sites_MA = read_csv('data/MiddleArm_Final_Appendix Darwin Harbour sediment study 2019_location.csv', trim_ws=TRUE,
                    col_types = cols(Sample=readr::col_factor(NULL),
                                     Date=readr::col_date(format='%d/%m/%y'),
                                     Area=readr::col_factor(NULL))) 
head(sites_MA)
save(sites_MA, file='data/processed/sites_MA.RData')
g=ggplot(middle_arm.fort) +
    geom_polygon(color='black',aes(y=lat, x=long))+
    geom_point(data=sites_MA, aes(y=Latitude, x=Longitude), color='red') +
    coord_map()
##g


chem_MA = readr::read_csv('data/primary/MiddleArm_Final_Appendix Darwin Harbour sediment study 2019.csv', trim_ws=TRUE,
                 col_types = cols(Sample=readr::col_factor(NULL)))
chem_MA
## some of the fields have < entries and so are considered character vectors.
## the > values indicate that the values are less than the limit of detection.
## I will replace them with half the next smallest value for that chemical.
chem_MA = chem_MA %>%
    dplyr::select(-Date, -Area) %>% 
    mutate_if(is.character, function(x) {
    xx = x[x!='<'];
    Min=min(as.numeric(xx), na.rm=TRUE);
    x[x=='<'] = 0.5*Min;
    as.numeric(x)
    })
sediment.data_MA = chem_MA
save(sediment.data_MA, file='data/processed/sediment.data_MA.RData')
## ----end

## Designated data ==========================================================
## ---- Read in the site data
designated_MA = read_csv('data/MIddleArm_DesignatedSite_Locations.csv', trim_ws=TRUE) 
head(designated_MA)

g=ggplot(middle_arm.fort) +
    geom_polygon(color='black',aes(y=lat, x=long))+
    geom_point(data=designated_MA, aes(y=Latitude, x=Longitude), color='red') +
    coord_map()
##g
priority_sites_MA = designated_MA
## ----end

## ---- Site maps
g=ggplot() +
    geom_sf(data=sf.poly_MA, fill='grey99') +
    coord_sf() + #(xlim=range(c(sediment.data_MA$Longitude,outerharbour.fort$long)),ylim=range(c(sediment.data_MA$Latitude,outerharbour.fort$lat))) +
    geom_sf(data=sf.poly_MA %>% filter(Name=='Middle Arm'), alpha=0.8) +
    geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude)) +
    geom_point(data=priority_sites_MA, aes(y=Latitude, x=Longitude), color='red', shape=16) +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
g

ggsave('output/Map2_MA.pdf', g, width=7, height=6)
ggsave('output/Map2_MA.png', g, width=7, height=6, units='in', dpi=300)

g=ggplot() +
    geom_sf(data=sf.poly_MA %>% filter(Name=='Middle Arm'), fill='#336A98') +
                                        #geom_sf(data=sf.poly_MA, fill='grey99') +
    geom_sf(data=st_as_sf(hard_seabed_MA.sp), fill='white', color='white') +
    geom_sf(data=st_as_sf(shipping_MA.sp), fill='white', color='white') +
    geom_sf(data=st_as_sf(exclusion_MA.sp), fill='white', color='white') +
    geom_sf(data=sf.poly_MA %>% filter(Name=='Middle Arm'), fill=NA, color='black') +
    coord_sf() + #(xlim=range(c(sediment.data_MA$Longitude,outerharbour.fort$long)),ylim=range(c(sediment.data_MA$Latitude,outerharbour.fort$lat))) +
    geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude)) +
    geom_point(data=priority_sites_MA, aes(y=Latitude, x=Longitude), color='red', shape=16) +
    annotation_scale(location='bl') +
    annotation_north_arrow(location='bl', which_north='true',
                           pad_x=unit(0.05,'npc'),pad_y=unit(0.05,'npc'),
                           style=north_arrow_fancy_orienteering) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank())
g
ggsave('output/Map3_MA.pdf', g, width=7, height=6)
ggsave('output/Map3_MA.png', g, width=7, height=6, units='in', dpi=300)

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

##outer_east.sp = rbind(east_arm.sp,outerharbour.sp,makeUniqueIDs = TRUE)
hydro.middle_arm.sp = crop(mask(hydro, middle_arm.sp), middle_arm.sp)
hydro.middle_arm.df = rasterToPoints(hydro.middle_arm.sp) %>% as.data.frame
save(hydro.middle_arm.sp, file='data/processed/hydro.middle_arm.sp.RData')
save(hydro.middle_arm.df, file='data/processed/hydro.middle_arm.df.RData')

plot_func = function (df, name) {
    ggplot() +
        geom_tile(data=df, aes(y=y, x=x, fill=Value)) +
        coord_sf(crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84') +
        scale_fill_gradientn(name,colors=tim.colors(12)) +
        theme_bw()+
        theme(axis.title.y=element_blank(), axis.title.x=element_blank())
}
dat = hydro.middle_arm.df %>% gather(key=Variable, value=Value,-x,-y)
dat = dat %>% group_by(Variable) %>% nest() %>%
    mutate(g=map2(data, Variable, plot_func))
grid.arrange(grobs=dat$g)
ggsave(filename='output/hydro_middle_arm.pdf', grid.arrange(grobs=dat$g), width=10, height=7)
ggsave(filename='output/hydro_middle_arm.png', grid.arrange(grobs=dat$g), width=10, height=7, units='in', dpi=300)

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
waves_ubot.middle_arm.sp = b = crop(mask(waves_ubot, middle_arm.sp), middle_arm.sp)
waves_ubot.middle_arm.df = rasterToPoints(waves_ubot.middle_arm.sp) %>% as.data.frame
save(waves_ubot.middle_arm.sp, file='data/processed/waves_ubot.middle_arm.sp.RData')
save(waves_ubot.middle_arm.df, file='data/processed/waves_ubot.middle_arm.df.RData')
waves.middle_arm.sp = b = crop(mask(waves, middle_arm.sp), middle_arm.sp)
waves.middle_arm.df = rasterToPoints(waves.middle_arm.sp) %>% as.data.frame
save(waves.middle_arm.sp, file='data/processed/waves.middle_arm.sp.RData')
save(waves.middle_arm.df, file='data/processed/waves.middle_arm.df.RData')

dat = waves_ubot.middle_arm.df %>% gather(key=Variable, value=Value,-x,-y)
dat = dat %>% group_by(Variable) %>% nest() %>%
    mutate(g=map2(data, Variable, plot_func))
grid.arrange(grobs=dat$g)
ggsave(filename='output/waves_ubot_middle_arm.pdf', grid.arrange(grobs=dat$g), width=11, height=7)
ggsave(filename='output/waves_ubot_middle_arm.png', grid.arrange(grobs=dat$g), width=11, height=7, units='in', dpi=300)

## ----end
