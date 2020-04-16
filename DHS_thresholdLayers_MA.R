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

## Middle Arm
g=hydro.middle_arm.df %>% gather(key=Variable, value=Value, -x,-y) %>%
    ggplot()+
    geom_histogram(aes(x=Value), fill='grey80', color='grey40') +
    facet_wrap(~Variable, scales='free') +
    scale_y_continuous('Frequency') +
    scale_x_continuous('') +
    theme_bw()
ggsave(filename='output/hist_hydro_middle_arm.pdf', g, width=6, height=6)
ggsave(filename='output/hist_hydro_middle_arm.png', g, width=6, height=6, units='in', dpi=300)

hist(hydro.middle_arm.sp)
plot(hydro.middle_arm.sp)

g=waves_ubot.middle_arm.df %>% gather(key=Variable, value=Value, -x,-y) %>%
    ggplot()+
    geom_histogram(aes(x=Value), fill='grey80', color='grey40') +
    facet_wrap(~Variable, scales='free', nrow=3) +
    scale_y_continuous('Frequency') +
    scale_x_continuous('') +
    theme_bw()
ggsave(filename='output/hist_waves_middle_arm.pdf', g, width=6, height=9)
ggsave(filename='output/hist_waves_middle_arm.png', g, width=6, height=9, units='in', dpi=300)

hist(waves_ubot.middle_arm.sp)
plot(waves_ubot.middle_arm.sp)


## East Arm
rasters = list('bedShear_50p'=list(Name='bedShear_50p', Threshold=0.2),
               'bedShear_75p'=list(Name='bedShear_75p', Threshold=0.3),
               'velocity_50p'=list(Name='velocity_50p', Threshold=0.2),
               'velocity_75p'=list(Name='velocity_75p', Threshold=0.3),
               
               'beagle_10_00_Ubot'=list(Name='beagle_10_00_Ubot', Threshold=0.2),
               'beagle_10_90_Ubot'=list(Name='beagle_10_90_Ubot', Threshold=0.2),
               'beagle_10_140_Ubot'=list(Name='beagle_10_140_Ubot', Threshold=0.2),
               'beagle_10_270_Ubot'=list(Name='beagle_10_270_Ubot', Threshold=0.2),
               'beagle_10_315_Ubot'=list(Name='beagle_10_315_Ubot', Threshold=0.2))
raster.masks = list()
for (r in rasters) {
    print(r$Name)
    if (grepl('(bed.*|velocity.*)',r$Name)) {
        x = subset(hydro.middle_arm.sp, r$Name)
    } else {
        x = subset(waves_ubot.middle_arm.sp, r$Name)
    }
    median(values(x), na.rm=TRUE)
    quantile(values(x), p=0.99, na.rm=TRUE)
    threshold=r$Threshold #median(values(bedShear_50p), na.rm=TRUE) #0.2
    reclassification.matrix=rbind(c(0,threshold,1),
                                  c(threshold,Inf,NA))
    raster.masks[[r$Name]] = reclassify(x, rcl=reclassification.matrix)
}
nms = names(raster.masks)
names(raster.masks) <- NULL
raster.masks = do.call('stack', raster.masks)
names(raster.masks) <- nms

raster.masks.df = rasterToPoints(raster.masks) %>% as.data.frame %>% gather(key=Variable, value=Value, -x, -y)
g=ggplot() +
    geom_tile(data=raster.masks.df, aes(y=y, x=x, fill=Value), show.legend=FALSE) +
    geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~Variable) +
    theme_bw() +
    coord_equal() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename='output/masks_MA.pdf', g, width=10, height=8)
ggsave(filename='output/masks_MA.png', g, width=10, height=8, units='in', dpi=300)

# Combine together bedShear50p and all of the beagle_Ubot layers into a single mask
raster.masks.comb_MA = calc(raster.masks[[c(1,5:9)]],sum)
raster.masks.comb_MA.df = rasterToPoints(raster.masks.comb_MA) %>% as.data.frame
save(raster.masks.comb_MA, file='data/processed/raster.masks.comb_MA.RData')
save(raster.masks.comb_MA.df, file='data/processed/raster.masks.comb_MA.df.RData')
g=ggplot() +
    geom_tile(data=raster.masks.comb_MA.df, aes(y=y, x=x, fill=layer), show.legend=FALSE) +
    geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    theme_bw() +
    coord_equal() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(filename='output/masks_MA1.pdf', g, width=5, height=4)
ggsave(filename='output/masks_MA1.png', g, width=5, height=4, units='in', dpi=300)
