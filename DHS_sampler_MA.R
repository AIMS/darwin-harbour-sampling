source('DHS_functions.R')
DHS_checkPackages()
library(BalancedSampling)

source('DHS_config.R')

load(file='data/processed/raster.masks.comb_MA.RData')
load(file='data/processed/hard_seabed_MA.sp.RData')
load(file='data/processed/shipping_MA.sp.RData')
load(file='data/processed/exclusion_MA.sp.RData')

l = list.files(path='data/processed', pattern='middle_arm.pred.rast')
nms = gsub('middle_arm.pred.rast\\_(.*)\\..*','\\1',l)

rasters=stack()
for (i in l) {
    load(file=paste0('data/processed/',i))
    rasters= stack(rasters, middle_arm.pred.rast)
}
names(rasters) <- nms

#Generate shapefile of 
combined.masks_MA=rasterToPolygons(raster.masks.comb_MA)
temp_data = combined.masks_MA@data %>% summarize(OBJECT='Middle Arm energy mask')
combined.masks_MA  <- unionSpatialPolygons(combined.masks_MA ,ID=rep(1, length(combined.masks_MA$layer)))
combined.masks_MA = SpatialPolygonsDataFrame(combined.masks_MA, data=temp_data)
tmaptools::write_shape(combined.masks_MA, file='data/processed/CombinedMask_MA.shp')


## Apply the masks
rasters_MA = mask(rasters, rasterToPolygons(raster.masks.comb_MA))
plot(rasters_MA)
rasters_MA = mask(rasters_MA, shipping_MA.sp, inverse=TRUE)
rasters_MA = mask(rasters_MA, exclusion_MA.sp, inverse=TRUE)
rasters_MA = mask(rasters_MA, hard_seabed_MA.sp, inverse=TRUE)



s = rasterToPoints(rasters_MA, spatial=TRUE)
A=NULL
clhs.mods = list()
clhs.data = list()
rand.data = list()
regular.data = list()
#mbh.data = list()
spatiallybalanced.data = list()
spatiallybalanced.data2 = list()

## In preparation for spatially balanced design exploration, we need to generate:
## 1. the full spatial domain to sample from (this needs to be based on the mask
## 2. probabilities that a cell can be sampled.  Although technically this is to
##    facilitate the inclusion of legacy sites (inclusion probability of 1), I am
##    going to use some property of the sediment sampling data.
##    Ideally, what we want is more samples where the underlying sediment chemical
##    communities are most heterogeneous and fewer where they are most homogeneous.
##    One way to explore this heterogeneity is to explore a multivariate distance
##    metric (like euclidean distance).  If each of the chemical species are first
##    standardized (to a mean of zero and a standard deviation of 1), the euclidean
##    distance can provide a way to express the multidimensional dissimilarities between
##    all pairs of cells.  The average dissimilarity associated with a cell can then
##    be used to represent the degree to which this cell differs from other cells.
##    The larger the value, the more heterogeneous this cell is to others

a = as.data.frame(s)
X1 = vegan::decostand(a,method='range',MARGIN=2) %>% as.matrix
X2 = vegan::decostand(cbind(a$x,a$y),method='range',MARGIN=2) %>% as.matrix

for (i in c(5,10,20,30,40,50,60,70,80,90,100,120,200,1000)) {
    print(i)
    temp.data = NULL
    clhs.mods[[i]] = list()
    temp.rand.data = NULL
    temp.regular.data = NULL
    #temp.mbh.data = NULL
    temp.spatiallybalanced.data = NULL
    temp.spatiallybalanced.data2 = NULL
    for (j in 1:5) {
      set.seed(j)
        s.clhs <- clhs(s, size = i, progress = FALSE, iter = 10000, simple = FALSE)
    
        clhs.mods[[i]][[j]] = s.clhs
        temp.data = rbind(temp.data, cbind(Rep=j, as.data.frame(s.clhs$sampled_data)))

      ## purely random
      set.seed(j)
        wch = sample(nrow(s), size=i, replace=FALSE)
        temp.rand.data = rbind(temp.rand.data, cbind(Rep=j, as.data.frame(s[wch,])))

      ## regular grid
      set.seed(j)
        ## wch=sampleRegular(rasters, size=i, na.rm=TRUE, xy=TRUE)
        ## raster.masks.comb.df = rasterToPolygons(raster.masks.comb)
        ## wch = spsample(raster.masks.comb.df, size=i, replace=FALSE, type='regular')
        
        ss = SpatialPixels(s)
        wch = spsample(ss, n=i, replace=FALSE, type='regular')
        dat = as.data.frame(wch) %>% dplyr::rename(x=x1, y=x2) %>% mutate(Rep=j)
        for (v in names(s)) {
            dat = dat %>% mutate(!!v := inlaPredict(v, wch))
        }
        temp.regular.data = rbind(temp.regular.data, dat)

      ## Spatially balanced design
      set.seed(j)
        p1=rep(i/nrow(a), nrow(a))
        wch=lpm2(prob=p1,x=X1)
        dat = as.data.frame(a) %>% mutate(Rep=j)
        dat = dat[wch,]
        #for (v in names(s)) {
        #    dat = dat %>% mutate(!!v := inlaPredict(v, wch))
        #}
        temp.spatiallybalanced.data = rbind(temp.spatiallybalanced.data, dat)

      ## Spatially balanced design (no environmental covariates)
      set.seed(j)
        p1=rep(i/nrow(a), nrow(a))
        wch=lpm2(prob=p1,x=X2)
        dat = as.data.frame(a) %>% mutate(Rep=j)
        dat = dat[wch,]
        temp.spatiallybalanced.data2 = rbind(temp.spatiallybalanced.data2, dat)
        
        ## ##qq=quasiSamp(n=i, dimension=2, potential.sites=a %>% dplyr::select(x,y) %>% as.matrix, inclusion.probs=sss$fit)
        ## qq=quasiSamp(n=i, dimension=2, potential.sites=a %>% dplyr::select(x,y) %>% as.matrix)
        ## wch = SpatialPoints(coords = qq[,1:2])
        ## dat = qq %>% dplyr::select(x,y) %>% mutate(Rep=j)
        ## for (v in names(s)) {
        ##     dat = dat %>% mutate(!!v := inlaPredict(v, wch))
        ## }
        ## temp.mbh.data = rbind(temp.mbh.data, dat)
        #ggplot(ss) + geom_tile(aes(y=y, x=x, fill=(fit))) + geom_point(data=qq, aes(y=y, x=x), color='red')
        #clhs.data[[paste(i)]] = as.data.frame(s.clhs$sampled_data)
        target=colMeans(s.clhs$initial_object@data)
#        A=rbind(A,cbind(Rep=j, cbind(colMeans(s.clhs$sampled_data@data))))
    }
    clhs.data[[paste(i)]]=temp.data
    rand.data[[paste(i)]]=temp.rand.data
    regular.data[[paste(i)]]=temp.regular.data
    #mbh.data[[paste(i)]]=temp.mbh.data
    spatiallybalanced.data[[paste(i)]]=temp.spatiallybalanced.data
    spatiallybalanced.data2[[paste(i)]]=temp.spatiallybalanced.data2
}

save(clhs.data, file='data/processed/clhs.data_MA.RData')
save(target, file='data/processed/target_MA.RData')
save(rand.data, file='data/processed/rand.data_MA.RData')
save(regular.data, file='data/processed/regular.data_MA.RData')
#save(mbh.data, file='data/processed/mbh.data.RData')
save(spatiallybalanced.data, file='data/processed/spatiallybalanced.data_MA.RData')
save(spatiallybalanced.data2, file='data/processed/spatiallybalanced.data2_MA.RData')

load(file='data/processed/clhs.data_MA.RData')
load(file='data/processed/target_MA.RData')
load(file='data/processed/rand.data_MA.RData')
load(file='data/processed/regular.data_MA.RData')
load(file='data/processed/spatiallybalanced.data_MA.RData')
load(file='data/processed/spatiallybalanced.data2_MA.RData')
load(file='data/processed/middle_arm.sp.RData')
load(file='data/processed/middle_arm.df.RData')


S=as.vector(s %>% as.data.frame %>% dplyr::select(-x,-y) %>% summarize_all(function(x) diff(range(x))) %>% as.matrix)

## clhs model
clhs.data.all = do.call('rbind', clhs.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))


a = clhs.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
                                        #ss=sweep(A,2,STAT=target, FUN=function(x,y) (x-y)^2/y)
                                        #ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
ss = t(apply(A,1, function(x) {abs(x-target)/S}))
dat = cbind(Error = apply(ss,1, mean),
            ErrorMax=apply(ss,1,function(x) max(x)),
            ErrorMin=apply(ss,1,function(x) min(x))
            ) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0),
                           BestMin=ifelse(ErrorMin==min(ErrorMin),1,0))
## dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##     group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))

clhs.data.best = clhs.data.all %>% right_join(dat %>% filter(Best==1))
clhs.data.bestMax = clhs.data.all %>% right_join(dat %>% filter(BestMax==1))
clhs.data.bestMin = clhs.data.all %>% right_join(dat %>% filter(BestMin==1))
g=ggplot() +
    geom_point(data=clhs.data.best, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhs_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhs_MA.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=clhs.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhsMax_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhsMax_MA.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=clhs.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhsMin_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhsMin_MA.png', g, width=9, height=9, units='in', dpi=300)

clhs.design_100_MA = clhs.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(clhs.design_100_MA, path='output/clhs.design_100_MA.csv')
clhs.design_50_MA = clhs.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(clhs.design_50_MA, path='output/clhs.design_50_MA.csv')


## random sample 
rand.data.all = do.call('rbind', rand.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = rand.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
ss = t(apply(A,1, function(x) {abs(x-target)/S}))
rand.dat = cbind(Error = apply(ss,1, mean),
                 ErrorMax=apply(ss,1,function(x) max(x)),
                 ErrorMin=apply(ss,1,function(x) min(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0),
                           BestMin=ifelse(ErrorMin==min(ErrorMin),1,0))

rand.data.best = rand.data.all %>% right_join(rand.dat %>% filter(Best==1))
rand.data.bestMax = rand.data.all %>% right_join(rand.dat %>% filter(BestMax==1))
rand.data.bestMin = rand.data.all %>% right_join(rand.dat %>% filter(BestMin==1))
g=ggplot() +
    geom_point(data=rand.data.best, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_rand_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_rand_MA.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=rand.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesignMax_rand_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesignMax_rand_MA.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=rand.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesignMin_rand_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesignMin_rand_MA.png', g, width=6, height=5, units='in', dpi=300)

rand.design_100_MA = rand.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(rand.design_100_MA, path='output/rand.design_100_MA.csv')
rand.design_50_MA = rand.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(rand.design_50_MA, path='output/rand.design_50_MA.csv')


## regular sample 
regular.data.all = do.call('rbind', regular.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = regular.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
ss = t(apply(A,1, function(x) {abs(x-target)/S}))
regular.dat = cbind(Error = apply(ss,1, mean),
                    ErrorMax=apply(ss,1,function(x) max(x)),
                    ErrorMin=apply(ss,1,function(x) min(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0),
                           BestMin=ifelse(ErrorMin==min(ErrorMin),1,0))
regular.data.best = regular.data.all %>% right_join(regular.dat %>% filter(Best==1))
regular.data.bestMax = regular.data.all %>% right_join(regular.dat %>% filter(BestMax==1))
regular.data.bestMin = regular.data.all %>% right_join(regular.dat %>% filter(BestMin==1))
g=ggplot() +
    geom_point(data=regular.data.best, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regular_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_regular_MA.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=regular.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regularMax_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_regularMax_MA.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=regular.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regularMin_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_regularMin_MA.png', g, width=6, height=5, units='in', dpi=300)

regular.design_100_MA = regular.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_100_MA, path='output/regular.design_100_MA.csv')
regular.design_50_MA = regular.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_50_MA, path='output/regular.design_50_MA.csv')


## Spatially balanced design
spatiallybalanced.data.all = do.call('rbind', spatiallybalanced.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = spatiallybalanced.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
ss = t(apply(A,1, function(x) {abs(x-target)/S}))
spatiallybalanced.dat = cbind(Error = apply(ss,1, mean),
                              ErrorMax=apply(ss,1,function(x) max(x)),
                              ErrorMin=apply(ss,1,function(x) min(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0),
                           BestMin=ifelse(ErrorMin==min(ErrorMin),1,0))
spatiallybalanced.data.best = spatiallybalanced.data.all %>% right_join(spatiallybalanced.dat %>% filter(Best==1))
spatiallybalanced.data.bestMax = spatiallybalanced.data.all %>% right_join(spatiallybalanced.dat %>% filter(BestMax==1))
spatiallybalanced.data.bestMin = spatiallybalanced.data.all %>% right_join(spatiallybalanced.dat %>% filter(BestMin==1))
g=ggplot() +
    geom_point(data=spatiallybalanced.data.best, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_spatiallybalanced_MA.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMax_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMax_MA.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMin_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMin_MA.png', g, width=6, height=5, units='in', dpi=300)

spatiallybalanced.design_100_MA = spatiallybalanced.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_100_MA, path='output/spatiallybalanced.design_100_MA.csv')
spatiallybalanced.design_50_MA = spatiallybalanced.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_50_MA, path='output/spatiallybalanced.design_50_MA.csv')


## Spatially balanced design (no environmental covariates)
spatiallybalanced.data2.all = do.call('rbind', spatiallybalanced.data2) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = spatiallybalanced.data2.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
ss = t(apply(A,1, function(x) {abs(x-target)/S}))
spatiallybalanced.dat2 = cbind(Error = apply(ss,1, mean),
                               ErrorMax=apply(ss,1,function(x) max(x)),
                               ErrorMin=apply(ss,1,function(x) min(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0),
                           BestMin=ifelse(ErrorMin==min(ErrorMin),1,0))
spatiallybalanced.data2.best = spatiallybalanced.data2.all %>% right_join(spatiallybalanced.dat2 %>% filter(Best==1))
spatiallybalanced.data2.bestMax = spatiallybalanced.data2.all %>% right_join(spatiallybalanced.dat2 %>% filter(BestMax==1))
spatiallybalanced.data2.bestMin = spatiallybalanced.data2.all %>% right_join(spatiallybalanced.dat2 %>% filter(BestMin==1))
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.best, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2_MA.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_spatiallybalanced2_MA.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2Max_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2Max_MA.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2Min_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2Min_MA.png', g, width=6, height=5, units='in', dpi=300)

spatiallybalanced2.design_100_MA = spatiallybalanced.data2.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_100_MA, path='output/spatiallybalanced2.design_100_MA.csv')
spatiallybalanced2.design_90_MA = spatiallybalanced.data2.best %>% filter(N==90) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_90_MA, path='output/spatiallybalanced2.design_90_MA.csv')
spatiallybalanced2.design_80_MA = spatiallybalanced.data2.best %>% filter(N==80) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_80_MA, path='output/spatiallybalanced2.design_80_MA.csv')
spatiallybalanced2.design_50_MA = spatiallybalanced.data2.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_50_MA, path='output/spatiallybalanced2.design_50_MA.csv')

g=ggplot() +
    geom_point(data=spatiallybalanced2.design_100, aes(y=Latitude, x=Longitude)) +
    geom_polygon(data=middle_arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2_final_MA.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2_final_MA.png', g, width=6, height=5, units='in', dpi=300)



## Based on mean error
g1=ggplot() +
    geom_smooth(data=dat, aes(y=Error, x=as.numeric(as.character(N))), color='blue', fill='blue', alpha=0.3,method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=dat, aes(y=Error, x=as.numeric(as.character(N))), color='blue', linetype='dashed', stat='summary') +
    geom_point(data=dat, aes(y=Error, x=as.numeric(as.character(N))), color='blue') +
    geom_smooth(data=rand.dat, aes(y=Error, x=as.numeric(as.character(N))), color='green', fill='green', alpha=0.3, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=rand.dat, aes(y=Error, x=as.numeric(as.character(N))), color='green', linetype='dashed', stat='summary') +
    geom_point(data=rand.dat, aes(y=Error, x=as.numeric(as.character(N))), color='green') +
    geom_smooth(data=regular.dat, aes(y=Error, x=as.numeric(as.character(N))), color='red', fill='red', alpha=0.3, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=regular.dat, aes(y=Error, x=as.numeric(as.character(N))), color='red', linetype='dashed', stat='summary') +
    geom_point(data=regular.dat, aes(y=Error, x=as.numeric(as.character(N))), color='red') +    
    geom_smooth(data=spatiallybalanced.dat, aes(y=Error, x=as.numeric(as.character(N))), color='purple', fill='orange', alpha=0.3, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat, aes(y=Error, x=as.numeric(as.character(N))), color='purple', linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat, aes(y=Error, x=as.numeric(as.character(N))), color='purple') +    
    geom_smooth(data=spatiallybalanced.dat2, aes(y=Error, x=as.numeric(as.character(N))), color='orange', fill='orange', alpha=0.3, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat2, aes(y=Error, x=as.numeric(as.character(N))), color='orange', linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat2, aes(y=Error, x=as.numeric(as.character(N))), color='orange') +    
    scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Mean Error') +
    #coord_trans(y='log2') +
    theme_bw() 

g2=ggplot() +
    geom_smooth(data=dat, aes(y=Error, x=as.numeric(as.character(N)), color='cLHS', fill='cLHS'), alpha=0.1,method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=dat, aes(y=Error, x=as.numeric(as.character(N)), color='cLHS'),linetype='dashed', stat='summary') +
    geom_point(data=dat, aes(y=Error, x=as.numeric(as.character(N)),color='cLHS')) +
    geom_smooth(data=rand.dat, aes(y=Error, x=as.numeric(as.character(N)), color='Random', fill='Random'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=rand.dat, aes(y=Error, x=as.numeric(as.character(N)), color='Random'), linetype='dashed', stat='summary') +
    geom_point(data=rand.dat, aes(y=Error, x=as.numeric(as.character(N)), color='Random')) +
    geom_smooth(data=regular.dat, aes(y=Error, x=as.numeric(as.character(N)), color='Regular', fill='Regular'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=regular.dat, aes(y=Error, x=as.numeric(as.character(N)), color='Regular'), linetype='dashed', stat='summary') +
    geom_point(data=regular.dat, aes(y=Error, x=as.numeric(as.character(N)), color='Regular')) +
    geom_smooth(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced', fill='Spatially balanced'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced'), linetype='dashed', stat='summary') +
  geom_point(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced')) +
  geom_smooth(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2', fill='Spatially balanced 2'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2'), linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2')) +
        scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Mean Error', breaks=scales::log_breaks(base=2)) +
    scale_color_manual('Method', labels=c('cLHS','Random','Regular','2D Spatially balanced','nD Spatially balanced'), breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    scale_fill_manual('Method', labels=c('cLHS','Random','Regular','2D Spatially balanced','nD Spatially balanced'), breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffort_MA.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffort_MA.png', g2, width=10, height=5, units='in', dpi=300)


# Based on max error
g1=ggplot() +
    geom_smooth(data=dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='blue', fill='blue', alpha=0.3,method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='blue', linetype='dashed', stat='summary') +
    geom_point(data=dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='blue') +
    geom_smooth(data=rand.dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='green', fill='green', alpha=0.3, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=rand.dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='green', linetype='dashed', stat='summary') +
    geom_point(data=rand.dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='green') +
    geom_smooth(data=regular.dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='red', fill='red', alpha=0.3, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=regular.dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='red', linetype='dashed', stat='summary') +
    geom_point(data=regular.dat, aes(y=ErrorMax, x=as.numeric(as.character(N))), color='red') +
    geom_smooth(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced', fill='Spatially balanced'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced'), linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced')) +
    geom_smooth(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2', fill='Spatially balanced 2'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2'), linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2')) +
    scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Maximum Error') +
    #coord_trans(y='log2') +
    theme_bw() 

g2=ggplot() +
    geom_smooth(data=dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='cLHS', fill='cLHS'), alpha=0.1,method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='cLHS'),linetype='dashed', stat='summary') +
    geom_point(data=dat, aes(y=ErrorMax, x=as.numeric(as.character(N)),color='cLHS')) +
    geom_smooth(data=rand.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Random', fill='Random'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=rand.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Random'), linetype='dashed', stat='summary') +
    geom_point(data=rand.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Random')) +
    geom_smooth(data=regular.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Regular', fill='Regular'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=regular.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Regular'), linetype='dashed', stat='summary') +
    geom_point(data=regular.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Regular')) +
    geom_smooth(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced', fill='Spatially balanced'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced'), linetype='dashed', stat='summary') +
  geom_point(data=spatiallybalanced.dat, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced')) +
  geom_smooth(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2', fill='Spatially balanced 2'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2'), linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat2, aes(y=ErrorMax, x=as.numeric(as.character(N)), color='Spatially balanced 2')) +
    scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Maximum Error', breaks=scales::log_breaks(base=2)) +
    scale_color_manual('Method', labels=c('cLHS','Random','Regular','2D Spatially balanced','nD Spatially balanced'), breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    scale_fill_manual('Method', labels=c('cLHS','Random','Regular','2D Spatially balanced','nD Spatially balanced'), breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffortMax_MA.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMax_MA.png', g2, width=10, height=5, units='in', dpi=300)


g2=ggplot() +
    geom_smooth(data=dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='cLHS', fill='cLHS'), alpha=0.1,method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='cLHS'),linetype='dashed', stat='summary') +
    geom_point(data=dat, aes(y=ErrorMin, x=as.numeric(as.character(N)),color='cLHS')) +
    geom_smooth(data=rand.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Random', fill='Random'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=rand.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Random'), linetype='dashed', stat='summary') +
    geom_point(data=rand.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Random')) +
    geom_smooth(data=regular.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Regular', fill='Regular'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=regular.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Regular'), linetype='dashed', stat='summary') +
    geom_point(data=regular.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Regular')) +
    geom_smooth(data=spatiallybalanced.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Spatially balanced', fill='Spatially balanced'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Spatially balanced'), linetype='dashed', stat='summary') +
  geom_point(data=spatiallybalanced.dat, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Spatially balanced')) +
  geom_smooth(data=spatiallybalanced.dat2, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Spatially balanced 2', fill='Spatially balanced 2'), alpha=0.1, method='gam', formula=y~s(x, bs='ps'), method.args=list(family=Gamma(link='log'))) +
    geom_line(data=spatiallybalanced.dat2, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Spatially balanced 2'), linetype='dashed', stat='summary') +
    geom_point(data=spatiallybalanced.dat2, aes(y=ErrorMin, x=as.numeric(as.character(N)), color='Spatially balanced 2')) +
    scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Minimum Error', breaks=scales::log_breaks(base=2)) +
    scale_color_manual('Method', labels=c('cLHS','Random','Regular','2D Spatially balanced','nD Spatially balanced'), breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    scale_fill_manual('Method', labels=c('cLHS','Random','Regular','2D Spatially balanced','nD Spatially balanced'), breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffortMin_MA.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMin_MA.png', g2, width=10, height=5, units='in', dpi=300)
