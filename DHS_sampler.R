source('DHS_functions.R')
DHS_checkPackages()
library(BalancedSampling)

source('DHS_config.R')

load(file='data/processed/raster.masks.comb.RData')
load(file='data/processed/hard_seabed_EA.sp.RData')
load(file='data/processed/shipping.sp.RData')
load(file='data/processed/port.sp.RData')
load(file='data/processed/boun.sp.RData')
load(file='data/processed/naval.sp.RData')
load(file='data/processed/hard_seabed_EA.sp.RData')

l = list.files(path='data/processed', pattern='east.arm.pred.rast')
nms = gsub('east.arm.pred.rast\\_(.*)\\..*','\\1',l)

rasters=stack()
for (i in l) {
    load(file=paste0('data/processed/',i))
    rasters= stack(rasters, east.arm.pred.rast)
}
names(rasters) <- nms

#Generate shapefile of 
combined.masks=rasterToPolygons(raster.masks.comb)
tmaptools::write_shape(combined.masks, file='data/processed/CombinedMask.shp')


## Apply the masks
rasters = mask(rasters, rasterToPolygons(raster.masks.comb))
plot(rasters)
rasters = mask(rasters, shipping.sp, inverse=TRUE)
rasters = mask(rasters, port.sp, inverse=TRUE)
rasters = mask(rasters, boun.sp, inverse=TRUE)
rasters = mask(rasters, naval.sp, inverse=TRUE)
rasters = mask(rasters, hard_seabed_EA.sp, inverse=TRUE)


## plot(east.arm.pred.rast)
## plot(raster.masks.comb, add=TRUE)
## plot(port.sp, add=TRUE)
## plot(shipping.sp, add=TRUE)
## plot(boun.sp, add=TRUE)
## plot(naval.sp, add=TRUE)

## plot(rasters)

s = rasterToPoints(rasters, spatial=TRUE)
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


## p1=rep(100/nrow(a), nrow(a))
## #a$Value=a$Cd
##                                         #X1 = vegan::decostand(cbind(a$x,a$y, a$Value),method='range',MARGIN=2)

## #X1 = vegan::decostand(cbind(a$Value),method='range',MARGIN=2)

## s1 = NULL
## for (i in 1:1) s1=c(s1,lpm2(prob=p1,x=X1))
## #X2 = data.frame(Longitude=a$x, Latitude=a$y, Value=a$Value)
## X2 = data.frame(Longitude=a$x, Latitude=a$y)
## X3 = X2[s1,]

## #Xa = vegdist(X1, method='euclidean')
## #Xb = as.matrix(Xa)[,1]
## #X4 = data.frame(Longitude=a$x, Latitude=a$y, Value=Xb)

## ggplot(X2) +
##                                         #geom_tile(aes(y=Latitude, x=Longitude, fill=Value)) +
##   #geom_tile(data=X4,aes(y=Latitude, x=Longitude, fill=Value)) +
##   #geom_contour(aes(y=Latitude, x=Longitude, z=Value)) +
##   scale_fill_gradientn(colors=tim.colors(12)) +
##   geom_point(data=X3, aes(y=Latitude, x=Longitude))

## plot(X1[,1],X1[,2]); # plot population
## points(X1[s1,1],X1[s1,2], pch=16,col='red'); # plot sample
## #s2=lpm2(prob=p1,x=X1)
## #points(X1[s2,1],X1[s2,2], pch=16,col='blue'); # plot sample

## a1 = a %>% dplyr::select(-x,-y)
## a2=vegan::decostand(a1,method='standardize',MARGIN=2)
## a3 = vegdist(a2, method='euclidean')
## a4 = as.matrix(a3)
## a4[lower.tri(a3)==FALSE] <- NA
## a5 = colMeans(a4, na.rm=TRUE)
## sss = a %>% dplyr::select(x,y) %>% mutate(fit=a5, fit=scales::rescale(fit, from=c(0,max(fit,na.rm=TRUE)),to=c(1,0)))
## p <- scales::rescale(sss$x-130, to=c(0,1)) #1-exp(-ss$x)
##                                         #standardise to get n samples
##                                         #p <- length(ss$x) * p / sum( p)
## sss$fit=1-p

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
save(clhs.data, file='data/processed/clhs.data.RData')
save(target, file='data/processed/target.RData')
save(rand.data, file='data/processed/rand.data.RData')
save(regular.data, file='data/processed/regular.data.RData')
#save(mbh.data, file='data/processed/mbh.data.RData')
save(spatiallybalanced.data, file='data/processed/spatiallybalanced.data.RData')
save(spatiallybalanced.data2, file='data/processed/spatiallybalanced.data2.RData')


load(file='data/processed/clhs.data.RData')
load(file='data/processed/target.RData')
load(file='data/processed/rand.data.RData')
load(file='data/processed/regular.data.RData')
#load(file='data/processed/mbh.data.RData')
load(file='data/processed/spatiallybalanced.data.RData')
load(file='data/processed/spatiallybalanced.data2.RData')
load(file='data/processed/East.arm.sp.RData')
load(file='data/processed/East.arm.df.RData')


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
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhs.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhs.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=clhs.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhsMax.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhsMax.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=clhs.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhsMin.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhsMin.png', g, width=9, height=9, units='in', dpi=300)


clhs.design_100 = clhs.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(clhs.design_100, path='output/clhs.design_100.csv')
clhs.design_50 = clhs.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(clhs.design_50, path='output/clhs.design_50.csv')

## random sample 
rand.data.all = do.call('rbind', rand.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = rand.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##rand.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_rand.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_rand.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=rand.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesignMax_rand.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesignMax_rand.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=rand.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesignMin_rand.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesignMin_rand.png', g, width=6, height=5, units='in', dpi=300)

rand.design_100 = rand.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(rand.design_100, path='output/rand.design_100.csv')
rand.design_50 = rand.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(rand.design_50, path='output/rand.design_50.csv')

## regular sample 
regular.data.all = do.call('rbind', regular.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = regular.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##regular.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##   group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regular.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_regular.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=regular.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regularMax.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_regularMax.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=regular.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regularMin.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_regularMin.png', g, width=6, height=5, units='in', dpi=300)

regular.design_100 = regular.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_100, path='output/regular.design_100.csv')
regular.design_50 = regular.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_50, path='output/regular.design_50.csv')

## Spatially balanced design
spatiallybalanced.data.all = do.call('rbind', spatiallybalanced.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = spatiallybalanced.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##regular.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##   group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_spatiallybalanced.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMax.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMax.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMin.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMin.png', g, width=6, height=5, units='in', dpi=300)

spatiallybalanced.design_100 = spatiallybalanced.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_100, path='output/spatiallybalanced.design_100.csv')
spatiallybalanced.design_50 = spatiallybalanced.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_50, path='output/spatiallybalanced.design_50.csv')


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
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_spatiallybalanced2.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2Max.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2Max.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2Min.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2Min.png', g, width=6, height=5, units='in', dpi=300)

spatiallybalanced2.design_100 = spatiallybalanced.data2.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_100, path='output/spatiallybalanced2.design_100.csv')
spatiallybalanced2.design_90 = spatiallybalanced.data2.best %>% filter(N==90) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_90, path='output/spatiallybalanced2.design_90.csv')
spatiallybalanced2.design_80 = spatiallybalanced.data2.best %>% filter(N==80) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_80, path='output/spatiallybalanced2.design_80.csv')
spatiallybalanced2.design_50 = spatiallybalanced.data2.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced2.design_50, path='output/spatiallybalanced2.design_50.csv')

g=ggplot() +
    geom_point(data=spatiallybalanced2.design_100, aes(y=Latitude, x=Longitude)) +
    geom_polygon(data=East.arm.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2_final.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2_final.png', g, width=6, height=5, units='in', dpi=300)




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

ggsave(filename='output/SamplingEffort.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffort.png', g2, width=10, height=5, units='in', dpi=300)

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

ggsave(filename='output/SamplingEffortMax.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMax.png', g2, width=10, height=5, units='in', dpi=300)

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

ggsave(filename='output/SamplingEffortMin.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMin.png', g2, width=10, height=5, units='in', dpi=300)



## Now for the Outer Harbour ==================*****************************++++++++++++++++++++===============================================
load(file='data/processed/raster.masks.comb_outer.RData')

l = list.files(path='data/processed', pattern='outer.pred.rast')
nms = gsub('outer.pred.rast\\_(.*)\\..*','\\1',l)
rasters=stack()
for (i in l) {
    load(file=paste0('data/processed/',i))
    rasters= stack(rasters, outer.pred.rast)
}
names(rasters) <- nms
rasters = mask(rasters, rasterToPolygons(raster.masks.comb))
plot(rasters)
rasters = mask(rasters, shipping.sp, inverse=TRUE)
rasters = mask(rasters, hard_seabed_OH.sp, inverse=TRUE)

s = rasterToPoints(rasters, spatial=TRUE)
A=NULL
clhs.mods = list()
clhs.data = list()
rand.data = list()
regular.data = list()
#mbh.data = list()
spatiallybalanced.data = list()
spatiallybalanced.data2 = list()

a = as.data.frame(s)
X1 = vegan::decostand(a,method='range',MARGIN=2) %>% as.matrix
X2 = vegan::decostand(cbind(a$x,a$y),method='range',MARGIN=2) %>% as.matrix

for (i in c(5,10,20,30,40,50,60,70,80,90,100,120,200,1000)) {
    print(i)
    temp.data = NULL
    clhs.mods[[i]] = list()
    temp.rand.data = NULL
    temp.regular.data = NULL
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
        ss = SpatialPixels(s)
        wch = spsample(ss, n=i, replace=FALSE, type='regular')
        dat = as.data.frame(wch) %>% dplyr::rename(x=x1, y=x2) %>% mutate(Rep=j)
        for (v in names(s)) {
            dat = dat %>% mutate(!!v := inlaPredict(v, wch, area='Outer'))
        }
        temp.regular.data = rbind(temp.regular.data, dat)

      ## Spatially balanced design
      set.seed(j)
        p1=rep(i/nrow(a), nrow(a))
        wch=lpm2(prob=p1,x=X1)
        dat = as.data.frame(a) %>% mutate(Rep=j)
        dat = dat[wch,]
        temp.spatiallybalanced.data = rbind(temp.spatiallybalanced.data, dat)

      ## Spatially balanced design
      set.seed(j)
        p1=rep(i/nrow(a), nrow(a))
        wch=lpm2(prob=p1,x=X2)
        dat = as.data.frame(a) %>% mutate(Rep=j)
        dat = dat[wch,]
        temp.spatiallybalanced.data2 = rbind(temp.spatiallybalanced.data2, dat)
        
        target=colMeans(s.clhs$initial_object@data)
   }
    clhs.data[[paste(i)]]=temp.data
    rand.data[[paste(i)]]=temp.rand.data
    regular.data[[paste(i)]]=temp.regular.data
    spatiallybalanced.data[[paste(i)]]=temp.spatiallybalanced.data
    spatiallybalanced.data2[[paste(i)]]=temp.spatiallybalanced.data2
}
save(clhs.data, file='data/processed/clhs.data_outer.RData')
save(target, file='data/processed/target_outer.RData')
save(rand.data, file='data/processed/rand.data_outer.RData')
save(regular.data, file='data/processed/regular.data_outer.RData')
#save(mbh.data, file='data/processed/mbh.data.RData')
save(spatiallybalanced.data, file='data/processed/spatiallybalanced.data_outer.RData')
save(spatiallybalanced.data2, file='data/processed/spatiallybalanced.data2_outer.RData')

load(file='data/processed/clhs.data_outer.RData')
load(file='data/processed/target_outer.RData')
load(file='data/processed/rand.data_outer.RData')
load(file='data/processed/regular.data_outer.RData')
#load(file='data/processed/mbh.data.RData')
load(file='data/processed/spatiallybalanced.data_outer.RData')
load(file='data/processed/spatiallybalanced.data2_outer.RData')
load(file='data/processed/Outer.df.RData')

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
            ErrorMin=apply(ss,1,function(x) min(x))) %>%
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
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhs_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhs_outer.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=clhs.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhsMax_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhsMax_outer.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=clhs.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_clhsMin_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_clhsMin_outer.png', g, width=9, height=9, units='in', dpi=300)


clhs.design_100 = clhs.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(clhs.design_100, path='output/clhs.design_100_outer.csv')
clhs.design_50 = clhs.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(clhs.design_50, path='output/clhs.design_50_outer.csv')


## random sample 
rand.data.all = do.call('rbind', rand.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = rand.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##rand.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_rand_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_rand_outer.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=rand.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesignMax_rand_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesignMax_rand_outer.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=rand.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesignMin_rand_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesignMin_rand_outer.png', g, width=6, height=5, units='in', dpi=300)

rand.design_100 = rand.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(rand.design_100, path='output/rand.design_100_outer.csv')
rand.design_50 = rand.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(rand.design_50, path='output/rand.design_50_outer.csv')

## regular sample 
regular.data.all = do.call('rbind', regular.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = regular.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##regular.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##   group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regular_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_regular_outer.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=regular.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regularMax_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_regularMax_outer.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=regular.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_regularMin_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_regularMin_outer.png', g, width=6, height=5, units='in', dpi=300)

regular.design_100 = regular.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_100, path='output/regular.design_100_outer.csv')
regular.design_50 = regular.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_50, path='output/regular.design_50_outer.csv')

## Spatially balanced design
spatiallybalanced.data.all = do.call('rbind', spatiallybalanced.data) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = spatiallybalanced.data.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##regular.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##   group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_spatiallybalanced_outer.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMax_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMax_outer.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMin_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMin_outer.png', g, width=6, height=5, units='in', dpi=300)

spatiallybalanced.design_100 = spatiallybalanced.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_100, path='output/spatiallybalanced.design_100_outer.csv')
spatiallybalanced.design_50 = spatiallybalanced.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_50, path='output/spatiallybalanced.design_50_outer.csv')


## Spatially balanced design 2
spatiallybalanced.data2.all = do.call('rbind', spatiallybalanced.data2) %>%
    mutate(N=gsub('([0-9]*)\\..*','\\1',rownames(.))) %>%
    mutate(N=factor(N,levels=unique(N)))

a = spatiallybalanced.data2.all %>% group_by(N,Rep) %>%
    summarize_all(mean) %>% as.data.frame
a.f = a %>% dplyr::select(N,Rep)
A=a %>% dplyr::select(-N,-Rep, -x, -y) %>% as.matrix
#ss=sweep(A,2,STAT=target, FUN=function(x,y) abs((x-y)/y))
##ss = t(apply(A,1, function(x) {abs(x-target)/S}))
##regular.dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##   group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))
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
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2_outer.pdf', g, width=9, height=9)
ggsave(filename='output/SamplingDesign_spatiallybalanced2_outer.png', g, width=9, height=9, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.bestMax, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMax2_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMax2_outer.png', g, width=6, height=5, units='in', dpi=300)
g=ggplot() +
    geom_point(data=spatiallybalanced.data2.bestMin, aes(y=y, x=x)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    facet_wrap(~N) +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalancedMin2_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalancedMin2_outer.png', g, width=6, height=5, units='in', dpi=300)

spatiallybalanced.design_100 = spatiallybalanced.data2.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_100, path='output/spatiallybalanced2.design_100_outer.csv')
spatiallybalanced.design_50 = spatiallybalanced.data2.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(spatiallybalanced.design_50, path='output/spatiallybalanced2.design_50_outer.csv')

g=ggplot() +
    geom_point(data=spatiallybalanced.design_100, aes(y=Latitude, x=Longitude)) +
    geom_polygon(data=Outer.df, aes(y=lat, x=long, group=group), fill=NA, color='black') +
    theme_bw() +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    coord_equal()
ggsave(filename='output/SamplingDesign_spatiallybalanced2_final_outer.pdf', g, width=6, height=5)
ggsave(filename='output/SamplingDesign_spatiallybalanced2_final_outer.png', g, width=6, height=5, units='in', dpi=300)


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
    scale_color_manual('Method', breaks=c('cLHS','Random','Regular','Spatially balanced', 'Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    scale_fill_manual('Method', breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffort_outer.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffort_outer.png', g2, width=10, height=5, units='in', dpi=300)

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
    scale_color_manual('Method', breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    scale_fill_manual('Method', breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffortMax_outer.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMax_outer.png', g2, width=10, height=5, units='in', dpi=300)

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
    scale_color_manual('Method', breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    scale_fill_manual('Method', breaks=c('cLHS','Random','Regular','Spatially balanced','Spatially balanced 2'), values=c('blue','green','red','purple','orange')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffortMin_outer.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMin_outer.png', g2, width=10, height=5, units='in', dpi=300)

