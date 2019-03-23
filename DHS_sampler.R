source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/raster.masks.comb.RData')

l = list.files(path='data/processed', pattern='east.arm.pred.rast')
nms = gsub('east.arm.pred.rast\\_(.*)\\..*','\\1',l)

rasters=stack()
for (i in l) {
    load(file=paste0('data/processed/',i))
    rasters= stack(rasters, east.arm.pred.rast)
}
names(rasters) <- nms

rasters = mask(rasters, rasterToPolygons(raster.masks.comb))
plot(rasters)



s = rasterToPoints(rasters, spatial=TRUE)
A=NULL
clhs.mods = list()
clhs.data = list()
rand.data = list()
regular.data = list()
for (i in c(5,10,20,30,40,50,60,70,80,90,100,120,200,1000)) {
    print(i)
    temp.data = NULL
    clhs.mods[[i]] = list()
    temp.rand.data = NULL
    temp.regular.data = NULL
    for (j in 1:5) {
        s.clhs <- clhs(s, size = i, progress = FALSE, iter = 10000, simple = FALSE)
    
        clhs.mods[[i]][[j]] = s.clhs
        temp.data = rbind(temp.data, cbind(Rep=j, as.data.frame(s.clhs$sampled_data)))

        ## purely random
        wch = sample(nrow(s), size=i, replace=FALSE)
        temp.rand.data = rbind(temp.rand.data, cbind(Rep=j, as.data.frame(s[wch,])))

        ## regular grid
        
        
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
        #clhs.data[[paste(i)]] = as.data.frame(s.clhs$sampled_data)
        target=colMeans(s.clhs$initial_object@data)
#        A=rbind(A,cbind(Rep=j, cbind(colMeans(s.clhs$sampled_data@data))))
    }
    clhs.data[[paste(i)]]=temp.data
    rand.data[[paste(i)]]=temp.rand.data
    regular.data[[paste(i)]]=temp.regular.data
}
save(clhs.data, file='data/processed/clhs.data.RData')
save(rand.data, file='data/processed/rand.data.RData')
save(regular.data, file='data/processed/regular.data.RData')


load(file='data/processed/clhs.data.RData')
load(file='data/processed/rand.data.RData')
load(file='data/processed/regular.data.RData')

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
dat = cbind(Error = apply(ss,1, mean), ErrorMax=apply(ss,1,function(x) max(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0))
## dat = data.frame(Error=rowMeans(ss)) %>% cbind(a.f) %>%
##     group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0))

clhs.data.best = clhs.data.all %>% right_join(dat %>% filter(Best==1))
clhs.data.bestMax = clhs.data.all %>% right_join(dat %>% filter(BestMax==1))
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
rand.dat = cbind(Error = apply(ss,1, mean), ErrorMax=apply(ss,1,function(x) max(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0))

rand.data.best = rand.data.all %>% right_join(rand.dat %>% filter(Best==1))
rand.data.bestMax = rand.data.all %>% right_join(rand.dat %>% filter(BestMax==1))
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
regular.dat = cbind(Error = apply(ss,1, mean), ErrorMax=apply(ss,1,function(x) max(x))) %>%
    as.data.frame %>% cbind(a.f) %>%
    group_by(N) %>% mutate(Best=ifelse(Error==min(Error),1,0),
                           BestMax=ifelse(ErrorMax==min(ErrorMax),1,0))
regular.data.best = regular.data.all %>% right_join(regular.dat %>% filter(Best==1))
regular.data.bestMax = regular.data.all %>% right_join(regular.dat %>% filter(BestMax==1))
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

regular.design_100 = regular.data.best %>% filter(N==100) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_100, path='output/regular.design_100.csv')
regular.design_50 = regular.data.best %>% filter(N==50) %>%
    dplyr::select(Longitude=x, Latitude=y)
write_csv(regular.design_50, path='output/regular.design_50.csv')

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
    scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Mean Error', breaks=scales::log_breaks(base=2)) +
    scale_color_manual('Method', breaks=c('cLHS','Random','Regular'), values=c('blue','green','red')) +
    scale_fill_manual('Method', breaks=c('cLHS','Random','Regular'), values=c('blue','green','red')) +
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
    scale_x_log10('Number of samples', breaks=as.numeric(as.character(unique(dat$N))))+
    scale_y_continuous('Maximum Error', breaks=scales::log_breaks(base=2)) +
    scale_color_manual('Method', breaks=c('cLHS','Random','Regular'), values=c('blue','green','red')) +
    scale_fill_manual('Method', breaks=c('cLHS','Random','Regular'), values=c('blue','green','red')) +
    coord_trans(y='log2') +
    theme_bw()  + theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))

ggsave(filename='output/SamplingEffortMax.pdf', g2, width=10, height=5)
ggsave(filename='output/SamplingEffortMax.png', g2, width=10, height=5, units='in', dpi=300)


