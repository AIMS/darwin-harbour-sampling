source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

## Read in the various data sources.
## - Darwin Harbour shapefiles
## - Munksgaard sediment chemical data
## - hydrodynamic model data (current velocity and tidal current)
## - waves model data 
source('DHS_readData.R')



## Process
## 1. isolate the sediment data for East Arm
## 2. use a barrier model to create a 2d surface over East Arm from the sediment data (that is a separate layer for each Chemical)
##    This will be used as a multivariate layer to inform conditioned latent hypercube sampling
## 3. isolate the hydro and wave data for East Arm
## 4. use the hydro/wave layers as a 'cost' layer for clhs - alternatively use a threshold to turn this into a binary for masking
## 5. run clhs for the Mg (and others) with 'cost' included to yield 50 samples.

## 1. Isolate the sediment data for East Arm
source('DHS_moreSpatialObjects.R')

## 2. This needs a recent version - I am currently doing this on a virtualbox
source('DHS_fitINLAmodels.R')
source('DHS_summariseINLAmodels.R')

## 3. Isolate hydro and wave data for East Arm
## This was performed in DHS_moreSpatialObjects.R

## 4. Threshold hydro and waves data
source('DHS_thresholdLayers.R')

## 5. Run clhs
source('DHS_sampler.R')




## bedShear_50p = subset(hydro_east.sp, 'bedShear_50p')
## median(values(bedShear_50p), na.rm=TRUE)
## threshold=0.2 #median(values(bedShear_50p), na.rm=TRUE) #0.2
## reclassification.matrix=rbind(c(0,threshold,1),
##                               c(threshold,Inf,NA))
## b = reclassify(bedShear_50p, rcl=reclassification.matrix)
## plot(b)
## points(Latitude~Longitude, data=sediment.data)




## east.arm.pred.stack = stack(east.arm.pred.rast)

## library(clhs)
## s = rasterToPoints(east.arm.pred.stack, spatial=TRUE)
## A=list()
## for (i in c(5,10,20,50,100,200,1000,10000)) {
##   print(i)
##   s.clhs <- clhs(s, size = i, progress = FALSE, iter = 10000, simple = FALSE)
  
##   target=mean(s.clhs$initial_object@data$z)
##   A[i]=mean(s.clhs$sampled_data@data$z)
## }
## AA=do.call('c',A)
## AA-target

## plot(density(s.clhs$initial_object@data$z))
## plot(density(s.clhs$sampled_data@data$z), new=TRUE)

## plot(s.clhs, mode=c("obj", 'hist'))
## subset.idx <- s.clhs$index_samples
## par(mar = c(1,1,1,1))
## plot(newdata.rast, axes=FALSE)
## points(newdata[subset.idx, ], bg = 'red', pch=21)



## ## Use one of the hydro layers to mask..
## plot(hydro.outer_east.sp)
## a=subset(hydro.outer_east.sp, 'velocity_75p')
## a1=subset(DH.zone.sp, Name=='East Arm')
## plot(a1)

## a2=crop(mask(a, a1), a1)



## aa = coordinates(a2)
## aa = as.data.frame(aa) %>% mutate(z=values(a2)) %>% filter(!is.na(z)) %>%
##     mutate(z=1-(z-min(z))/max(z))
## plot(y~x, aa)

## aaa = sample(nrow(aa) ,size=1000, prob=aa$z)
## points(y~x, aa[aaa,], col='red')

## hist(a2)
## median(values(a), na.rm=TRUE)
## threshold=0.3
## reclassification.matrix=rbind(c(0,threshold,1),
##                               c(threshold,Inf,NA))
## b = reclassify(a2, rcl=reclassification.matrix)
## plot(b)
## points(Latitude~Longitude, data=sediment.data)
## plot(a1, add=TRUE)
## d=rasterToPolygons(b, dissolve=TRUE)

## points(spsample(d, 50, type='random'), col='red', pch=16)

## values(a) %>% head

## ## Try the library(clhs)
## a1=subset(DH.zone.sp, Name=='East Arm')
## a2=crop(mask(hydro.outer_east.sp, a1), a1)

## library(clhs)
## s = rasterToPoints(a2, spatial=TRUE)
## s.clhs <- clhs(s, size = 100, progress = FALSE, iter = 1000, simple = FALSE)
## plot(s.clhs, mode=c("obj", "box"))
## subset.idx <- s.clhs$index_samples
## par(mar = c(1,1,1,1))
## plot(a2, axes=FALSE)
## points(s[subset.idx, ], bg = 'red', pch=21)


## reclassification.matrix=rbind(c(0,threshold,0),
##                               c(threshold,Inf,1000))
## b = reclassify(a2, rcl=reclassification.matrix)
## bb = rasterToPoints(b, spatial=TRUE)
## s = rasterToPoints(a2, spatial=TRUE)
## #s = sampleRegular(a2, size=10000, sp=TRUE)
## s@data$velocity_75p = bb@data[,4]#s@data$velocity_75p^10
## s@data = s@data[,c(2,4)]
## s = s[which(!is.na(s$bedShear_75p)),]
## s.clhs <- clhs(s, size = 100, progress = FALSE, iter = 1000, simple = FALSE, cost='velocity_75p')
## plot(s.clhs, mode=c("obj", "box"))
## subset.idx <- s.clhs$index_samples
## a3=subset(hydro.outer_east.sp, 'velocity_75p')
## a4=crop(mask(a3, a1), a1)
## #a4=crop(mask(b, a1), a1)

## par(mar = c(1,1,1,1))
## plot(a4, axes=FALSE)
## points(s[subset.idx, ], bg = 'red', pch=21)


## ## INLA attempt 1 ================================================================================

## ## INLA barrier model to predict the spatial patterns throughout a grid
## max.edge=0.02
## bound.outer=0.01
## mesh = inla.mesh.2d(boundary = DH.wh.sp,
##                     loc=cbind(sediment.data$Longitude, sediment.data$Latitude),
##                      max.edge = c(1,2)*max.edge,
##                     cutoff = 0.005,
##                     offset = c(max.edge, bound.outer))
## plot(mesh, main="Our mesh", lwd=0.5)
## points(sediment.data$Longitude,sediment.data$Latitude, col="red")
## lines(DH.wh.sp, col='black',cex=3)

## A.i.s = inla.spde.make.A(mesh, loc=cbind(sediment.data$Longitude, sediment.data$Latitude))
## stk = inla.stack(data=list(y=sediment.data$Mg), 
##                     effects=list(s=1:mesh$n,
##                                  m = rep(1, nrow(sediment.data))),
##                    A=list(A.i.s, 1),
##                     remove.unused = FALSE, tag='est')
		    
## ## The stationary model
## prior.range = c(1, .5)
## prior.sigma = c(3, 0.01)
## spde = inla.spde2.pcmatern(mesh, prior.range=prior.range, prior.sigma=prior.sigma)
## hyper.iid = list(prec = list(prior='pc.prec', param=prior.sigma))

## mesh = dt.mesh.addon.posTri(mesh)
## ## - compute the triangle positions
## posTri = SpatialPoints(mesh$posTri)
## proj4string(posTri) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
## normal = over(DH.wh.sp, posTri, returnList=T)
## # - checking which mesh triangles are inside the normal area
## normal = unlist(normal)
## Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
## Omega.SP = dt.polygon.omega(mesh, Omega)
## plot(Omega.SP[[2]], col="grey", main="The barrier region (in grey)")

## Q.barrier = dt.create.Q(mesh, Omega, 
##                         fixed.ranges = c(NA, 0.5))
## # - We fix the barrier range to a different value than we 
## #   used for simulations
## # - - Why? It does not matter, as long as it is 'small' 
## #     the models are very
## #     similar
## # - - This shows that you do not need to know the 
## #     true 'barrier range'!
## # - time: Ca 1 min

## log.prior = dt.create.prior.log.exp(
##   prior.param = c(-log(prior.sigma[2])/prior.sigma[1], -log(prior.range[2])/prior.range[1])) 
##     #c(-log(0.01)/3, -log(0.5)*6))
## # - The prior parameters are the lambdas in the exponential 
## #   priors for standard deviation and inverse-range
## # - the first is log(prob)/exceed, the second log(prob)*exceed
## # - the second is exponential for inverse range, therefore multiplication!

## barrier.model = dt.inla.model(
##   Q = Q.barrier, log.prior=log.prior)

## mod = list(shortname="barrier-model")
## mod$formula = y~ -1+m + f(s, model=barrier.model)

## stack=stk

## ## Running all the models
## ## Initial values
## # - speeds up computations
## # - improves accuracy of computations
## # - set these to NULL the first time you run a model
## mod$init = c(NULL, NULL, NULL)

## mod$res = inla(mod$formula,
##                   data=inla.stack.data(stack),
##                   control.predictor=list(A=inla.stack.A(stack), link=1, compute=TRUE),
##                   family="gamma", 
##                   control.family = list(hyper = hyper.iid, link='log'),
##                   control.inla= list(int.strategy = "eb"),
##                   control.mode=list(restart=T, theta=mod$init))
## summary(mod$res)
## print(paste(round(mod$res$internal.summary.hyperpar$mode, 3), collapse = ','))

## field = mod$res$summary.random$s$mean + mod$res$summary.fixed['m', 'mean']
## xlim = DH.wh.sp@bbox[1, ] 
## ylim = DH.wh.sp@bbox[2, ]
## proj = inla.mesh.projector(mesh, xlim = xlim, 
##                            ylim = ylim, dims=c(300, 300))
## field.proj = inla.mesh.project(proj, field)
## zlim = range(field.proj, na.rm=TRUE)
## image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
##              xlim = xlim, ylim = ylim, asp=1) 
## contour(x = proj$x, y=proj$y,
##         z = field.proj,
##         levels=seq(zlim[1], zlim[2],length.out = 10),
##         add=TRUE, drawlabels=F, col="white")
## plot(Omega.SP[[2]], add=T, border="black", col="white")
## points(sediment.data$Longitude, sediment.data$Latitude)




## ## INLA attempt 2 ================================================================================
## ## mesh = inla.mesh.2d(boundary = DH.wh.sp,
## ##                      loc=cbind(sediment.data$Longitude, sediment.data$Latitude),
## ##                     max.edge = c(1,2)*max.edge,
## ##                     cutoff = 0.005,
## ##                     offset = c(max.edge, bound.outer))
## ## plot(mesh)
## ## points(sediment.data$Longitude,sediment.data$Latitude, col="red")
## ## lines(DH.wh.sp, col='black',cex=3)

## ## ##make the stack
## ## phi = inla.spde.make.A(mesh=mesh,
## ##                        loc=cbind(sediment.data$Longitude, sediment.data$Latitude))
## ## stack_est = inla.stack(data=list(y=sediment.data$Mg),
## ##                        A=list(phi,1,1),
## ##                        effects = list(s=1:mesh$n,
## ##                                       data.frame(Intercept=1,
## ##                                                  lon = sediment.data$Longitude,
## ##                                                  lat = sediment.data$Latitude),
## ##                                       iidx=1:nrow(sediment.data)),
## ##                        remove.unused=FALSE,tag='est')
## ## spde = inla.spde2.pcmatern(mesh, prior.range = c(6, .5), prior.sigma = c(3, 0.01))
## ## hyper.iid = list(prec = list(prior = 'pc.prec', param = c(3, 0.01)))

## ## ##Barrier model
## ## tl =length(mesh$graph$tv[,1]) ##number of triangles in the mesh
## ## posTri = matrix(0, tl, 2)
## ## for (t in 1:tl) {
## ##     temp=mesh$loc[mesh$graph$tv[t,],]
## ##     posTri[t,] = colMeans(temp)[c(1,2)]
## ## }
## ## posTri = SpatialPoints(posTri) ##compute the triangle positions
## ## proj4string(posTri) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
## ## normal = over(DH.wh.sp, posTri, returnList = TRUE) ## check which mesh triangles are inside the normal area
## ## normal = unlist(normal)
## ## barrier.triangles=setdiff(1:tl, normal)
## ## poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)
## ## plot(poly.barrier, col='grey')

## ## barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles,
## ##                                       prior.range=c(6,0.5),
## ##                                       prior.sigma=c(3,0.01))
## ## formula = y~ -1+Intercept + lon + lat +
## ##     f(s, model=barrier.model) +
## ##     f(iidx, model="iid", hyper=hyper.iid)

## ## output = inla(formula,
## ##               data=inla.stack.data(stack_est),
## ##               control.predictor=list(A=inla.stack.A(stack_est), compute=TRUE),
## ##               verbose=TRUE
## ## #              control.inla = list(int.strategy = 'eb')
## ## )

















## ## ## Get the spatial knots
## ## dat.st=sediment.data %>% dplyr::select(Date,Latitude,Longitude,Mg) %>% distinct %>% #filter(Measure=='Mg', !is.na(Value)) %>%
## ##     mutate(Dt.num=decimal_date(Date), Dt = as.numeric(as.Date(Date) - min(as.Date(Date))))
## ## a = point.in.polygon(point.x=dat.st$Longitude, point.y=dat.st$Latitude, pol.x=DH.wh.df$long, pol.y=DH.wh.df$lat)
## ## dat.st = dat.st[a==1,]
## ## nrow(dat.st)

## ## ggplot(dat.st, aes(y=Latitude, x=Longitude)) +
## ##     geom_point(aes(color=Mg)) + scale_color_gradientn(colors=terrain.colors(10))


## ## getKnots = function(seed=5, water, n=25) {
## ##     water.coords=water@polygons[[1]]@Polygons[[1]]@coords
## ##     bound = list(list(Longitude=water.coords[,1], Latitude=water.coords[,2]))
## ##     knots.sp = spsample(water, n=n, type='nonaligned')
## ##     knots=as.data.frame(knots.sp@coords)
## ##     names(knots) <- c('Longitude','Latitude')
## ##     list(knots=knots, knots.sp=knots.sp, bound=bound)
## ## }

## ## knots = getKnots(seed=5, water=DH.wh.sp, n=50)
## ## g1=ggplot(DH.wh.df) +
## ##     geom_polygon(aes(y=lat, x=long,group=group))+
## ##     geom_point(data=dat.st, aes(y=Latitude, x=Longitude), colour='black', shape='+')+
## ##     geom_point(data=as.data.frame(knots$knots.sp), aes(y=x2, x=x1), color='red')+
## ##     coord_map()+
## ##     theme_bw()+theme(axis.title=element_blank())
## ## g1

## ## kn=(data.frame(Longitude=dat.st$Longitude,Latitude=dat.st$Latitude) %>% distinct)#[-1,]
## ## kn=kn[c(1:15,24:26,28:37,40,42:43,46:52,54:57,60:67,71:77, 79:83, 86:88,90:102,104:106, 108:110,112,114,
## ##         117:121,123:128, 130, 133, 136,138,140:153, 155:160, 163:164, 166:167, 170:173, 177:184,
## ##         186:188, 190:192, 194, 196),]
## ## kn=kn[sample(1:nrow(kn), size=20),]
## ## ggplot(dat.st, aes(y=Latitude, x=Longitude)) + geom_point(aes(color=Mg)) + scale_color_gradientn(colors=terrain.colors(10)) +
## ##     geom_point(data=kn, shape=21)

## ## bound=knots$bound

## ## kn =spsample(DH.wh.sp, n=100, type='regular') %>% as.data.frame %>% rename(Latitude=x2, Longitude=x1)

## ## mod = gam(Mg~
## ##               te(Longitude, Latitude, bs=c('sf'), k=c(nrow(kn)), d=c(2),
## ##                  xt=list(bnd=bound)
## ##                  ),
## ##               family=Gamma(link='log'),
## ##           data=dat.st,
## ##           knots=c(kn),
## ##           method='REML')

## ## #plot(mod, scheme=1)
## ## plot(mod, scheme=2)



## ## #kn=kn[c(-16,-18,-19),]

## ## ## plot(bound[[1]][['Latitude']] ~ bound[[1]][['Longitude']], type='l')
## ## ## points(Latitude~Longitude, kn)
## ## ## library(mgcv)
## ## ## mod = gam(Mg~
## ## ##               te(Longitude, Latitude, Dt.num, bs=c('sf','cr'), k=c(nrow(kn), k), d=c(2,1),
## ## ##                  xt=list(list(bnd=bound), NULL)
## ## ##                  ) +
## ## ##               s(Dt.num, bs='gp', k=k),
## ## ##           family=Gamma(link='log'),
## ## ##           data=dat.st,
## ## ##           knots=c(kn, k),
## ## ##           method='REML')

## ## ## mod = gam(Mg~
## ## ##               te(Longitude, Latitude, bs=c('sf'), k=c(nrow(kn)), d=c(2),
## ## ##                  xt=list(bnd=bound)
## ## ##                  ),
## ## ##               family=Gamma(link='log'),
## ## ##           data=dat.st,
## ## ##           knots=c(kn),
## ## ##           method='REML')

## ## p =spsample(DH.wh.sp, n=1000, type='regular')
## ## plot(p)
## ## pp = as.data.frame(p) %>% dplyr::rename(Latitude=x2, Longitude=x1)

## ## pp = pp %>% mutate(fit=as.vector(predict(mod, newdata=pp, type='response')))

## ## ggplot(pp, aes(y=Latitude, x=Longitude)) + geom_point(aes(color=fit)) + scale_color_gradientn(colors=terrain.colors(10)) +
## ##     geom_point(data=kn)




## ## library(INLA)
## ## INLA:::inla.dynload.workaround()
## ## max.edge=0.02
## ## bound.outer=0.01
## ## mesh = inla.mesh.2d(boundary = DH.wh.sp,
## ##                     loc=cbind(sediment.data$Longitude, sediment.data$Latitude),
## ##                      max.edge = c(1,2)*max.edge,
## ##                     cutoff = 0.005,
## ##                     offset = c(max.edge, bound.outer))
## ## plot(mesh, main="Our mesh", lwd=0.5)
## ## points(sediment.data$Longitude,sediment.data$Latitude, col="red")
## ## lines(DH.wh.sp, col='black',cex=3)

## ## A.i.s = inla.spde.make.A(mesh, loc=cbind(sediment.data$Longitude, sediment.data$Latitude))
## ## stk = inla.stack(data=list(y=sediment.data$Mg), 
## ##                     effects=list(s=1:mesh$n,
## ##                                  m = rep(1, nrow(sediment.data))),
## ##                    A=list(A.i.s, 1),
## ##                     remove.unused = FALSE, tag='est') 
## ## ## Stationary model
## ## prior.range = c(1, .5)
## ## prior.sigma = c(3, 0.01)
## ## spde = inla.spde2.pcmatern(mesh, prior.range=prior.range, prior.sigma=prior.sigma)

## ## M = list()
## ## M[[1]] = list(shortname="stationary-model")
## ## M[[1]]$formula = y~ -1+m + f(s, model=spde)

## ## ## Barrier models
## ## source('https://haakonbakka.bitbucket.io/functions-barriers-dt-models-march2017.R')
## ## mesh = dt.mesh.addon.posTri(mesh)
## ## ## - compute the triangle positions
## ## posTri = SpatialPoints(mesh$posTri)
## ## proj4string(posTri) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
## ## normal = over(DH.wh.sp, posTri, returnList=T)
## ## # - checking which mesh triangles are inside the normal area
## ## normal = unlist(normal)
## ## Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
## ## Omega.SP = dt.polygon.omega(mesh, Omega)
## ## plot(Omega.SP[[2]], col="grey", main="The barrier region (in grey)")

## ## Q.barrier = dt.create.Q(mesh, Omega, 
## ##                         fixed.ranges = c(NA, 0.5))
## ## # - We fix the barrier range to a different value than we 
## ## #   used for simulations
## ## # - - Why? It does not matter, as long as it is 'small' 
## ## #     the models are very
## ## #     similar
## ## # - - This shows that you do not need to know the 
## ## #     true 'barrier range'!
## ## # - time: Ca 1 min

## ## log.prior = dt.create.prior.log.exp(
## ##   prior.param = c(-log(prior.sigma[2])/prior.sigma[1], -log(prior.range[2])/prior.range[1])) 
## ##     #c(-log(0.01)/3, -log(0.5)*6))
## ## # - The prior parameters are the lambdas in the exponential 
## ## #   priors for standard deviation and inverse-range
## ## # - the first is log(prob)/exceed, the second log(prob)*exceed
## ## # - the second is exponential for inverse range, therefore multiplication!

## ## barrier.model = dt.inla.model(
## ##   Q = Q.barrier, log.prior=log.prior)

## ## M[[2]] = list(shortname="barrier-model")
## ## M[[2]]$formula = y~ -1+m + f(s, model=barrier.model)

## ## ## Stationary model where the mesh stops at the boundary
## ## mesh2 = inla.mesh.2d(boundary=DH.wh.sp,
## ##                     max.edge = max.edge,
## ##                     #cutoff = 0.1,
## ##                     cutoff = 0.01)

## ## plot(mesh2, main="The second mesh", lwd=0.5)

## ## ## This Stack and A matrix
## ## A.i.s2 = inla.spde.make.A(mesh2, loc=cbind(sediment.data$Longitude, sediment.data$Latitude))
## ## stk2 = inla.stack(data=list(y=sediment.data$Mg), 
## ##                     effects=list(s=1:mesh2$n,
## ##                                  m = rep(1, nrow(sediment.data))),
## ##                    A=list(A.i.s2, 1),
## ##                     remove.unused = FALSE, tag='est') 
## ## spde2 = inla.spde2.pcmatern(mesh2, prior.range=prior.range, prior.sigma=prior.sigma)

## ## M[[3]] = list(shortname="neumann-model")
## ## M[[3]]$formula = y~ -1+m + f(s, model=spde2)
## ## M[[3]]$stack = stk2

## ## ## Running all the models
## ## ## Initial values
## ## # - speeds up computations
## ## # - improves accuracy of computations
## ## # - set these to NULL the first time you run a model
## ## M[[1]]$init = c(7.142,0.314,-0.648)
## ## M[[2]]$init = c(6.986,-1.135,-0.603)
## ## M[[3]]$init = c(7.221,-0.953,-1.383)

## ## hyper.iid = list(prec = list(prior='pc.prec', param=prior.sigma))
## ## i = 1
## ## stack=stk
## ## if (!is.null(M[[i]]$stack)) stack = M[[i]]$stack
## ## M[[i]]$res = inla(M[[i]]$formula,
## ##                   data=inla.stack.data(stack),
## ##                   control.predictor=list(A=inla.stack.A(stack)),
## ##                   family="gaussian", 
## ##                   control.family = list(hyper = hyper.iid),
## ##                                         #control.family = list(hyper = hyper.fixed),
## ##                   control.inla= list(int.strategy = "eb"),
## ##                                         #verbose=T,
## ##                   control.mode=list(restart=T, theta=M[[i]]$init))

## ## print(paste(round(M[[i]]$res$internal.summary.hyperpar$mode, 3), collapse = ','))

## ## field = M[[i]]$res$summary.random$s$mean + M[[i]]$res$summary.fixed['m', 'mean']
## ## xlim = DH.wh.sp@bbox[1, ] 
## ## ylim = DH.wh.sp@bbox[2, ]
## ## proj = inla.mesh.projector(mesh, xlim = xlim, 
## ##                            ylim = ylim, dims=c(300, 300))
## ## field.proj = inla.mesh.project(proj, field)
## ## zlim = range(field.proj, na.rm=TRUE)



## ## local.plot.field(field, mesh, main=paste(), zlim=global.zlim)

## ## ##make the stack
## ## phi = inla.spde.make.A(mesh=mesh,
## ##                        loc=cbind(sediment.data$Longitude, sediment.data$Latitude))
## ## stack_est = inla.stack(data=list(y=sediment.data$Mg),
## ##                        A=list(phi,1,1),
## ##                        effects = list(s=1:mesh$n,
## ##                                       data.frame(Intercept=1,
## ##                                                  lon = sediment.data$Longitude,
## ##                                                  lat = sediment.data$Latitude),
## ##                                       iidx=1:nrow(sediment.data)),
## ##                        remove.unused=FALSE,tag='est')
## ## spde = inla.spde2.pcmatern(mesh, prior.range = c(6, .5), prior.sigma = c(3, 0.01))
## ## hyper.iid = list(prec = list(prior = 'pc.prec', param = c(3, 0.01)))
## ## ##Barrier model
## ## tl =length(mesh$graph$tv[,1]) ##number of triangles in the mesh
## ## posTri = matrix(0, tl, 2)
## ## for (t in 1:tl) {
## ##     temp=mesh$loc[mesh$graph$tv[t,],]
## ##     posTri[t,] = colMeans(temp)[c(1,2)]
## ## }
## ## posTri = SpatialPoints(posTri) ##compute the triangle positions
## ## proj4string(posTri) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
## ## normal = over(DH.wh.sp, posTri, returnList = TRUE) ## check which mesh triangles are inside the normal area
## ## normal = unlist(normal)
## ## barrier.triangles=setdiff(1:tl, normal)
## ## poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)
## ## plot(poly.barrier)

## ## barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles,
## ##                                       prior.range=c(6,0.5),
## ##                                       prior.sigma=c(3,0.01))
## ## formula = y~ -1+Intercept + lon + lat +
## ##     f(s, model=barrier.model) +
## ##     f(iidx, model="iid", hyper=hyper.iid)

## ## output = inla(formula,
## ##               data=inla.stack.data(stack_est),
## ##               control.predictor=list(A=inla.stack.A(stack_est), compute=TRUE),
## ##               verbose=TRUE
## ## #              control.inla = list(int.strategy = 'eb')
## ## )



## ## spde = inla.spde2.matern(mesh=mesh, alpha=2)
## ## s_index = inla.spde.make.index(name='Spatial.field', n.spde=mesh$n)
## ## phi = inla.spde.make.A(mesh=mesh,
## ##                        loc=cbind(sediment.data$Longitude, sediment.data$Latitude))
## ## stack_est = inla.stack(data=list(y=sediment.data$Mg),
## ##                        A=list(phi,1,1,1),
## ##                        effects = list(s_index,
## ##                                       list(Intercept=rep(1,nrow(sediment.data))),
## ##                                       list(lon = sediment.data$Longitude),
## ##                                       list(lat = sediment.data$Latitude)),
## ##                        tag='est')

## ## pts=as.data.frame(spsample(DH.wh.sp,n=500, type='regular'))
## ## colnames(pts) <- c('Longitude','Latitude')
## ## pts = data.frame(Longitude=mesh$loc[,1], Latitude=mesh$loc[,2])
## ## points(pts)

## ## phi_pred = Diagonal(n=nrow(pts))
## ## stack_pred = inla.stack(data=list(y=NA),
## ##                         A=list(phi_pred,1,1,1),
## ##                         effects = list(s_index,
## ##                                        list(Intercept=rep(1,nrow(pts))),
## ##                                        list(lon=pts$Longitude),
## ##                                        list(lat=pts$Latitude)),
## ##                         tag='pred')

## ## stack = inla.stack(stack_est, stack_pred)
## ## formula = y ~ -1 + lat + long + Intercept + f(Spatial.field, model=spde)
## ## output = inla(formula,
## ##               data=inla.stack.data(stack, spde=spde),
## ##               control.predictor=list(A=inla.stack.A(stack), compute=TRUE))







## ## outer <-read.xls('data/HabitatMappingOffset_ShallowOuterDarwinHarbour.xlsx',
## ##                  sheet='All_Element_GS_Data (ppm)',header=TRUE)
## ## outer = outer %>% dplyr::select(-matches('X.[0-9]')) %>% dplyr::select(-X)

## ## outer.s = outer %>%
## ##     mutate(Missing = apply(., 1, function(x) sum(is.na(x)))) %>% filter(Missing < 3) %>% #exclude rows (sites) with more than 2 NA's
## ##     dplyr::select_if(~-sum(is.na(.))==0) # remove columns with any missing values

## ## outer.x = outer.s %>%
## ##     dplyr::select(-Survey,-GA.Sample.number,-GA.Sample.ID, -Sample.type, -Longitude, -Latitude) %>%
## ##     as.matrix %>% # convert to matrix
## ##     decostand(method='standardize', na.rm=TRUE) %>% #standardize
## ##     vegdist(method='euclidean') #euclidean distance

## ## outer.hclust = hclust(outer.x, labels=outer.s)
## ## cl=kmeans(outer.x, 5)
## ## plot(Latitude~Longitude, outer.s, col=cl$cluster)

## ## ggplot(outer, aes(y=Latitude, x=Longitude)) + geom_point()
