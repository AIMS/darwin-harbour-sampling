###################################################################
## The following function checks to ensure that all the required ##
## packages are available on the system.                         ##
###################################################################
DHS_checkPackages <- function() {
    requiredPackages <- c('gdata', 'vegan','sp','maptools','tidyverse','raster',
                          'lubridate','INLA','mgcv','inlabru','gridExtra','clhs',
                          'ggspatial','sf','MBHdesign','vegan')
    lapply(requiredPackages, require, character.only = TRUE)    
    p=installed.packages()[,'Package']
    if (!all(requiredPackages %in% p)) {
        wch = requiredPackages[!requiredPackages %in% p]
        install.packages(wch, repos='http://cran.csiro.au')
        return(paste('The following packages are missing:',wch,' I have attempted to install them.  Please try again..'))
    }

}


#########################################################################
## The following function appends a log statement into a log file      ##
## parameters:                                                         ##
##     status:    a string indicating either 'FAILURE', 'SUCCESS' or   ##
##                'WARNING'                                            ##
##     logFile:   a character string representation of the log file    ##
##                name (including path relative to current working     ##
##                directory)                                           ##
##     Category:  a character string with a category to appear         ##
##                verbatim in the log                                  ##
##     msg1:      a character string with a message to appear verbatim ##
##                in the log                                           ##
#########################################################################
DHS_log <- function (status, logFile='data/logs/env.log',Category, msg1) {
    ## Check if the log file exists, and if it does not, create it
    d=dirname(logFile)
    files <- list.files(d)
    ## If the file does not exist - create it
    if (!all(files %in% list.files(path=d, pattern='*.log'))) file.create(logFile)
    #if(!any(grepl(paste0('^',logFile,'$'),files))) file.create(logFile) #system(paste0('touch ',logFile))
    now <- Sys.time()

    msg <- paste0(now, '|',status, ': ', Category, ' ',msg1)
    if( !is.null(msg)){ write(msg,file=paste0(logFile),append=TRUE)}

}

DHS_tryCatch <- function(expr, logFile,Category, expectedClass=NULL, msg=NULL, return=NULL, showWarnings=FALSE) {
    #msg <- paste0(now, '| ', msg)
    options(digits.secs=2)              ## switch to subsecond display
    max.warnings<-4
    warnings<-0
    W <- NULL
    w.handler <- function(w){ # warning handler
        m<-w$message
        if ((warnings < max.warnings) && (grepl ('DHS_WARNING', m)>0)) {
            DHS_log('WARNING', logFile,Category, paste(warnings, msg, m))
            warnings<<-warnings+1
        }
        invokeRestart("muffleWarning")
    }
    ret <- list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                    warning = w.handler),warning = W)
    if(!is.atomic(ret$value) && !is.null(ret$value$message)){
        ## An error occurred
        class(ret) <- "try-error"
        DHS_log('FAILED', logFile,Category, paste(msg, ret$value$message))
        if(!is.null(return)) {
            FALSE
        }#else return()
    } else {    #no error check for warning
        DHS_log('SUCCESS', logFile, Category, msg)
        if(!is.null(return)) {
            TRUE
        }
    }
}


##############################################################################################
## The following functions are taken from:
## source('https://haakonbakka.bitbucket.io/functions-barriers-dt-models-march2017.R')
## They are reproduced here incase the aforementioned dissappears..
##############################################################################################
source('functions-barriers-dt-models-march2017.R')


## The following function fits a barrier model in INLA

fitINLA.barriermodel = function(bndry, data, var, mesh.type='points',prior.range = c(1, .5), prior.sigma = c(3, 0.01), max.edge=0.02) {

    data = data %>% mutate_(.dots=setNames(var, 'Value'))
    #max.edge=2000 #max.edge=0.02
    bound.outer=0.01
    if (mesh.type=='boundary') {
        mesh = inla.mesh.2d(boundary = bndry,
                            loc=cbind(data$Longitude, data$Latitude),
                            max.edge = c(1,2)*max.edge,
                            cutoff = 0.005,
                            offset = c(max.edge, bound.outer))
    } else if (mesh.type=='points') {
        mesh = inla.mesh.2d(loc=cbind(data$Longitude, data$Latitude),
                            max.edge = c(0.5,1)*max.edge)#,
                            #cutoff = 0.0005,
                            #offset = c(max.edge, bound.outer*4))
    }
    ## plot(mesh, main="Our mesh", lwd=0.5)
    ## points(data$Longitude,data$Latitude, col="red")
    ## lines(bndry, col='black',cex=3)

    A.i.s = inla.spde.make.A(mesh, loc=cbind(data$Longitude, data$Latitude))
    stk = inla.stack(data=list(y=data$Value), 
                     effects=list(s=1:mesh$n,
                                  m = rep(1, nrow(data))),
                     A=list(A.i.s, 1),
                     remove.unused = FALSE, tag='est')

    ## The stationary model
    #prior.range = c(1, .5)
    #prior.sigma = c(3, 0.01)
    spde = inla.spde2.pcmatern(mesh, prior.range=prior.range, prior.sigma=prior.sigma)
    #spde = inla.spde2.matern(mesh)
    hyper.iid = list(prec = list(prior='pc.prec', param=prior.sigma))

    mesh = dt.mesh.addon.posTri(mesh)
    ## - compute the triangle positions
    posTri = SpatialPoints(mesh$posTri)
    proj4string(posTri) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
    #proj4string(posTri) <- '+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
    normal = over(bndry, posTri, returnList=T)
                                        # - checking which mesh triangles are inside the normal area
    normal = unlist(normal)
    Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
    Omega.SP = dt.polygon.omega(mesh, Omega)
    #plot(Omega.SP[[2]], col="grey", main="The barrier region (in grey)")

    Q.barrier = dt.create.Q(mesh, Omega, 
                            fixed.ranges = c(NA, 0.5))
                                        # - We fix the barrier range to a different value than we 
                                        #   used for simulations
                                        # - - Why? It does not matter, as long as it is 'small' 
                                        #     the models are very
                                        #     similar
                                        # - - This shows that you do not need to know the 
                                        #     true 'barrier range'!
                                        # - time: Ca 1 min

    log.prior = dt.create.prior.log.exp(
        prior.param = c(-log(prior.sigma[2])/prior.sigma[1], -log(prior.range[2])/prior.range[1])) 
                                        #c(-log(0.01)/3, -log(0.5)*6))
                                        # - The prior parameters are the lambdas in the exponential 
                                        #   priors for standard deviation and inverse-range
                                        # - the first is log(prob)/exceed, the second log(prob)*exceed
                                        # - the second is exponential for inverse range, therefore multiplication!

    barrier.model = dt.inla.model(
        Q = Q.barrier, log.prior=log.prior)

    mod = list(shortname="barrier-model")
    mod$formula = y~ -1+m + f(s, model=barrier.model)

    stack=stk

    ## Running all the models
    ## Initial values
                                        # - speeds up computations
                                        # - improves accuracy of computations
                                        # - set these to NULL the first time you run a model
    mod$init = c(NULL, NULL, NULL)

    mod$res = inla(mod$formula,
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack), link=1, compute=TRUE),
                   family="gamma", 
                   control.family = list(hyper = hyper.iid, link='log'),
                   control.inla= list(int.strategy = "eb"),
                   control.mode=list(restart=T, theta=mod$init))
    #summary(mod$res)

    #print(paste(round(mod$res$internal.summary.hyperpar$mode, 3), collapse = ','))

    field = mod$res$summary.random$s$mean + mod$res$summary.fixed['m', 'mean']
    xlim = bndry@bbox[1, ] 
    ylim = bndry@bbox[2, ]
    proj = inla.mesh.projector(mesh, xlim = xlim, 
                               ylim = ylim, dims=c(300, 300))
    field.proj = inla.mesh.project(proj, field)
    zlim = range(field.proj, na.rm=TRUE)

    ## image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
    ##            xlim = xlim, ylim = ylim, asp=1) 
    ## contour(x = proj$x, y=proj$y,
    ##         z = field.proj,
    ##         levels=seq(zlim[1], zlim[2],length.out = 10),
    ##         add=TRUE, drawlabels=F, col="white")
    ## plot(Omega.SP[[2]], add=T, border="black", col="white")
    ## points(data$Longitude, data$Latitude)
    coords.grid = as.matrix(expand.grid(Longitude=proj$x, Latitude=proj$y))
    field.proj = inla.mesh.projector(mesh, loc=coords.grid)
    newdata = data.frame(coords.grid, fit = exp(inla.mesh.project(field.proj, field)))
    return(list(fit=newdata, mesh=mesh, Omega.SP=Omega.SP, mod=mod))
}



inlaPredict = function(v, wch, area='East Arm') {
    if (area=='East Arm') {
        load(file=paste0('data/processed/east.arm.mod_',v,'.RData'))
        mod=east.arm.mod
    } else {
        load(file=paste0('data/processed/outer.mod_',v,'.RData'))
            mod=outer.mod
    }

    field = mod$mod$res$summary.random$s$mean + mod$mod$res$summary.fixed['m', 'mean']
    proj = inla.mesh.projector(mod$mesh, loc=wch)
    if (area=='East Arm') rm('east.arm.mod','mod')
    if (area!='East Arm') rm('outer.mod','mod')
    gc()
    
    exp(inla.mesh.project(proj, field))
}
