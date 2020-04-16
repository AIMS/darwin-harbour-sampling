source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/hydro.middle_arm.sp.RData')
load(file='data/processed/hydro.middle_arm.df.RData')
load(file='data/processed/sediment.data_MA.RData')
load(file='data/processed/sites_MA.RData')
load(file='data/processed/DH.zone_MA.sp.RData')
load(file='data/processed/DH.zone_MA.df.RData')
load(file='data/processed/middle_arm.sp.RData')
load(file='data/processed/middle_arm.df.RData')

vars = names(sediment.data_MA)
vars = vars[!vars %in% c("Sample","Comments","Zone","Date","Latitude","Longitude","Datum","Residue","depth")]

for (v in vars) {  
    print(v)
    (d=diff(range(sediment.data_MA$Latitude))/2)
    (f=log(median(sediment.data_MA[[v]], na.rm=TRUE)))
    if (f<0.1) f=3
    prior.range=c(d,0.5) #prior.range=c(1,0.5) #priors[[v]]$prior.range
    prior.sigma=c(f,0.01) #priors[[v]]$prior.sigma

    ## Incase we need to fit the model on UTM projected data
    ## East.arm.sp = spTransform(East.arm.sp, CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
    ## coordinates(sediment.data) <- ~Longitude+Latitude
    ## proj4string(sediment.data) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
    ## sediment.data = spTransform(sediment.data, CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
    ## sediment.data = as.data.frame(sediment.data) 

  ## fit the Middle Arm model
  middle_arm.mod = fitINLA.barriermodel(bndry=middle_arm.sp, data=sediment.data_MA, var=v, mesh.type='boundary', prior.range = prior.range, prior.sigma = prior.sigma, max.edge=0.02)
  save(middle_arm.mod, file=paste0('data/processed/middle_arm.mod_',v,'.RData'))

  ## fit the Whole harbour model
  #harbour.mod = fitINLA.barriermodel(bndry=DH.wh.sp, data=sediment.data, var=v, prior.range = prior.range, prior.sigma = prior.sigma)
  #save(harbour.mod, file=paste0('data/processed/harbour.mod_',v,'.RData'))
}
