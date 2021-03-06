source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/hydro.outer_east.sp.RData')
load(file='data/processed/hydro.outer_east.df.RData')
load(file='data/processed/sediment.data.RData')
load(file='data/processed/outer_sites.RData')
load(file='data/processed/DH.zone.sp.RData')
load(file='data/processed/DH.zone.df.RData')
load(file='data/processed/East.arm.sp.RData')
load(file='data/processed/East.arm.df.RData')

## Fit the INLA barrier models (East Arm)
vars = names(sediment.data)
vars = vars[!vars %in% c("Sample","Comments","Zone","Date","Latitude","Longitude","Datum","Residue")]
##for (v in c('Mg','Al','P','S','Ca')) {
for (v in vars) {  
    print(v)
    (d=diff(range(sediment.data$Latitude))/2)
    (f=log(median(sediment.data[[v]])))
    if (f<0.1) f=3
    prior.range=c(d,0.5) #prior.range=c(1,0.5) #priors[[v]]$prior.range
    prior.sigma=c(f,0.01) #priors[[v]]$prior.sigma

    ## Incase we need to fit the model on UTM projected data
    ## East.arm.sp = spTransform(East.arm.sp, CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
    ## coordinates(sediment.data) <- ~Longitude+Latitude
    ## proj4string(sediment.data) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
    ## sediment.data = spTransform(sediment.data, CRS('+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
    ## sediment.data = as.data.frame(sediment.data) 

  ## fit the East Arm model
  east.arm.mod = fitINLA.barriermodel(bndry=East.arm.sp, data=sediment.data, var=v, prior.range = prior.range, prior.sigma = prior.sigma, max.edge=0.02)
  save(east.arm.mod, file=paste0('data/processed/east.arm.mod_',v,'.RData'))

  ## fit the Whole harbour model
  #harbour.mod = fitINLA.barriermodel(bndry=DH.wh.sp, data=sediment.data, var=v, prior.range = prior.range, prior.sigma = prior.sigma)
  #save(harbour.mod, file=paste0('data/processed/harbour.mod_',v,'.RData'))
}


## Fit the INLA barrier models (Outer Harbour)
outer_sites = outer_sites %>%
    dplyr::select(-Survey,-Location,-`GA Sample number`,-`GA Sample ID`, `Sample type`,-`X9`,-`Water depth (m)`,`Sand %`,
                         `Mud %`, `Gravel %`) %>%
    setNames(gsub('([a-zA-Z]*)\ .*','\\1',names(.)))
vars = names(outer_sites)
vars = vars[!vars %in% c("Sample","Survey","Location","GA Sample number","GA Sample ID", "Sample type","Latitude","Longitude","X9","Water depth (m)","Sand",
                         "Mud", "Gravel","Total")]
for (v in vars) {  
    print(v)
    (d=diff(range(outer_sites$Latitude))/2)
    (f=log(median(outer_sites[[v]], na.rm=TRUE)))
    if (f<0.1) f=3
    prior.range=c(d,0.5) #prior.range=c(1,0.5) #priors[[v]]$prior.range
    prior.sigma=c(f,0.01) #priors[[v]]$prior.sigma

  ## fit the Outer Harbour model
  outer.mod = fitINLA.barriermodel(bndry=Outer.sp, data=outer_sites, var=v, prior.range = prior.range, prior.sigma = prior.sigma, max.edge=0.07)
  save(outer.mod, file=paste0('data/processed/outer.mod_',v,'.RData'))
}
