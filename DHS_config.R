source('DHS_functions.R')
DHS_checkPackages()

DHS_tryCatch(
    {
        config = readLines('parameters/DHS.conf')
        config = gsub('(.*Date)=(.*)','\\1=as.Date(\'\\2\')',config)
        eval(parse(text=config))
        return=NULL
    }, 'logs/all.log','--Config--',msg='load general configurations', return=NULL)



DHS_tryCatch(
    {
        files <- list.files('data')
        if(!any(grepl('^primary$',files))) system('mkdir data/primary')
        if(!any(grepl('^processed$',files))) system('mkdir data/processed')

        files <- list.files()
        if(!any(grepl('^output$',files))) system('mkdir output')
        if(!any(grepl('^output$',files))) system('cd output; mkdir figures;')

        files <- list.files()
        if(!any(grepl('^logs$',files))) system('mkdir logs')
        return=NULL
    }, 'logs/all.log','--Config--',msg='configure necessary folders', return=NULL)

