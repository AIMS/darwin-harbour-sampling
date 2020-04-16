source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/sediment.data_MA.RData')
load(file='data/processed/DH.wh_MA.sp.RData')
load(file='data/processed/DH.wh_MA.df.RData')
load(file='data/processed/middle_arm.sp.RData')
load(file='data/processed/middle_arm.df.RData')
load(file='data/processed/sediment.data_MA.RData')
load(file='data/processed/sites_MA.RData')


## East Arm models
vars = names(sediment.data_MA)
vars = vars[!vars %in% c("Sample","Comments","Zone","Date","Latitude","Longitude","Datum","Residue","depth")]
for (v in vars) {
    load(file=paste0('data/processed/middle_arm.mod_',v,'.RData'))
    ## East arm model
    g1=ggplot() + inlabru:::gg.inla.mesh(middle_arm.mod$mesh) +
        geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude), size=1) +
        geom_polygon(data=DH.wh_MA.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_polygon(data=middle_arm.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        coord_equal() +
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='a)')

    g2=ggplot() + inlabru:::gg.inla.mesh(middle_arm.mod$mesh) +
        geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude), size=1) +
        geom_polygon(data=middle_arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        #geom_polygon(data=broom::tidy(east.arm.mod$Omega.SP[[2]]), aes(y=lat, x=long, group=id), fill='red', color=NA, show.legend=FALSE, alpha=0.5) +
        geom_sf(data=st_as_sf(middle_arm.mod$Omega.SP[[2]]), fill='red',color=NA, show.legend=FALSE, alpha=0.5) +
                                        #coord_equal() +
        coord_sf() +
        theme_bw()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='b)')

    ## It would be good to use point.in.polygon, but it does not work properly when there are multiple disconnected polygons
    ##wch = point.in.polygon(east.arm.mod$fit$Longitude, east.arm.mod$fit$Latitude,East.arm.df$long, East.arm.df$lat)
    wch=matrix(nrow=nrow(middle_arm.mod$fit), ncol=length(unique(middle_arm.df$piece)))
    i =0
    for (p in unique(middle_arm.df$group)) { # for each polygon, assess whether the point is in the polygon or not
        i = i+1
        e.a.d=middle_arm.df %>% filter(group==p)
        wch[,i] = point.in.polygon(middle_arm.mod$fit$Longitude, middle_arm.mod$fit$Latitude,e.a.d$long, e.a.d$lat)
    }
    wch=rowSums(wch)
    wch=ifelse(wch==0,0,1)
    newdata = middle_arm.mod$fit[wch==1,]
    sediment.data_MA = sediment.data_MA %>% mutate_(.dots=setNames(v,'Value'))
    Range = range(sediment.data_MA$Value, newdata$fit, na.rm=TRUE)

    g3=ggplot() +
        geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
        scale_fill_gradientn('Predicted', colors=tim.colors(64), limits=Range) +
        geom_polygon(data=middle_arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
        geom_point(data=sediment.data_MA, aes(y=Latitude, x=Longitude, color=Value)) +
        scale_color_gradientn('Observed',colors=tim.colors(64), limits=Range) +
        coord_equal()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='c)')

    g4=ggplot() +
        geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
        scale_fill_gradientn('Predicted', colors=tim.colors(64)) +
        geom_polygon(data=middle_arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        #geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
        #geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude, color=Value)) +
        #scale_color_gradientn('Observed',colors=tim.colors(64), limits=Range) +
        coord_equal()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='d)')

    ggsave(filename=paste0('output/inla_model_middle_arm',v,'.pdf'),
           grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
           width=10, height=10)
    ggsave(filename=paste0('output/inla_model_middle_arm_',v,'.png'),
           grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
           width=10, height=10, units='in', dpi=300)

    ggsave(filename=paste0('output/inla_model_middle_arm1',v,'.pdf'),
           grid.arrange(g1,g3, g2, g4, layout_matrix=rbind(c(1,2), c(3,4)), ncol=2, widths=c(0.7,1)),
           width=10, height=10)
    ggsave(filename=paste0('output/inla_model_middle_arm1_',v,'.png'),
           grid.arrange(g1,g3, g2, g4, layout_matrix=rbind(c(1,2), c(3,4)), ncol=2, widths=c(0.7,1)),
           width=10, height=10, units='in', dpi=300)

    middle_arm.pred = newdata
    middle_arm.pred.rast = rasterFromXYZ(middle_arm.pred %>% dplyr::select(x=Longitude,y=Latitude,z=fit), crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
    save(middle_arm.pred, file=paste0('data/processed/middle_arm.pred_',v,'.RData'))
    save(middle_arm.pred.rast, file=paste0('data/processed/middle_arm.pred.rast_',v,'.RData'))
    
    ## ## Harbour model
    ## g1=ggplot() + inlabru:::gg.inla.mesh(harbour.mod$mesh) +
    ##     geom_point(data=sediment.data, aes(y=Latitude, x=Longitude), size=1) +
    ##     geom_polygon(data=DH.wh.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
    ##     geom_polygon(data=East.arm.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
    ##     coord_equal() +
    ##     theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    ##     annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='a)')

    ## g2=ggplot() + inlabru:::gg.inla.mesh(harbour.mod$mesh) +
    ##     geom_point(data=sediment.data, aes(y=Latitude, x=Longitude), size=1) +
    ##     geom_polygon(data=East.arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
    ##     geom_polygon(data=broom::tidy(harbour.mod$Omega.SP[[2]]), aes(y=lat, x=long, group=id), fill='red', color=NA, show.legend=FALSE, alpha=0.5) +
    ##     coord_equal() +
    ##     theme_bw()+
    ##     theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
    ##     annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='b)')


    ## Harbour.df = broom::tidy(DH.wh.sp)
    ## wch = point.in.polygon(harbour.mod$fit$Longitude, harbour.mod$fit$Latitude,
    ##                        DH.wh.df$long, DH.wh.df$lat)
    ## newdata = harbour.mod$fit[wch==1,]
    ## sediment.data = sediment.data %>% mutate_(.dots=setNames(v,'Value'))
    ## Range = range(sediment.data$Value, newdata$fit, na.rm=TRUE)
    ##                                     #Range1 = range(sediment.data.east$Value)
    ##                                     #Range2 = range(newdata$fit, na.rm=TRUE)
    ## g3=ggplot() +
    ##     geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
    ##     scale_fill_gradientn(colors=tim.colors(64), limits=Range) +
    ##     geom_polygon(data=DH.wh.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
    ##     geom_point(data=sediment.data, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
    ##     geom_point(data=sediment.data, aes(y=Latitude, x=Longitude, color=Value)) +
    ##     scale_color_gradientn(colors=tim.colors(64), limits=Range) +
    ##     coord_equal()+
    ##     theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
    ##     annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='c)')


    ## ggsave(filename=paste0('output/inla_model_harbour_',v,'.pdf'),
    ##        grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
    ##        width=10, height=10)
    ## ggsave(filename=paste0('output/inla_model_harbour_',v,'.png'),
    ##        grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
    ##        width=10, height=10, units='in', dpi=300)

}
