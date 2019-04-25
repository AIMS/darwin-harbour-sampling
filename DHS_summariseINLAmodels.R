source('DHS_functions.R')
DHS_checkPackages()

source('DHS_config.R')

load(file='data/processed/sediment.data.RData')
load(file='data/processed/DH.wh.sp.RData')
load(file='data/processed/DH.wh.df.RData')
load(file='data/processed/East.arm.sp.RData')
load(file='data/processed/East.arm.df.RData')
load(file='data/processed/sediment.data.east.RData')
load(file='data/processed/outer_sites.RData')
load(file='data/processed/Outer.sp.RData')
load(file='data/processed/Outer.df.RData')

## East Arm models
vars = names(sediment.data)
vars = vars[!vars %in% c("Sample","Comments","Zone","Date","Latitude","Longitude","Datum","Residue")]
for (v in vars) {
    load(file=paste0('data/processed/east.arm.mod_',v,'.RData'))
    ## East arm model
    g1=ggplot() + inlabru:::gg.inla.mesh(east.arm.mod$mesh) +
        geom_point(data=sediment.data, aes(y=Latitude, x=Longitude), size=1) +
        geom_polygon(data=DH.wh.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_polygon(data=East.arm.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        coord_equal() +
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='a)')

    g2=ggplot() + inlabru:::gg.inla.mesh(east.arm.mod$mesh) +
        geom_point(data=sediment.data, aes(y=Latitude, x=Longitude), size=1) +
        geom_polygon(data=East.arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        #geom_polygon(data=broom::tidy(east.arm.mod$Omega.SP[[2]]), aes(y=lat, x=long, group=id), fill='red', color=NA, show.legend=FALSE, alpha=0.5) +
        geom_sf(data=st_as_sf(east.arm.mod$Omega.SP[[2]]), fill='red',color=NA, show.legend=FALSE, alpha=0.5) +
                                        #coord_equal() +
        coord_sf() +
        theme_bw()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='b)')

    ## It would be good to use point.in.polygon, but it does not work properly when there are multiple disconnected polygons
    ##wch = point.in.polygon(east.arm.mod$fit$Longitude, east.arm.mod$fit$Latitude,East.arm.df$long, East.arm.df$lat)
    wch=matrix(nrow=nrow(east.arm.mod$fit), ncol=length(unique(East.arm.df$piece)))
    i =0
    for (p in unique(East.arm.df$group)) { # for each polygon, assess whether the point is in the polygon or not
        i = i+1
        e.a.d=East.arm.df %>% filter(group==p)
        wch[,i] = point.in.polygon(east.arm.mod$fit$Longitude, east.arm.mod$fit$Latitude,e.a.d$long, e.a.d$lat)
    }
    wch=rowSums(wch)
    wch=ifelse(wch==0,0,1)
    newdata = east.arm.mod$fit[wch==1,]
    sediment.data.east = sediment.data.east %>% mutate_(.dots=setNames(v,'Value'))
    Range = range(sediment.data.east$Value, newdata$fit, na.rm=TRUE)

    g3=ggplot() +
        geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
        scale_fill_gradientn('Predicted', colors=tim.colors(64), limits=Range) +
        geom_polygon(data=East.arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
        geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude, color=Value)) +
        scale_color_gradientn('Observed',colors=tim.colors(64), limits=Range) +
        coord_equal()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='c)')

    g4=ggplot() +
        geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
        scale_fill_gradientn('Predicted', colors=tim.colors(64)) +
        geom_polygon(data=East.arm.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        #geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
        #geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude, color=Value)) +
        #scale_color_gradientn('Observed',colors=tim.colors(64), limits=Range) +
        coord_equal()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='d)')

    ggsave(filename=paste0('output/inla_model_east.arm',v,'.pdf'),
           grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
           width=10, height=10)
    ggsave(filename=paste0('output/inla_model_east.arm_',v,'.png'),
           grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
           width=10, height=10, units='in', dpi=300)

    ggsave(filename=paste0('output/inla_model_east.arm1',v,'.pdf'),
           grid.arrange(g1,g3, g2, g4, layout_matrix=rbind(c(1,2), c(3,4)), ncol=2, widths=c(0.7,1)),
           width=10, height=10)
    ggsave(filename=paste0('output/inla_model_east.arm1_',v,'.png'),
           grid.arrange(g1,g3, g2, g4, layout_matrix=rbind(c(1,2), c(3,4)), ncol=2, widths=c(0.7,1)),
           width=10, height=10, units='in', dpi=300)

    east.arm.pred = newdata
    east.arm.pred.rast = rasterFromXYZ(east.arm.pred %>% dplyr::select(x=Longitude,y=Latitude,z=fit), crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
    save(east.arm.pred, file=paste0('data/processed/east.arm.pred_',v,'.RData'))
    save(east.arm.pred.rast, file=paste0('data/processed/east.arm.pred.rast_',v,'.RData'))
    
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


## Outer Harbour models
outer_sites = outer_sites %>%
    dplyr::select(-Survey,-Location,-`GA Sample number`,-`GA Sample ID`, `Sample type`,-`X9`,-`Water depth (m)`,`Sand %`,
                         `Mud %`, `Gravel %`) %>%
    setNames(gsub('([a-zA-Z]*)\ .*','\\1',names(.)))
vars = names(outer_sites)
vars = vars[!vars %in% c("Sample","Survey","Location","GA Sample number","GA Sample ID", "Sample type","Latitude","Longitude","X9","Water depth (m)","Sand",
                         "Mud", "Gravel","Total")]
for (v in vars) {
    load(file=paste0('data/processed/outer.mod_',v,'.RData'))
    ## East arm model
    g1=ggplot() + inlabru:::gg.inla.mesh(outer.mod$mesh) +
        geom_point(data=outer_sites, aes(y=Latitude, x=Longitude), size=1) +
        geom_polygon(data=DH.wh.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_polygon(data=Outer.sp, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        coord_equal() +
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='a)')

    g2=ggplot() + inlabru:::gg.inla.mesh(outer.mod$mesh) +
        geom_point(data=outer_sites, aes(y=Latitude, x=Longitude), size=1) +
        geom_polygon(data=Outer.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_sf(data=st_as_sf(outer.mod$Omega.SP[[2]]), fill='red',color=NA, show.legend=FALSE, alpha=0.5) +
        #geom_polygon(data=broom::tidy(outer.mod$Omega.SP[[2]]), aes(y=lat, x=long, group=id), fill='red', color=NA, show.legend=FALSE, alpha=0.5) +
        coord_sf() +
        theme_bw()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='b)')

    wch = point.in.polygon(outer.mod$fit$Longitude, outer.mod$fit$Latitude,
                           Outer.df$long, Outer.df$lat)
    newdata = outer.mod$fit[wch==1,]
    outer_sites = outer_sites %>% mutate_(.dots=setNames(v,'Value'))
    Range = range(outer_sites$Value, newdata$fit, na.rm=TRUE)

    g3=ggplot() +
        geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
        scale_fill_gradientn('Predicted', colors=tim.colors(64), limits=Range) +
        geom_polygon(data=Outer.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
        geom_point(data=outer_sites, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
        geom_point(data=outer_sites, aes(y=Latitude, x=Longitude, color=Value)) +
        scale_color_gradientn('Observed',colors=tim.colors(64), limits=Range) +
        coord_equal()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='c)')

    g4=ggplot() +
        geom_tile(data=newdata, aes(y=Latitude, x=Longitude, fill=fit)) +
        scale_fill_gradientn('Predicted', colors=tim.colors(64)) +
        geom_polygon(data=Outer.df, aes(y=lat,x=long, group=group), fill=NA, color='black') +
                                        #geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude), color='black', size=1.5) +
                                        #geom_point(data=sediment.data.east, aes(y=Latitude, x=Longitude, color=Value)) +
                                        #scale_color_gradientn('Observed',colors=tim.colors(64), limits=Range) +
        coord_equal()+
        theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
        annotate(geom='text', y=Inf, x=-Inf, hjust=-0.1, vjust=1.1, label='d)')

    ggsave(filename=paste0('output/inla_model_outer',v,'.pdf'),
           grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
           width=10, height=10)
    ggsave(filename=paste0('output/inla_model_outer_',v,'.png'),
           grid.arrange(g1,g2, g3, layout_matrix=rbind(c(1,2), c(3,3)), ncol=2),
           width=10, height=10, units='in', dpi=300)

    ggsave(filename=paste0('output/inla_model_outer1',v,'.pdf'),
           grid.arrange(g1,g3, g2, g4, layout_matrix=rbind(c(1,2), c(3,4)), ncol=2, widths=c(0.7,1)),
           width=10, height=10)
    ggsave(filename=paste0('output/inla_model_outer1_',v,'.png'),
           grid.arrange(g1,g3, g2, g4, layout_matrix=rbind(c(1,2), c(3,4)), ncol=2, widths=c(0.7,1)),
           width=10, height=10, units='in', dpi=300)

    outer.pred = newdata
    outer.pred.rast = rasterFromXYZ(outer.pred %>% dplyr::select(x=Longitude,y=Latitude,z=fit), crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'))
    save(outer.pred, file=paste0('data/processed/outer.pred_',v,'.RData'))
    save(outer.pred.rast, file=paste0('data/processed/outer.pred.rast_',v,'.RData'))
    
}
