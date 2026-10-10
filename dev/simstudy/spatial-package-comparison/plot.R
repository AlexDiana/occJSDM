# Scientific figures from the independently verified compact bundle only.
plot_spatial_comparison <- function(bundle,output_dir) {
  verify_teaching_bundle(bundle)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  library(ggplot2)
  cfg <- bundle$config; files <- character()
  save_plot <- function(p,name,width=11,height=7) {
    path <- file.path(output_dir,paste0("lesson-5-spatial-",name,".png"))
    ggsave(path,p,width=width,height=height,dpi=160,bg="white",limitsize=FALSE)
    files <<- c(files,path)
  }
  theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title.position="plot"))
  labels <- c(env="Environment only",trend="Quadratic coordinates",mem="Moran eigenvectors",gp="GP field",dnn="Coordinate network",spatial2="Two spatial factors")
  definitions <- data.frame(x=1,y=4:1,
    name=c("GP field","Quadratic trend","Moran eigenvectors","Coordinate network"),
    detail=c("Nearby sites share a field through a covariance kernel.\nThe fitted kernel controls smoothness and spatial distance.",
      "Five fixed columns: x, y, x², xy, y².\nEach species estimates its own coefficients.",
      "A graph of training neighbours defines fixed spatial columns.\nEach species estimates its own coefficients.",
      "A neural network learns coordinate features from observations.\nPenalties constrain the learned spatial contribution."))
  p <- ggplot(definitions,aes(x,y))+geom_label(aes(label=name),x=.22,hjust=0,size=4,label.size=NA,fill="#e7eff5")+
    geom_text(aes(label=detail),x=.9,hjust=0,size=3.5,lineheight=1.2)+xlim(0,3.1)+ylim(.5,4.5)+theme_void()+
    labs(title="Four ways to represent spatial effects",subtitle="Every representation is judged on the same training observations and held-out sites.")
  save_plot(p,"representations",12,5)
  # Fixed examples: community 1, species 2 and 6, all four arrangements.
  for(package in names(cfg$arms)) for(sp in cfg$species[c(2,6)]) {
    maps <- list(); missing <- list()
    for(arr in cfg$arrangements) {
      ds <- paste("rep01",arr,sep="-"); ti <- bundle$truth_index[[ds]]; land <- bundle$landscapes$rep01
      j <- match(sp,land$species); ix <- ti$lattice_index; xy <- land$points[ix,,drop=FALSE]
      baseline <- plogis(land$B0[j]+land$environment[ix]*land$B[j]+mean(land$field[ti$training_index,j]))
      maps[[length(maps)+1L]] <- data.frame(x=xy[,1],y=xy[,2],arrangement=arr,model="Truth",value=land$psi[ix,j]-baseline)
      for(arm in cfg$arms[[package]]) {
        key <- paste(package,arm,ds,sep="-"); pred <- bundle$predictions[[key]]$lattice
        if(is.null(pred)) missing[[length(missing)+1L]] <- data.frame(x=.5,y=.5,arrangement=arr,model=labels[[arm]],label="Prediction unavailable") else
          maps[[length(maps)+1L]] <- data.frame(x=xy[,1],y=xy[,2],arrangement=arr,model=labels[[arm]],value=pred$delta_p[,match(sp,pred$species_id)])
      }
    }
    d <- do.call(rbind,maps); d$model <- factor(d$model,levels=c("Truth",unname(labels[cfg$arms[[package]]])))
    d$arrangement <- factor(d$arrangement,levels=cfg$arrangements)
    p <- ggplot(d,aes(x,y,fill=value))+geom_tile()+facet_grid(model~arrangement,drop=FALSE)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
      scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",midpoint=0,limits=c(-1,1),name="Spatial probability\ncontrast")+
      labs(title=paste(package,sp,"on community 1"),subtitle="Truth and fitted delta_p use training-centred spatial effects; identical probability scale.",x="x coordinate",y="y coordinate")+
      theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())
    if(length(missing)) {
      z <- do.call(rbind,missing); z$model <- factor(z$model,levels=levels(d$model)); z$arrangement <- factor(z$arrangement,levels=cfg$arrangements)
      p <- p+geom_text(data=z,aes(x,y,label=label),inherit.aes=FALSE,size=2.7)
    }
    save_plot(p,paste("maps",package,sp,sep="-"),12,12)
  }
  s <- bundle$scores$species; s <- s[s$package!="oracle",,drop=FALSE]
  s$model <- paste(s$package,s$arm,sep=" / ")
  lattice <- s[s$scope=="lattice",,drop=FALSE]
  for(metric in c("delta_rmse","delta_rms","delta_correlation","field_rmse","field_attenuation","rmse")) {
    good <- lattice[is.finite(lattice[[metric]]),,drop=FALSE]
    if(!nrow(good)) next
    d <- aggregate(good[[metric]],good[c("model","community","arrangement","prevalence")],mean); names(d)[5] <- "value"
    p <- ggplot(d,aes(value,model,color=community))+geom_point(position=position_dodge(width=.55),size=1.8)+facet_grid(prevalence~arrangement,scales="free_x")+
      labs(title=switch(metric,delta_rmse="Spatial-contrast recovery",delta_rms="Spatial-contrast amplitude",delta_correlation="Spatial-contrast pattern recovery",field_rmse="Native logit-field recovery",field_attenuation="Native logit-field attenuation",rmse="Held-out occupancy-probability error"),
        subtitle=if(grepl("^field_",metric)) "Only matching logit fields; probit fields have no fabricated recovery scores." else "One point per community and prevalence group; failed predictions are counted separately.",x=switch(metric,delta_rmse="RMSE of delta_p",delta_rms="RMS of delta_p (black crosses: truth)",delta_correlation="Correlation of fitted and true delta_p",field_rmse="RMSE in logit field units",field_attenuation="Slope of fitted field on true field",rmse="RMSE of occurrence probability"),y=NULL,color="Community")
    if(metric=="delta_rms") {
      truth <- aggregate(good$truth_delta_rms,good[c("model","community","arrangement","prevalence")],mean); names(truth)[5] <- "value"
      p <- p+geom_point(data=truth,aes(group=community),shape=4,color="black",size=1.8,position=position_dodge(width=.55))
    }
    if(metric=="field_attenuation") p <- p+geom_vline(xintercept=1,linetype=2,color="grey50")
    save_plot(p,paste0("recovery-",metric),14,12)
  }
  bins <- bundle$scores$strata
  if(!is.null(bins) && nrow(bins)) {
    d <- bins[bins$scope=="lattice" & bins$stratification=="distance" & bins$package!="oracle" & is.finite(bins$rmse),,drop=FALSE]
    if(nrow(d)) {
      d <- aggregate(rmse~package+arm+community+prevalence+bin+target,d,mean)
      d$bin <- factor(d$bin,levels=c("<=0.02","(0.02,0.05]","(0.05,0.1]",">0.1"))
      p <- ggplot(d,aes(bin,rmse,color=target,group=interaction(community,target)))+geom_line(alpha=.4)+geom_point(size=1)+facet_grid(package+arm~prevalence)+
        labs(title="Held-out error by distance from a training site",subtitle="Lines preserve separate communities; empty bins remain in the data with zero count.",x="Nearest training-site distance",y="Probability RMSE",color="Prediction")+
        theme(axis.text.x=element_text(angle=30,hjust=1))
      save_plot(p,"distance-errors",13,23)
    }
  }
  a <- bundle$scores$accounting; a <- a[a$kind=="comparison",,drop=FALSE]
  counts <- as.data.frame(table(paste(a$package,a$arm,sep=" / "),a$status),stringsAsFactors=FALSE); names(counts) <- c("model","status","cells")
  save_plot(ggplot(counts,aes(cells,model,fill=status))+geom_col()+labs(title="Every comparison cell remains accounted for",x="Number of cells (12 per model)",y=NULL,fill="Fit status"),"diagnostics",12,9)
  ci <- lattice[is.finite(lattice$probability_width),,drop=FALSE]
  if(nrow(ci)) save_plot(ggplot(ci,aes(probability_width,probability_inclusion,color=package,shape=factor(prevalence)))+
    geom_point(alpha=.65)+facet_wrap(~arrangement)+geom_hline(yintercept=.95,linetype=2,color="grey50")+
    labs(title="Interval width and descriptive truth inclusion",subtitle="gllvm intervals condition on fitted coefficients, loadings and ranges; point-only arms have no fabricated intervals.",x="Mean width of 95% probability interval",y="Fraction of sites including truth",shape="Prevalence"),"uncertainty",12,8)
  # Species-level expected attribution is computed directly from exported matrices.
  attribution <- list()
  for(i in which(a$package %in% names(cfg$arms))) {
    job <- a[i,]; p <- bundle$predictions[[job$key]]$training
    if(is.null(p)) next
    idx <- bundle$truth_index[[job$dataset_id]]; land <- bundle$landscapes[[idx$community]]
    truep <- land$psi[idx$training_index,,drop=FALSE]
    logscore <- function(q) { q <- pmin(pmax(q,1e-12),1-1e-12); truep*log(q)+(1-truep)*log1p(-q) }
    spatial <- colMeans(.5*((logscore(p$p_env_zero)-logscore(p$p_both_zero))+(logscore(p$p_full)-logscore(p$p_space_zero))))
    matchrow <- which(lattice$key==job$key)
    if(length(matchrow)) {
      r <- lattice[matchrow,]; r <- r[match(p$species_id,r$species_id),]
      r$training_spatial_gain <- spatial; attribution[[length(attribution)+1L]] <- r
    }
  }
  if(length(attribution)) {
    d <- do.call(rbind,attribution); d <- d[is.finite(d$space_zero_benefit),]
    if(nrow(d)) save_plot(ggplot(d,aes(training_spatial_gain,space_zero_benefit,color=package,shape=community))+
      geom_hline(yintercept=0,color="grey60")+geom_vline(xintercept=0,color="grey60")+geom_point(alpha=.65)+facet_grid(prevalence~arrangement)+
      labs(title="Does fitted spatial structure improve held-out prediction?",subtitle="Negative gains remain visible; the top row contains the 5% species.",x="Expected spatial log-score gain at training sites",y="Held-out RMSE benefit versus zeroed spatial effect"),"rare-species-attribution",13,10)
  }
  invisible(files)
}
