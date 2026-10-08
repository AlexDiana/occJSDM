native_hyperparameters <- function(fit_record,generating_range=NULL) {
  f <- fit_record; rows <- list()
  add <- function(name,values,units,definition,comparable=FALSE,true_value=NA_real_) {
    q <- quantile(values,c(.025,.5,.975),names=FALSE)
    posterior <- f$package %in% c("occJSDM","Hmsc","oracle")
    rows[[length(rows)+1L]] <<- data.frame(parameter=name,lower=if(posterior) q[1] else NA_real_,median=q[2],upper=if(posterior) q[3] else NA_real_,units=units,
      definition=definition,range_recovery_supported=comparable,true_value=true_value,
      interval_type=if(posterior) "posterior interval" else "point estimate; interval unavailable")
  }
  if(f$package=="occJSDM" && f$arm=="gp") {
    js <- f$native$results_output$jsdm_output; ell <- f$native$infos$l_s_grid[js$idx_ls_output]
    scaling <- f$native$infos$list_Xs_mat$sd_df
    comparable <- length(scaling)==2L && max(abs(scaling-scaling[1]))<1e-12
    add("native_length_scale",ell,"standardised coordinate units","Base kernel exp(-d^2/(2*l^2)); fitted field is subset of regressors")
    for(k in seq_along(scaling)) add(paste0("half_correlation_distance_axis_",k),sqrt(2*log(2))*ell*scaling[k],"raw unit-square coordinate distance","Base squared-exponential kernel along one raw axis")
    if(comparable) add("raw_isotropic_length_scale",ell*scaling[1],"raw unit-square coordinate distance","Equal axis scales; matching squared-exponential definition",TRUE,if(is.null(generating_range)) NA_real_ else generating_range)
    add("spatial_prior_sd",js$sigmabs_output,"logit field units","Shared prior SD of spatial coefficients; distinct from realised field SD")
  } else if(f$package=="Hmsc" && f$arm=="spatial2") {
    for(k in 1:2) {
      a <- unlist(lapply(f$native$postList,function(ch) vapply(ch,function(d) sort(f$native$ranLevels$site$alphapw[d$Alpha[[1]],1])[k],numeric(1))))
      add(paste0("ordered_factor_half_correlation_distance_",k),a*log(2),"raw unit-square coordinate distance","Exponential factor kernel; factors ordered within each draw; zero denotes spatial independence")
    }
  } else if(f$package=="gllvm" && f$arm=="spatial2") {
    for(k in 1:2) add(paste0("factor_half_correlation_distance_",k),f$native$params$rho.lv[k]*log(2),"raw unit-square coordinate distance","Exponential factor kernel; not a squared-exponential generating range")
  } else return(data.frame(parameter="range",lower=NA_real_,median=NA_real_,upper=NA_real_,units="unavailable",definition="No fitted range in this arm",range_recovery_supported=FALSE,true_value=NA_real_,interval_type="unavailable"))
  do.call(rbind,rows)
}
