validate_spatial_diagnostics <- function(a) {
  S<-nrow(a$species);stopifnot(length(S)==1L,S>0L)
  expected<-c('amplitude',paste0('field_mean_',seq_len(S)),
    paste0('field_rms_',seq_len(S)),paste0('field_projection_',seq_len(S)))
  d<-a$spatial$diagnostics
  stopifnot(is.data.frame(d),!anyDuplicated(d$quantity),setequal(d$quantity,expected),
    all(c('rhat','ess_mean') %in% names(d)))
  invisible(TRUE)
}

pair_needs_long <- function(ig,hc) length(ig$reasons)>0L || length(hc$reasons)>0L

extension_gate <- function(pairs,selected,probe) {
  stopifnot(nrow(pairs)==9L,!anyDuplicated(pairs$community),
    identical(as.integer(table(pairs$grid_index)),c(3L,3L,3L)),
    setequal(pairs$grid_index,c(4,6,8)),length(selected)>0L,!is.null(probe))
  expected_names<-as.vector(outer(pairs$key,c('inverse_gamma','half_cauchy'),paste,sep=':'))
  stopifnot(length(selected)==18L,!anyDuplicated(names(selected)),
    setequal(names(selected),expected_names))
  invisible(lapply(selected,validate_spatial_diagnostics));validate_spatial_diagnostics(probe)
  change <- function(metric)mean(pairs[[paste0(metric,'_hc')]]-pairs[[paste0(metric,'_ig')]])
  raw <- change('raw_rmse');centred <- change('centred_rmse')
  correlation <- change('centred_correlation');range <- change('range_mae')
  occupancy <- change('occupancy_mae')
  improved <- sum(pairs$centred_rmse_hc<pairs$centred_rmse_ig)
  spatial <- sum(vapply(selected,function(a)length(spatial_flags(a$spatial$diagnostics))>0L,logical(1)))
  start <- length(spatial_flags(probe$spatial$diagnostics))
  stopifnot(all(is.finite(c(raw,centred,correlation,range,occupancy))))
  data.frame(criterion=c('Lower mean raw field RMSE','Lower mean centred field RMSE',
    'Higher mean centred field correlation','Centred RMSE improved in at least six communities',
    'No unresolved spatial diagnostics in either prior','Initialization sensitivity diagnostics pass',
    'Range MAE increase at most 0.01','Occupancy MAE increase at most 0.01'),
    value=c(raw,centred,correlation,improved,spatial,start,range,occupancy),
    pass=c(raw<0,centred<0,correlation>0,improved>=6,spatial==0,start==0,range<=.01,occupancy<=.01))
}

result_row <- function(a,prior,phase) {
  group <- function(name)a$groups[a$groups$metric==name & a$groups$group=='all',,drop=FALSE]
  o <- group('occupancy');r <- group('range');s <- group('spatial_sd')
  stopifnot(nrow(o)==1L,nrow(r)==1L,nrow(s)==1L)
  data.frame(key=a$job$key,community=a$job$community,grid_index=a$job$grid_index,
    replicate=a$job$replicate,arm=a$job$arm,prior=prior,phase=phase,
    a$spatial$score,amplitude_mean=s$estimate,amplitude_bias=s$bias,
    amplitude_coverage=s$coverage,range_mean=r$estimate,range_mae=r$mae,
    range_coverage=r$coverage,occupancy_bias=o$bias,occupancy_mae=o$mae,
    occupancy_rmse=o$rmse,occupancy_coverage=o$coverage,
    spatial_flag_count=length(spatial_flags(a$spatial$diagnostics)),
    flag_count=length(a$reasons),reasons=paste(a$reasons,collapse='; '))
}

pair_rows <- function(rows) {
  by <- c('key','community','grid_index','replicate','arm')
  ig <- rows[rows$prior=='inverse_gamma',];hc <- rows[rows$prior=='half_cauchy',]
  stopifnot(!anyDuplicated(ig$key),!anyDuplicated(hc$key),setequal(ig$key,hc$key))
  merge(ig,hc,by=by,suffixes=c('_ig','_hc'),sort=TRUE)
}

paired_summary <- function(pairs) {
  metrics <- c('raw_rmse','centred_rmse','centred_correlation','centred_slope',
    'posterior_mean_rms','amplitude_mean','range_mae','range_coverage',
    'occupancy_bias','occupancy_mae','occupancy_rmse','occupancy_coverage')
  do.call(rbind,lapply(split(pairs,pairs$arm),function(a)do.call(rbind,lapply(metrics,function(m){
    ig <- a[[paste0(m,'_ig')]];hc <- a[[paste0(m,'_hc')]]
    data.frame(arm=a$arm[1],metric=m,inverse_gamma=mean(ig),half_cauchy=mean(hc),
      as.list(stratified_mean(hc-ig,a$grid_index)),decreased=sum(hc<ig),increased=sum(hc>ig))
  }))))
}

read_checked_result <- function(path) {
  if(!file.exists(path))stop('Required result missing: ',path)
  a <- readRDS(path)
  validate_spatial_diagnostics(a)
  stopifnot(file.exists(a$source_fit),
    identical(unname(tools::md5sum(a$source_fit)),a$source_fit_md5))
  if(!is.null(a$source_result))stopifnot(file.exists(a$source_result),
    identical(unname(tools::md5sum(a$source_result)),a$source_result_md5))
  # Recompute decision reasons so obsolete cached flags cannot control selection.
  expected <- unique(c(diagnostic_reasons(a),spatial_flags(a$spatial$diagnostics)))
  stopifnot(identical(a$reasons,expected))
  a
}

validate_choice <- function(a,initial,prior,phase,new_fit_hashes) {
  fields <- c('key','community','grid_index','replicate','arm','knots','input_md5')
  stopifnot(identical(a$job[fields],initial$job[fields]),a$job$knots==100L,
    identical(a$fit_hashes,if(prior=='inverse_gamma' && !is.null(a$source_result))
      initial$fit_hashes else new_fit_hashes))
  expected_priors <- if(prior=='half_cauchy')list(sigma_bs_prior=prior,sigma_bs_scale=1) else list()
  if(is.null(a$source_result))stopifnot(identical(a$priors,expected_priors))
  expected <- switch(phase,initial=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L),
    long=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L),
    'starts-initial'=list(nchain=4L,nburn=3000L,niter=5000L,nthin=1L),
    'starts-long'=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L))
  stopifnot(!is.null(expected),identical(a$mcmc,expected))
  if(startsWith(phase,'starts-'))stopifnot(identical(a$starts,c(.1,.3,1,3)),
    a$job$key=='range6-rep01-binary-k100') else stopifnot(is.null(a$starts))
  invisible(TRUE)
}
