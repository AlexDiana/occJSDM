# The three spatial ranges are fixed design strata. Replicate communities,
# rather than species or site/species cells, determine Monte Carlo uncertainty.
validate_result_protocol <- function(result,manifest_row,phase) {
  expected <- switch(phase,initial=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L),
    long=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L))
  stopifnot(!is.null(expected),identical(result$mcmc,expected),nrow(manifest_row)==1L)
  fields <- c('key','community','grid_index','replicate','arm','knots','input_md5')
  stopifnot(all(vapply(fields,function(f)identical(result$job[[f]],manifest_row[[f]]),logical(1))))
  invisible(TRUE)
}

community_sensitivity <- function(initial,selected) {
  keys <- c('community','grid_index','arm','knots','metric','group')
  stopifnot(!anyDuplicated(initial[keys]),!anyDuplicated(selected[keys]))
  joined <- merge(initial,selected,by=keys,suffixes=c('_initial','_selected'))
  stopifnot(nrow(joined)==nrow(initial),nrow(joined)==nrow(selected))
  for(q in c('bias','mae','rmse','coverage','interval_width'))
    joined[[paste0(q,'_change')]] <- joined[[paste0(q,'_selected')]]-joined[[paste0(q,'_initial')]]
  joined
}

stratified_mean <- function(value,grid) {
  stopifnot(length(value)==length(grid),all(is.finite(value)),
    setequal(unique(grid),c(4,6,8)))
  cells <- split(value,grid); sizes <- lengths(cells)
  stopifnot(all(sizes>=2L))
  estimate <- mean(vapply(cells,mean,numeric(1)))
  variance_components <- vapply(cells,var,numeric(1))/sizes/length(cells)^2
  variance <- sum(variance_components)
  df <- if(variance>0) variance^2/sum(variance_components^2/(sizes-1)) else Inf
  se <- sqrt(variance); half <- if(se>0) qt(.975,df)*se else 0
  c(mean=estimate,se=se,lower=estimate-half,upper=estimate+half,
    df=df,n=length(value),minimum=min(value),maximum=max(value))
}

diagnostic_reasons <- function(a) {
  reasons <- character()
  if(any(a$groups$rhat>1.05,na.rm=TRUE)) reasons <- c(reasons,'group Rhat > 1.05')
  if(any(a$species$rhat>1.05,na.rm=TRUE)) reasons <- c(reasons,'species Rhat > 1.05')
  if(any(a$elements$rhat>1.05,na.rm=TRUE)) reasons <- c(reasons,'element Rhat > 1.05')
  primary <- a$groups$metric=='occupancy'
  if(any(a$groups$ess_mean[primary]<100,na.rm=TRUE)) reasons <- c(reasons,'occupancy group ESS < 100')
  if(any(a$species$ess_mean<100,na.rm=TRUE)) reasons <- c(reasons,'occupancy species ESS < 100')
  if(any(!is.finite(a$elements$rhat[a$elements$metric!='range'])))
    reasons <- c(reasons,'non-range element Rhat unavailable')
  if(any(grepl('Convergence|Low ESS',a$warnings))) reasons <- c(reasons,'native convergence warning')
  reasons
}

read_study_selection <- function(study,final=FALSE) {
  manifest <- readRDS(file.path(study,'input-manifest.rds'))
  stopifnot(nrow(manifest)==81L)
  selected <- initial <- list(); rows <- list()
  for(key in manifest$key) {
    initial_file <- file.path(study,'initial',paste0(key,'-result.rds'))
    if(!file.exists(initial_file)) stop('Initial batch incomplete: ',key)
    a <- readRDS(initial_file)
    validate_result_protocol(a,manifest[manifest$key==key,,drop=FALSE],'initial')
    initial[[key]] <- a
    reasons <- diagnostic_reasons(a)
    needs_long <- length(reasons)>0L
    long_file <- file.path(study,'long',paste0(key,'-result.rds'))
    b <- a;selected_file <- initial_file;phase <- 'initial'
    if(final && needs_long) {
      if(!file.exists(long_file)) stop('Prespecified longer check missing: ',key)
      b <- readRDS(long_file);selected_file <- long_file;phase <- 'long'
      validate_result_protocol(b,manifest[manifest$key==key,,drop=FALSE],'long')
      stopifnot(identical(b$job,a$job),identical(b$fit_hashes,a$fit_hashes),
        identical(b$score_hash,a$score_hash))
    }
    selected[[key]] <- b
    rows[[key]] <- data.frame(key=key,community=a$job$community,arm=a$job$arm,
      knots=a$job$knots,grid_index=a$job$grid_index,replicate=a$job$replicate,
      needs_long=needs_long,initial_reasons=paste(reasons,collapse='; '),
      phase=phase,selected_file=selected_file,
      selected_md5=unname(tools::md5sum(selected_file)),
      selected_reasons=paste(diagnostic_reasons(b),collapse='; '),
      max_group_rhat=max(b$groups$rhat,na.rm=TRUE),
      max_element_rhat=max(b$elements$rhat,na.rm=TRUE),
      min_occupancy_ess=min(b$groups$ess_mean[b$groups$metric=='occupancy'],na.rm=TRUE),
      warning_count=length(b$warnings),reconstruction=b$reconstruction)
  }
  list(selected=selected,initial=initial,manifest=do.call(rbind,rows))
}

collect_study_rows <- function(results,table) {
  do.call(rbind,lapply(results,function(a) {
    prefix <- data.frame(key=a$job$key,community=a$job$community,
      grid_index=a$job$grid_index,replicate=a$job$replicate,arm=a$job$arm,knots=a$job$knots)
    cbind(prefix[rep(1,nrow(a[[table]])),,drop=FALSE],a[[table]])
  }))
}

aggregate_scores <- function(groups) {
  values <- c('bias','mae','rmse','coverage','interval_width')
  by <- split(groups,interaction(groups$arm,groups$knots,groups$metric,groups$group,drop=TRUE))
  do.call(rbind,lapply(by,function(a) {
    stopifnot(nrow(a)==9L,!anyDuplicated(a$community))
    do.call(rbind,lapply(values,function(value) data.frame(
      a[1,c('arm','knots','metric','group')],quantity=value,
      as.list(stratified_mean(a[[value]],a$grid_index)))))
  }))
}

paired_effects <- function(groups) {
  a <- groups[groups$metric=='occupancy',]
  comparisons <- list()
  for(arm in c('binary','low','high')) for(k in c(50L,100L))
    comparisons[[paste0(arm,'_k',k,'_minus_k20')]] <- list(
      a[a$arm==arm & a$knots==k,],a[a$arm==arm & a$knots==20L,])
  for(k in c(20L,50L,100L)) for(arm in c('low','high'))
    comparisons[[paste0(arm,'_minus_binary_k',k)]] <- list(
      a[a$arm==arm & a$knots==k,],a[a$arm=='binary' & a$knots==k,])
  out <- list()
  for(name in names(comparisons)) {
    pair <- comparisons[[name]]
    join <- merge(pair[[1]],pair[[2]],by=c('community','grid_index','replicate','group'),suffixes=c('_a','_b'))
    stopifnot(nrow(join)==nrow(pair[[1]]),nrow(join)==nrow(pair[[2]]))
    for(g in unique(join$group)) for(value in c('bias','mae','rmse','coverage')) {
      d <- join[join$group==g,]
      stopifnot(nrow(d)==9L,!anyDuplicated(d$community))
      change <- d[[paste0(value,'_a')]]-d[[paste0(value,'_b')]]
      out[[length(out)+1L]] <- data.frame(comparison=name,group=g,quantity=value,
        as.list(stratified_mean(change,d$grid_index)),
        decreased=sum(change<0),increased=sum(change>0))
    }
  }
  do.call(rbind,out)
}
