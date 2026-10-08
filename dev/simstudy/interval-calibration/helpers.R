# Research scoring only. Production fitting functions and priors are unchanged.
interval_rows <- function(arr, truth, block, scenario, replicate) {
  dd <- dim(arr); ni <- dd[length(dd)-1L]; nc <- dd[length(dd)]
  stopifnot(length(dd)>=2L,prod(head(dd,-2L))==length(truth),all(is.finite(arr)),all(is.finite(truth)))
  xx <- array(arr,c(length(truth),ni,nc))
  do.call(rbind,lapply(seq_along(truth),function(k) {
    z <- matrix(xx[k,,],ni,nc); qq <- quantile(z,c(.025,.975),names=FALSE)
    data.frame(scenario=scenario,replicate=replicate,block=block,element=k,
      truth=as.vector(truth)[k],post_mean=mean(z),lower=qq[1],upper=qq[2],
      post_sd=sd(as.vector(z)),rhat=posterior::rhat(z),ess_bulk=posterior::ess_bulk(z),
      ess_tail=posterior::ess_tail(z),mcse=posterior::mcse_mean(z))
  }))
}

summarise_intervals <- function(rows) {
  rows$covered <- as.integer(rows$truth>=rows$lower & rows$truth<=rows$upper)
  rows$bias <- rows$post_mean-rows$truth
  rows$width <- rows$upper-rows$lower
  rows$squared_error <- rows$bias^2
  rows$lower_miss <- as.integer(rows$truth<rows$lower)
  rows$upper_miss <- as.integer(rows$truth>rows$upper)
  if(!'group' %in% names(rows)) rows$group <- 'all'
  variables <- c('covered','bias','width','squared_error','lower_miss','upper_miss')
  rr <- aggregate(rows[variables],rows[c('scenario','block','group','replicate')],mean)
  groups <- split(rr,interaction(rr$scenario,rr$block,rr$group,drop=TRUE))
  result <- do.call(rbind,lapply(groups,function(z) {
    n <- nrow(z); se <- sd(z$covered)/sqrt(n); bs <- sd(z$bias)/sqrt(n)
    crit <- if(n>1L) qt(.975,n-1L) else NA_real_
    data.frame(scenario=z$scenario[1],block=z$block[1],group=z$group[1],communities=n,
      coverage=mean(z$covered),coverage_mcse=se,
      coverage_mc_lower=max(0,mean(z$covered)-crit*se),coverage_mc_upper=min(1,mean(z$covered)+crit*se),
      bias=mean(z$bias),bias_mcse=bs,bias_mc_lower=mean(z$bias)-crit*bs,bias_mc_upper=mean(z$bias)+crit*bs,
      width=mean(z$width),rmse=sqrt(mean(z$squared_error)),
      lower_miss=mean(z$lower_miss),upper_miss=mean(z$upper_miss))
  }))
  rownames(result)<-NULL
  list(rows=rows,replicates=rr,summary=result)
}

write_interval_tables <- function(rows,prefix,out) {
  result<-summarise_intervals(rows)
  for(nm in names(result))write.csv(result[[nm]],file.path(out,paste0(prefix,'-',nm,'.csv')),row.names=FALSE)
  invisible(result)
}
