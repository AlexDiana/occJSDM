# Research-only sensitivity; no production namespace or default is changed.
weighted_quantile <- function(x,w,probs=c(.025,.975)) {
  stopifnot(length(x)==length(w),all(is.finite(x)),all(is.finite(w)),all(w>=0),sum(w)>0)
  ord<-order(x);cum<-cumsum(w[ord])/sum(w)
  vapply(probs,function(p)x[ord[which(cum>=p)[1L]]],numeric(1))
}

draw_matrix <- function(x) {
  dd<-dim(x);stopifnot(length(dd)>=2L)
  matrix(x,nrow=prod(head(dd,-2L)))
}

joint_log_ratio <- function(beta,theta,collection_mean=0,theta_b=20) {
  stopifnot(length(dim(beta))==4L,dim(beta)[1L]==2L,length(dim(theta))==3L,
    identical(dim(beta)[-1L],dim(theta)),all(theta>0 & theta<1),
    collection_mean %in% c(0,1),theta_b %in% c(20,30))
  slopes<-draw_matrix(beta[2,,,,drop=FALSE]);t<-draw_matrix(theta)
  # Sum over every species at each joint draw; intercept prior is unchanged.
  collection_mean*colSums(slopes)/2 - nrow(slopes)*collection_mean^2/4 +
    (theta_b-20)*colSums(log1p(-t)) + nrow(t)*log(theta_b/20)
}

weight_info <- function(logw,ni,nc) {
  w<-exp(logw-max(logw));w<-w/sum(w)
  if(diff(range(logw))<1e-12) {
    k<-NA_real_; reff<-1; psis_ess<-length(w)
  }else {
    inv<-exp(-logw-max(-logw))
    reff<-as.numeric(loo::relative_eff(inv,chain_id=rep(seq_len(nc),each=ni)))
    stopifnot(is.finite(reff),reff>0)
    ps<-loo::psis(logw,r_eff=reff)
    k<-as.numeric(loo::pareto_k_values(ps));psis_ess<-as.numeric(loo::psis_n_eff_values(ps))
  }
  list(weights=w,pareto_k=k,weight_ess=1/sum(w^2),psis_ess=psis_ess,r_eff=reff,
    chain_weight=colSums(matrix(w,ni,nc)),weight_flag=(!is.na(k) && k>.7) || psis_ess<400)
}

block_diagnostics <- function(arr,block) {
  dd<-dim(arr);ni<-dd[length(dd)-1L];nc<-dd[length(dd)]
  xx<-array(arr,c(prod(head(dd,-2L)),ni,nc))
  do.call(rbind,lapply(seq_len(dim(xx)[1L]),function(k){
    z<-matrix(xx[k,,],ni,nc)
    data.frame(block=block,element=k,rhat=posterior::rhat(z),
      ess_bulk=posterior::ess_bulk(z),ess_tail=posterior::ess_tail(z),mcse=posterior::mcse_mean(z))
  }))
}

source_diagnostics <- function(fit,warnings) {
  ro<-fit$results_output
  z<-do.call(rbind,lapply(c('theta0','beta_theta','p','q'),function(nm)
    block_diagnostics(ro[[paste0(nm,'_output')]],nm)))
  z$flag<-!is.finite(z$rhat) | z$rhat>1.05 | !is.finite(z$ess_bulk) | z$ess_bulk<400 |
    !is.finite(z$ess_tail) | z$ess_tail<400
  list(elements=z,flag=any(z$flag) || length(warnings)>0L,
    max_rhat=max(z$rhat),min_ess=min(z$ess_bulk,z$ess_tail),warnings=length(warnings))
}

score_theta <- function(theta,truth,w) {
  z<-draw_matrix(theta);stopifnot(nrow(z)==length(truth),ncol(z)==length(w))
  do.call(rbind,lapply(seq_len(nrow(z)),function(j){
    q<-weighted_quantile(z[j,],w)
    data.frame(species=j,truth=truth[j],post_mean=sum(z[j,]*w),lower=q[1],upper=q[2],
      covered=as.integer(truth[j]>=q[1] && truth[j]<=q[2]))
  }))
}

mean_fitter <- function(collection_mean) {
  stopifnot(collection_mean %in% c(0,1));original<-occJSDM::runOccJSDM
  target<-quote(b_betatheta <- rep(0,ncov_theta)); count<-0L
  replace<-function(e) {
    if(missing(e))return(quote(expr=))
    if(identical(e,target)) {count<<-count+1L;return(substitute(b_betatheta<-rep(M,ncov_theta),list(M=collection_mean)))}
    if(is.call(e))for(i in seq_along(e))e[i]<-lapply(e[i],replace)
    e
  }
  f<-original;body(f)<-replace(body(original));stopifnot(count==1L,identical(formals(f),formals(original)),
    identical(environment(f),environment(original)))
  if(collection_mean==0)stopifnot(identical(body(f),body(original)))
  f
}

fit_input <- function(input,fitter,mcmc,theta_b=20) {
  s<-input$scenario
  fitter(input$sim$data_list,listParams=list(n_factors=s$d,n_lattrait=s$gt),
    listPriors=list(a_theta0=1,b_theta0=theta_b,a_p=5,b_p=1,a_q=1,b_q=20,b_betatheta_slope_var=2),
    threshold=1,occCovariates=paste0('X_psi.EnvCov.',seq_len(s$ncov_psi)),
    collCovariates='X_theta',spatCovariates=NULL,MCMCparams=mcmc)
}

summaries <- function(rows) {
  rows$bias<-rows$post_mean-rows$truth;rows$abs_error<-abs(rows$bias);rows$width<-rows$upper-rows$lower
  rows$lower_miss<-as.integer(rows$truth<rows$lower);rows$upper_miss<-as.integer(rows$truth>rows$upper)
  vals<-c('covered','bias','abs_error','width','lower_miss','upper_miss')
  rr<-aggregate(rows[vals],rows[c('contamination','design','replicate','collection_mean','theta_b')],mean)
  gg<-split(rr,interaction(rr$contamination,rr$design,rr$collection_mean,rr$theta_b,drop=TRUE))
  ss<-do.call(rbind,lapply(gg,function(z){
    n<-nrow(z); data.frame(z[1,c('contamination','design','collection_mean','theta_b')],
      communities=n,coverage=mean(z$covered),coverage_mcse=sd(z$covered)/sqrt(n),
      bias=mean(z$bias),bias_mcse=sd(z$bias)/sqrt(n),abs_error=mean(z$abs_error),
      width=mean(z$width),lower_miss=mean(z$lower_miss),upper_miss=mean(z$upper_miss))
  }))
  rownames(rr)<-rownames(ss)<-NULL
  list(rows=rows,communities=rr,summary=ss)
}

save_tables <- function(rows,prefix,out) {
  z<-summaries(rows)
  for(nm in names(z))write.csv(z[[nm]],file.path(out,paste0(prefix,'-',nm,'.csv')),row.names=FALSE)
  invisible(z)
}
