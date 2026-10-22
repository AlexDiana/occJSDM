# Version 1: see ESTIMAND-AMENDMENT.md. Never interpret finite-chain spatial
# averages as finite posterior expectations under the unbounded half-Cauchy.
robust_trace_diagnostics <- function(x) {
  x<-as.matrix(x)
  stopifnot(all(is.finite(x)),nrow(x)>=4L,ncol(x)>=2L)
  constant<-any(apply(x,2,function(z)all(z==z[1])))
  call<-function(f)if(constant)NA_real_ else suppressWarnings(as.numeric(f(x)))
  quantile_ess<-function(p)call(function(z)posterior::ess_quantile(z,probs=p))
  quantile_mcse<-function(p)call(function(z)posterior::mcse_quantile(z,probs=p))
  c(rhat=call(posterior::rhat),ess_bulk=call(posterior::ess_bulk),
    ess_median=quantile_ess(.5),ess_q025=quantile_ess(.025),ess_q975=quantile_ess(.975),
    mcse_median=quantile_mcse(.5),mcse_q025=quantile_mcse(.025),mcse_q975=quantile_mcse(.975),
    chain_median_gap=diff(range(apply(x,2,median))),constant_chain=as.numeric(constant))
}

robust_spatial_flags <- function(d) {
  required<-c('quantity','rhat','ess_bulk','ess_median','ess_q025','ess_q975')
  stopifnot(is.data.frame(d),nrow(d)>0L,all(required %in% names(d)),!anyDuplicated(d$quantity))
  reasons<-character()
  for(i in seq_len(nrow(d))) {
    bad<-required[-1][!is.finite(unlist(d[i,required[-1]]))]
    if(length(bad))reasons<-c(reasons,paste(d$quantity[i],paste(bad,collapse='/'),'unavailable'))
    if(is.finite(d$rhat[i]) && d$rhat[i]>1.05)reasons<-c(reasons,paste(d$quantity[i],'Rhat > 1.05'))
    for(nm in required[-(1:2)])if(is.finite(d[[nm]][i]) && d[[nm]][i]<100)
      reasons<-c(reasons,paste(d$quantity[i],nm,'< 100'))
  }
  reasons
}

robust_flag_rows <- function(d) {
  names<-c('rhat','ess_bulk','ess_median','ess_q025','ess_q975')
  stopifnot(nrow(d)>0L,all(names %in% colnames(d)))
  x<-as.matrix(d[names])
  !apply(is.finite(x),1,all) | (!is.na(d$rhat) & d$rhat>1.05) |
    apply(x[,-1,drop=FALSE],1,function(z)any(z<100,na.rm=TRUE))
}

robust_field_draws <- function(coefficients,index,bases,truth,amplitude) {
  d<-dim(coefficients);stopifnot(length(d)==4L,length(dim(truth))==2L)
  ps<-d[1];S<-d[2];ni<-d[3];nc<-d[4];n<-nrow(truth)
  stopifnot(ncol(truth)==S,identical(dim(index),c(ni,nc)),identical(dim(amplitude),c(ni,nc)),
    all(is.finite(coefficients)),all(is.finite(amplitude)),all(index %in% seq_along(bases)),
    all(vapply(bases,function(b)identical(dim(b),c(n,ps)),logical(1))))
  tc<-sweep(truth,2,colMeans(truth),'-');denom<-colSums(tc^2);stopifnot(all(denom>0))
  labels<-c('amplitude',paste0('field_mean_',seq_len(S)),paste0('field_rms_',seq_len(S)),
    paste0('field_projection_',seq_len(S)))
  traces<-array(NA_real_,c(length(labels),ni,nc),dimnames=list(labels,NULL,NULL));traces[1,,]<-amplitude
  med<-lo<-hi<-matrix(NA_real_,n,S);field_diagnostics<-vector('list',n*S)
  # One species at a time bounds storage even for the longer schedule.
  for(s in seq_len(S)) {
    draws<-matrix(NA_real_,n,ni*nc)
    for(ch in seq_len(nc))for(g in sort(unique(index[,ch]))) {
      ii<-which(index[,ch]==g)
      for(first in seq.int(1L,length(ii),by=100L)) {
        its<-ii[first:min(length(ii),first+99L)];k<-length(its)
        field<-bases[[g]]%*%matrix(coefficients[,s,its,ch,drop=FALSE],ps,k)
        stopifnot(all(is.finite(field)))
        draws[,(ch-1L)*ni+its]<-field
        traces[1L+s,its,ch]<-colMeans(field)
        traces[1L+S+s,its,ch]<-sqrt(colMeans(field^2))
        traces[1L+2L*S+s,its,ch]<-colSums(field*tc[,s])/denom[s]
      }
    }
    for(site in seq_len(n)) {
      x<-matrix(draws[site,],ni,nc);q<-quantile(x,c(.025,.5,.975),names=FALSE)
      lo[site,s]<-q[1];med[site,s]<-q[2];hi[site,s]<-q[3]
      id<-(s-1L)*n+site
      field_diagnostics[[id]]<-data.frame(quantity=paste0('field_',id),site=site,species=s,
        truth=truth[site,s],median=q[2],bias=q[2]-truth[site,s],lower=q[1],upper=q[3],
        covered=q[1]<=truth[site,s] && truth[site,s]<=q[3],as.list(robust_trace_diagnostics(x)))
    }
  }
  fc<-sweep(med,2,colMeans(med),'-')
  score<-data.frame(raw_rmse=sqrt(mean((med-truth)^2)),raw_mae=mean(abs(med-truth)),
    centred_rmse=sqrt(mean((fc-tc)^2)),centred_mae=mean(abs(fc-tc)),
    centred_correlation=cor(c(fc),c(tc)),centred_slope=sum(fc*tc)/sum(tc^2),
    median_field_rms=sqrt(mean(med^2)),true_rms=sqrt(mean(truth^2)),
    raw_bias=mean(med-truth),raw_coverage=mean(lo<=truth & truth<=hi),
    raw_interval_width=mean(hi-lo))
  diagnostics<-do.call(rbind,lapply(seq_along(labels),function(i)
    data.frame(quantity=labels[i],as.list(robust_trace_diagnostics(matrix(traces[i,,],ni,nc))))))
  aq<-quantile(amplitude,c(.025,.5,.975),names=FALSE)
  list(field_median=med,field_lower=lo,field_upper=hi,score=score,traces=traces,
    diagnostics=diagnostics,field_diagnostics=do.call(rbind,field_diagnostics),
    amplitude=data.frame(lower=aq[1],median=aq[2],upper=aq[3]))
}

validate_robust_spatial <- function(a) {
  S<-nrow(a$species);n<-nrow(a$spatial$field_median)
  stopifnot(S>0L,n>0L,ncol(a$spatial$field_median)==S)
  expected<-c('amplitude',paste0('field_mean_',seq_len(S)),paste0('field_rms_',seq_len(S)),
    paste0('field_projection_',seq_len(S)))
  stopifnot(setequal(a$spatial$diagnostics$quantity,expected),
    nrow(a$spatial$diagnostics)==length(expected),
    identical(a$spatial$field_diagnostics$quantity,paste0('field_',seq_len(n*S))))
  robust_spatial_flags(a$spatial$diagnostics);robust_spatial_flags(a$spatial$field_diagnostics)
  invisible(TRUE)
}

robust_reasons <- function(a) {
  # All Rhat values in the archived scorer are already rank-normalized.
  # Retain its bounded occupancy ESS rules, excluding generic native warnings.
  copy<-a;copy$warnings<-character();reasons<-diagnostic_reasons(copy)
  for(nm in c('diagnostics','field_diagnostics')) {
    n<-sum(robust_flag_rows(a$spatial[[nm]]))
    if(n)reasons<-c(reasons,paste(n,if(nm=='diagnostics')'spatial summary traces' else
      'pointwise spatial fields','fail rank/quantile checks'))
  }
  unique(reasons)
}

robust_extension_gate <- function(pairs,selected,probe) {
  stopifnot(nrow(pairs)==9L,!anyDuplicated(pairs$community),
    identical(as.integer(table(pairs$grid_index)),c(3L,3L,3L)),setequal(pairs$grid_index,c(4,6,8)))
  expected<-as.vector(outer(pairs$key,c('inverse_gamma','half_cauchy'),paste,sep=':'))
  stopifnot(length(selected)==18L,!anyDuplicated(names(selected)),setequal(names(selected),expected))
  invisible(lapply(selected,validate_robust_spatial));validate_robust_spatial(probe)
  failed<-function(a)any(robust_flag_rows(a$spatial$diagnostics)) ||
    any(robust_flag_rows(a$spatial$field_diagnostics))
  change<-function(m)mean(pairs[[paste0(m,'_hc')]]-pairs[[paste0(m,'_ig')]])
  raw<-change('raw_rmse');centred<-change('centred_rmse');correlation<-change('centred_correlation')
  range<-change('range_mae');occupancy<-change('occupancy_mae')
  improved<-sum(pairs$centred_rmse_hc<pairs$centred_rmse_ig)
  spatial<-sum(vapply(selected,failed,logical(1)));start<-as.integer(failed(probe))
  stopifnot(all(is.finite(c(raw,centred,correlation,range,occupancy))))
  data.frame(criterion=c('Lower mean raw median-field RMSE','Lower mean centred median-field RMSE',
    'Higher mean centred median-field correlation','Centred RMSE improved in at least six communities',
    'No unresolved spatial rank/quantile diagnostics in either prior','Initialization rank/quantile diagnostics pass',
    'Range MAE increase at most 0.01','Occupancy MAE increase at most 0.01'),
    value=c(raw,centred,correlation,improved,spatial,start,range,occupancy),
    pass=c(raw<0,centred<0,correlation>0,improved>=6,spatial==0,start==0,range<=.01,occupancy<=.01))
}

robust_result_row <- function(a,prior,phase) {
  group<-function(name)a$groups[a$groups$metric==name & a$groups$group=='all',,drop=FALSE]
  o<-group('occupancy');r<-group('range');stopifnot(nrow(o)==1L,nrow(r)==1L)
  amp<-a$spatial$amplitude
  trace_flags<-sum(robust_flag_rows(a$spatial$diagnostics));field_flags<-sum(robust_flag_rows(a$spatial$field_diagnostics))
  data.frame(key=a$job$key,community=a$job$community,grid_index=a$job$grid_index,
    replicate=a$job$replicate,arm=a$job$arm,prior=prior,phase=phase,a$spatial$score,
    amplitude_median=amp$median,amplitude_bias=amp$median-1,amplitude_lower=amp$lower,
    amplitude_upper=amp$upper,amplitude_coverage=amp$lower<=1 && amp$upper>=1,
    range_mean=r$estimate,range_mae=r$mae,range_coverage=r$coverage,
    occupancy_bias=o$bias,occupancy_mae=o$mae,occupancy_rmse=o$rmse,occupancy_coverage=o$coverage,
    trace_flag_count=trace_flags,field_flag_count=field_flags,spatial_flag_count=trace_flags+field_flags,
    flag_count=length(a$reasons),reasons=paste(a$reasons,collapse='; '),native_warning_count=length(a$warnings))
}

robust_paired_summary <- function(pairs) {
  metrics<-c('raw_rmse','centred_rmse','centred_correlation','centred_slope','median_field_rms',
    'raw_bias','raw_coverage','raw_interval_width','amplitude_median','amplitude_coverage',
    'range_mae','range_coverage','occupancy_bias','occupancy_mae','occupancy_rmse','occupancy_coverage')
  do.call(rbind,lapply(split(pairs,pairs$arm),function(a)do.call(rbind,lapply(metrics,function(m){
    ig<-a[[paste0(m,'_ig')]];hc<-a[[paste0(m,'_hc')]]
    data.frame(arm=a$arm[1],metric=m,inverse_gamma=mean(ig),half_cauchy=mean(hc),
      as.list(stratified_mean(hc-ig,a$grid_index)),decreased=sum(hc<ig),increased=sum(hc>ig))
  }))))
}
