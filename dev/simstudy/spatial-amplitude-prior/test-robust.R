library(testthat)
source('../spatial-targeted-recheck/analysis.R')
source('metrics.R')
source('analysis.R')
source('robust.R')

test_that('field medians, quantiles and centring use every correctly indexed draw', {
  set.seed(824)
  n<-3L;ps<-2L;S<-2L;ni<-200L;nc<-2L
  co<-array(rexp(ps*S*ni*nc),c(ps,S,ni,nc))
  idx<-matrix(sample(1:2,ni*nc,TRUE),ni,nc)
  bases<-list(matrix(c(1,0,2,0,1,-1),n,ps),matrix(c(.3,1,-1,2,0,1),n,ps))
  truth<-matrix(c(-1,0,1,2,-1,0),n,S)
  amp<-matrix(rexp(ni*nc),ni,nc)
  got<-robust_field_draws(co,idx,bases,truth,amp)
  reference<-array(NA_real_,c(n,S,ni,nc))
  for(ch in 1:nc)for(it in 1:ni)reference[,,it,ch]<-bases[[idx[it,ch]]]%*%co[,,it,ch]
  med<-apply(reference,1:2,median)
  lo<-apply(reference,1:2,quantile,.025)
  hi<-apply(reference,1:2,quantile,.975)
  expect_equal(got$field_median,med,tolerance=1e-12)
  expect_equal(got$field_lower,lo,tolerance=1e-12)
  expect_equal(got$field_upper,hi,tolerance=1e-12)
  centred<-sweep(med,2,colMeans(med),'-')
  tc<-sweep(truth,2,colMeans(truth),'-')
  expect_equal(got$score$centred_rmse,sqrt(mean((centred-tc)^2)))
  expect_equal(got$score$raw_rmse,sqrt(mean((med-truth)^2)))
  expect_equal(got$score$raw_coverage,mean(lo<=truth & truth<=hi))
  expect_equal(got$amplitude$median,median(amp))
  expect_gt(max(abs(med-apply(reference,1:2,mean))),.1)
  x<-matrix(reference[2,1,,],ni,nc)
  expect_equal(got$field_diagnostics$ess_median[2],unname(posterior::ess_quantile(x,probs=.5)))
  expect_false(any(c('ess_mean','mcse_mean') %in% names(got$diagnostics)))
})

test_that('rank and quantile screening cannot silently accept missing or poor mixing', {
  d<-data.frame(quantity='amplitude',rhat=1.01,ess_bulk=200,ess_median=200,ess_q025=200,ess_q975=200)
  expect_length(robust_spatial_flags(d),0L)
  for(name in c('ess_bulk','ess_median','ess_q025','ess_q975')) {
    b<-d;b[[name]]<-99
    expect_true(any(grepl(name,robust_spatial_flags(b),fixed=TRUE)))
  }
  b<-d;b$rhat<-1.051;expect_true(any(grepl('Rhat',robust_spatial_flags(b))))
  b<-d;b$ess_median<-NA_real_;expect_true(any(grepl('unavailable',robust_spatial_flags(b))))
  expect_error(robust_spatial_flags(d[FALSE,]))
})

test_that('the amended gate uses pointwise and trace diagnostics', {
  pairs<-data.frame(key=paste0('c',1:9),community=paste0('c',1:9),grid_index=rep(c(4,6,8),each=3),
    raw_rmse_ig=1,raw_rmse_hc=.9,centred_rmse_ig=1,centred_rmse_hc=.8,
    centred_correlation_ig=.2,centred_correlation_hc=.3,
    range_mae_ig=.02,range_mae_hc=.02,occupancy_mae_ig=.1,occupancy_mae_hc=.1)
  d<-data.frame(quantity=c('amplitude',paste0('field_mean_',1:2),paste0('field_rms_',1:2),
    paste0('field_projection_',1:2)),rhat=1.01,ess_bulk=200,ess_median=200,ess_q025=200,ess_q975=200)
  fd<-d[rep(1,6),];fd$quantity<-paste0('field_',1:6)
  a<-list(spatial=list(diagnostics=d,field_diagnostics=fd,field_median=matrix(0,3,2)),
    species=data.frame(species=1:2))
  selected<-rep(list(a),18)
  names(selected)<-as.vector(outer(pairs$key,c('inverse_gamma','half_cauchy'),paste,sep=':'))
  expect_true(all(robust_extension_gate(pairs,selected,a)$pass))
  selected[[1]]$spatial$field_diagnostics$ess_median[1]<-90
  expect_false(all(robust_extension_gate(pairs,selected,a)$pass))
  expect_error(robust_extension_gate(pairs,selected[-1],a))
})
