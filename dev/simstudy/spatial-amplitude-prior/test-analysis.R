library(testthat)
source('../spatial-targeted-recheck/analysis.R')
source('metrics.R')
source('analysis.R')

test_that('a problem in either prior selects matched longer runs', {
  a <- list(reasons=character());b<-list(reasons='amplitude ESS < 100')
  expect_false(pair_needs_long(a,a))
  expect_true(pair_needs_long(a,b));expect_true(pair_needs_long(b,a))
})

test_that('extension requires field improvement, diagnostics and ancillary safeguards', {
  pairs <- data.frame(key=paste0('c',1:9),community=paste0('c',1:9),grid_index=rep(c(4,6,8),each=3),
    raw_rmse_ig=1,raw_rmse_hc=.9,centred_rmse_ig=.9,centred_rmse_hc=.8,
    centred_correlation_ig=.4,centred_correlation_hc=.5,
    range_mae_ig=.03,range_mae_hc=.03,occupancy_mae_ig=.1,occupancy_mae_hc=.1)
  good <- data.frame(quantity=c('amplitude',paste0('field_mean_',1:8),
    paste0('field_rms_',1:8),paste0('field_projection_',1:8)),rhat=1.01,ess_mean=200)
  records <- rep(list(list(spatial=list(diagnostics=good),species=data.frame(species=1:8))),18)
  names(records)<-as.vector(outer(pairs$key,c('inverse_gamma','half_cauchy'),paste,sep=':'))
  probe <- list(spatial=list(diagnostics=good),species=data.frame(species=1:8))
  expect_true(all(extension_gate(pairs,records,probe)$pass))
  pairs$raw_rmse_hc <- 1.1
  expect_false(all(extension_gate(pairs,records,probe)$pass))
  pairs$raw_rmse_hc <- .9
  pairs$centred_rmse_hc[1:4] <- .91
  expect_false(all(extension_gate(pairs,records,probe)$pass))
  pairs$centred_rmse_hc <- .8
  pairs$range_mae_hc <- .041
  expect_false(all(extension_gate(pairs,records,probe)$pass))
  pairs$range_mae_hc <- .03;pairs$occupancy_mae_hc <- .111
  expect_false(all(extension_gate(pairs,records,probe)$pass))
  pairs$occupancy_mae_hc <- .1
  probe$spatial$diagnostics$ess_mean <- 50
  expect_false(all(extension_gate(pairs,records,probe)$pass))
  probe$spatial$diagnostics <- good;records[[1]]$spatial$diagnostics$rhat <- 1.06
  expect_false(all(extension_gate(pairs,records,probe)$pass))
  expect_error(extension_gate(pairs[-1,],records,probe))
  expect_error(extension_gate(pairs,records[-1],probe))
  records[[1]]$spatial$diagnostics<-good[FALSE,]
  expect_error(extension_gate(pairs,records,probe))
})

test_that('selected longer results cannot silently change prior scale or fitting code', {
  job<-list(key='range4-rep01-binary-k100',community='range4-rep01',grid_index=4L,
    replicate=1L,arm='binary',knots=100L,input_md5='same-input')
  a<-list(job=job,fit_hashes=c(code='frozen'),priors=list(sigma_bs_prior='half_cauchy',sigma_bs_scale=1),
    mcmc=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L),starts=NULL)
  initial<-a
  expect_true(validate_choice(a,initial,'half_cauchy','long',c(code='frozen')))
  a$priors$sigma_bs_scale<-2
  expect_error(validate_choice(a,initial,'half_cauchy','long',c(code='frozen')))
  a$priors$sigma_bs_scale<-1;a$fit_hashes<-c(code='different')
  expect_error(validate_choice(a,initial,'half_cauchy','long',c(code='frozen')))
  a$fit_hashes<-c(code='frozen');a$job$key<-'wrong'
  expect_error(validate_choice(a,initial,'half_cauchy','long',c(code='frozen')))
})
