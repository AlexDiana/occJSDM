library(testthat)
source('analysis.R')

test_that('an initial-budget result cannot masquerade as a longer check', {
  row <- data.frame(key='case',community='range4-rep01',grid_index=4L,
    replicate=1L,arm='low',knots=20L,input_md5='abc')
  a <- list(job=as.list(row),mcmc=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L))
  expect_no_error(validate_result_protocol(a,row,'initial'))
  expect_error(validate_result_protocol(a,row,'long'))
  a$mcmc <- list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L)
  expect_no_error(validate_result_protocol(a,row,'long'))
  a$job$arm <- 'binary'
  expect_error(validate_result_protocol(a,row,'long'))
})

test_that('community changes that cancel remain visible in substitution sensitivity', {
  a <- data.frame(community=c('a','b'),grid_index=c(4,4),arm='low',knots=20L,
    metric='occupancy',group='low',bias=c(.1,.1),mae=c(.1,.1),rmse=c(.1,.1),coverage=c(.8,.8),interval_width=c(.2,.2))
  b <- a;b$bias <- c(.2,0)
  d <- community_sensitivity(a,b)
  expect_equal(d$bias_change,c(.1,-.1))
  expect_equal(mean(d$bias_change),0,tolerance=1e-12)
  expect_equal(max(abs(d$bias_change)),.1)
})

test_that('range-stratified uncertainty treats communities as independent units', {
  d <- data.frame(grid_index=rep(c(4,6,8),each=3),value=c(1,2,3,11,12,13,21,22,23))
  a <- stratified_mean(d$value,d$grid_index)
  expect_equal(a['mean'],c(mean=12))
  expect_equal(a['se'],c(se=1/3))
  expect_equal(a['df'],c(df=6))
  expect_equal(a['n'],c(n=9))
  # A missing range cannot silently turn an incomplete batch into a full summary.
  expect_error(stratified_mean(d$value[1:6],d$grid_index[1:6]))
})

test_that('long-run selection follows diagnostics and ignores scientific direction', {
  clean <- list(job=list(key='case'),warnings=character(),
    groups=data.frame(metric=c('occupancy','range'),rhat=c(1,NA),ess_mean=c(500,NA)),
    species=data.frame(rhat=1,ess_mean=500),
    elements=data.frame(metric=c('occupancy','range'),rhat=c(1,NA)))
  expect_length(diagnostic_reasons(clean),0)
  a<-clean;a$groups$rhat[1]<-1.051
  expect_true('group Rhat > 1.05' %in% diagnostic_reasons(a))
  a<-clean;a$species$ess_mean<-99
  expect_true('occupancy species ESS < 100' %in% diagnostic_reasons(a))
  a<-clean;a$elements$rhat[1]<-NA_real_
  expect_true('non-range element Rhat unavailable' %in% diagnostic_reasons(a))
  a<-clean;a$warnings<-'Low ESS in Bs'
  expect_true('native convergence warning' %in% diagnostic_reasons(a))
  a<-clean;a$groups$rhat[1]<-1.05;a$species$ess_mean<-100
  expect_length(diagnostic_reasons(a),0)
})

test_that('paired contrasts join shuffled communities and retain uncertainty of differences', {
  d <- expand.grid(grid_index=c(4L,6L,8L),replicate=1:3,
    arm=c('binary','low','high'),knots=c(20L,50L,100L),stringsAsFactors=FALSE)
  d$community <- paste(d$grid_index,d$replicate,sep='-')
  d$metric <- 'occupancy';d$group <- 'all'
  d$bias <- .01*d$replicate+(.01*d$replicate)*match(d$knots,c(20L,50L,100L))+
    c(binary=0,low=.1,high=.2)[d$arm]
  d$mae <- d$bias;d$rmse <- d$bias;d$coverage <- .9
  set.seed(2);d <- d[sample(nrow(d)),]
  a <- paired_effects(d)
  selected <- a[a$comparison=='low_k50_minus_k20' & a$quantity=='mae',]
  expect_equal(selected$mean,.02)
  expect_equal(selected$se,1/300)
  expect_equal(selected$df,6)
  selected <- a[a$comparison=='high_minus_binary_k100' & a$quantity=='bias',]
  expect_equal(selected$mean,.2)
  expect_equal(selected$increased,9)
  expect_error(paired_effects(d[-1,]))
})
