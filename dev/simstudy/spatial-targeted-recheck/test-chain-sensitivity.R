library(testthat)
source('chain-sensitivity.R')

test_that('chain scoring preserves site/species order and rare prevalence groups', {
  target <- rep(c(.01,.05,.25,.75),each=2)
  truth <- outer(c(.5,1,1.5),target)
  truth[,7:8] <- c(.5,.85,.9)
  delta <- rep(seq_len(8)*.0001,each=3)
  result <- list(probability=truth,chain_probability=cbind(as.vector(truth)-delta,as.vector(truth)+delta),
    elements=data.frame(metric='occupancy',element=seq_len(24),truth=as.vector(truth)),
    species=data.frame(target=target))
  d <- score_observed_probability_chains(result)
  expect_equal(nrow(d),16L)
  expect_equal(d$bias[d$group=='all'],c(-.00045,.00045))
  expect_equal(d$mae[d$group=='prevalence_1pct'],rep(.00015,2))
  expect_equal(d$mae[d$group=='prevalence_5pct'],rep(.00035,2))
  result$elements <- result$elements[24:1,]
  expect_error(score_observed_probability_chains(result))
})

test_that('chain-only ranges retain every community and allow pooled MAE below them', {
  d <- expand.grid(grid_index=c(4,6,8),replicate=1:3,chain=1:2,
    arm='binary',knots=20L,group='all',stringsAsFactors=FALSE)
  d$community <- paste(d$grid_index,d$replicate,sep='-')
  d$bias <- ifelse(d$chain==1,-.1,.1);d$mae <- .1;d$rmse <- .1
  pooled <- unique(d[c('community','grid_index','replicate','arm','knots','group')])
  pooled$bias <- pooled$mae <- pooled$rmse <- 0
  a <- summarise_observed_chains(d,pooled)
  expect_equal(a$pooled,rep(0,3))
  expect_equal(a$chain_min_mean[a$quantity=='bias'],-.1)
  expect_equal(a$chain_max_mean[a$quantity=='bias'],.1)
  expect_equal(a$chain_min_mean[a$quantity=='mae'],.1)
  expect_equal(a$n,rep(9L,3))
  expect_error(summarise_observed_chains(d[d$community!='4-1',],pooled))
})

test_that('paired envelopes use matched communities and every observed chain', {
  d <- expand.grid(grid_index=c(4,6,8),replicate=1:3,chain=1:2,
    arm=c('binary','low','high'),knots=c(20L,50L,100L),group='all',stringsAsFactors=FALSE)
  d$community <- paste(d$grid_index,d$replicate,sep='-')
  d$bias <- d$replicate+.01*d$chain+ifelse(d$knots==20,0,.2)+ifelse(d$arm=='binary',0,.3)
  d$mae <- d$rmse <- d$bias
  pooled <- aggregate(d[c('bias','mae','rmse')],d[c('community','grid_index','replicate','arm','knots','group')],mean)
  set.seed(3);d <- d[sample(nrow(d)),]
  a <- paired_observed_chains(d,pooled)
  support <- a[a$comparison=='binary_k100_minus_k20' & a$quantity=='bias',]
  expect_equal(support$pooled,.2)
  expect_equal(support$chain_min_mean,.19)
  expect_equal(support$chain_max_mean,.21)
  detection <- a[a$comparison=='high_minus_binary_k100' & a$quantity=='mae',]
  expect_equal(detection$pooled,.3)
  expect_equal(detection$chain_min_mean,.29)
  expect_equal(detection$chain_max_mean,.31)
  expect_error(paired_observed_chains(d[!(d$community=='4-1' & d$arm=='low' & d$knots==50),],pooled))
})

test_that('range audit distinguishes shared and distinct constant chains', {
  d <- expand.grid(chain=1:2,range=c(.1,.2,.3))
  d$frequency <- as.numeric(d$range==.2)
  a <- diagnose_range_chains(d)
  expect_equal(a$state,'same_point_mass')
  expect_equal(a$max_total_variation,0)
  d$frequency <- as.numeric((d$chain==1 & d$range==.1)|(d$chain==2 & d$range==.3))
  a <- diagnose_range_chains(d)
  expect_equal(a$state,'distinct_point_masses')
  expect_equal(a$max_total_variation,1)
  expect_equal(a$chain_mean_gap,.2)
  d$frequency <- c(.5,0,.5,1,0,0)
  expect_equal(diagnose_range_chains(d)$state,'some_constant_chains')
})
