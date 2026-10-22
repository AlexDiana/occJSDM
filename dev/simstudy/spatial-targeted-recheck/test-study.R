library(testthat)
source('generator.R')

test_that('rare-species targets refer to species prevalence and seeds reproduce whole communities', {
  a <- make_spatial_input(4L, 1L)
  expect_identical(a, make_spatial_input(4L, 1L))
  expect_equal(unname(colMeans(a$truth$psi)), c(.01,.01,.05,.05,.25,.25,.75,.75), tolerance=1e-10)
  expect_false(identical(a$truth$psi, make_spatial_input(4L, 2L)$truth$psi))
  expect_equal(plogis(sweep(a$truth$X %o% a$truth$B+a$truth$field, 2, a$truth$B0, '+')), a$truth$psi)
  expect_equal(as.vector(a$truth$Xs), as.vector(scale(as.matrix(a$data$binary$info[c('longitude','latitude')]))))
})

source('score.R')
test_that('missing saved binary probabilities cannot pass a numerical reconstruction check', {
  estimate <- matrix(c(.1,.3),2,1)
  expect_true(is.na(verify_saved_probability(estimate,NULL)))
  expect_error(verify_saved_probability(estimate,matrix(.1,1,1)))
  expect_error(verify_saved_probability(estimate,matrix(c(.1,NA),2,1)))
  expect_error(verify_saved_probability(estimate,matrix(c(.1,.4),2,1)))
  expect_equal(verify_saved_probability(estimate,estimate),0)
})
test_that('posterior probabilities average transformed draws and keep one-species dimensions', {
  xs <- matrix(c(0,1,0,1),2,2)
  fit <- list(Xs=xs,X_psi=matrix(c(-1,1),2,1),
    infos=list(ps=2L,n_factors=0L,l_s_grid=c(.1,.2),
      list_Xs=list(X_s=xs,X_tilde=xs,Xs_index=1:2)),
    results_output=list(jsdm_output=list(B0_output=array(rep(c(-4,0,4,4),2),c(1,4,2)),
      B_output=array(0,c(1,1,4,2)),Bs_output=array(0,c(2,1,4,2)),
      idx_ls_output=matrix(1L,4,2))))
  out <- reconstruct_spatial_draws(fit)
  expected <- (plogis(-4)+.5+2*plogis(4))/4
  expect_equal(dim(out$probability),c(2,4,2))
  expect_equal(rowMeans(matrix(out$probability,2)),rep(expected,2))
  expect_gt(abs(expected-plogis(1)),.1)
  expect_equal(out$field_mean,matrix(0,2,1))
})

test_that('rare-species groups use whole-species prevalence rather than individual low-probability cells', {
  truth <- matrix(c(.001,.019,.01,.99),2,2)
  g <- spatial_groups(truth,c(.01,.5))
  expect_identical(which(g$prevalence_1pct),1:2)
  expect_identical(which(g$low),1:3)
  expect_identical(which(g$high),4L)
  expect_false(any(g$prevalence_5pct))
})

test_that('point scores expose cancellation and score intervals against probability truth', {
  x <- array(c(.1,.9,.1,.9,.1,.9,.1,.9),c(2,2,2))
  y <- c(.2,.8)
  a <- score_draw_block(x,y,'probe')
  expect_equal(a$elements$bias,c(-.1,.1))
  expect_equal(a$summary$bias,0,tolerance=1e-12)
  expect_equal(a$summary$mae,.1)
  expect_equal(a$summary$rmse,.1)
  expect_equal(a$summary$coverage,0)
})

test_that('binary and two-stage arms share biology and PCR blocks map to their actual samples', {
  a <- make_spatial_input(6L, 1L)
  expect_identical(a$data$binary$OTU, a$truth$z)
  lo <- a$data$low; hi <- a$data$high
  expect_identical(lo$info, hi$info)
  expect_equal(dim(lo$OTU), c(2400,8))
  expect_equal(as.numeric(table(lo$info$Site)), rep(24,100))
  expect_equal(as.numeric(table(lo$info$Sample,lo$info$Primer)), rep(6,400))
  expect_true(all(lo$OTU %in% 0:1))
  wrows <- a$truth$w[lo$info$Sample,]
  expect_identical(lo$OTU[wrows==1], hi$OTU[wrows==1])
  expect_true(all(hi$OTU[wrows==0] >= lo$OTU[wrows==0]))
  expect_equal(a$truth$theta, plogis(a$truth$Xt %*% a$truth$beta_theta))
  expect_equal(a$truth$p, a$truth$p_positive)
  expect_equal(a$truth$q$low, a$truth$q_positive$low)
})

test_that('zero-occupancy species are retained and study fitting controls are bounded', {
  a <- make_spatial_input(8L, 1L)
  # This scenario deliberately includes expected occupied-site counts of one.
  allz <- lapply(1:9, function(r) make_spatial_input(6L,r)$truth$z)
  expect_true(any(vapply(allz,function(z) any(colSums(z)==0),logical(1))))
  expect_true(all(vapply(allz,ncol,integer(1))==8L))
  expect_error(make_spatial_input(0L,1L))
  expect_error(make_spatial_input(4L,0L))
})
