spatial_range_fixture <- function() {
  basis <- array(c(.2,.5,.7, .8,.3,.1,
                   .4,.8,.9, .7,.6,.4,
                   .8,.9,1, .9,.8,.7), c(3,2,3))
  list(Bs=matrix(c(1,-.5,.3,1.2),2,2),
       offset=matrix(c(.2,-.1,.3,-.5,.8,.1),3,2),
       Omega=matrix(c(.7,1.1,.4,.8,.5,.9),3,2),
       centers=matrix(rep(1:2,each=3),3,2),
       summaries=list(Ks_all=basis,l_s_grid=c(.05,.15,.25)))
}

test_that("range weights match the conditional Gaussian observation likelihood", {
  x <- spatial_range_fixture()
  y <- matrix(c(.7,-.5,1.3,.4,-.2,.8),3,2)
  expected <- vapply(1:3, function(j) {
    eta <- x$offset + x$summaries$Ks_all[,,j] %*% x$Bs
    sum(dnorm(y, eta, 1/sqrt(x$Omega), log=TRUE)) +
      dgamma(x$summaries$l_s_grid[j],2,3,log=TRUE)
  }, numeric(1))
  actual <- spatial_range_logweights(x$Bs,x$offset,y*x$Omega,x$Omega,
                                      x$centers,x$summaries,2,3)
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-12)
})

test_that("range weights use the Polya-Gamma likelihood for binary states", {
  x <- spatial_range_fixture()
  kappa <- matrix(c(-.5,.5,.5,-.5,-.5,.5),3,2)
  expected <- vapply(1:3, function(j) {
    eta <- x$offset + x$summaries$Ks_all[,,j] %*% x$Bs
    sum(kappa*eta - x$Omega*eta^2/2) +
      dgamma(x$summaries$l_s_grid[j],1,1,log=TRUE)
  }, numeric(1))
  actual <- spatial_range_logweights(x$Bs,x$offset,kappa,x$Omega,
                                      x$centers,x$summaries,1,1)
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-12)
  # Changing the data with the same current spatial coefficients must alter
  # the range posterior; a density of the old field alone cannot do this.
  changed <- spatial_range_logweights(x$Bs,x$offset,-kappa,x$Omega,
                                       x$centers,x$summaries,1,1)
  expect_false(isTRUE(all.equal(actual-actual[1],changed-changed[1])))
})

test_that("the range sampler uses the grid prior when the field is zero", {
  x <- spatial_range_fixture()
  x$Bs[] <- 0
  kappa <- matrix(0,3,2)
  actual <- spatial_range_logweights(x$Bs,x$offset,kappa,x$Omega,
                                      x$centers,x$summaries,2,3)
  expected <- dgamma(x$summaries$l_s_grid,2,3,log=TRUE)
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-12)
  set.seed(711)
  expected_draws <- sample.int(3,300,replace=TRUE,prob=exp(expected-max(expected)))
  set.seed(711)
  draws <- replicate(300,sample_ls(1L,x$Bs,x$offset,kappa,x$Omega,
                                    x$centers,x$summaries,2,3))
  expect_identical(as.integer(draws),as.integer(expected_draws))
})

test_that("spatial coefficient updates retain the trait-predicted prior mean", {
  # With no observation precision the exact posterior is the supplied prior.
  # Its nonzero spatial mean must survive every coefficient sampler.
  args <- list(k=matrix(0,4,2),X=matrix(0,4,0),Tr=matrix(c(1,2),2,1),
               U=matrix(0,4,0),G=matrix(0,1,0),A=matrix(0,2,0),C=matrix(0,0,0),
               sigma_b=.3,Gs=matrix(c(1.5,-.5),1,2),As=matrix(0,2,0),
               Cs=matrix(0,0,2),sigma_bs=.2,Ks=matrix(0,4,2),
               Xs_centers=matrix(rep(1:2,each=4),4,2),
               Omega=matrix(0,4,2),model="continuous")
  expected <- t(args$Tr %*% args$Gs)
  for (sampler in list(sample_BBsL_cpp,sample_BBsL_parallel,sample_BBsL)) {
    set.seed(919)
    setOccJSDMSeed(919)
    draws <- replicate(256,do.call(sampler,args)$Bs)
    expect_equal(apply(draws,c(1,2),mean),expected,tolerance=.075)
  }
})

test_that("a spatial fit supports one knot without dropping matrix dimensions", {
  sim <- simulate_fixture(model="continuous")
  fit <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    occCovariates=fixture_occ_covariates(),spatCovariates=fixture_spat_covariates(),
    listParams=list(n_supportpoints=1L,n_factors=0L,n_lattrait=0L),
    MCMCparams=FIXTURE_MCMC)))
  expect_equal(dim(fit$results_output$jsdm_output$Bs_output)[1],1)
  pred <- predictNewSites(fit,X_psi=sim$data_list$info[,fixture_occ_covariates()],
                          X_s=sim$data_list$info[,fixture_spat_covariates()],verbose=FALSE)
  expect_true(all(is.finite(pred)))
})

test_that("routine spatial precomputation omits unused full-GP matrices", {
  set.seed(817)
  locations <- computeSpatialSummaries(matrix(rnorm(40),20,2),6)
  out <- precomputeSORmatrices(c(.1,.2),locations)
  expect_null(out$Lm1_grid)
  expect_null(out$logDetKuu_grid)
  diagnostic <- precomputeSORmatrices(c(.1,.2),locations,full_gp=TRUE)
  expect_equal(out$Ks_all,diagnostic$Ks_all)
  expect_equal(dim(diagnostic$Lm1_grid),c(20,20,2))
})

test_that("spatial fitting and prediction share the full SoR covariance", {
  raw <- as.matrix(expand.grid(x=seq(-1,1,length.out=4),y=seq(-1,1,length.out=4)))
  raw <- rbind(raw, raw[1,,drop=FALSE]) # duplicate sites keep their row mapping
  transformed <- create_covariates_matrix(as.data.frame(raw),spline_vars=FALSE,
                                          remove_intercept=TRUE)
  set.seed(818)
  locations <- computeSpatialSummaries(transformed$X,9,maxPoints=5)
  grid <- c(.15,.3)
  fit_basis <- precomputeSORmatrices(grid,locations)$Ks_all
  pred_basis <- createSpatialPredMatrix(as.data.frame(raw),grid,
                                        locations$X_tilde,transformed$list_matrix)
  # Compare implied covariance, independently of the internal index layout.
  for (j in seq_along(grid)) {
    H <- KsBproduct(fit_basis[,,j],diag(9),locations$Xs_centers)
    expect_equal(H,pred_basis[,,j],tolerance=1e-12)
    K <- K2(transformed$X,locations$X_tilde,1,grid[j])
    expected <- K %*% solve(K2(locations$X_tilde,locations$X_tilde,1,grid[j]) +
                              diag(1e-5,9),t(K))
    expect_equal(tcrossprod(H),expected,tolerance=1e-10)
    reversed <- locations
    reversed$X_tilde <- locations$X_tilde[9:1,]
    other <- precomputeSORmatrices(grid[j],reversed)$Ks_all[,,1]
    H_reversed <- KsBproduct(other,diag(9),reversed$Xs_centers)
    expect_equal(tcrossprod(H_reversed),expected,tolerance=1e-10)
  }
})
