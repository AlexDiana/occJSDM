spatial_range_fixture <- function() {
  basis <- array(c(.2,.5,.7, .8,.3,.1,
                   .4,.8,.9, .7,.6,.4,
                   .8,.9,1, .9,.8,.7), c(3,2,3))
  list(X=matrix(c(-1,0,1),3,1), U=matrix(c(.2,-.5,.8),3,1),
       M_B=matrix(c(.4,-.3),1,2), M_Bs=matrix(c(1,-.5,.3,1.2),2,2),
       sigma_b=.6,sigma_bs=1.2,
       Omega=matrix(c(.7,1.1,.4,.8,.5,.9),3,2),
       Xs_centers=matrix(rep(1:2,each=3),3,2),
       list_SoRSummaries=list(Ks_all=basis,l_s_grid=c(.05,.15,.25)),
       a_l_s=2,b_l_s=3)
}

# Independent reference integrates the coefficients in observation space,
# giving an n x n marginal covariance rather than the fitted precision matrix.
range_marginal_reference <- function(x,y) {
  vapply(1:3,function(j) {
    Z <- cbind(1,x$X,x$U,x$list_SoRSummaries$Ks_all[,,j])
    V <- diag(c(1,x$sigma_b^2,1,rep(x$sigma_bs^2,2)))
    sum(vapply(1:2,function(s) {
      mean <- Z %*% c(0,x$M_B[,s],0,x$M_Bs[,s])
      cov <- diag(1/x$Omega[,s]) + Z %*% V %*% t(Z)
      L <- chol(cov)
      r <- forwardsolve(t(L),y[,s]-mean)
      -.5*sum(r^2)-sum(log(diag(L)))
    },numeric(1))) + dgamma(x$list_SoRSummaries$l_s_grid[j],2,3,log=TRUE)
  },numeric(1))
}

test_that("blocked range weights match an independently integrated Gaussian model", {
  x <- spatial_range_fixture()
  y <- matrix(c(.7,-.5,1.3,.4,-.2,.8),3,2)
  expected <- range_marginal_reference(x,y)
  actual <- do.call(spatial_range_logweights,c(x,list(kappa=y*x$Omega)))
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-12)
})

test_that("blocked range weights integrate the augmented binary likelihood", {
  x <- spatial_range_fixture()
  kappa <- matrix(c(-.5,.5,.5,-.5,-.5,.5),3,2)
  # Completing the PG kernel's square yields the pseudo-response kappa/Omega.
  expected <- range_marginal_reference(x,kappa/x$Omega)
  actual <- do.call(spatial_range_logweights,c(x,list(kappa=kappa)))
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-12)
  changed <- do.call(spatial_range_logweights,c(x,list(kappa=-kappa)))
  expect_false(isTRUE(all.equal(actual-actual[1],changed-changed[1])))
})

test_that("blocked range updates use only the grid prior without observations", {
  x <- spatial_range_fixture()
  x$Omega[] <- 0
  actual <- do.call(spatial_range_logweights,c(x,list(kappa=matrix(0,3,2))))
  expected <- dgamma(x$list_SoRSummaries$l_s_grid,2,3,log=TRUE)
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-12)
  set.seed(711)
  expected_draws <- sample.int(3,300,replace=TRUE,prob=exp(expected-max(expected)))
  set.seed(711)
  draws <- replicate(300,sample_ls(actual))
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

test_that("the integrated binary kernel agrees with numerical quadrature", {
  basis <- array(c(.2,.6,1.1),c(1,1,3))
  grid <- c(.05,.15,.25)
  actual <- spatial_range_logweights(
    matrix(0,1,0),matrix(0,1,0),matrix(0,0,1),matrix(.6,1,1),1,.8,
    matrix(.5,1,1),matrix(.7,1,1),matrix(1,1,1),
    list(Ks_all=basis,l_s_grid=grid),2,3)
  expected <- vapply(seq_along(grid),function(j) {
    h <- basis[1,1,j]
    integral <- integrate(function(eta) {
      exp(.5*eta-.7*eta^2/2)*dnorm(eta,h*.6,sqrt(1+(h*.8)^2))
    },-Inf,Inf,rel.tol=1e-10)$value
    log(integral)+dgamma(grid[j],2,3,log=TRUE)
  },numeric(1))
  expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-9)
})

test_that("the full fitter redraws coefficients immediately after choosing range", {
  original_draw <- sample_BBsL_cpp
  events <- character()
  aligned <- logical()
  update <- update_jSDMcoef
  update_env <- new.env(parent=environment(update))
  update_env$sample_ls <- function(logweights) {
    events <<- c(events,"range")
    10L
  }
  update_env$sample_BBsL_cpp <- function(...) {
    events <<- c(events,"coefficients")
    state <- parent.frame()
    aligned <<- c(aligned,isTRUE(all.equal(state$Ks,
      matrix(state$list_SoRSummaries$Ks_all[,,10],nrow=nrow(state$z)))))
    original_draw(...)
  }
  environment(update) <- update_env
  fitter <- runOccJSDM
  fit_env <- new.env(parent=environment(fitter))
  fit_env$update_jSDMcoef <- update
  environment(fitter) <- fit_env
  sim <- simulate_fixture(model="binary")
  suppressMessages(suppressWarnings(fitter(sim$data_list,
    occCovariates=fixture_occ_covariates(),spatCovariates=fixture_spat_covariates(),
    listParams=list(n_supportpoints=FIXTURE_KNOTS),
    MCMCparams=list(nchain=2L,nburn=2L,niter=2L,nthin=1L))))
  expect_gt(length(aligned),0)
  expect_true(all(aligned))
  expect_identical(events,rep(c("range","coefficients"),length(aligned)))
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
