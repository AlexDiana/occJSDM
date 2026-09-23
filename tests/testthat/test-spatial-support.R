test_that("requesting every spatial location retains the full covariance", {
  unique_locations <- rbind(c(0,0),c(.04,0),c(.4,.1),c(.7,.6),c(1,1),c(.2,.8))
  order <- c(3,1,6,3,2,5,4,1)
  observed <- unique_locations[order,,drop=FALSE]
  locations <- computeSpatialSummaries(observed,6)
  expect_equal(nrow(locations$X_tilde),6)
  expect_equal(locations$X_tilde,unique(observed))
  basis <- precomputeSORmatrices(.12,locations)$Ks_all[,,1]
  H <- KsBproduct(basis,diag(nrow(locations$X_tilde)),locations$Xs_centers)
  # Independent dense covariance at all unique observed coordinates.
  K <- exp(-as.matrix(dist(unique_locations))^2/(2*.12^2))
  expected <- K %*% solve(K+diag(1e-5,6),K)
  expect_equal(tcrossprod(H),unname(expected[order,order]),tolerance=1e-10)
})

test_that("one unique location can supply one requested support point", {
  observed <- matrix(rep(c(.3,.8),each=4),4,2)
  locations <- computeSpatialSummaries(observed,5)
  expect_equal(locations$X_tilde,matrix(c(.3,.8),1,2))
  H <- precomputeSORmatrices(.12,locations)$Ks_all[,,1,drop=FALSE]
  expect_equal(as.vector(H),rep(1/sqrt(1+1e-5),4),tolerance=1e-12)
})

test_that("the public fitter honours one support point per observed location", {
  sim <- simulate_fixture(model="continuous")
  knots <- nrow(unique(sim$data_list$info[,fixture_spat_covariates(),drop=FALSE]))
  fit <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    occCovariates=fixture_occ_covariates(),spatCovariates=fixture_spat_covariates(),
    listParams=list(n_supportpoints=knots,n_factors=0L,n_lattrait=0L),
    MCMCparams=list(nchain=2L,nburn=2L,niter=2L,nthin=1L))))
  expect_equal(fit$infos$ps,knots)
  expect_equal(dim(fit$results_output$jsdm_output$Bs_output)[1],knots)
  pred <- predictNewSites(fit,X_psi=sim$data_list$info[,fixture_occ_covariates()],
                         X_s=sim$data_list$info[,fixture_spat_covariates()],verbose=FALSE)
  expect_true(all(is.finite(pred)))
})
