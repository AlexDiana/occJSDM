test_that("spatial half-Cauchy SD draws match a directly integrated conditional density", {
  # Removing the half-Cauchy factor, using variance instead of SD, or counting
  # only species (rather than all coefficients) changes this reference.
  p <- 4L; S <- 5L; m <- p*S
  Tr <- matrix(seq_len(S)/S,S,1)
  G <- matrix(c(.3,-.2,.7,-.4),1,p)
  A <- matrix(seq_len(S)/10,S,1); C <- matrix(c(.1,.2,-.1,.3),1,p)
  mean <- t(Tr %*% G + A %*% C)
  for (rms in c(1e-4,.3,2)) {
    residual <- matrix(rep(c(-1,1),length.out=m)*rms,p,S)
    scale <- .7
    log_density <- function(u) -(m-1)*u-sum(residual^2)/(2*exp(2*u))-
      log1p((exp(u)/scale)^2)
    bounds <- log(rms)+c(-8,8)
    peak <- optimize(log_density,bounds,maximum=TRUE)$objective
    mass <- integrate(function(u) exp(log_density(u)-peak),bounds[1],bounds[2])$value
    expected <- integrate(function(u) exp(u+log_density(u)-peak),bounds[1],bounds[2])$value/mass
    sigma <- rms
    set.seed(729);setOccJSDMSeed(729)
    draws <- numeric(14000)
    for(i in seq_along(draws)) {
      sigma <- sample_spatial_sd_half_cauchy(mean+residual,Tr,G,A,C,sigma,scale)
      draws[i] <- sigma
    }
    expect_lte(abs(mean(draws[-seq_len(2000)])/expected-1),.015)
    expect_true(all(is.finite(draws) & draws>0))
  }
})

test_that("spatial amplitude updates subtract trait means and respect scale units", {
  B <- matrix(c(.3,-.2,.7,.8,.2,-.1),2,3)
  Tr <- matrix(c(1,-1,2),3,1);G <- matrix(c(.5,-.1),1,2)
  A <- matrix(c(.2,.1,.4),3,1);C <- matrix(c(.3,.6),1,2)
  zeroTr <- matrix(0,3,0);zeroG <- matrix(0,0,2)
  draw <- function(B,Tr,G,A,C,sigma,scale) {
    set.seed(826);setOccJSDMSeed(826)
    sample_spatial_sd_half_cauchy(B,Tr,G,A,C,sigma,scale)
  }
  a <- draw(B+t(Tr%*%G+A%*%C),Tr,G,A,C,.4,.6)
  expect_equal(a,draw(B,zeroTr,zeroG,zeroTr,zeroG,.4,.6),tolerance=1e-12)
  expect_equal(draw(7*B,zeroTr,zeroG,zeroTr,zeroG,7*.4,7*.6),7*a,tolerance=1e-12)
})

test_that("public spatial fits apply an explicit prior and preserve the default RNG path", {
  sim <- simulate_fixture(model="binary",useSpatField=TRUE)
  sim$data_list$traits <- NULL
  run <- function(priors=list(),spatial=TRUE) {
    set.seed(932)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates=fixture_occ_covariates(),
      spatCovariates=if(spatial)fixture_spat_covariates() else NULL,
      listParams=list(n_factors=0L,n_lattrait=0L,n_supportpoints=8L),
      listPriors=priors,MCMCparams=list(nchain=2L,nburn=12L,niter=15L,nthin=1L))))
  }
  legacy <- run()
  explicit <- run(list(sigma_bs_prior="inverse_gamma"))
  expect_identical(explicit$results_output,legacy$results_output)
  expect_identical(legacy$infos$spatial_sd_prior,list(type="inverse_gamma",shape=10,rate=1))
  small <- run(list(sigma_bs_prior="half_cauchy",sigma_bs_scale=.1))
  large <- run(list(sigma_bs_prior="half_cauchy",sigma_bs_scale=10))
  expect_false(isTRUE(all.equal(small$results_output$jsdm_output$sigmabs_output,
                               large$results_output$jsdm_output$sigmabs_output)))
  expect_identical(small$infos$spatial_sd_prior,list(type="half_cauchy",scale=.1))
  expect_null(small$infos$noise_prior)
  expect_identical(run(list(sigma_bs_prior=c(named="half_cauchy"),sigma_bs_scale=.1))$results_output,
                   small$results_output)
  for(priors in list(list(sigma_bs_prior="typo"),list(sigma_bs_prior=NA_character_),
    list(sigma_bs_prior=c("half_cauchy","inverse_gamma")),
    list(sigma_bs_prior="half_cauchy",sigma_bs_scale=0),
    list(sigma_bs_prior="half_cauchy",sigma_bs_scale=Inf),
    list(sigma_bs_prior="half_cauchy",sigma_bs_scale=c(1,2)),
    list(sigma_bs_prior="half_cauchy",sigma_bs_scale="1")))
    expect_error(run(priors),"sigma_bs")
  expect_error(run(list(sigma_bs_prior="half_cauchy"),FALSE),"spatial")
  expect_error(run(list(sigma_bs_scale=1)),"sigma_bs_scale")
  expect_null(run(spatial=FALSE)$infos$spatial_sd_prior)
})
