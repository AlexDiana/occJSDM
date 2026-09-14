test_that("half-Cauchy noise draws match direct integration over the noise SD", {
  # This reference integrates likelihood times prior directly, without the
  # inverse-gamma auxiliary variables used by the proposed sampler.
  n <- 40L
  for (rms in c(.1,1,3)) {
    scale <- .7
    log_density <- function(u) -(n-1)*u-n*rms^2/(2*exp(2*u))-
      log1p((exp(u)/scale)^2)
    bounds <- log(rms)+c(-8,8)
    peak <- optimize(log_density,bounds,maximum=TRUE)$objective
    mass <- integrate(function(u) exp(log_density(u)-peak),bounds[1],bounds[2])$value
    expected <- integrate(function(u) exp(u+log_density(u)-peak),bounds[1],bounds[2])$value/mass
    z <- matrix(rms,n,1)
    tau <- rms
    set.seed(726)
    draws <- numeric(12000)
    for (i in seq_along(draws)) {
      tau <- sample_tau_half_cauchy(z,z*0,tau,scale)
      draws[i] <- tau
    }
    expect_equal(mean(draws[-seq_len(1000)]),expected,tolerance=.01*expected)
  }
})

test_that("noise SD units follow the supplied half-Cauchy scale", {
  z <- matrix(c(.3,-.2,.7,.8),2,2)
  eta <- matrix(c(.1,-.1,.2,.5),2,2)
  set.seed(827)
  original <- sample_tau_half_cauchy(z,eta,c(.4,.9),.6)
  set.seed(827)
  rescaled <- sample_tau_half_cauchy(7*z,7*eta,7*c(.4,.9),7*.6)
  expect_equal(rescaled,7*original,tolerance=1e-12)
})

test_that("the public continuous fit uses and records the requested noise prior", {
  sim <- simulate_fixture(model="continuous",useSpatField=FALSE)
  sim$data_list$traits <- NULL
  run <- function(priors) {
    set.seed(328)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates=fixture_occ_covariates(),
      listParams=list(n_factors=0L,n_lattrait=0L),listPriors=priors,
      MCMCparams=list(nchain=2L,nburn=10L,niter=10L,nthin=1L))))
  }
  small <- run(list(tau_prior="half_cauchy",tau_scale=.1))
  large <- run(list(tau_prior="half_cauchy",tau_scale=10))
  expect_false(isTRUE(all.equal(small$results_output$jsdm_output$tau_output,
                               large$results_output$jsdm_output$tau_output)))
  expect_equal(small$infos$noise_prior,list(type="half_cauchy",scale=.1))
  legacy <- run(list(tau_prior="inverse_gamma",a_tau=5,b_tau=5))
  expect_equal(legacy$infos$noise_prior,list(type="inverse_gamma",shape=5,rate=5))
  other <- run(list(tau_prior="inverse_gamma",a_tau=5,b_tau=50))
  expect_gt(mean(other$results_output$jsdm_output$tau_output),
            mean(legacy$results_output$jsdm_output$tau_output))
  named_half <- run(list(tau_prior=c(choice="half_cauchy"),tau_scale=.1))
  expect_equal(named_half$results_output,small$results_output)
  named_ig <- run(list(tau_prior=c(choice="inverse_gamma"),a_tau=5,b_tau=50))
  expect_equal(named_ig$results_output,other$results_output)
  for (priors in list(list(tau_prior="typo"),list(tau_prior=NA_character_),
    list(tau_prior="half_cauchy",tau_scale=0),list(tau_prior="half_cauchy",tau_scale=Inf),
    list(tau_prior="inverse_gamma",a_tau=-1),list(tau_prior="inverse_gamma",b_tau=c(1,2))))
    expect_error(run(priors),"tau")
})
