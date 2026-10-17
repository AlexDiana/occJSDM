# sigma_b is the standard deviation of each species' environmental responses
# B around Tr G + A C. Its prior was a hard-coded InvGamma(10, 1) on sigma_b^2;
# listPriors can now place a half-Cauchy on sigma_b (sigma_b_prior,
# sigma_b_scale) or set the inverse-gamma shape and rate (a_sigmab, b_sigmab).
# The default must leave fits bit-identical, and the choice must reach the
# one place sigma_b is drawn, update_jSDMcoef(), without touching the spatial
# SD prior that shares the half-Cauchy sampler.

test_that("read_sigma_b_prior validates and defaults", {
  rd <- occJSDM:::read_sigma_b_prior
  default <- list(type = "inverse_gamma", shape = 10, rate = 1)
  expect_identical(rd(list()), default)
  expect_identical(rd(NULL), default)
  expect_identical(rd(list(sigma_b_prior = "inverse_gamma")), default)
  expect_identical(rd(list(a_sigmab = 10, b_sigmab = 1)), default)
  expect_identical(rd(list(a_sigmab = 2, b_sigmab = 1 / 9)),
                   list(type = "inverse_gamma", shape = 2, rate = 1 / 9))
  expect_identical(rd(list(sigma_b_prior = c(named = "half_cauchy"))),
                   list(type = "half_cauchy", scale = 1))
  expect_identical(rd(list(sigma_b_prior = "half_cauchy", sigma_b_scale = c(s = 2.5))),
                   list(type = "half_cauchy", scale = 2.5))
  for (bad in list(list(sigma_b_prior = "typo"), list(sigma_b_prior = NA_character_),
                   list(sigma_b_prior = c("half_cauchy", "inverse_gamma")),
                   list(sigma_b_prior = 1)))
    expect_error(rd(bad), "sigma_b_prior must be 'half_cauchy' or 'inverse_gamma'")
  for (bad in list(0, -1, NA_real_, Inf, NaN, c(1, 2), numeric(0), "1", TRUE)) {
    expect_error(rd(list(sigma_b_prior = "half_cauchy", sigma_b_scale = bad)),
                 "sigma_b_scale must be a finite positive number")
    expect_error(rd(list(a_sigmab = bad)), "a_sigmab must be a finite positive number")
    expect_error(rd(list(b_sigmab = bad)), "b_sigmab must be a finite positive number")
  }
  expect_error(rd(list(sigma_b_scale = 1)), "sigma_b_scale requires sigma_b_prior = 'half_cauchy'")
  expect_error(rd(list(sigma_b_prior = "inverse_gamma", sigma_b_scale = 1)), "sigma_b_scale requires")
  expect_error(rd(list(sigma_b_prior = "half_cauchy", a_sigmab = 2)),
               "a_sigmab and b_sigmab require sigma_b_prior = 'inverse_gamma'")
  expect_error(rd(list(sigma_b_prior = "half_cauchy", b_sigmab = 2)), "require sigma_b_prior")
})

test_that("sample_sigmab draws sigma_b^2 from InvGamma(shape + pS/2, rate + SS/2)", {
  # Documents the parameterisation the user-set a_sigmab and b_sigmab enter:
  # b_sigmab is the rate (the numerator of the inverse-gamma's exponential
  # term), so a larger b_sigmab pulls sigma_b up. The posterior mean of
  # sigma_b^2 is (b + SS/2) / (a + pS/2 - 1).
  p <- 2L; S <- 5L
  Tr <- matrix(0, S, 0); G <- matrix(0, 0, p); A <- matrix(0, S, 0); C <- matrix(0, 0, p)
  B <- matrix(c(.3, -.4, .5, -.2, .1, .6, -.3, .2, -.5, .4), p, S)
  SS <- sum(B^2); M <- 20000L
  for (ab in list(c(10, 1), c(2, 1 / 9), c(.5, 3))) {
    a <- ab[1]; b <- ab[2]
    shape <- a + p * S / 2; rate <- b + SS / 2
    want_mean <- rate / (shape - 1)
    want_var <- rate^2 / ((shape - 1)^2 * (shape - 2))
    set.seed(61); setOccJSDMSeed(61)
    v <- replicate(M, sample_sigmab(B, Tr, G, A, C, a, b))^2
    expect_lt(abs(mean(v) - want_mean), 4 * sqrt(want_var / M))
  }
})

sigma_b_fixture_run <- function(sim, priors, spatial = FALSE, seed = 57,
                                MCMCparams = list(nchain = 2L, nburn = 12L, niter = 15L, nthin = 1L)) {
  set.seed(seed)
  suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    occCovariates = fixture_occ_covariates(),
    spatCovariates = if (spatial) fixture_spat_covariates() else NULL,
    listParams = c(list(n_factors = 0L, n_lattrait = 0L),
                   if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),
    listPriors = priors, MCMCparams = MCMCparams)))
}

test_that("default fits are identical with and without the explicit default sigma_b prior", {
  sim <- simulate_fixture(model = "binary")
  default <- list(type = "inverse_gamma", shape = 10, rate = 1)
  for (spatial in c(FALSE, TRUE)) {
    base <- sigma_b_fixture_run(sim, list(), spatial)
    expect_identical(base$infos$sigma_b_prior, default)
    expect_identical(sigma_b_fixture_run(sim, list(sigma_b_prior = "inverse_gamma"), spatial)$results_output,
                     base$results_output)
    expect_identical(sigma_b_fixture_run(sim, list(a_sigmab = 10, b_sigmab = 1), spatial)$results_output,
                     base$results_output)
    weak <- sigma_b_fixture_run(sim, list(a_sigmab = 2, b_sigmab = 1 / 9), spatial)
    expect_false(identical(weak$results_output$jsdm_output$sigmab_output,
                           base$results_output$jsdm_output$sigmab_output))
    expect_identical(weak$infos$sigma_b_prior, list(type = "inverse_gamma", shape = 2, rate = 1 / 9))
    small <- sigma_b_fixture_run(sim, list(sigma_b_prior = "half_cauchy", sigma_b_scale = .1), spatial)
    large <- sigma_b_fixture_run(sim, list(sigma_b_prior = "half_cauchy", sigma_b_scale = 10), spatial)
    expect_false(isTRUE(all.equal(small$results_output$jsdm_output$sigmab_output,
                                  large$results_output$jsdm_output$sigmab_output)))
    expect_identical(small$infos$sigma_b_prior, list(type = "half_cauchy", scale = .1))
    expect_identical(sigma_b_fixture_run(sim, list(sigma_b_prior = c(named = "half_cauchy"),
                                                   sigma_b_scale = .1), spatial)$results_output,
                     small$results_output)
    # The sigma_b choice leaves the other recorded priors alone.
    expect_identical(small$infos$spatial_sd_prior, base$infos$spatial_sd_prior)
    expect_identical(small$infos$intercept_prior, base$infos$intercept_prior)
    expect_null(small$infos$noise_prior)
  }
  expect_error(sigma_b_fixture_run(sim, list(sigma_b_prior = "typo")), "sigma_b_prior")
  expect_error(sigma_b_fixture_run(sim, list(sigma_b_scale = 1)), "sigma_b_scale")
  expect_error(sigma_b_fixture_run(sim, list(sigma_b_prior = "half_cauchy", a_sigmab = 2)), "a_sigmab")
})

test_that("the coefficient update draws sigma_b under the requested prior and sigma_bs under its own", {
  # Record what the two SD samplers receive inside the MCMC. The spatial SD
  # shares the half-Cauchy sampler, so the spatial fits set it to a different
  # scale and the inverse-gamma fits leave it at its (10, 1) default: each
  # call is then attributable to sigma_b or sigma_bs by its arguments, and the
  # counts show the sigma_b draw happens once per iteration either way.
  original_ig <- sample_sigmab
  original_hc <- sample_spatial_sd_half_cauchy
  sim <- simulate_fixture(model = "binary")
  fit_with <- function(priors, spatial) {
    seen <- list(ig = NULL, hc = NULL)
    update <- update_jSDMcoef
    update_env <- new.env(parent = environment(update))
    update_env$sample_sigmab <- function(B, Tr, G, A, C, a_sigmab, b_sigmab) {
      seen$ig <<- rbind(seen$ig, data.frame(p = nrow(B), a = a_sigmab, b = b_sigmab))
      original_ig(B, Tr, G, A, C, a_sigmab, b_sigmab)
    }
    update_env$sample_spatial_sd_half_cauchy <- function(B, Tr, G, A, C, sigma, scale) {
      seen$hc <<- rbind(seen$hc, data.frame(p = nrow(B), scale = scale))
      original_hc(B, Tr, G, A, C, sigma, scale)
    }
    environment(update) <- update_env
    # Mocked in the namespace, where run_one_chain() looks it up.
    local_mocked_bindings(update_jSDMcoef = update, .package = "occJSDM")
    sigma_b_fixture_run(sim, priors, spatial,
                        MCMCparams = list(nchain = 1L, nburn = 3L, niter = 4L, nthin = 1L))
    seen
  }
  # B has one row per occupancy covariate. Bs has one row per spatial basis
  # column, which is not the number of spatial covariates, so the spatial
  # calls are recognised as the ones whose row count is not B's.
  n_iter <- 7L; p <- length(fixture_occ_covariates())
  spatial_rows <- function(calls) {
    ps <- unique(calls$p[calls$p != p])
    expect_length(ps, 1L); expect_gt(ps, 0)
    ps
  }

  plain <- fit_with(list(), FALSE)
  expect_null(plain$hc)
  expect_identical(nrow(plain$ig), n_iter)
  expect_true(all(plain$ig$p == p & plain$ig$a == 10 & plain$ig$b == 1))

  weak <- fit_with(list(a_sigmab = 2, b_sigmab = 1 / 9), TRUE)
  expect_null(weak$hc)
  expect_identical(nrow(weak$ig), 2L * n_iter)
  ps <- spatial_rows(weak$ig)
  expect_identical(sum(weak$ig$p == p & weak$ig$a == 2 & weak$ig$b == 1 / 9), n_iter)
  expect_identical(sum(weak$ig$p == ps & weak$ig$a == 10 & weak$ig$b == 1), n_iter)

  hc <- fit_with(list(sigma_b_prior = "half_cauchy", sigma_b_scale = 3), FALSE)
  expect_null(hc$ig)
  expect_identical(nrow(hc$hc), n_iter)
  expect_true(all(hc$hc$p == p & hc$hc$scale == 3))

  both <- fit_with(list(sigma_b_prior = "half_cauchy", sigma_b_scale = 3,
                        sigma_bs_prior = "half_cauchy", sigma_bs_scale = .5), TRUE)
  expect_null(both$ig)
  expect_identical(nrow(both$hc), 2L * n_iter)
  expect_identical(spatial_rows(both$hc), ps)
  expect_identical(sum(both$hc$p == p & both$hc$scale == 3), n_iter)
  expect_identical(sum(both$hc$p == ps & both$hc$scale == .5), n_iter)

  mixed <- fit_with(list(sigma_b_prior = "half_cauchy", sigma_b_scale = 3), TRUE)
  expect_identical(nrow(mixed$hc), n_iter)
  expect_true(all(mixed$hc$p == p & mixed$hc$scale == 3))
  expect_identical(nrow(mixed$ig), n_iter)
  expect_true(all(mixed$ig$p == ps & mixed$ig$a == 10 & mixed$ig$b == 1))
})

test_that("a weak prior lets sigma_b leave the default's 0.33 when the responses are spread", {
  # With B's residual spread well above 0.33, the InvGamma(10, 1) prior holds
  # sigma_b down and the weaker priors let it rise. A fit with a dataset whose
  # true sigma_b is 1 shows the ordering; the test asks only for the direction.
  sim <- simulate_fixture(model = "binary", useSpatField = FALSE, seed = 7L)
  run <- function(priors) sigma_b_fixture_run(sim, priors, FALSE, seed = 91,
    MCMCparams = list(nchain = 2L, nburn = 150L, niter = 150L, nthin = 1L))
  sb <- function(fit) mean(fit$results_output$jsdm_output$sigmab_output)
  default <- sb(run(list()))
  weak <- sb(run(list(a_sigmab = 2, b_sigmab = 1 / 9)))
  cauchy <- sb(run(list(sigma_b_prior = "half_cauchy", sigma_b_scale = 1)))
  expect_lt(default, .5)
  expect_gt(weak, default)
  expect_gt(cauchy, default)
})
