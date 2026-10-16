# Experimental community mean for the collection intercepts
# (listPriors$theta_intercept_prior = "community"), the collection-side
# companion of b0_prior = "community" (test-b0-community-prior.R). Each
# species' collection intercept beta_theta[1, s] gets a Normal(mu_theta, 1)
# prior and mu_theta ~ Normal(0, mu_theta_sd^2) is drawn from its conjugate
# conditional every iteration, after the collection coefficients. The prior
# reaches the collection sampler through b_betatheta[1], the only place the
# collection intercepts' prior enters; both C++ samplers read it there. The
# default, theta_intercept_prior = "fixed", must leave fits bit for bit as
# they were.

test_that("read_collection_intercept_prior keeps the default and validates the option", {
  read <- occJSDM:::read_collection_intercept_prior
  for (model in c("two_stage", "occupancy", "binary", "continuous")) {
    expect_null(read(list(), model))
    expect_null(read(list(theta_intercept_prior = "fixed"), model))
  }
  for (model in c("two_stage", "occupancy")) {
    expect_identical(read(list(theta_intercept_prior = "community"), model),
                     list(type = "community", sd = 1, mu_theta_mean = 0, mu_theta_sd = 2))
    expect_identical(read(list(theta_intercept_prior = "community",
                               mu_theta_sd = c(named = 3)), model),
                     list(type = "community", sd = 1, mu_theta_mean = 0, mu_theta_sd = 3))
  }
  for (model in c("binary", "continuous")) {
    expect_error(read(list(theta_intercept_prior = "community"), model),
                 "requires an occupancy or two-stage model")
  }
  for (bad in list("Community", "estimated", NA_character_, c("fixed", "community"),
                   1, TRUE, character(0))) {
    expect_error(read(list(theta_intercept_prior = bad), "two_stage"),
                 "theta_intercept_prior must be 'fixed' or 'community'")
  }
  for (bad in list(0, -1, NA_real_, Inf, NaN, c(1, 2), numeric(0), "a", TRUE)) {
    expect_error(read(list(theta_intercept_prior = "community", mu_theta_sd = bad),
                      "two_stage"),
                 "mu_theta_sd must be a finite positive number")
  }
  expect_error(read(list(mu_theta_sd = 2), "two_stage"),
               "mu_theta_sd requires theta_intercept_prior = 'community'")
  expect_error(read(list(theta_intercept_prior = "fixed", mu_theta_sd = 2), "occupancy"),
               "mu_theta_sd requires theta_intercept_prior = 'community'")
})

test_that("the collection community mean is drawn from its conjugate conditional", {
  # beta_theta[1, s] ~ N(mu, s0^2) and mu ~ N(0, h^2) give mu | . ~ N(m, v)
  # with v = 1 / (1/h^2 + S/s0^2) and m = v * sum(beta_theta[1, ]) / s0^2.
  intercepts <- c(-1.1, -0.4, -2.0, -0.7, -1.6, 0.3)
  for (case in list(c(1, 2), c(.5, 2), c(1, .3), c(1, 5))) {
    s0 <- case[1]; h <- case[2]
    v <- 1 / (1 / h^2 + length(intercepts) / s0^2)
    m <- v * sum(intercepts) / s0^2
    set.seed(7)
    one <- sample_mu_theta(intercepts, s0, h)
    set.seed(7)
    expect_equal(one, rnorm(1, m, sqrt(v)), tolerance = 1e-12)
    M <- 20000L
    set.seed(8)
    draws <- replicate(M, sample_mu_theta(intercepts, s0, h))
    info <- paste("sd =", s0, "mu_theta_sd =", h)
    expect_lt(abs(mean(draws) - m), 5 * sqrt(v / M), label = info)
    expect_lt(abs(var(draws) / v - 1), 5 * sqrt(2 / (M - 1)), label = info)
  }
})

test_that("both collection samplers centre the intercept on the prior mean they are given", {
  # One species, one occupied site, two field samples of which one caught it:
  # the intercept's posterior is N(mu, 1) times logistic(b) * (1 - logistic(b)),
  # whose mean is a one-dimensional integral. Repeated Polya-Gamma steps form
  # a chain with that posterior as its stationary distribution. A sampler
  # ignoring the prior mean could not pass: the posterior means at mu = 0 and
  # mu = -2.2 differ by about 1.5, against a tolerance of about 0.1.
  w <- matrix(c(1, 0), 2, 1); z <- matrix(1, 1, 1); idx_z <- c(1L, 1L)
  X_theta <- matrix(1, 2, 1)
  M <- 6000L
  posterior_mean <- function(mu) {
    kernel <- function(b) dnorm(b, mu, 1) * plogis(b) * plogis(-b)
    integrate(function(b) b * kernel(b), -Inf, Inf)$value /
      integrate(kernel, -Inf, Inf)$value
  }
  for (sampler_name in c("sample_betatheta_cpp", "sample_betatheta_cpp_parallel")) {
    sampler <- get(sampler_name)
    for (mu in c(0, -2.2, 1.4)) {
      set.seed(31); setOccJSDMSeed(31)
      beta <- matrix(0, 1, 1)
      draws <- numeric(M)
      for (i in seq_len(M)) {
        beta <- sampler(w, z, beta, idx_z, X_theta, mu, matrix(1, 1, 1))
        draws[i] <- beta[1, 1]
      }
      draws <- draws[-(1:100)]
      # Lag-one autocorrelation is small here; a quarter of the draws as the
      # effective sample size is conservative.
      tolerance <- 5 * sd(draws) / sqrt(length(draws) / 4)
      expect_lt(abs(mean(draws) - posterior_mean(mu)), tolerance,
                label = paste(sampler_name, "mu =", mu))
    }
  }
})

test_that("unoccupied species' collection coefficients are drawn from their prior", {
  z <- cbind(c(1, 0, 1), c(0, 0, 0), c(0, 1, 0))
  idx_z <- c(1L, 1L, 2L, 3L, 3L)
  beta <- matrix(c(-1, .5, 2, -3, .1, .2), 2, 3)
  b <- c(-.8, 0); B <- diag(c(1, 2))
  set.seed(4)
  out <- sample_betatheta_unoccupied(beta, z, idx_z, b, B)
  expect_identical(out[, c(1, 3)], beta[, c(1, 3)])
  expect_false(any(out[, 2] == beta[, 2]))
  M <- 20000L
  set.seed(5)
  draws <- replicate(M, sample_betatheta_unoccupied(beta, z, idx_z, b, B)[, 2])
  expect_true(all(abs(rowMeans(draws) - b) < 5 * sqrt(diag(B) / M)))
  expect_true(all(abs(apply(draws, 1, var) / diag(B) - 1) < 5 * sqrt(2 / (M - 1))))
  # With every species occupied somewhere nothing changes and no random
  # number is drawn.
  z_all <- cbind(c(1, 0, 1), c(0, 1, 0), c(0, 1, 0))
  set.seed(6); before <- .Random.seed
  expect_identical(sample_betatheta_unoccupied(beta, z_all, idx_z, b, B), beta)
  expect_identical(.Random.seed, before)
})

theta_community_fit <- function(priors, model = "two_stage", spatial = FALSE,
                                mcmc = list(nchain = 2L, nburn = 12L, niter = 15L,
                                            nthin = 1L), seed = 41) {
  sim <- simulate_fixture(model = model)
  set.seed(seed)
  suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    occCovariates = fixture_occ_covariates(),
    collCovariates = if (model %in% c("occupancy", "two_stage")) "X_theta" else NULL,
    spatCovariates = if (spatial) fixture_spat_covariates() else NULL,
    listParams = c(list(n_factors = 1L, n_lattrait = 1L),
                   if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),
    listPriors = priors, MCMCparams = mcmc)))
}

test_that("each iteration's collection draw uses the community mean drawn in the iteration before", {
  # Record the prior mean the collection sampler receives and the intercepts
  # mu_theta is drawn from. The first iteration starts at 0; afterwards each
  # must use the value drawn at the end of the one before, which for retained
  # iterations is the previous saved draw, and mu_theta must be drawn from
  # that iteration's intercepts with the sampler's own prior sd.
  original_sampler <- sample_betatheta_cpp_parallel
  original_mu <- sample_mu_theta
  fit_with <- function(priors, model) {
    seen <- list(mean = numeric(), slope_mean = numeric(), intercepts = list(),
                 sd = numeric(), hyper = numeric())
    local_mocked_bindings(
      sample_betatheta_cpp_parallel = function(w, z, beta_theta, idx_z, X_theta,
                                               b_betatheta, B_betatheta) {
        seen$mean <<- c(seen$mean, b_betatheta[1])
        seen$slope_mean <<- c(seen$slope_mean, b_betatheta[-1])
        original_sampler(w, z, beta_theta, idx_z, X_theta, b_betatheta, B_betatheta)
      },
      sample_mu_theta = function(beta_theta0, sd, mu_theta_sd) {
        seen$intercepts[[length(seen$intercepts) + 1L]] <<- beta_theta0
        seen$sd <<- c(seen$sd, sd); seen$hyper <<- c(seen$hyper, mu_theta_sd)
        original_mu(beta_theta0, sd, mu_theta_sd)
      },
      .package = "occJSDM")
    fit <- theta_community_fit(priors, model,
                               mcmc = list(nchain = 1L, nburn = 2L, niter = 4L, nthin = 1L))
    c(seen, list(fit = fit))
  }
  for (model in c("two_stage", "occupancy")) {
    community <- fit_with(list(theta_intercept_prior = "community", mu_theta_sd = 1.5), model)
    draws <- community$fit$results_output$mu_theta_output
    expect_identical(dim(draws), c(4L, 1L))
    expect_length(community$mean, 6)
    expect_identical(community$mean[1], 0)
    expect_identical(community$mean[4:6], draws[1:3, 1])
    expect_true(all(community$mean[-1] != 0))
    # Only the intercept's prior mean moves; the slopes keep theirs.
    expect_identical(community$slope_mean, rep(0, 6))
    # mu_theta is drawn from the retained intercepts with the sampler's prior sd.
    saved <- community$fit$results_output$beta_theta_output
    for (t in 1:4) expect_identical(unname(community$intercepts[[t + 2L]]),
                                    unname(saved[1, , t, 1]))
    expect_identical(community$sd, rep(1, 6))
    expect_identical(community$hyper, rep(1.5, 6))
    # By default the sampler receives a zero mean throughout and mu_theta is
    # never drawn or saved.
    default <- fit_with(list(), model)
    expect_identical(default$mean, rep(0, 6))
    expect_length(default$intercepts, 0)
    expect_null(default$fit$results_output$mu_theta_output)
  }
})

test_that("default fits are identical with and without theta_intercept_prior = 'fixed'", {
  for (case in list(list("two_stage", FALSE), list("two_stage", TRUE),
                    list("occupancy", FALSE))) {
    model <- case[[1]]; spatial <- case[[2]]
    base <- theta_community_fit(list(), model, spatial)
    fixed <- theta_community_fit(list(theta_intercept_prior = "fixed"), model, spatial)
    expect_identical(fixed$results_output, base$results_output)
    expect_false("mu_theta_output" %in% names(base$results_output))
    expect_false("collection_intercept_prior" %in% names(base$infos))
    community <- theta_community_fit(list(theta_intercept_prior = "community"), model, spatial)
    expect_false(identical(community$results_output$beta_theta_output,
                           base$results_output$beta_theta_output))
    expect_identical(community$infos$collection_intercept_prior,
                     list(type = "community", sd = 1, mu_theta_mean = 0, mu_theta_sd = 2))
    # The occupancy intercepts keep their fixed prior.
    expect_identical(community$infos$intercept_prior, base$infos$intercept_prior)
    expect_null(community$results_output$jsdm_output$mu_b0_output)
    draws <- community$results_output$mu_theta_output
    expect_identical(dim(draws), c(15L, 2L))
    expect_true(all(is.finite(draws)))
    # The other blocks keep their shapes.
    expect_identical(lapply(community$results_output[names(base$results_output)], dim),
                     lapply(base$results_output, dim))
  }
  expect_error(theta_community_fit(list(theta_intercept_prior = "community"), "binary"),
               "requires an occupancy or two-stage model")
  expect_error(theta_community_fit(list(theta_intercept_prior = "hierarchical")),
               "theta_intercept_prior must be 'fixed' or 'community'")
  expect_error(theta_community_fit(list(mu_theta_sd = 1)),
               "mu_theta_sd requires theta_intercept_prior = 'community'")
})

test_that("both community means can be estimated together", {
  both <- theta_community_fit(list(b0_prior = "community", theta_intercept_prior = "community"))
  b0_only <- theta_community_fit(list(b0_prior = "community"))
  expect_identical(dim(both$results_output$mu_theta_output), c(15L, 2L))
  expect_identical(dim(both$results_output$jsdm_output$mu_b0_output), c(15L, 2L))
  expect_identical(both$infos$intercept_prior, b0_only$infos$intercept_prior)
  expect_identical(both$infos$collection_intercept_prior$type, "community")
  expect_false(identical(both$results_output$jsdm_output$mu_b0_output,
                         b0_only$results_output$jsdm_output$mu_b0_output))
})

test_that("mu_theta draws are thinned and reported with the other diagnostics", {
  fit <- theta_community_fit(list(theta_intercept_prior = "community", mu_theta_sd = 1.5),
                             mcmc = list(nchain = 2L, nburn = 10L, niter = 12L, nthin = 1L),
                             seed = 42)
  expect_identical(fit$infos$collection_intercept_prior$mu_theta_sd, 1.5)
  thinned <- occJSDM:::thinOutput(fit, thin = 3)
  expect_identical(thinned$results_output$mu_theta_output,
                   fit$results_output$mu_theta_output[c(1, 4, 7, 10), , drop = FALSE])
  diag <- returnConvergenceDiagnostics(fit)
  row <- diag[diag$param == "mu_theta", ]
  expect_identical(nrow(row), 1L)
  expect_equal(row$mean, mean(fit$results_output$mu_theta_output))
})

test_that("the community mean recovers a community's average collection intercept", {
  skip_on_cran()
  # One-stage data, 30 species whose collection rates are 10 to 25% (true
  # intercepts about -2.2 to -1.1) rather than near the fixed prior's 50%,
  # with 4 field samples per site and no collection covariate, so the true
  # intercepts need no rescaling. The estimated community mean should cover
  # the realised average, and the species' intercepts should be pulled
  # towards it rather than towards 0, so lie below the fixed prior's.
  #
  # Not asserted: that the community fit's average intercept error is the
  # smaller. It is not, here (-0.14 against +0.14, measured on 16 October
  # 2026), because B0's fixed Normal(0, 1) prior still pulls each rare
  # species' occupancy up and its collection intercept down, which the fixed
  # collection prior at 0 partly offset. Over 8 further seeds at 2 and at 4
  # samples per site the community fit's average intercept was below the
  # fixed fit's in all 16 and its interval for mu_theta covered the realised
  # mean in all 16.
  set.seed(20261016)
  S <- 30L; n <- 80L; M <- 4L
  sim <- simulateOccJSDMData(
    list(n = n, S = S, g = 1L, M = rep(M, n), P = 1L, K = rep(1L, n * M),
         ncov_psi = 1L, ncov_theta = 0L),
    list(p = matrix(.5, 1, S), q = matrix(.02, 1, S), theta0 = rep(.02, S),
         theta_baseline = runif(S, .1, .25)),
    list(gt = 0L, d = 0L, ds = 0, sigma_b = .5, sigma_bs = .5, sigma_ts = .5,
         sigma_h = 1, sigma_s = .5, l_s = .3, tau = rep(1, S), useSpatField = FALSE),
    model = "occupancy")
  truth <- as.vector(sim$true_params$beta_theta_true[1, ])
  run <- function(priors) {
    set.seed(5)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates = "X_psi",
      listParams = list(n_factors = 0L, n_lattrait = 0L),
      listPriors = priors,
      MCMCparams = list(nchain = 2L, nburn = 500L, niter = 500L, nthin = 1L))))
  }
  community <- run(list(theta_intercept_prior = "community"))
  fixed <- run(list())
  mu <- as.vector(community$results_output$mu_theta_output)
  expect_lt(abs(mean(mu) - mean(truth)), .5)
  expect_true(quantile(mu, .025) < mean(truth) && mean(truth) < quantile(mu, .975))
  intercept_error <- function(fit) {
    mean(apply(fit$results_output$beta_theta_output[1, , , , drop = FALSE], 2, mean)) -
      mean(truth)
  }
  expect_lt(intercept_error(community), intercept_error(fixed))
  expect_gt(intercept_error(fixed), 0)
})
