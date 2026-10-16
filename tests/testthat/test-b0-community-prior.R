# Experimental community mean for the occupancy intercepts (listPriors$b0_prior
# = "community"). Each species' B0 gets a Normal(mu_b0, sigma_b0^2) prior and
# mu_b0 ~ Normal(0, mu_b0_sd^2) is drawn from its conjugate conditional every
# iteration. Like sigma_b0 (test-intercept-prior.R), the prior mean is encoded
# in five places: both C++ samplers, the R reference sampler and both range
# log-weight functions. If the range update integrated B0 under a different
# mean from the one used to draw it, the range moves would target the wrong
# posterior, so each site is tested against an independent reference here. The
# default, b0_prior = "fixed", must leave fits bit for bit as they were.

test_that("read_intercept_prior keeps the default and validates the community option", {
  expect_identical(occJSDM:::read_intercept_prior(list()), list(mean = 0, sd = 1))
  expect_identical(occJSDM:::read_intercept_prior(list(b0_prior = "fixed")),
                   list(mean = 0, sd = 1))
  expect_identical(occJSDM:::read_intercept_prior(list(b0_prior = "fixed", sigma_b0 = 2)),
                   list(mean = 0, sd = 2))
  expect_identical(occJSDM:::read_intercept_prior(list(b0_prior = "community")),
                   list(type = "community", sd = 1, mu_b0_mean = 0, mu_b0_sd = 2))
  expect_identical(
    occJSDM:::read_intercept_prior(list(b0_prior = "community", sigma_b0 = 1.5,
                                        mu_b0_sd = c(named = 3))),
    list(type = "community", sd = 1.5, mu_b0_mean = 0, mu_b0_sd = 3))
  for (bad in list("Community", "estimated", NA_character_, c("fixed", "community"),
                   1, TRUE, character(0))) {
    expect_error(occJSDM:::read_intercept_prior(list(b0_prior = bad)),
                 "b0_prior must be 'fixed' or 'community'")
  }
  for (bad in list(0, -1, NA_real_, Inf, NaN, c(1, 2), numeric(0), "a", TRUE)) {
    expect_error(occJSDM:::read_intercept_prior(list(b0_prior = "community", mu_b0_sd = bad)),
                 "mu_b0_sd must be a finite positive number")
  }
  expect_error(occJSDM:::read_intercept_prior(list(mu_b0_sd = 2)),
               "mu_b0_sd requires b0_prior = 'community'")
  expect_error(occJSDM:::read_intercept_prior(list(b0_prior = "fixed", mu_b0_sd = 2)),
               "mu_b0_sd requires b0_prior = 'community'")
})

test_that("the community mean is drawn from its conjugate conditional", {
  # B0[s] ~ N(mu, s0^2) and mu ~ N(0, h^2) give mu | B0 ~ N(m, v) with
  # v = 1 / (1/h^2 + S/s0^2) and m = v * sum(B0) / s0^2.
  B0 <- c(-3.1, -2.4, -0.7, -2.9, -1.6)
  for (case in list(c(1, 2), c(.5, 2), c(1, .3))) {
    s0 <- case[1]; h <- case[2]
    v <- 1 / (1 / h^2 + length(B0) / s0^2)
    m <- v * sum(B0) / s0^2
    set.seed(7)
    one <- sample_mu_b0(B0, s0, h)
    set.seed(7)
    expect_equal(one, rnorm(1, m, sqrt(v)), tolerance = 1e-12)
    M <- 20000L
    set.seed(8)
    draws <- replicate(M, sample_mu_b0(B0, s0, h))
    info <- paste("sigma_b0 =", s0, "mu_b0_sd =", h)
    expect_lt(abs(mean(draws) - m), 5 * sqrt(v / M), label = info)
    expect_lt(abs(var(draws) / v - 1), 5 * sqrt(2 / (M - 1)), label = info)
  }
})

# Intercept-only block, as in test-intercept-prior.R.
community_intercept_only_args <- function(k, Omega, model) {
  n <- nrow(k); S <- ncol(k)
  list(k = k, X = matrix(0, n, 0), Tr = matrix(0, S, 0), U = matrix(0, n, 0),
       G = matrix(0, 0, 0), A = matrix(0, S, 0), C = matrix(0, 0, 0), sigma_b = 1,
       Gs = matrix(0, 0, 0), As = matrix(0, S, 0), Cs = matrix(0, 0, 0),
       sigma_bs = 1, Ks = matrix(0, n, 0), Xs_centers = matrix(0, n, 0),
       Omega = Omega, model = model)
}

test_that("B0 conditional draws centre on mu_b0 in every sampler", {
  # Posterior N((linear + mu/s0^2) / precision, 1 / precision) with
  # precision = sum(Omega) + 1/s0^2; linear as in test-intercept-prior.R.
  k <- cbind(c(.5, .5, .5), c(-.5, -.5, -.5))
  Omega <- cbind(c(.5, .8, 1.2), c(.4, .6, .5))
  M <- 4000L
  analytic <- function(model, s0, mu) {
    linear <- if (model == "binary") colSums(k) else colSums(k * Omega)
    precision <- colSums(Omega) + 1 / s0^2
    list(mean = (linear + mu / s0^2) / precision, var = 1 / precision)
  }
  for (model in c("binary", "continuous")) {
    # A sampler that ignored mu_b0 must not pass: the means at mu = 0 and
    # mu = -2.5 differ by more than 10 Monte Carlo SE.
    zero <- analytic(model, 1, 0); moved <- analytic(model, 1, -2.5)
    expect_true(all(abs(zero$mean - moved$mean) > 10 * sqrt(moved$var / M)))
    for (sampler_name in c("sample_BBsL_cpp", "sample_BBsL_parallel", "sample_BBsL")) {
      sampler <- get(sampler_name)
      for (case in list(c(1, -2.5), c(.5, 1.7), c(3, -4))) {
        s0 <- case[1]; mu <- case[2]
        args <- c(community_intercept_only_args(k, Omega, model),
                  list(sigma_b0 = s0, mu_b0 = mu))
        set.seed(12); setOccJSDMSeed(12)
        draws <- replicate(M, as.vector(do.call(sampler, args)$B0))
        want <- analytic(model, s0, mu)
        info <- paste(sampler_name, model, "sigma_b0 =", s0, "mu_b0 =", mu)
        expect_true(all(abs(rowMeans(draws) - want$mean) < 5 * sqrt(want$var / M)),
                    info = info)
        expect_true(all(abs(apply(draws, 1, var) / want$var - 1) <
                          5 * sqrt(2 / (M - 1))), info = info)
      }
    }
  }
})

test_that("mu_b0 moves only the intercept's prior mean", {
  # No observation precision, so each draw coordinate is its prior: B0 must
  # centre on mu_b0 while B, L and Bs keep mean zero (G and Gs are zero here),
  # ruling out an off-by-one write into b_current.
  n <- 4L; S <- 2L
  args <- list(k = matrix(0, n, S), X = matrix(0, n, 1), Tr = matrix(c(1, 2), S, 1),
               U = matrix(0, n, 1), G = matrix(0, 1, 1), A = matrix(0, S, 0),
               C = matrix(0, 0, 1), sigma_b = .4, Gs = matrix(0, 1, 2),
               As = matrix(0, S, 0), Cs = matrix(0, 0, 2), sigma_bs = .7,
               Ks = matrix(0, n, 2), Xs_centers = matrix(rep(1:2, each = n), n, 2),
               Omega = matrix(0, n, S), model = "continuous")
  M <- 3000L
  sds <- c(1.3, .4, 1, .7, .7)
  for (sampler_name in c("sample_BBsL_cpp", "sample_BBsL_parallel", "sample_BBsL")) {
    sampler <- get(sampler_name)
    set.seed(24); setOccJSDMSeed(24)
    draws <- replicate(M, {
      r <- do.call(sampler, c(args, list(sigma_b0 = 1.3, mu_b0 = -2.2)))
      c(r$B0[1], r$B[, 1], r$L[, 1], r$Bs[, 1])
    })
    want <- c(-2.2, 0, 0, 0, 0)
    expect_true(all(abs(rowMeans(draws) - want) < 5 * sds / sqrt(M)), info = sampler_name)
    expect_true(all(abs(apply(draws, 1, var) / sds^2 - 1) < 5 * sqrt(2 / (M - 1))),
                info = sampler_name)
  }
})

test_that("the samplers default to a zero intercept prior mean", {
  k <- cbind(c(.5, .5, .5), c(-.5, -.5, -.5))
  Omega <- cbind(c(.5, .8, 1.2), c(.4, .6, .5))
  args <- community_intercept_only_args(k, Omega, "binary")
  for (sampler in list(sample_BBsL_cpp, sample_BBsL_parallel, sample_BBsL)) {
    set.seed(3); setOccJSDMSeed(3)
    implicit <- do.call(sampler, args)
    set.seed(3); setOccJSDMSeed(3)
    expect_identical(do.call(sampler, c(args, list(mu_b0 = 0))), implicit)
    set.seed(3); setOccJSDMSeed(3)
    expect_false(identical(implicit, do.call(sampler, c(args, list(mu_b0 = -1)))))
  }
})

test_that("the integrated binary kernel uses the intercept prior mean", {
  # One site, one species: eta = B0 + h * Bs with B0 ~ N(mu, s0^2) and
  # Bs ~ N(.6, .8^2), so eta ~ N(mu + h * .6, s0^2 + (h * .8)^2).
  basis <- array(c(.2, .6, 1.1), c(1, 1, 3))
  grid <- c(.05, .15, .25)
  weights <- function(...) spatial_range_logweights(
    matrix(0, 1, 0), matrix(0, 1, 0), matrix(0, 0, 1), matrix(.6, 1, 1), 1, .8,
    matrix(.5, 1, 1), matrix(.7, 1, 1), matrix(1, 1, 1),
    list(Ks_all = basis, l_s_grid = grid), 2, 3, ...)
  for (case in list(c(1, -2.5), c(2, 1.4))) {
    s0 <- case[1]; mu <- case[2]
    actual <- weights(sigma_b0 = s0, mu_b0 = mu)
    expected <- vapply(seq_along(grid), function(j) {
      h <- basis[1, 1, j]
      integral <- integrate(function(eta) {
        exp(.5 * eta - .7 * eta^2 / 2) * dnorm(eta, mu + h * .6, sqrt(s0^2 + (h * .8)^2))
      }, -Inf, Inf, rel.tol = 1e-10)$value
      log(integral) + dgamma(grid[j], 2, 3, log = TRUE)
    }, numeric(1))
    expect_equal(actual - actual[1], expected - expected[1], tolerance = 1e-9)
  }
  default <- weights()
  expect_identical(default, weights(mu_b0 = 0))
  moved <- weights(mu_b0 = -2.5)
  expect_gt(max(abs((moved - moved[1]) - (default - default[1]))), 1e-4)
})

# The mixed fixture and observation-space reference of test-intercept-prior.R,
# with a prior mean for B0.
community_range_fixture <- function() {
  set.seed(2026101601)
  n <- 9L; m <- 4L; S <- 3L; p <- 2L; d <- 1L; ps <- 2L
  location <- c("s1", "s2", "s1", "s3", "s4", "s2", "s1", "s3", "s4")
  group <- match(location, unique(location))
  grid <- c(.04, .1, .2)
  basis <- array(0, c(n, ps, length(grid)))
  for (j in seq_along(grid)) basis[, , j] <- matrix(rnorm(m * ps, sd = j / 2), m, ps)[group, ]
  list(X = matrix(rnorm(n * p), n, p), U = matrix(rnorm(n * d), n, d),
       M_B = matrix(rnorm(p * S, mean = .4), p, S),
       M_Bs = matrix(rnorm(ps * S, mean = -.3), ps, S),
       sigma_b = .6, sigma_bs = 1.2,
       kappa = matrix(rnorm(n * S), n, S), Omega = matrix(exp(rnorm(n * S, sd = .8)), n, S),
       Xs_centers = matrix(rep(seq_len(ps), each = n), n, ps),
       list_SoRSummaries = list(Ks_all = basis, l_s_grid = grid),
       a_l_s = 2.3, b_l_s = 3.4, location = location)
}

community_range_reference <- function(a, s0, mu) {
  n <- nrow(a$X); p <- ncol(a$X); d <- ncol(a$U); ps <- nrow(a$M_Bs); S <- ncol(a$Omega)
  v <- c(s0^2, rep(a$sigma_b^2, p), rep(1, d), rep(a$sigma_bs^2, ps))
  m <- rbind(rep(mu, S), a$M_B, matrix(0, d, S), a$M_Bs)
  vapply(seq_along(a$list_SoRSummaries$l_s_grid), function(j) {
    H <- matrix(0, n, ps)
    Ks <- matrix(a$list_SoRSummaries$Ks_all[, , j], nrow = n)
    for (i in seq_len(n)) H[i, a$Xs_centers[i, ]] <- Ks[i, ]
    Z <- cbind(1, a$X, a$U, H)
    total <- 0
    for (s in seq_len(S)) {
      V <- diag(1 / a$Omega[, s], n) + tcrossprod(sweep(Z, 2, sqrt(v), "*"))
      L <- chol(V)
      r <- forwardsolve(t(L), a$kappa[, s] / a$Omega[, s] - drop(Z %*% m[, s]))
      total <- total - sum(log(diag(L))) - .5 * sum(r^2)
    }
    total + dgamma(a$list_SoRSummaries$l_s_grid[j], a$a_l_s, a$b_l_s, log = TRUE)
  }, numeric(1))
}

test_that("dense and grouped range log-weights integrate B0 under mu_b0", {
  a <- community_range_fixture()
  dense_args <- a[names(a) != "location"]
  for (case in list(c(1, -2.7), c(.4, 1.9))) {
    s0 <- case[1]; mu <- case[2]
    expected <- community_range_reference(a, s0, mu)
    extra <- list(sigma_b0 = s0, mu_b0 = mu)
    dense <- do.call(spatial_range_logweights, c(dense_args, extra))
    # Repeated labels and non-constant Omega route the public function to the
    # grouped implementation, so this also checks that mu_b0 is forwarded.
    routed <- do.call(spatial_range_logweights, c(a, extra))
    grouped <- do.call(spatial_range_logweights_grouped, c(a, extra))
    expect_equal(dense - dense[1], expected - expected[1], tolerance = 1e-10)
    expect_equal(grouped - grouped[1], expected - expected[1], tolerance = 1e-10)
    expect_equal(routed, grouped, tolerance = 1e-10)
    expect_equal(grouped, dense, tolerance = 1e-10)
  }
  # The default is a zero mean, and a non-zero mean moves the weights.
  zero <- community_range_reference(a, 1, 0)
  for (case in list(list(spatial_range_logweights, dense_args),
                    list(spatial_range_logweights_grouped, a))) {
    implicit <- do.call(case[[1]], case[[2]])
    expect_identical(do.call(case[[1]], c(case[[2]], list(mu_b0 = 0))), implicit)
    expect_equal(implicit - implicit[1], zero - zero[1], tolerance = 1e-10)
  }
  moved <- community_range_reference(a, 1, -2.7)
  expect_gt(max(abs((moved - moved[1]) - (zero - zero[1]))), 1e-4)
})

test_that("each iteration's draws use the community mean drawn in the iteration before", {
  # Record the mu_b0 that both sites receive. The first iteration starts at 0;
  # afterwards each iteration must use the value drawn at the end of the one
  # before, which for retained iterations is the previous saved draw. Values
  # are read by name, so the update must also pass mu_b0 by name.
  original_cpp <- sample_BBsL_cpp
  original_range <- spatial_range_logweights
  fit_with <- function(priors, spatial) {
    seen <- list(coefficients = numeric(), range = numeric())
    update <- update_jSDMcoef
    update_env <- new.env(parent = environment(update))
    update_env$sample_BBsL_cpp <- function(...) {
      seen$coefficients <<- c(seen$coefficients, list(...)$mu_b0)
      original_cpp(...)
    }
    update_env$spatial_range_logweights <- function(...) {
      seen$range <<- c(seen$range, list(...)$mu_b0)
      original_range(...)
    }
    environment(update) <- update_env
    # Mocked in the namespace, where run_one_chain() looks it up.
    local_mocked_bindings(update_jSDMcoef = update, .package = "occJSDM")
    sim <- simulate_fixture(model = "binary")
    fit <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates = fixture_occ_covariates(),
      spatCovariates = if (spatial) fixture_spat_covariates() else NULL,
      listParams = c(list(n_factors = 0L, n_lattrait = 0L),
                     if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),
      listPriors = priors,
      MCMCparams = list(nchain = 1L, nburn = 2L, niter = 4L, nthin = 1L))))
    c(seen, list(fit = fit))
  }
  spatial <- fit_with(list(b0_prior = "community"), TRUE)
  expect_length(spatial$coefficients, 6)
  expect_identical(spatial$range, spatial$coefficients)
  expect_identical(spatial$coefficients[1], 0)
  draws <- spatial$fit$results_output$jsdm_output$mu_b0_output
  expect_identical(dim(draws), c(4L, 1L))
  expect_identical(spatial$coefficients[4:6], draws[1:3, 1])
  expect_true(all(spatial$coefficients[-1] != 0))
  # Non-spatial fits never choose a range but still draw with mu_b0.
  plain <- fit_with(list(b0_prior = "community"), FALSE)
  expect_length(plain$range, 0)
  expect_identical(plain$coefficients[4:6],
                   plain$fit$results_output$jsdm_output$mu_b0_output[1:3, 1])
  # By default both sites receive a zero mean throughout and nothing is saved.
  default <- fit_with(list(), TRUE)
  expect_identical(default$coefficients, rep(0, 6))
  expect_identical(default$range, rep(0, 6))
  expect_null(default$fit$results_output$jsdm_output$mu_b0_output)
})

test_that("default fits are identical with and without b0_prior = 'fixed'", {
  sim <- simulate_fixture(model = "binary")
  run <- function(priors, spatial = FALSE) {
    set.seed(41)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates = fixture_occ_covariates(),
      spatCovariates = if (spatial) fixture_spat_covariates() else NULL,
      listParams = c(list(n_factors = 0L, n_lattrait = 0L),
                     if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),
      listPriors = priors,
      MCMCparams = list(nchain = 2L, nburn = 12L, niter = 15L, nthin = 1L))))
  }
  for (spatial in c(FALSE, TRUE)) {
    base <- run(list(), spatial)
    fixed <- run(list(b0_prior = "fixed"), spatial)
    expect_identical(fixed$results_output, base$results_output)
    expect_identical(base$infos$intercept_prior, list(mean = 0, sd = 1))
    expect_false("mu_b0_output" %in% names(base$results_output$jsdm_output))
    community <- run(list(b0_prior = "community"), spatial)
    expect_false(identical(community$results_output$jsdm_output$B0_output,
                           base$results_output$jsdm_output$B0_output))
    expect_identical(community$infos$intercept_prior,
                     list(type = "community", sd = 1, mu_b0_mean = 0, mu_b0_sd = 2))
    draws <- community$results_output$jsdm_output$mu_b0_output
    expect_identical(dim(draws), c(15L, 2L))
    expect_true(all(is.finite(draws)))
    # The other blocks keep their shapes.
    expect_identical(lapply(community$results_output$jsdm_output[names(base$results_output$jsdm_output)], dim),
                     lapply(base$results_output$jsdm_output, dim))
  }
  expect_error(run(list(b0_prior = "hierarchical")), "b0_prior must be 'fixed' or 'community'")
  expect_error(run(list(mu_b0_sd = 1)), "mu_b0_sd requires b0_prior = 'community'")
})

test_that("mu_b0 draws are thinned and reported with the other diagnostics", {
  sim <- simulate_fixture(model = "two_stage")
  set.seed(42)
  fit <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    occCovariates = fixture_occ_covariates(), collCovariates = "X_theta",
    listParams = list(n_factors = 1L, n_lattrait = 1L),
    listPriors = list(b0_prior = "community", mu_b0_sd = 1.5),
    MCMCparams = list(nchain = 2L, nburn = 10L, niter = 12L, nthin = 1L))))
  expect_identical(fit$infos$intercept_prior$mu_b0_sd, 1.5)
  thinned <- occJSDM:::thinOutput(fit, thin = 3)
  expect_identical(thinned$results_output$jsdm_output$mu_b0_output,
                   fit$results_output$jsdm_output$mu_b0_output[c(1, 4, 7, 10), , drop = FALSE])
  diag <- returnConvergenceDiagnostics(fit)
  row <- diag[diag$param == "mu_b0", ]
  expect_identical(nrow(row), 1L)
  expect_equal(row$mean, mean(fit$results_output$jsdm_output$mu_b0_output))
})

test_that("the community mean recovers a rare community's average intercept", {
  skip_on_cran()
  # Binary data, 30 species whose true intercepts are drawn around -2 rather
  # than the fixed prior's 0. The estimated community mean should cover the
  # realised average, and the species' intercepts should be pulled towards it
  # rather than towards 0. simstudy_shifted_simulator() is in helper-simstudy.R.
  set.seed(20261016)
  S <- 30L; n <- 100L
  sim <- simstudy_shifted_simulator(-2)(
    list(n = n, S = S, g = 1L, M = rep(1L, n), P = 1L, K = rep(1L, n),
         ncov_psi = 1L, ncov_theta = 1L),
    list(p = matrix(.5, 1, S), q = matrix(.02, 1, S), theta0 = rep(.05, S),
         theta_baseline = rep(.3, S)),
    list(gt = 0L, d = 0L, ds = 0, sigma_b = .5, sigma_bs = .5, sigma_ts = .5,
         sigma_h = 1, sigma_s = .5, l_s = .3, tau = rep(1, S), useSpatField = FALSE),
    model = "binary")
  # The fit centres the covariate, so its intercept is the true one plus the
  # slope times the covariate's mean.
  jp <- sim$true_params$jsdmParams_true
  B0_true <- jp$B0 + mean(sim$data_list$info$EnvCov.1) * as.vector(jp$B)
  run <- function(priors) {
    set.seed(5)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates = "EnvCov.1",
      listParams = list(n_factors = 0L, n_lattrait = 0L),
      listPriors = priors,
      MCMCparams = list(nchain = 2L, nburn = 500L, niter = 500L, nthin = 1L))))
  }
  community <- run(list(b0_prior = "community"))
  fixed <- run(list())
  mu <- as.vector(community$results_output$jsdm_output$mu_b0_output)
  expect_lt(abs(mean(mu) - mean(B0_true)), .5)
  expect_true(quantile(mu, .025) < mean(B0_true) && mean(B0_true) < quantile(mu, .975))
  # The fixed prior pulls the species towards 0, the community mean much less.
  # Measured on 16 October 2026: average error 0.05 with the community mean and
  # 0.22 with the fixed prior; posterior mean of mu_b0 -1.75 against -1.82.
  b0_error <- function(fit) mean(apply(fit$results_output$jsdm_output$B0_output, 1, mean)) -
    mean(B0_true)
  expect_lt(abs(b0_error(community)), .12)
  expect_gt(b0_error(fixed), .15)
  expect_lt(abs(b0_error(community)), abs(b0_error(fixed)))
})
