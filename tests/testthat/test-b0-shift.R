# The collapsed shift of mu_b0 and the species' occupancy intercepts, for the
# experimental b0_prior = "community" in occupancy and two-stage fits. Each
# iteration, before the latent occupancy states z are drawn, a common shift
# delta ~ Normal(0, sd^2) of every B0 and of mu_b0 is proposed and accepted by
# Metropolis under the first-stage likelihood with z integrated out, site by
# site and species by species, times the hyperprior ratio of mu_b0 (the
# Normal(mu_b0, sigma_b0^2) prior terms of the B0 cancel because the
# deviations B0 - mu_b0 are unchanged). The z draw that follows, from its
# full conditional at the shifted predictor, is what makes the step valid.
# Design: dev/simstudy/mu-b0-drift-20261018/REPORT.md. The settings live in
# MCMCparams (b0_shift, b0_shift_sd, b0_shift_proposals); the proposal sd is
# adapted during burn-in unless set, and the acceptance rates are recorded in
# infos$b0_shift. Fits with the default prior, and community fits with the
# move switched off, must leave the sampler untouched.

community_prior <- function() list(type = "community", sd = 1, mu_b0_mean = 0, mu_b0_sd = 2)
fixed_prior <- function() list(mean = 0, sd = 1)

test_that("read_b0_shift is NULL unless the community prior is on in an occupancy or two-stage fit", {
  for (params in list(list(), NULL, list(nchain = 2))) {
    expect_null(occJSDM:::read_b0_shift(params, fixed_prior(), "two_stage"))
    expect_null(occJSDM:::read_b0_shift(params, fixed_prior(), "occupancy"))
    expect_null(occJSDM:::read_b0_shift(params, community_prior(), "binary"))
    expect_null(occJSDM:::read_b0_shift(params, community_prior(), "continuous"))
  }
  for (model in c("occupancy", "two_stage")) {
    expect_identical(occJSDM:::read_b0_shift(list(nchain = 2), community_prior(), model),
                     list(enabled = TRUE, sd = NULL, proposals = 5L))
    expect_identical(occJSDM:::read_b0_shift(list(b0_shift = TRUE), community_prior(), model),
                     list(enabled = TRUE, sd = NULL, proposals = 5L))
    expect_identical(occJSDM:::read_b0_shift(list(b0_shift = FALSE), community_prior(), model),
                     list(enabled = FALSE))
  }
  expect_identical(
    occJSDM:::read_b0_shift(list(b0_shift_sd = c(named = .3), b0_shift_proposals = 2),
                            community_prior(), "occupancy"),
    list(enabled = TRUE, sd = .3, proposals = 2L))
})

test_that("read_b0_shift validates its settings", {
  for (key in c("b0_shift", "b0_shift_sd", "b0_shift_proposals")) {
    params <- stats::setNames(list(if (key == "b0_shift") TRUE else 1), key)
    expect_error(occJSDM:::read_b0_shift(params, fixed_prior(), "two_stage"),
                 "requires b0_prior = 'community' in an occupancy or two-stage fit")
    expect_error(occJSDM:::read_b0_shift(params, community_prior(), "binary"),
                 "requires b0_prior = 'community' in an occupancy or two-stage fit")
  }
  for (bad in list(NA, c(TRUE, FALSE), "TRUE", 1, logical(0))) {
    expect_error(occJSDM:::read_b0_shift(list(b0_shift = bad), community_prior(), "occupancy"),
                 "b0_shift must be TRUE or FALSE")
  }
  for (bad in list(0, -1, NA_real_, Inf, NaN, c(1, 2), numeric(0), "a", TRUE)) {
    expect_error(occJSDM:::read_b0_shift(list(b0_shift_sd = bad), community_prior(), "occupancy"),
                 "b0_shift_sd must be a finite positive number")
  }
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), numeric(0), "a", TRUE)) {
    expect_error(occJSDM:::read_b0_shift(list(b0_shift_proposals = bad), community_prior(), "occupancy"),
                 "b0_shift_proposals must be a single whole number of at least 1")
  }
  for (extra in list(list(b0_shift_sd = .2), list(b0_shift_proposals = 3))) {
    expect_error(occJSDM:::read_b0_shift(c(list(b0_shift = FALSE), extra), community_prior(), "occupancy"),
                 "b0_shift_sd and b0_shift_proposals require b0_shift = TRUE")
  }
})

# A tiny first stage: 3 sites with 2, 1 and 3 field samples, 2 species.
tiny_first_stage <- function() {
  set.seed(2026101801)
  n <- 3L; S <- 2L; M <- c(2L, 1L, 3L); N <- sum(M)
  idx_z_w <- rep(seq_len(n), M)
  list(n = n, S = S, N = N, idx_z_w = idx_z_w,
       w = matrix(rbinom(N * S, 1, .5), N, S),
       theta = matrix(runif(N * S, .2, .8), N, S),
       theta0 = runif(S, .02, .2),
       eta = matrix(rnorm(n * S, -.5, 1), n, S))
}

# The first-stage likelihood of w with z summed out by enumerating every one
# of the 2^(n S) occupancy configurations.
brute_force_loglik <- function(eta, w, theta, theta0, idx_z_w) {
  n <- nrow(eta); S <- ncol(eta)
  psi <- plogis(eta)
  configs <- as.matrix(expand.grid(rep(list(0:1), n * S)))
  lp <- apply(configs, 1, function(zv) {
    z <- matrix(zv, n, S)
    total <- sum(dbinom(z, 1, psi, log = TRUE))
    for (r in seq_len(nrow(w))) for (s in seq_len(S)) {
      p <- if (z[idx_z_w[r], s] == 1) theta[r, s] else theta0[s]
      total <- total + dbinom(w[r, s], 1, p, log = TRUE)
    }
    total
  })
  m <- max(lp)
  m + log(sum(exp(lp - m)))
}

test_that("the per-pair log-likelihoods are those of the site's field samples given z", {
  a <- tiny_first_stage()
  pair <- b0_shift_pair_loglik(a$w, a$theta, a$theta0, a$idx_z_w, a$n)
  expect_identical(dim(pair$l1), c(a$n, a$S))
  expect_identical(dim(pair$l0), c(a$n, a$S))
  for (i in seq_len(a$n)) for (s in seq_len(a$S)) {
    rows <- which(a$idx_z_w == i)
    expect_equal(pair$l1[i, s], sum(dbinom(a$w[rows, s], 1, a$theta[rows, s], log = TRUE)),
                 tolerance = 1e-12)
    expect_equal(pair$l0[i, s], sum(dbinom(a$w[rows, s], 1, a$theta0[s], log = TRUE)),
                 tolerance = 1e-12)
  }
})

test_that("the per-pair log-likelihoods are the terms the chain's own first-stage likelihood forms", {
  # computeModelLoglikFirstStage_cpp() scores the first stage in the chain
  # (R/chains.R), and ZProbWorker (src/functions.cpp) forms the same per-site
  # sums for the z draw. The step integrates z out of exactly these terms, so
  # its l1 and l0 must equal their sums over each site's field samples at
  # z = 1 and z = 0, to rounding.
  a <- tiny_first_stage()
  pair <- b0_shift_pair_loglik(a$w, a$theta, a$theta0, a$idx_z_w, a$n)
  by_site <- function(zval) {
    z <- matrix(zval, a$n, a$S)
    ll <- computeModelLoglikFirstStage_cpp(a$w, z, a$theta, a$theta0, a$idx_z_w)
    # one value per field sample and species, species fastest
    rowsum(matrix(ll, nrow(a$w), a$S, byrow = TRUE), a$idx_z_w, reorder = TRUE)
  }
  expect_equal(unname(pair$l1), unname(by_site(1)), tolerance = 1e-12)
  expect_equal(unname(pair$l0), unname(by_site(0)), tolerance = 1e-12)
})

test_that("the z-integrated log-likelihood equals the sum over every occupancy configuration", {
  a <- tiny_first_stage()
  pair <- b0_shift_pair_loglik(a$w, a$theta, a$theta0, a$idx_z_w, a$n)
  for (shift in c(0, -1.3, 2.1)) {
    expect_equal(b0_shift_loglik(a$eta + shift, pair$l1, pair$l0),
                 brute_force_loglik(a$eta + shift, a$w, a$theta, a$theta0, a$idx_z_w),
                 tolerance = 1e-10)
  }
  # Extreme predictors must not overflow or lose the small term.
  big <- matrix(c(40, -40, 40, -40, 0, 0), a$n, a$S)
  expect_true(is.finite(b0_shift_loglik(big, pair$l1, pair$l0)))
  expect_equal(b0_shift_loglik(big, pair$l1, pair$l0),
               brute_force_loglik(big, a$w, a$theta, a$theta0, a$idx_z_w), tolerance = 1e-10)
})

test_that("the acceptance ratio is the brute-force likelihood ratio times the hyperprior ratio", {
  a <- tiny_first_stage()
  pair <- b0_shift_pair_loglik(a$w, a$theta, a$theta0, a$idx_z_w, a$n)
  mu_b0 <- -.8
  for (hyper in list(list(mu_b0_mean = 0, mu_b0_sd = 2), list(mu_b0_mean = .5, mu_b0_sd = .7))) {
    for (delta in c(-.7, .3, 1.9)) {
      want <- brute_force_loglik(a$eta + delta, a$w, a$theta, a$theta0, a$idx_z_w) -
        brute_force_loglik(a$eta, a$w, a$theta, a$theta0, a$idx_z_w) +
        dnorm(mu_b0 + delta, hyper$mu_b0_mean, hyper$mu_b0_sd, log = TRUE) -
        dnorm(mu_b0, hyper$mu_b0_mean, hyper$mu_b0_sd, log = TRUE)
      got <- b0_shift_logratio(delta, a$eta, pair$l1, pair$l0, mu_b0, hyper)
      expect_equal(got$logratio, want, tolerance = 1e-10)
      expect_equal(got$proposed, brute_force_loglik(a$eta + delta, a$w, a$theta, a$theta0, a$idx_z_w),
                   tolerance = 1e-10)
      # A cached current value gives the same ratio.
      cached <- b0_shift_logratio(delta, a$eta, pair$l1, pair$l0, mu_b0, hyper,
                                  current = b0_shift_loglik(a$eta, pair$l1, pair$l0))
      expect_equal(cached$logratio, want, tolerance = 1e-10)
    }
  }
  # The hyperprior centre is used, not assumed to be 0.
  centred <- b0_shift_logratio(.3, a$eta, pair$l1, pair$l0, mu_b0, list(mu_b0_mean = 0, mu_b0_sd = 2))
  moved <- b0_shift_logratio(.3, a$eta, pair$l1, pair$l0, mu_b0, list(mu_b0_mean = 1, mu_b0_sd = 2))
  expect_false(isTRUE(all.equal(centred$logratio, moved$logratio)))
})

test_that("the step proposes a common shift from R's RNG and accepts it by Metropolis", {
  a <- tiny_first_stage()
  pair <- b0_shift_pair_loglik(a$w, a$theta, a$theta0, a$idx_z_w, a$n)
  hyper <- list(mu_b0_mean = 0, mu_b0_sd = 2)
  B0 <- c(-.4, .9); mu_b0 <- .2; sd <- .6
  # The step draws, per proposal, delta then u; replay that here.
  replay <- function(seed, proposals, adapt, t) {
    set.seed(seed)
    b0 <- B0; mu <- mu_b0; eta <- a$eta; cur_sd <- sd; accepted <- 0L
    for (k in seq_len(proposals)) {
      delta <- rnorm(1, 0, cur_sd); u <- runif(1)
      lr <- b0_shift_logratio(delta, eta, pair$l1, pair$l0, mu, hyper)$logratio
      if (log(u) < lr) { b0 <- b0 + delta; mu <- mu + delta; eta <- eta + delta; accepted <- accepted + 1L }
      if (adapt) {
        t <- t + 1L
        cur_sd <- exp(log(cur_sd) + (min(1, exp(lr)) - .44) / sqrt(t))
      }
    }
    list(B0 = b0, mu_b0 = mu, eta = eta, accepted = accepted, sd = cur_sd, t = t)
  }
  for (case in list(list(seed = 3, proposals = 1L, adapt = FALSE, t = 0L),
                    list(seed = 4, proposals = 3L, adapt = FALSE, t = 0L),
                    list(seed = 5, proposals = 4L, adapt = TRUE, t = 0L),
                    list(seed = 6, proposals = 2L, adapt = TRUE, t = 37L))) {
    set.seed(case$seed)
    got <- b0_shift_step(B0, mu_b0, a$eta, pair$l1, pair$l0, hyper, sd = sd,
                         proposals = case$proposals, adapt = case$adapt, t = case$t)
    want <- replay(case$seed, case$proposals, case$adapt, case$t)
    expect_equal(got, want, tolerance = 1e-12)
    # The deviations from the community mean never change.
    expect_equal(got$B0 - got$mu_b0, B0 - mu_b0, tolerance = 1e-12)
    expect_equal(got$eta - a$eta, matrix(got$mu_b0 - mu_b0, a$n, a$S), tolerance = 1e-12)
    expect_true(got$accepted >= 0L && got$accepted <= case$proposals)
    if (!case$adapt) expect_identical(got$sd, sd)
  }
  # Over many proposals at a sensible scale some are accepted and some rejected.
  set.seed(8)
  many <- b0_shift_step(B0, mu_b0, a$eta, pair$l1, pair$l0, hyper, sd = sd,
                        proposals = 200L, adapt = FALSE, t = 0L)
  expect_gt(many$accepted, 20L)
  expect_lt(many$accepted, 180L)
})

# The occupancy fixture with 3 field samples per site; no spatial term or
# factors, so a fit takes well under a second.
shift_fixture <- function() cached_fixture("shift_occupancy", function() {
  simulate_fixture(model = "occupancy", M = 3L)
})

fit_shift <- function(priors, mcmc, seed = 41, model = NULL, spatial = FALSE) {
  sim <- if (is.null(model)) shift_fixture() else simulate_fixture(model = model, M = 3L)
  set.seed(seed)
  suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    occCovariates = fixture_occ_covariates(),
    collCovariates = if (is.null(model) || model %in% c("occupancy", "two_stage")) "X_theta" else NULL,
    spatCovariates = if (spatial) fixture_spat_covariates() else NULL,
    listParams = c(list(n_factors = 0L, n_lattrait = 0L),
                   if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),
    listPriors = priors, MCMCparams = mcmc)))
}

test_that("the step never runs with the fixed prior, with the move off, or without latent states", {
  local_mocked_bindings(b0_shift_step = function(...) stop("the shift step ran"), .package = "occJSDM")
  mcmc <- list(nchain = 2L, nburn = 5L, niter = 6L, nthin = 1L)
  fixed <- fit_shift(list(), mcmc)
  expect_null(fixed$infos$b0_shift)
  expect_false("b0_shift" %in% names(fixed$infos))
  off <- fit_shift(list(b0_prior = "community"), c(mcmc, list(b0_shift = FALSE)))
  expect_identical(off$infos$b0_shift, list(enabled = FALSE))
  binary <- fit_shift(list(b0_prior = "community"), mcmc, model = "binary")
  expect_null(binary$infos$b0_shift)
  expect_true(all(is.finite(binary$results_output$jsdm_output$mu_b0_output)))
})

test_that("the step runs once per iteration with the community prior and is recorded", {
  calls <- list()
  original <- b0_shift_step
  local_mocked_bindings(b0_shift_step = function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args[c("sd", "proposals", "adapt", "t")]
    original(...)
  }, .package = "occJSDM")
  S <- FIXTURE_S
  mcmc <- list(nchain = 2L, nburn = 20L, niter = 10L, nthin = 2L)
  fit <- fit_shift(list(b0_prior = "community"), mcmc)
  iterations <- mcmc$nburn + mcmc$niter * mcmc$nthin
  expect_length(calls, mcmc$nchain * iterations)
  expect_true(all(vapply(calls, function(x) x$proposals == 5L, TRUE)))
  # Adapted during burn-in only, from 2 / sqrt(S); the proposal counter runs
  # on through burn-in and stops with it.
  per_chain <- split(calls, rep(seq_len(mcmc$nchain), each = iterations))
  for (ch in per_chain) {
    expect_identical(ch[[1]]$sd, 2 / sqrt(S))
    expect_identical(ch[[1]]$t, 0L)
    expect_identical(vapply(ch, function(x) x$adapt, TRUE),
                     c(rep(TRUE, mcmc$nburn), rep(FALSE, mcmc$niter * mcmc$nthin)))
    expect_identical(ch[[mcmc$nburn + 1L]]$t, 5L * mcmc$nburn)
    kept_sd <- vapply(ch[-seq_len(mcmc$nburn)], function(x) x$sd, 0)
    expect_true(all(kept_sd == kept_sd[1]))
    expect_false(kept_sd[1] == 2 / sqrt(S))
  }
  rec <- fit$infos$b0_shift
  expect_identical(names(rec), c("enabled", "proposals", "adapted", "sd", "acceptance", "acceptance_burnin"))
  expect_identical(rec$enabled, TRUE)
  expect_identical(rec$proposals, 5L)
  expect_identical(rec$adapted, TRUE)
  expect_identical(rec$sd, vapply(per_chain, function(ch) ch[[iterations]]$sd, 0, USE.NAMES = FALSE))
  for (field in c("acceptance", "acceptance_burnin")) {
    expect_length(rec[[field]], mcmc$nchain)
    expect_true(all(rec[[field]] >= 0 & rec[[field]] <= 1))
  }
  # The move changes the draws: the same seed with the move off gives others.
  off <- fit_shift(list(b0_prior = "community"), c(mcmc, list(b0_shift = FALSE)))
  expect_false(identical(off$results_output$jsdm_output$B0_output,
                         fit$results_output$jsdm_output$B0_output))
  expect_false(identical(off$results_output$jsdm_output$mu_b0_output,
                         fit$results_output$jsdm_output$mu_b0_output))
  # Everything else keeps its shape.
  expect_identical(lapply(fit$results_output$jsdm_output, dim), lapply(off$results_output$jsdm_output, dim))
  expect_identical(fit$infos$intercept_prior, off$infos$intercept_prior)
})

test_that("a user-set proposal sd is used as it is, with the given number of proposals", {
  calls <- list()
  original <- b0_shift_step
  local_mocked_bindings(b0_shift_step = function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args[c("sd", "proposals", "adapt")]
    original(...)
  }, .package = "occJSDM")
  mcmc <- list(nchain = 2L, nburn = 6L, niter = 4L, nthin = 1L, b0_shift_sd = .3, b0_shift_proposals = 2)
  fit <- fit_shift(list(b0_prior = "community"), mcmc)
  expect_length(calls, 20L)
  expect_true(all(vapply(calls, function(x) identical(x$sd, .3) && x$proposals == 2L && !x$adapt, TRUE)))
  rec <- fit$infos$b0_shift
  expect_identical(rec[c("enabled", "proposals", "adapted", "sd")],
                   list(enabled = TRUE, proposals = 2L, adapted = FALSE, sd = c(.3, .3)))
  # Without burn-in nothing is adapted and the burn-in acceptance is undefined.
  none <- fit_shift(list(b0_prior = "community"), list(nchain = 1L, nburn = 0L, niter = 4L, nthin = 1L))
  expect_identical(none$infos$b0_shift$sd, 2 / sqrt(FIXTURE_S))
  expect_identical(none$infos$b0_shift$adapted, TRUE)
  expect_true(is.na(none$infos$b0_shift$acceptance_burnin))
  expect_false(is.na(none$infos$b0_shift$acceptance))
})

test_that("the settings are rejected where the move cannot apply", {
  mcmc <- list(nchain = 1L, nburn = 2L, niter = 2L, nthin = 1L)
  expect_error(fit_shift(list(), c(mcmc, list(b0_shift = TRUE))),
               "b0_shift requires b0_prior = 'community' in an occupancy or two-stage fit")
  expect_error(fit_shift(list(b0_prior = "community"), c(mcmc, list(b0_shift_sd = .2)), model = "binary"),
               "b0_shift_sd requires b0_prior = 'community' in an occupancy or two-stage fit")
})

test_that("a fit with the move runs on workers and records every chain", {
  skip_on_cran()
  mcmc <- list(nchain = 2L, nburn = 8L, niter = 6L, nthin = 1L)
  parallel <- fit_shift(list(b0_prior = "community"), c(mcmc, list(ncores = 2L)), model = "two_stage")
  rec <- parallel$infos$b0_shift
  expect_identical(rec$enabled, TRUE)
  expect_length(rec$sd, 2L)
  expect_true(all(is.finite(rec$sd)) && all(rec$sd > 0))
  expect_length(rec$acceptance, 2L)
  expect_true(all(rec$acceptance >= 0 & rec$acceptance <= 1))
})

# Whether a fit is identical in process and on workers depends on the BLAS,
# not on this step, which does no BLAS work. The other cross-worker identity
# tests (test-parallel-chains.R) are gated by the declaration set after
# checking that a single-threaded reference BLAS is in use; the gate is
# repeated here because that file's helper is local to it.
test_that("a fit with the move is identical in process and on workers", {
  skip_on_cran()
  if (!nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")))
    skip("fit identity across worker counts needs a declared single-threaded reference BLAS")
  mcmc <- list(nchain = 2L, nburn = 8L, niter = 6L, nthin = 1L)
  serial <- fit_shift(list(b0_prior = "community"), mcmc, model = "two_stage")
  parallel <- fit_shift(list(b0_prior = "community"), c(mcmc, list(ncores = 2L)), model = "two_stage")
  expect_identical(occJSDM:::drop_execution(parallel), occJSDM:::drop_execution(serial))
})

test_that("long chains with and without the move agree on the posterior of mu_b0 and B0", {
  skip_on_cran()
  # An occupancy fit with 8 species at 60 sites, 3 field samples each: the
  # move changes the sampler, not its target, so the posterior means and sds
  # of mu_b0 and of every B0 must agree within Monte Carlo error, taken from
  # each chain pair's own bulk ESS. The chains without the move mix mu_b0
  # slowly, which widens their error and so the tolerance.
  sim <- simulate_fixture(model = "occupancy", n = 60L, S = 8L, M = 3L, g = 1L, gt = 0L, d = 0L)
  run <- function(extra) {
    set.seed(2026101802)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates = fixture_occ_covariates(), collCovariates = "X_theta",
      listParams = list(n_factors = 0L, n_lattrait = 0L),
      listPriors = list(b0_prior = "community"),
      MCMCparams = c(list(nchain = 2L, nburn = 300L, niter = 2000L, nthin = 1L), extra))))
  }
  with_move <- run(list())
  without <- run(list(b0_shift = FALSE))
  draws <- function(fit) {
    jo <- fit$results_output$jsdm_output
    c(list(mu_b0 = jo$mu_b0_output),
      lapply(seq_len(dim(jo$B0_output)[1]), function(s) jo$B0_output[s, , ]))
  }
  a <- draws(with_move); b <- draws(without)
  expect_gt(with_move$infos$b0_shift$acceptance[1], .1)
  for (k in seq_along(a)) {
    ess_a <- posterior::ess_bulk(a[[k]]); ess_b <- posterior::ess_bulk(b[[k]])
    sd_a <- sd(as.vector(a[[k]])); sd_b <- sd(as.vector(b[[k]]))
    mean_tol <- 4 * sqrt(sd_a^2 / ess_a + sd_b^2 / ess_b)
    sd_tol <- 4 * sqrt(1 / (2 * ess_a) + 1 / (2 * ess_b))
    info <- sprintf("block %d: means %.3f and %.3f (tolerance %.3f), sds %.3f and %.3f, ESS %.0f and %.0f",
                    k, mean(a[[k]]), mean(b[[k]]), mean_tol, sd_a, sd_b, ess_a, ess_b)
    expect_lt(abs(mean(a[[k]]) - mean(b[[k]])), mean_tol, label = info)
    expect_lt(abs(log(sd_a / sd_b)), sd_tol, label = info)
  }
})
