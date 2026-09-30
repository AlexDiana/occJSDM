# The species occupancy intercept B0 has a Normal(0, sigma_b0^2) prior. The
# default sigma_b0 = 1 is the previous hard-coded prior, so default fits must
# not change. The prior is encoded in five places (both C++ samplers, the R
# reference sampler, and both range log-weight functions); a single missed site
# would make the range update integrate B0 under a different prior from the one
# used to draw it, so each site is tested against an independent reference here.

test_that("read_intercept_prior validates and defaults", {
  expect_equal(occJSDM:::read_intercept_prior(list()), list(mean = 0, sd = 1))
  expect_equal(occJSDM:::read_intercept_prior(NULL), list(mean = 0, sd = 1))
  expect_equal(occJSDM:::read_intercept_prior(list(sigma_b0 = 3))$sd, 3)
  expect_identical(occJSDM:::read_intercept_prior(list(sigma_b0 = c(named = 2)))$sd, 2)
  for (bad in list(0, -1, NA_real_, NA_integer_, Inf, -Inf, NaN, c(1, 2), numeric(0),
                   "a", TRUE, list(1))) {
    expect_error(occJSDM:::read_intercept_prior(list(sigma_b0 = bad)),
                 "sigma_b0 must be a finite positive number")
  }
})

# Intercept-only block: n sites, S species, no covariates, factors or spatial
# terms, so sampleB_SoR() sees a single column of ones.
intercept_only_args <- function(k, Omega, model) {
  n <- nrow(k); S <- ncol(k)
  list(k = k, X = matrix(0, n, 0), Tr = matrix(0, S, 0), U = matrix(0, n, 0),
       G = matrix(0, 0, 0), A = matrix(0, S, 0), C = matrix(0, 0, 0), sigma_b = 1,
       Gs = matrix(0, 0, 0), As = matrix(0, S, 0), Cs = matrix(0, 0, 0),
       sigma_bs = 1, Ks = matrix(0, n, 0), Xs_centers = matrix(0, n, 0),
       Omega = Omega, model = model)
}

test_that("B0 conditional draws match the analytic posterior at every sigma_b0", {
  # sampleB_SoR() draws N(precision^-1 * linear, precision^-1) with
  # precision = sum(Omega) + 1/sigma_b0^2. The linear term is sum(k) for the
  # binary model (k is already the PG-augmented kappa) and sum(k * Omega) for
  # the continuous model. Three sites and modest Omega keep the likelihood weak
  # enough that sigma_b0 = 1, 3 and 0.5 give clearly different posteriors.
  k <- cbind(c(.5, .5, .5), c(-.5, -.5, -.5))
  Omega <- cbind(c(.5, .8, 1.2), c(.4, .6, .5))
  M <- 4000L
  analytic <- function(model, s0) {
    linear <- if (model == "binary") colSums(k) else colSums(k * Omega)
    precision <- colSums(Omega) + 1 / s0^2
    list(mean = linear / precision, var = 1 / precision)
  }
  for (model in c("binary", "continuous")) {
    # Monte Carlo error of M independent draws: SE(mean) = sd / sqrt(M) and
    # SE(var) = var * sqrt(2 / (M - 1)). Tolerances below are 5 SE. The design
    # must separate sigma_b0 = 1 from 3 by more than 10 SE, so that a sampler which
    # ignored sigma_b0 could not pass.
    one <- analytic(model, 1); three <- analytic(model, 3)
    expect_true(all(abs(one$mean - three$mean) > 10 * sqrt(three$var / M)))
    expect_true(all(abs(one$var / three$var - 1) > 10 * sqrt(2 / (M - 1))))
    for (sampler_name in c("sample_BBsL_cpp", "sample_BBsL_parallel", "sample_BBsL")) {
      sampler <- get(sampler_name)
      for (s0 in c(1, 3, .5)) {
        args <- c(intercept_only_args(k, Omega, model), list(sigma_b0 = s0))
        set.seed(11); setOccJSDMSeed(11)
        # The C++ samplers return B0 as an S x 1 matrix, the R sampler a vector.
        draws <- replicate(M, as.vector(do.call(sampler, args)$B0))
        want <- analytic(model, s0)
        info <- paste(sampler_name, model, "sigma_b0 =", s0)
        expect_true(all(abs(rowMeans(draws) - want$mean) < 5 * sqrt(want$var / M)),
                    info = info)
        expect_true(all(abs(apply(draws, 1, var) / want$var - 1) <
                          5 * sqrt(2 / (M - 1))), info = info)
      }
    }
  }
})

test_that("sigma_b0 sets only the intercept prior variance", {
  # With no observation precision the posterior is the prior, so each draw
  # coordinate has the variance of its own prior. B0 must follow sigma_b0 while
  # B, L and Bs keep sigma_b^2, 1 and sigma_bs^2, ruling out an off-by-one fill.
  n <- 4L; S <- 2L
  args <- list(k = matrix(0, n, S), X = matrix(0, n, 1), Tr = matrix(c(1, 2), S, 1),
               U = matrix(0, n, 1), G = matrix(0, 1, 1), A = matrix(0, S, 0),
               C = matrix(0, 0, 1), sigma_b = .4, Gs = matrix(0, 1, 2),
               As = matrix(0, S, 0), Cs = matrix(0, 0, 2), sigma_bs = .7,
               Ks = matrix(0, n, 2), Xs_centers = matrix(rep(1:2, each = n), n, 2),
               Omega = matrix(0, n, S), model = "continuous")
  M <- 3000L
  for (sampler_name in c("sample_BBsL_cpp", "sample_BBsL_parallel", "sample_BBsL")) {
    sampler <- get(sampler_name)
    for (s0 in c(1, 3)) {
      set.seed(23); setOccJSDMSeed(23)
      draws <- replicate(M, {
        r <- do.call(sampler, c(args, list(sigma_b0 = s0)))
        c(r$B0[1], r$B[, 1], r$L[, 1], r$Bs[, 1])
      })
      want <- c(s0^2, .4^2, 1, .7^2, .7^2)
      expect_true(all(abs(apply(draws, 1, var) / want - 1) < 5 * sqrt(2 / (M - 1))),
                  info = paste(sampler_name, "sigma_b0 =", s0))
    }
  }
})

test_that("the samplers default to the previous unit intercept prior", {
  k <- cbind(c(.5, .5, .5), c(-.5, -.5, -.5))
  Omega <- cbind(c(.5, .8, 1.2), c(.4, .6, .5))
  args <- intercept_only_args(k, Omega, "binary")
  for (sampler in list(sample_BBsL_cpp, sample_BBsL_parallel, sample_BBsL)) {
    set.seed(3); setOccJSDMSeed(3)
    implicit <- do.call(sampler, args)
    set.seed(3); setOccJSDMSeed(3)
    explicit <- do.call(sampler, c(args, list(sigma_b0 = 1)))
    expect_identical(implicit, explicit)
    set.seed(3); setOccJSDMSeed(3)
    expect_false(identical(implicit, do.call(sampler, c(args, list(sigma_b0 = 3)))))
  }
})

test_that("the integrated binary kernel uses the intercept prior sigma_b0", {
  # One site, one species, no covariates: eta = B0 + h * Bs with B0 ~ N(0, s0^2)
  # and Bs ~ N(.6, .8^2), so eta ~ N(h * .6, s0^2 + (h * .8)^2). Here s0 = 3.
  basis <- array(c(.2, .6, 1.1), c(1, 1, 3))
  grid <- c(.05, .15, .25)
  weights <- function(...) spatial_range_logweights(
    matrix(0, 1, 0), matrix(0, 1, 0), matrix(0, 0, 1), matrix(.6, 1, 1), 1, .8,
    matrix(.5, 1, 1), matrix(.7, 1, 1), matrix(1, 1, 1),
    list(Ks_all = basis, l_s_grid = grid), 2, 3, ...)
  actual <- weights(sigma_b0 = 3)
  expected <- vapply(seq_along(grid), function(j) {
    h <- basis[1, 1, j]
    integral <- integrate(function(eta) {
      exp(.5 * eta - .7 * eta^2 / 2) * dnorm(eta, h * .6, sqrt(9 + (h * .8)^2))
    }, -Inf, Inf, rel.tol = 1e-10)$value
    log(integral) + dgamma(grid[j], 2, 3, log = TRUE)
  }, numeric(1))
  expect_equal(actual - actual[1], expected - expected[1], tolerance = 1e-9)
  # The default is the unit prior, and sigma_b0 = 3 is distinguishable from it.
  default <- weights()
  expect_identical(default, weights(sigma_b0 = 1))
  expect_gt(max(abs((actual - actual[1]) - (default - default[1]))), 1e-4)
})

# Small mixed fixture with repeated site labels: every coefficient block
# (intercept, covariate, factor and spatial terms) is present.
intercept_range_fixture <- function() {
  set.seed(2026092901)
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

# Observation-space integration of the pseudo-data kappa / Omega under
# Normal(Z * m, diag(1 / Omega) + Z * V * Z'). Independent of the precision-form
# crossproducts used by both range log-weight functions.
intercept_range_reference <- function(a, s0) {
  n <- nrow(a$X); p <- ncol(a$X); d <- ncol(a$U); ps <- nrow(a$M_Bs); S <- ncol(a$Omega)
  v <- c(s0^2, rep(a$sigma_b^2, p), rep(1, d), rep(a$sigma_bs^2, ps))
  m <- rbind(rep(0, S), a$M_B, matrix(0, d, S), a$M_Bs)
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

test_that("dense and grouped range log-weights integrate B0 under sigma_b0", {
  a <- intercept_range_fixture()
  dense_args <- a[names(a) != "location"]
  for (s0 in c(3, .4)) {
    expected <- intercept_range_reference(a, s0)
    dense <- do.call(spatial_range_logweights, c(dense_args, list(sigma_b0 = s0)))
    # Repeated labels and non-constant Omega route the public function to the
    # grouped implementation, so this also checks that sigma_b0 is forwarded.
    routed <- do.call(spatial_range_logweights, c(a, list(sigma_b0 = s0)))
    grouped <- do.call(spatial_range_logweights_grouped, c(a, list(sigma_b0 = s0)))
    expect_equal(dense - dense[1], expected - expected[1], tolerance = 1e-10)
    expect_equal(grouped - grouped[1], expected - expected[1], tolerance = 1e-10)
    expect_equal(routed, grouped, tolerance = 1e-10)
    expect_equal(grouped, dense, tolerance = 1e-10)
  }
  # The default is the unit prior, and a non-default value moves the weights.
  unit <- intercept_range_reference(a, 1)
  for (case in list(list(spatial_range_logweights, dense_args),
                    list(spatial_range_logweights_grouped, a))) {
    implicit <- do.call(case[[1]], case[[2]])
    explicit <- do.call(case[[1]], c(case[[2]], list(sigma_b0 = 1)))
    expect_identical(implicit, explicit)
    expect_equal(implicit - implicit[1], unit - unit[1], tolerance = 1e-10)
  }
  moved <- intercept_range_reference(a, 3)
  expect_gt(max(abs((moved - moved[1]) - (unit - unit[1]))), 1e-4)
})

test_that("default fits are identical with and without an explicit sigma_b0 = 1", {
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
    expect_identical(run(list(sigma_b0 = 1), spatial)$results_output, base$results_output)
    expect_identical(base$infos$intercept_prior, list(mean = 0, sd = 1))
    wide <- run(list(sigma_b0 = 4), spatial)
    expect_false(identical(wide$results_output, base$results_output))
    expect_identical(wide$infos$intercept_prior, list(mean = 0, sd = 4))
  }
  expect_error(run(list(sigma_b0 = 0)), "sigma_b0 must be a finite positive number")
  expect_error(run(list(sigma_b0 = c(1, 2)), TRUE), "sigma_b0")
})

test_that("the MCMC coefficient update forwards sigma_b0 to every intercept site", {
  # Defaults would silently hide a missed forward, so run tiny fits with a
  # non-default value and record what each sampler actually received. The values
  # are read by name, which also requires the update to pass sigma_b0 by name.
  original_cpp <- sample_BBsL_cpp
  original_range <- spatial_range_logweights
  fit_with <- function(priors, spatial) {
    seen <- list(coefficients = numeric(), range = numeric())
    update <- update_jSDMcoef
    update_env <- new.env(parent = environment(update))
    update_env$sample_BBsL_cpp <- function(...) {
      seen$coefficients <<- c(seen$coefficients, list(...)$sigma_b0)
      original_cpp(...)
    }
    update_env$spatial_range_logweights <- function(...) {
      seen$range <<- c(seen$range, list(...)$sigma_b0)
      original_range(...)
    }
    environment(update) <- update_env
    fitter <- runOccJSDM
    fit_env <- new.env(parent = environment(fitter))
    fit_env$update_jSDMcoef <- update
    environment(fitter) <- fit_env
    sim <- simulate_fixture(model = "binary")
    suppressMessages(suppressWarnings(fitter(sim$data_list,
      occCovariates = fixture_occ_covariates(),
      spatCovariates = if (spatial) fixture_spat_covariates() else NULL,
      listParams = c(list(n_factors = 0L, n_lattrait = 0L),
                     if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),
      listPriors = priors,
      MCMCparams = list(nchain = 2L, nburn = 2L, niter = 2L, nthin = 1L))))
    seen
  }
  spatial <- fit_with(list(sigma_b0 = 4), TRUE)
  expect_gt(length(spatial$coefficients), 0)
  expect_true(all(spatial$coefficients == 4))
  expect_identical(length(spatial$range), length(spatial$coefficients))
  expect_true(all(spatial$range == 4))
  # Non-spatial fits never choose a range but still draw the intercept.
  plain <- fit_with(list(sigma_b0 = 4), FALSE)
  expect_gt(length(plain$coefficients), 0)
  expect_true(all(plain$coefficients == 4))
  expect_length(plain$range, 0)
  # Without a supplied value both sites receive the previous unit prior.
  default <- fit_with(list(), TRUE)
  expect_gt(length(default$coefficients), 0)
  expect_true(all(default$coefficients == 1))
  expect_identical(length(default$range), length(default$coefficients))
  expect_true(all(default$range == 1))
})

test_that("the coefficient update fails clearly without a stored intercept prior", {
  # A hand-built list_priors that lacks intercept_prior must not fall back to a
  # silent default or reach the samplers with a NULL prior SD.
  original_update <- update_jSDMcoef
  fitter <- runOccJSDM
  fit_env <- new.env(parent = environment(fitter))
  fit_env$update_jSDMcoef <- function(list_data, list_params, list_priors, ...) {
    list_priors$intercept_prior <- NULL
    original_update(list_data, list_params, list_priors, ...)
  }
  environment(fitter) <- fit_env
  sim <- simulate_fixture(model = "binary")
  expect_error(suppressMessages(suppressWarnings(fitter(sim$data_list,
    occCovariates = fixture_occ_covariates(),
    listParams = list(n_factors = 0L, n_lattrait = 0L),
    MCMCparams = list(nchain = 1L, nburn = 1L, niter = 1L, nthin = 1L)))),
    "list_priors$intercept_prior$sd is missing", fixed = TRUE)
})
