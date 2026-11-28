# Stage 1 of the count model: likelihood pieces tested alone. Nothing here is
# called by the sampler loop yet. Design:
# dev/superpowers/specs/2026-11-12-count-model-design.md.

test_that("the count model has one name, counts", {
  eta <- matrix(c(0, log(2), log(3), 0, 0, 0), 3, 2)
  expect_equal(occJSDM:::SD(eta, "counts"), c(stats::sd(c(1, 2, 3)), 0))
  expect_error(occJSDM:::SD(eta, "count"), "Using unavailable model")
})

test_that("Polya-Gamma moments match the defining series", {
  d <- function(k, z) 1 / (2 * pi^2 * ((k - 0.5)^2 + z^2 / (4 * pi^2)))
  for (z in c(0, 1e-4, 0.005, 0.5, 2, 12, 30)) {
    m <- pgMoments_cpp(c(1, 2.5), c(z, -z))
    expect_equal(m[, 1], c(1, 2.5) * sum(d(1:2e6, z)), tolerance = 1e-5)
    expect_equal(m[, 2], c(1, 2.5) * sum(d(1:2e6, z)^2), tolerance = 1e-5)
  }
  big <- pgMoments_cpp(c(1, 1), c(700, -5000))
  expect_true(all(is.finite(big)) && all(big > 0))
  expect_equal(big[1, ], c(1 / 1400, 1 / (2 * 700^3)), tolerance = 1e-8)
})

test_that("the negative-binomial log density matches dnbinom", {
  y <- matrix(c(0, 1, 7, 250, 1e5, 0, 3, NA), 4, 2)
  log_mean <- matrix(c(-30, -2, 1.5, 40, 11, -700, 0, 1), 4, 2)
  r <- c(0.3, 50)
  got <- nbLogPmf_cpp(y, log_mean, r)
  want <- stats::dnbinom(y, size = rep(r, each = 4), mu = exp(log_mean), log = TRUE)
  expect_equal(got, matrix(want, 4, 2), tolerance = 1e-10)
  expect_true(is.na(got[4, 2]))
  expect_identical(nbLogPmf_cpp(matrix(c(-1, 2.5), 2, 1), matrix(0, 2, 1), 1),
                   matrix(-Inf, 2, 1))
  expect_error(nbLogPmf_cpp(y, log_mean, c(0.3, 0)), "size must be positive")
  expect_error(nbLogPmf_cpp(y, log_mean, 0.3), "one size per species")
  expect_error(nbLogPmf_cpp(y, log_mean[1:2, ], r), "same dimensions")
})

test_that("at whole-number shapes the count sampler is the existing sampler, summed", {
  z <- c(-11, -2, 0, 0.3, 5, 12)
  setOccJSDMSeed(4242L)
  unit <- samplePGvariables_parallel(matrix(z, 3, 6, byrow = TRUE))
  setOccJSDMSeed(4242L)
  got <- samplePGcount_cpp(matrix(3, 1, 6), matrix(z, 1, 6))
  expect_equal(as.vector(got), colSums(unit), tolerance = 1e-12)
})

test_that("count sampler draws have the Polya-Gamma mean and variance at every kind of shape", {
  grid <- expand.grid(b = c(0.1, 0.5, 0.9, 1, 2.3, 13.5, 169.5, 171, 5000.25),
                      z = c(0, 1, -5, 12, -20, 40))
  setOccJSDMSeed(99L)
  for (i in seq_len(nrow(grid))) {
    b <- grid$b[i]; z <- grid$z[i]
    n <- if (b > 3) 2e4 else 2e5
    x <- as.vector(samplePGcount_cpp(matrix(b, n, 1), matrix(z, n, 1)))
    m <- pgMoments_cpp(b, z)
    expect_true(all(is.finite(x)) && all(x > 0), info = paste(b, z))
    expect_lt(abs(mean(x) - m[1]) / sqrt(m[2] / n), 4.5)
    expect_equal(stats::var(x), m[2], tolerance = 0.08, info = paste(b, z))
  }
})

test_that("the fractional and extreme-predictor routes have the right distribution, not only the right moments", {
  series <- function(n, b, z, K = 1000) {
    d <- 1 / (2 * pi^2 * ((seq_len(K) - 0.5)^2 + z^2 / (4 * pi^2)))
    as.vector(matrix(stats::rgamma(n * K, b), n, K) %*% d)
  }
  set.seed(20261125)
  setOccJSDMSeed(20261125L)
  for (case in list(c(0.5, 1), c(0.1, -8), c(0.3, 13), c(1, -20), c(5, 12.5))) {
    b <- case[1]; z <- case[2]
    x <- as.vector(samplePGcount_cpp(matrix(b, 5000, 1), matrix(z, 5000, 1)))
    p <- suppressWarnings(stats::ks.test(x, series(5000, b, z))$p.value)
    expect_gt(p, 0.001)
  }
})

test_that("count sampler draws are reproducible and validate their inputs", {
  shape <- matrix(c(0.3, 2, 7.5, 200), 2, 2)
  phi <- matrix(c(-1, 0, 3, -14), 2, 2)
  setOccJSDMSeed(7L); a <- samplePGcount_cpp(shape, phi)
  setOccJSDMSeed(7L); b <- samplePGcount_cpp(shape, phi)
  expect_identical(a, b)
  expect_error(samplePGcount_cpp(shape, phi[, 1, drop = FALSE]), "same dimensions")
  expect_error(samplePGcount_cpp(replace(shape, 1, 0), phi), "shape must be positive")
  expect_error(samplePGcount_cpp(shape, replace(phi, 1, NA)), "phi must be finite")
})

test_that("extreme-predictor draws stay positive and inverse Gaussian at tiny shapes", {
  # Added after the stage 1 review; not in the plan. pinvgauss() is the
  # distribution function of an inverse Gaussian with mean m and shape
  # parameter lambda.
  pinvgauss <- function(x, m, lambda) {
    stats::pnorm(sqrt(lambda / x) * (x / m - 1)) +
      exp(2 * lambda / m + stats::pnorm(-sqrt(lambda / x) * (x / m + 1), log.p = TRUE))
  }
  setOccJSDMSeed(1L)
  for (b in c(1e-6, 1e-8, 1e-9, 1e-12)) {
    x <- as.vector(samplePGcount_cpp(matrix(b, 2e5, 1), matrix(-13, 2e5, 1)))
    expect_true(all(is.finite(x)) && all(x > 0), info = paste(b))
    p <- suppressWarnings(stats::ks.test(x, pinvgauss, m = b / 26, lambda = b^2 / 4)$p.value)
    expect_gt(p, 0.001)
  }
  tiny <- samplePGcount_cpp(matrix(1e-200, 100, 1), matrix(40, 100, 1))
  expect_true(all(is.finite(tiny)) && all(tiny > 0))
})

test_that("count_log_mean adds each sample's log effort to its site's predictor", {
  eta <- matrix(1:6, 3, 2)
  expect_equal(occJSDM:::count_log_mean(eta, c(1L, 1L, 3L), log(c(1, 2, 4))),
               eta[c(1, 1, 3), ] + log(c(1, 2, 4)))
  expect_error(occJSDM:::count_log_mean(eta, c(1L, 4L), c(0, 0)), "idx_site")
  expect_error(occJSDM:::count_log_mean(eta, c(1L, 2L), c(0, NA)), "log_effort")
})

test_that("the size prior is a half-normal on the inverse square root of the size", {
  log_r <- seq(-3, 5, length.out = 9)
  r <- exp(log_r)
  direct <- log(2 * stats::dnorm(r^(-1/2), 0, 0.7) * 0.5 * r^(-1/2))
  gap <- direct - occJSDM:::count_size_logprior(log_r, 0.7)
  expect_equal(gap, rep(gap[1], 9), tolerance = 1e-10)
})

size_run <- function(y, log_mean, include = matrix(TRUE, nrow(y), ncol(y)),
                     iterations = 3000, burn = 1000) {
  S <- ncol(y)
  r <- rep(1, S); step <- rep(occJSDM:::count_size_initial_step, S); t <- 0L
  keep <- matrix(NA_real_, iterations - burn, S)
  accepted <- keep
  for (it in seq_len(iterations)) {
    o <- occJSDM:::count_size_step(r, y, log_mean, include, 1, step, it <= burn, t)
    r <- o$r; step <- o$step; t <- o$t
    if (it > burn) { keep[it - burn, ] <- r; accepted[it - burn, ] <- o$accepted }
  }
  list(draws = keep, acceptance = colMeans(accepted), step = step, t = t)
}

test_that("the size update recovers the size from counts with a known mean", {
  set.seed(20261125)
  N <- 400
  centre <- c(2, 2, 0, 0, -2, -2)
  truth <- c(0.5, 5, 0.5, 5, 0.5, 5)
  log_mean <- sapply(centre, function(m) stats::rnorm(N, m, 1))
  y <- sapply(1:6, function(s) stats::rnbinom(N, size = truth[s], mu = exp(log_mean[, s])))
  fit <- size_run(y, log_mean)
  q <- apply(fit$draws, 2, stats::quantile, c(0.005, 0.5, 0.995))
  expect_true(all(truth >= q[1, ] & truth <= q[3, ]))
  expect_true(all(abs(log(q[2, 1:4] / truth[1:4])) < log(1.5)))
  expect_true(all(fit$acceptance > 0.25 & fit$acceptance < 0.65))
  expect_identical(fit$t, 1000L)
})

test_that("the size update ignores missing and excluded cells and survives a species never counted", {
  set.seed(1)
  N <- 60
  log_mean <- matrix(stats::rnorm(N * 3), N, 3)
  y <- matrix(stats::rnbinom(N * 3, size = 2, mu = exp(log_mean)), N, 3)
  y[, 3] <- 0
  y[1:10, 1] <- NA
  include <- matrix(TRUE, N, 3)
  include[11:20, 2] <- FALSE
  changed <- y
  changed[11:20, 2] <- 999
  set.seed(5); a <- occJSDM:::count_size_step(rep(1, 3), y, log_mean, include, 1, rep(0.5, 3), TRUE, 0L)
  set.seed(5); b <- occJSDM:::count_size_step(rep(1, 3), changed, log_mean, include, 1, rep(0.5, 3), TRUE, 0L)
  expect_identical(a, b)
  expect_true(all(is.finite(a$r)) && all(a$r > 0))
  expect_identical(a$t, 1L)
  fixed <- occJSDM:::count_size_step(rep(1, 3), y, log_mean, include, 1, rep(0.5, 3), FALSE, 7L)
  expect_identical(fixed$step, rep(0.5, 3))
  expect_identical(fixed$t, 7L)
  expect_error(occJSDM:::count_size_step(rep(1, 2), y, log_mean, include, 1, rep(0.5, 3), TRUE, 0L),
               "one size per species")
})
