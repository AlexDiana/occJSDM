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
