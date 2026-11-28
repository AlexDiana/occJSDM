# The unit Polya-Gamma sampler samplepg_fast() switches away from the exact
# Devroye sampler when the predictor exceeds 12 in absolute value. These tests
# check the draws on both sides of that switch against the exact moments and
# the exact Laplace transform of PG(1, z).
#
# samplepg_fast() has two copies, and R reaches them by different routes:
#   - src/jsdm.cpp, through samplePGvariables_parallel(), which returns the
#     draws themselves;
#   - src/functions.cpp, only through sample_betatheta_cpp() and
#     sample_betatheta_cpp_parallel(), which return a coefficient drawn given
#     the Polya-Gamma value and never the value itself. pg_functions_copy()
#     below recovers it exactly.

pg_mean <- function(z) tanh(z / 2) / (2 * z)

pg_var <- function(z) (sinh(z) - z) / (4 * z^3 * cosh(z / 2)^2)

# E[exp(-t * omega)] for omega ~ PG(1, z) (Polson, Scott and Windle 2013).
# Written with exponentials so that it does not overflow at large t.
pg_laplace <- function(z, t) {
  a <- abs(z) / 2
  b <- sqrt(a^2 + t / 2)
  exp(a - b) * (1 + exp(-2 * a)) / (1 + exp(-2 * b))
}

# n draws of PG(1, z) from the copy in src/jsdm.cpp.
pg_jsdm_copy <- function(z, n, seed) {
  setOccJSDMSeed(seed)
  as.vector(samplePGvariables_parallel(matrix(z, n, 1)))
}

# n draws of PG(1, z) from the copy in src/functions.cpp.
#
# One sample, an intercept-only design and n species, each occupied and each
# with its coefficient set to z, so that every species has one cell with
# predictor z. With prior mean 0 the coefficient returned is
# k / p + e / sqrt(p), where p = omega + 1 / prior_var, k = w - 1/2 and e is a
# standard normal from R's generator. The Polya-Gamma draws come from the
# package's own generator instead, so two calls under the same pair of seeds,
# one with w = 1 and one with w = 0, share omega and e and differ by exactly
# 1 / p. The prior is made flat so that p is omega to rounding at any z.
pg_functions_copy <- function(z, n, seed,
                              sampler = sample_betatheta_cpp_parallel) {
  prior_var <- 1e300
  draw <- function(w) {
    set.seed(seed)
    setOccJSDMSeed(seed)
    as.vector(sampler(w = matrix(w, 1, n), z = matrix(1, 1, n),
                      beta_theta = matrix(z, 1, n), idx_z = 1L,
                      X_theta = matrix(1, 1, 1), b_betatheta = 0,
                      B_betatheta = matrix(prior_var, 1, 1)))
  }
  1 / (draw(1) - draw(0)) - 1 / prior_var
}

# z = 11.9 uses the exact sampler and the rest the large-predictor draw, so
# 11.9 and 12.1 show that the two sides of the switch agree.
PG_TEST_Z <- c(11.9, 12.1, 13, 20, -15)
PG_TEST_N <- 200000L

# Tolerances. With n = 200,000 the standard errors are about 0.1% of the mean,
# 0.5% of the variance, and 0.1% and 0.3% of the two Laplace transform values,
# so each tolerance is roughly eight to ten standard errors. The defect these
# tests were written for halved the variance, which moves the two Laplace
# transform values by about 3% and 17% at z = 13.
expect_pg_moments <- function(omega, z) {
  label <- paste0("z = ", z)
  expect_true(all(is.finite(omega) & omega > 0), label = label)
  expect_equal(mean(omega) / pg_mean(z), 1, tolerance = 0.01, label = label)
  expect_equal(var(omega) / pg_var(z), 1, tolerance = 0.04, label = label)
  # The Laplace transform fixes the whole distribution, not only two moments.
  expect_equal(mean(exp(-2 * abs(z) * omega)) / pg_laplace(z, 2 * abs(z)), 1,
               tolerance = 0.01, label = paste0(label, ", t = 2|z|"))
  expect_equal(mean(exp(-8 * abs(z) * omega)) / pg_laplace(z, 8 * abs(z)), 1,
               tolerance = 0.02, label = paste0(label, ", t = 8|z|"))
}

test_that("the moment formulas agree with the defining series of PG(1, z)", {
  # PG(1, z) is sum_k g_k / (2 pi^2 ((k - 1/2)^2 + z^2 / (4 pi^2))) with g_k
  # independent standard exponentials, which gives its mean and variance as
  # series. This guards the closed forms the other tests rely on.
  k <- 1:2e6
  for (z in c(0.5, 3, 13, 20)) {
    d <- 2 * pi^2 * ((k - 0.5)^2 + z^2 / (4 * pi^2))
    expect_equal(sum(1 / d), pg_mean(z), tolerance = 1e-5)
    expect_equal(sum(1 / d^2), pg_var(z), tolerance = 1e-9)
  }
})

test_that("jsdm.cpp PG draws have the exact moments beyond a predictor of 12", {
  for (z in PG_TEST_Z) {
    expect_pg_moments(pg_jsdm_copy(z, PG_TEST_N, seed = 20261127L), z)
  }
})

test_that("functions.cpp PG draws have the exact moments beyond a predictor of 12", {
  for (z in PG_TEST_Z) {
    expect_pg_moments(pg_functions_copy(z, PG_TEST_N, seed = 20261127L), z)
  }
})

test_that("the two copies of the PG sampler make the same draws from one seed", {
  for (z in c(0, 2.5, PG_TEST_Z)) {
    from_jsdm <- pg_jsdm_copy(z, 2000L, seed = 77L)
    from_functions <- pg_functions_copy(z, 2000L, seed = 77L)
    # pg_functions_copy() recovers the draw through 1 / (difference) - 1, so
    # agreement is to rounding rather than bit for bit.
    expect_equal(from_functions, from_jsdm, tolerance = 1e-9,
                 label = paste0("functions.cpp draws at z = ", z))
    serial <- pg_functions_copy(z, 2000L, seed = 77L,
                                sampler = sample_betatheta_cpp)
    expect_identical(serial, from_functions)
  }
})

test_that("PG draws are finite and positive at extreme predictors", {
  for (z in c(12 + 1e-9, -12.5, 50, 700, -1e3, 1e6, 1e12, 1e100, -1e200, 1e300)) {
    from_jsdm <- pg_jsdm_copy(z, 20000L, seed = 5L)
    from_functions <- pg_functions_copy(z, 20000L, seed = 5L)
    expect_true(all(is.finite(from_jsdm) & from_jsdm > 0),
                label = paste0("jsdm.cpp draws at z = ", z))
    expect_true(all(is.finite(from_functions) & from_functions > 0),
                label = paste0("functions.cpp draws at z = ", z))
  }
})
