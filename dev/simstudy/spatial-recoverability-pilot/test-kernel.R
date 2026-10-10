# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-kernel.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/kernel.R")

test_that("half-correlation distance matches the spec", {
  expect_equal(half_correlation(0.03), 0.03 * sqrt(2 * log(2)))
  expect_equal(round(half_correlation(c(0.03, 0.08, 0.15)), 3), c(0.035, 0.094, 0.177))
})

test_that("stable_factor takes the Cholesky route on a well-conditioned kernel", {
  set.seed(1); xy <- matrix(runif(200), 100, 2)
  f <- stable_factor(sq_exp(xy, xy, 0.03))
  expect_identical(f$route, "cholesky"); expect_equal(dim(f$lower), c(100L, 100L))
  expect_equal(f$lower %*% t(f$lower), sq_exp(xy, xy, 0.03) + diag(1e-8, 100), tolerance = 1e-10)
  expect_true(is.finite(f$condition))
})

test_that("stable_factor falls through to the eigen route and pads to square", {
  set.seed(2); xy <- matrix(runif(1200), 600, 2)
  K <- sq_exp(xy, xy, 0.5)   # extremely smooth: Cholesky fails or condition exceeds the ceiling
  f <- stable_factor(K)
  expect_identical(f$route, "eigen"); expect_equal(dim(f$lower), c(600L, 600L))
  expect_lt(f$error, 1e-5); expect_lt(f$rank, 600L)
  expect_equal(f$lower %*% t(f$lower), K, tolerance = 1e-5)
})

test_that("lattice_conditional reproduces the conditional covariance within tolerance", {
  set.seed(7); xy <- matrix(runif(1800), 900, 2)
  g <- (1:40 - 0.5) / 40; lat <- as.matrix(expand.grid(x = g, y = g))
  cf <- lattice_conditional(xy, lat, 0.15)
  expect_lt(cf$error, 1e-5); expect_lt(cf$rank, 1600L)
  expect_equal(dim(cf$A), c(1600L, 900L)); expect_equal(nrow(cf$lower), 1600L)
  expect_true(cf$route %in% c("cholesky", "eigen"))
})

test_that("lattice_conditional predicts a smooth field well at short range from dense sites", {
  set.seed(8); xy <- matrix(runif(1800), 900, 2)
  g <- (1:40 - 0.5) / 40; lat <- as.matrix(expand.grid(x = g, y = g))
  all <- rbind(xy, lat); L <- stable_factor(sq_exp(all, all, 0.15))$lower
  f <- as.vector(L %*% rnorm(nrow(all)))
  cf <- lattice_conditional(xy, lat, 0.15)
  pred <- as.vector(cf$A %*% f[1:900])
  expect_gt(cor(pred, f[-(1:900)]), 0.99)
})

test_that("Gauss-Hermite integrates a logistic-normal mean to 1e-6", {
  gh <- gauss_hermite(40L)
  expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-12)
  gh_value <- function(eta, sd) sum(gh$weights * plogis(eta + sqrt(2) * sd * gh$nodes)) / sqrt(pi)
  for (eta in c(-3, -0.5, 0, 1.2)) {
    exact <- integrate(function(f) plogis(eta + f) * dnorm(f, 0, 1), -Inf, Inf)$value
    expect_equal(gh_value(eta, 1), exact, tolerance = 1e-6)
  }
})
cat("Kernel tests passed.\n")
