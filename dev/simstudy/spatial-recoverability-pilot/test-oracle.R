# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-oracle.R
library(testthat)
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R")) source(file.path("dev/simstudy/spatial-recoverability-pilot", f))
compile_pilot_samplers(".")

test_that("the sweep sampler source is unmodified", {
  expect_identical(unname(tools::md5sum("dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp")), ELLIPSE_MD5)
})

test_that("the joint sampler reduces exactly to the single-field sampler", {
  set.seed(5); n <- 30L; xy <- matrix(runif(2 * n), n, 2)
  L <- stable_factor(sq_exp(xy, xy, 0.1))$lower
  y <- as.numeric(runif(n) < 0.5); offset <- rnorm(n)
  initial <- cbind(rep(0, n), L %*% matrix(rnorm(n * 3), n, 3))
  set.seed(99); a <- ellipse_draws_cpp(L, y, 1L, offset, FALSE, initial, 20L, 50L)
  set.seed(99); b <- joint_ellipse_draws_cpp(L, matrix(y, n, 1), matrix(offset, n, 1), matrix(1, 1, 1),
                                             1L, 0, matrix(0, n, 1), initial, 20L, 50L, 1L)
  expect_equal(a$draws, b$draws, tolerance = 1e-12)
  expect_equal(a$mean_proposals, b$mean_proposals)
})

test_that("the joint sampler recovers its prior when there is no data", {
  set.seed(6); n <- 20L; S <- 3L; xy <- matrix(runif(2 * n), n, 2)
  L <- stable_factor(sq_exp(xy, xy, 0.1))$lower
  w <- matrix(rnorm(2 * S), S, 2); w <- sqrt(0.9) * w / sqrt(rowSums(w^2))
  Y <- matrix(0, n, S); offset <- matrix(-50, n, S)   # likelihood flat in the state
  d <- (2L + S) * n
  initial <- matrix(rnorm(d * 4), d, 4)
  fit <- joint_ellipse_draws_cpp(L, Y, offset, w, 2L, sqrt(0.1), matrix(0, n, S), initial, 200L, 4000L, 1L)
  draws <- matrix(fit$draws, d, 4000 * 4)
  K <- L %*% t(L)
  for (b in c(1L, 3L)) expect_equal(cov(t(draws[(b - 1L) * n + 1:n, ])), K, tolerance = 0.15)
  expect_lt(max(abs(cov(t(draws[1:n, ]), t(draws[n + 1:n, ])))), 0.15)   # blocks independent
})

test_that("the joint sampler matches quadrature on a one-site logistic-normal case", {
  L <- matrix(1.2, 1, 1); y <- matrix(1, 1, 1); offset <- matrix(-0.5, 1, 1)
  set.seed(1); initial <- cbind(0, 1.2 * matrix(rnorm(3), 1, 3))
  fit <- joint_ellipse_draws_cpp(L, y, offset, matrix(1, 1, 1), 1L, 0, matrix(0, 1, 1), initial, 2000L, 20000L, 1L)
  num <- integrate(function(f) f * plogis(-0.5 + f) * dnorm(f, 0, 1.2), -Inf, Inf)$value
  den <- integrate(function(f) plogis(-0.5 + f) * dnorm(f, 0, 1.2), -Inf, Inf)$value
  expect_equal(mean(fit$draws), num / den, tolerance = 0.02)
})

test_that("joint traces match quantities recomputed from unthinned draws", {
  cm <- make_community("t1", 100L, 1L, range = 0.08)
  o <- run_joint_oracle(cm, "uniform", nburn = 5L, niter = 12L, thin = 1L, seed = 3L)
  idx <- cm$index$uniform; w <- joint_weights(cm); c_indep <- sqrt(1 - cm$shared_fraction)
  h <- function(s, it, ch) w[s, 1] * o$site_draws[[1]][, (ch - 1L) * 12L + it] + w[s, 2] * o$site_draws[[2]][, (ch - 1L) * 12L + it] +
    c_indep * o$site_draws[[2L + s]][, (ch - 1L) * 12L + it]
  truth <- cm$spatial[idx, 4]; tc <- truth - mean(truth)
  hs <- h(4L, 7L, 2L)
  expect_equal(o$traces[7, 2, 3 * 3 + 1], mean(hs))
  expect_equal(o$traces[7, 2, 3 * 3 + 2], sum(hs * tc) / sum(tc^2))
  expect_equal(o$traces[7, 2, 3 * 3 + 3], mean(plogis(o$offset[, 4] + hs)))
})

test_that("joint draws reproduce under a fixed seed and thinning keeps every thin-th draw", {
  cm <- make_community("t1", 100L, 1L, range = 0.08)
  a <- run_joint_oracle(cm, "uniform", nburn = 10L, niter = 40L, thin = 4L, seed = 3L)
  b <- run_joint_oracle(cm, "uniform", nburn = 10L, niter = 40L, thin = 4L, seed = 3L)
  expect_identical(a$site_draws, b$site_draws)
  expect_length(a$site_draws, 14L); expect_equal(dim(a$site_draws[[1]]), c(100L, 10L * 4L))
  expect_equal(dim(a$traces), c(40L, 4L, 36L))
  full <- run_joint_oracle(cm, "uniform", nburn = 10L, niter = 40L, thin = 1L, seed = 3L)
  expect_equal(a$site_draws[[1]][, 1:10], full$site_draws[[1]][, seq(1L, 40L, by = 4L)])
})

test_that("species oracle returns the agreed shapes and a screen table", {
  cm <- make_community("t2", 100L, 1L)
  r <- run_species_oracle(cm, "clustered", species = 5L, nburn = 10L, niter = 40L, thin = 4L, seed = 2L)
  expect_length(r$site_draws, 1L); expect_equal(dim(r$site_draws[[1]]), c(100L, 40L))
  expect_equal(dim(r$traces), c(40L, 4L, 3L))
  expect_equal(r$factor_record$range, 0.08)
  sc <- screen_traces(r$traces)
  expect_equal(nrow(sc), 3L); expect_true(all(c("rhat", "ess_bulk", "ess_tail", "pass") %in% names(sc)))
})
cat("Oracle tests passed.\n")
