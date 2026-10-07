# Run from the repository root: Rscript dev/simstudy/spatial-design-sweep/test-oracle.R
library(testthat)
source("dev/simstudy/spatial-design-sweep/generator.R")
source("dev/simstudy/spatial-design-sweep/oracle.R")
ELLIPSE_MD5 <- "8154a611d1b7ab926eb32bf441e9af3e"

test_that("the oracle sampler source is the unmodified diagnosis file", {
  h <- oracle_source_hashes(".")
  expect_identical(unname(h[basename(names(h)) == "ellipse.cpp"]), ELLIPSE_MD5)
})

test_that("the oracle covariance is the true isotropic kernel plus jitter", {
  xy <- make_lattice(5L)
  Q <- oracle_covariance(xy, 0.3)
  expect_equal(Q, sq_exp(xy, xy, 0.3) + diag(1e-8, 25), tolerance = 1e-12)
  expect_true(min(eigen(Q, symmetric = TRUE, only.values = TRUE)$values) > 0)
})

test_that("the sampler reproduces a one-site logistic-normal posterior mean", {
  compile_oracle(".")
  y <- 1; offset <- -0.5; sd <- 1.2
  lower <- matrix(sd, 1, 1)
  set.seed(1)
  initial <- cbind(0, lower %*% matrix(rnorm(3), 1, 3))
  fit <- ellipse_draws_cpp(lower, y, 1L, offset, FALSE, initial, 2000L, 20000L)
  draws <- as.vector(fit$draws)
  numerator <- integrate(function(f) f * plogis(offset + f) * dnorm(f, 0, sd), -Inf, Inf)$value
  denominator <- integrate(function(f) plogis(offset + f) * dnorm(f, 0, sd), -Inf, Inf)$value
  expect_equal(mean(draws), numerator / denominator, tolerance = 0.02)
})

test_that("run_oracle and score_oracle return the agreed shapes", {
  compile_oracle(".")
  input <- make_sweep_input(1L)
  res <- run_oracle(input, "grid", 6L, seed = 5L, nburn = 50L, niter = 100L)
  expect_identical(dim(res$draws), c(100L, 100L, 4L))
  tr <- input$surveys$grid$truth
  sc <- score_oracle(res$draws, tr$field[, 6], tr$psi[, 6], res$offset)
  expect_true(all(c("centred_rmse", "centred_correlation", "centred_slope", "zero_field_rmse",
                    "occupancy_bias", "occupancy_mae", "flag_count") %in% names(sc$metrics)))
  expect_length(sc$field_median, 100L)
  lat <- oracle_lattice(res$draws, tr$xy, input$landscape$points[input$landscape$index$lattice, ],
                        input$landscape$range,
                        input$landscape$B0[6] + input$landscape$B[6] * input$landscape$environment[input$landscape$index$lattice],
                        thin = 10L, seed = 7L)
  expect_length(lat$probability_mean, 1600L)
  expect_true(all(lat$probability_mean > 0 & lat$probability_mean < 1))
})
cat("Oracle tests passed.\n")
