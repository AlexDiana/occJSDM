# Run from the repository root with the installed occJSDM available.
library(testthat)
suppressPackageStartupMessages(library(occJSDM))
source("dev/simstudy/spatial-design-sweep/generator.R")
source("dev/simstudy/spatial-design-sweep/score.R")

test_that("sweep groups partition sites by band and species by prevalence", {
  psi <- matrix(c(.1, .5, .9, .3), 2, 2); prevalence <- c(.05, .75)
  g <- sweep_groups(psi, prevalence)
  expect_true(all(g$low + g$medium + g$high == 1))
  expect_identical(unname(g$prevalence_5pct[, 1]), c(TRUE, TRUE))
  expect_identical(unname(g$prevalence_5pct[, 2]), c(FALSE, FALSE))
})

test_that("a tiny binary fit scores and its lattice prediction matches the native median", {
  input <- make_sweep_input(1L)
  data <- input$surveys$grid$binary
  set.seed(3)
  fit <- suppressMessages(suppressWarnings(runOccJSDM(data,
    listParams = list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = 100L),
    occCovariates = "environment", spatCovariates = c("longitude", "latitude"),
    MCMCparams = list(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L))))
  expect_identical(fit$infos$ps, 100L)
  sc <- score_sweep_fit(fit, input, "grid", "binary")
  expect_true(all(c("occupancy", "intercept", "environment_slope", "range", "spatial_sd") %in% sc$groups$metric))
  expect_identical(dim(sc$field_median), c(100L, 8L))
  expect_identical(dim(sc$lattice_mean), c(1600L, 8L))
  expect_true(sc$verification$lattice < 1e-8)
  expect_true(sc$verification$basis < 1e-10)
  expect_true(all(c("bin", "mae", "bias") %in% names(sc$lattice)))
  expect_identical(sort(unique(sc$lattice$bin)), c("0.02 to 0.05", "0.05 to 0.1", "above 0.1", "all", "up to 0.02"))
})
cat("Score tests passed.\n")
