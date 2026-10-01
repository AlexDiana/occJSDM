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

  land <- input$landscape; lattice_index <- land$index$lattice; tr <- input$surveys$grid$truth
  # The distance bins partition the lattice in every arrangement, with no cell left unbinned.
  for (arr in names(input$surveys)) {
    bins <- distance_bins(land$points[lattice_index, ], input$surveys[[arr]]$truth$xy)
    expect_false(anyNA(bins$bin))
    expect_identical(sum(table(bins$bin)), 1600L)
  }

  # One group MAE by hand: the "all" occupancy row is the mean absolute error of the posterior mean probability.
  all_row <- sc$groups[sc$groups$metric == "occupancy" & sc$groups$group == "all", ]
  expect_identical(nrow(all_row), 1L)
  expect_equal(all_row$mae, mean(abs(sc$probability - tr$psi)), tolerance = 1e-12)

  # The no-spatial prediction, computed here from the fit's draws, matches predictNewSites with useSpatial = FALSE
  # (median of every draw) and the scored lattice_mean_nospatial (mean of every fourth draw) on the 50-cell subset.
  set.seed(1L); subset <- sort(sample(length(lattice_index), 50L))
  env_new <- land$environment[lattice_index][subset]
  native <- suppressMessages(predictNewSites(fit, X_psi = data.frame(environment = env_new),
    X_s = data.frame(longitude = land$points[lattice_index[subset], 1], latitude = land$points[lattice_index[subset], 2]),
    useEnvCov = TRUE, useSpatial = FALSE, useBiotic = FALSE, confidence = .95, verbose = FALSE))
  js <- fit$results_output$jsdm_output; S <- dim(js$B0_output)[1]; ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  X_new <- as.matrix(occJSDM:::transform_new_covariates(data.frame(environment = env_new),
                                                        fit$infos$list_X_psi_mat, remove_intercept = TRUE))
  draws <- array(NA_real_, c(length(subset), S, ni, nc))
  for (ch in seq_len(nc)) for (it in seq_len(ni))
    draws[, , it, ch] <- plogis(sweep(X_new %*% matrix(js$B_output[, , it, ch], ncol(X_new), S), 2, js$B0_output[, it, ch], "+"))
  expect_lt(max(abs(apply(draws, c(1, 2), median) - native[2, , ])), 1e-8)
  thinned <- draws[, , seq(1L, ni, by = 4L), , drop = FALSE]
  expect_lt(max(abs(apply(thinned, c(1, 2), mean) - sc$lattice_mean_nospatial[subset, ])), 1e-12)
})
cat("Score tests passed.\n")
