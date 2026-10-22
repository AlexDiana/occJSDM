# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-score.R
library(testthat)
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R", "predict.R", "score.R")) source(file.path("dev/simstudy/spatial-recoverability-pilot", f))
compile_pilot_samplers(".")

test_that("uninformed reference agrees with Monte Carlo integration and is not the plug-in", {
  set.seed(1); eta <- c(-2, 0, 1.5)
  mc <- vapply(eta, function(e) mean(plogis(e + rnorm(2e6))), numeric(1))
  expect_equal(uninformed_reference(eta), mc, tolerance = 2e-3)
  expect_true(all(abs(uninformed_reference(eta) - plogis(eta))[c(1, 3)] > 1e-3))
  expect_equal(uninformed_reference(0), 0.5, tolerance = 1e-12)
})

test_that("score_species reproduces known correlation, slope and RMSE on a synthetic case", {
  cm <- make_community("t2", 100L, 2L); lat <- cm$index$lattice; S <- 12L
  truth <- cm$spatial[lat, ]
  prediction <- list(field_mean = 0.5 * truth, probability_mean = cm$psi[lat, ],
                     site_field_mean = 0.5 * cm$spatial[cm$index$uniform, ], site_probability_mean = cm$psi[cm$index$uniform, ])
  screen <- data.frame(species = 1:S, status = c(rep("pass", 11), "unresolved"))
  sc <- score_species(cm, "uniform", prediction, screen)
  expect_equal(nrow(sc), S)
  expect_equal(sc$field_slope, rep(0.5, S), tolerance = 1e-10)
  expect_equal(sc$field_correlation, rep(1, S), tolerance = 1e-10)
  expect_equal(sc$probability_rmse_oracle, rep(0, S), tolerance = 1e-12)
  expect_equal(sc$contrast_correlation, rep(1, S), tolerance = 1e-10)
  expect_equal(sc$headroom_fraction, rep(1, S), tolerance = 1e-10)
  expect_equal(sc$headroom_pp, 100 * sc$probability_rmse_uninformed)
  expect_equal(sc$status, screen$status)
  expect_equal(sc$train_field_slope, rep(0.5, S), tolerance = 1e-10)
  expect_equal(sc$class, cm$species$class)
})

test_that("predict_lattice returns probabilities in (0,1), records the conditional and reuses it per range", {
  cm <- make_community("t2", 100L, 1L)
  o5 <- run_species_oracle(cm, "uniform", 5L, nburn = 20L, niter = 80L, thin = 4L, seed = 4L)
  o8 <- run_species_oracle(cm, "uniform", 8L, nburn = 20L, niter = 80L, thin = 4L, seed = 5L)   # same class as 5
  p <- predict_lattice(cm, "uniform", list(species = c(5L, 8L), oracles = list(o5, o8)), seed = 9L)
  expect_equal(dim(p$probability_mean), c(1600L, 2L)); expect_true(all(p$probability_mean > 0 & p$probability_mean < 1))
  expect_equal(nrow(p$conditional_records), 1L); expect_lt(max(p$conditional_records$error), 1e-5)
  expect_equal(colnames(p$field_mean), c("species05", "species08"))
})

test_that("predict_lattice for truth 1 combines blocks with the loadings", {
  cm <- make_community("t1", 100L, 1L, range = 0.08)
  o <- run_joint_oracle(cm, "clustered", nburn = 10L, niter = 20L, thin = 2L, seed = 3L)
  p <- predict_lattice(cm, "clustered", list(joint = o), seed = 9L)
  expect_equal(dim(p$field_mean), c(1600L, 12L))
  w <- joint_weights(cm); c_indep <- sqrt(1 - cm$shared_fraction)
  manual <- rowMeans(w[3, 1] * o$site_draws[[1]] + w[3, 2] * o$site_draws[[2]] + c_indep * o$site_draws[[5]])
  expect_equal(p$site_field_mean[, 3], manual)
})

test_that("layout_geometry reports nearest distances and coverage fractions", {
  cm <- make_community("t2", 100L, 1L)
  g <- layout_geometry(cm, "uniform")
  expect_true(g$within_half_short < g$within_half_medium && g$within_half_medium < g$within_half_long)
  expect_true(g$within_half_short > 0.2 && g$within_half_short < 0.45)   # spec: about 32% at n = 100
  expect_true(is.na(g$centre_separation)); expect_equal(layout_geometry(cm, "clustered")$centre_separation, 0.2)
})

test_that("gate logic counts unresolved cells as unresolved and applies every threshold", {
  th <- default_thresholds()
  make <- function(pattern, slope, headroom, status = "pass") data.frame(
    truth = "t2", cell = "t2-n300", group = "medium", n = 300L, range = NA, replicate = 1:10, arrangement = "uniform",
    contrast_correlation = pattern, field_slope = slope, headroom_pp = headroom, status = status)
  set.seed(3)
  expect_identical(gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 3, 0.2)), th)$decision, "pass")
  expect_identical(gate_cells(make(c(rep(0.7, 7), rep(0.3, 3)), rep(0.8, 10), rnorm(10, 3, 0.2)), th)$decision, "fail")
  expect_identical(gate_cells(make(rep(0.7, 10), rep(0.3, 10), rnorm(10, 3, 0.2)), th)$decision, "fail")
  expect_identical(gate_cells(make(rep(0.7, 10), rep(1.8, 10), rnorm(10, 3, 0.2)), th)$decision, "fail")
  expect_identical(gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 0.2, 0.5)), th)$decision, "fail")
  expect_identical(gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 1.05, 0.02)), th)$decision, "pass")
  expect_identical(gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 3, 0.2), status = c(rep("pass", 9), "unresolved")), th)$decision, "unresolved")
})

test_that("aggregate_scores averages within occupancy then pools, and keeps truth 2 classes apart", {
  cm <- make_community("t2", 100L, 2L); lat <- cm$index$lattice
  prediction <- list(field_mean = 0.5 * cm$spatial[lat, ], probability_mean = cm$psi[lat, ],
                     site_field_mean = 0.5 * cm$spatial[cm$index$uniform, ], site_probability_mean = cm$psi[cm$index$uniform, ])
  sc <- score_species(cm, "uniform", prediction, data.frame(species = 1:12, status = "pass"))
  agg <- aggregate_scores(sc)
  expect_equal(nrow(agg), 3L); expect_setequal(agg$group, c("short", "medium", "long"))
  expect_true(all(is.na(agg$range)))
  expect_equal(agg$field_slope, rep(0.5, 3), tolerance = 1e-10)
  sc$truth <- "t1"; sc$range <- 0.08; sc$class <- NA
  agg1 <- aggregate_scores(sc); expect_equal(nrow(agg1), 1L); expect_identical(agg1$group, "all"); expect_equal(agg1$range, 0.08)
})

test_that("read_map covers all four branches", {
  g <- function(t1_pass_n, t2_pass_n) {
    rows <- rbind(
      expand.grid(truth = "t1", range = c(0.03, 0.08, 0.15), class = NA_character_, n = c(100L, 300L, 900L), stringsAsFactors = FALSE),
      expand.grid(truth = "t2", range = NA_real_, class = c("short", "medium", "long"), n = c(100L, 300L, 900L), stringsAsFactors = FALSE))
    rows$decision <- "fail"
    rows$decision[rows$truth == "t1" & rows$range %in% 0.08 & rows$n %in% t1_pass_n] <- "pass"
    rows$decision[rows$truth == "t2" & rows$class %in% "medium" & rows$n %in% t2_pass_n] <- "pass"
    rows
  }
  expect_identical(read_map(g(c(300L, 900L), 900L))$branch, "both")
  expect_identical(read_map(g(c(300L, 900L), 900L))$recommended_n, 900L)
  expect_identical(read_map(g(c(300L, 900L), integer(0)))$branch, "truth1_only")
  expect_identical(read_map(g(c(300L, 900L), integer(0)))$recommended_n, 300L)
  expect_identical(read_map(g(integer(0), 900L))$branch, "truth2_only")
  expect_identical(read_map(g(integer(0), integer(0)))$branch, "stop")
})
cat("Score tests passed.\n")
