# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-score.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("config.R", "score.R")) source(file.path(D, f))
set.seed(1)
truth <- list(psi_lattice = matrix(plogis(rnorm(1600 * 3)), 1600, 3), field_lattice = matrix(rnorm(1600 * 3), 1600, 3),
              eta_lattice = matrix(0, 1600, 3), occupancy = c(0.25, 0.5, 0.75), species = paste0("s", 1:3), spatial_sd = rep(1, 3),
              psi_training = matrix(0.5, 10, 3), field_training = matrix(0, 10, 3))
pred <- function(p, field = NULL, link = "logit") list(p_full = p, delta_p = p - 0.5, field_mean = if (is.null(field)) p * 0 else field, link = link,
                                                        site_id = paste0("l", 1:1600), species_id = truth$species)
env <- pred(matrix(0.5, 1600, 3)); oracle <- list(probability_mean = truth$psi_lattice, status = "pass")

test_that("headroom share is 1 at the oracle, 0 at the env arm, negative when worse, NA at zero headroom", {
  s1 <- score_arm(pred(truth$psi_lattice), env, oracle, truth, "d", "occJSDM", "gp", "pass")
  expect_equal(s1$headroom_share, rep(1, 3), tolerance = 1e-12)
  s0 <- score_arm(env, env, oracle, truth, "d", "occJSDM", "env", "pass"); expect_equal(s0$headroom_share, rep(0, 3))
  worse <- pred(matrix(0.9, 1600, 3)); sw <- score_arm(worse, env, oracle, truth, "d", "Hmsc", "spatial2", "pass")
  expect_true(all(sw$headroom_share < 0))
  flat_oracle <- list(probability_mean = matrix(0.5, 1600, 3), status = "pass")
  sz <- score_arm(pred(truth$psi_lattice), env, flat_oracle, truth, "d", "occJSDM", "gp", "pass"); expect_true(all(is.na(sz$headroom_share)))
})

test_that("slope and correlation are computed only for logit arms and pattern for every arm", {
  s <- score_arm(pred(truth$psi_lattice, 0.5 * truth$field_lattice, "logit"), env, oracle, truth, "d", "occJSDM", "gp", "pass")
  expect_equal(s$field_slope, rep(0.5, 3), tolerance = 1e-10); expect_equal(s$field_correlation, rep(1, 3), tolerance = 1e-10)
  h <- score_arm(pred(truth$psi_lattice, 0.5 * truth$field_lattice, "probit"), env, oracle, truth, "d", "Hmsc", "spatial2", "pass")
  expect_true(all(is.na(h$field_slope))); expect_false(any(is.na(h$pattern)))
})

test_that("unsupported fits give NA scores with their status and stay in the table", {
  s <- score_arm(NULL, env, oracle, truth, "d", "gllvm", "spatial2", "unsupported")
  expect_equal(nrow(s), 3L); expect_true(all(is.na(s$headroom_share))); expect_true(all(s$status == "unsupported"))
})

test_that("a missing env prediction gives NA shares and keeps every row with its status", {
  s <- score_arm(pred(truth$psi_lattice), NULL, oracle, truth, "d", "occJSDM", "gp", "unsupported")
  expect_equal(nrow(s), 3L); expect_true(all(is.na(s$headroom_share))); expect_true(all(is.na(s$rmse_env)))
  expect_true(all(s$status == "unsupported")); expect_equal(s$rmse, rep(0, 3))
})

test_that("aggregation pools gate occupancies and records how many species were skipped", {
  s <- rbind(score_arm(pred(truth$psi_lattice), env, oracle, truth, "d1", "occJSDM", "gp", "pass"),
             score_arm(pred(truth$psi_lattice), env, oracle, truth, "d2", "occJSDM", "gp", "pass"))
  s$headroom_share[1] <- NA
  a <- aggregate_arm(s, c(0.25, 0.5, 0.75))
  expect_equal(nrow(a), 2L); expect_equal(a$skipped_species, c(1L, 0L)); expect_equal(a$headroom_share, c(1, 1))
})

test_that("paired comparisons apply the delta rule and the minimum pairing", {
  cfg <- stage2_config()
  mk <- function(arm, h, status = "pass") data.frame(cell = "t1-r0.08-n900", dataset_id = sprintf("d%d", seq_along(h)), package = sub("/.*", "", arm),
                                                      arm = sub(".*/", "", arm), key = arm, headroom_share = h, status = status, stringsAsFactors = FALSE)
  agg <- rbind(mk("occJSDM/gp", rep(0.5, 10)), mk("Hmsc/spatial2", rep(0.5, 10) + 0.15 + rnorm(10, 0, 0.02)),
               mk("Hmsc/spatial6", rep(0.5, 10) + 0.03), mk("gllvm/spatial2", rep(0.5, 10) + 0.15, status = c(rep("pass", 5), rep("unstable", 5))))
  p <- paired_comparisons(agg, cfg)
  r <- p[p$a == "occJSDM/gp" & p$b == "Hmsc/spatial2", ]; expect_identical(r$verdict, "b clearly better"); expect_equal(r$paired, 10L)
  r6 <- p[p$a == "occJSDM/gp" & p$b == "Hmsc/spatial6", ]; expect_identical(r6$verdict, "no clear difference")
  rg <- p[p$a == "occJSDM/gp" & p$b == "gllvm/spatial2", ]; expect_identical(rg$verdict, "inconclusive"); expect_equal(rg$paired, 5L)
})

test_that("the development reading covers the four branches", {
  v <- function(t1_gp_vs_factor, t2_gp_vs_factor) data.frame(
    cell = c("t1-r0.08-n900", "t2-n900"), a = "occJSDM/gp", b = "Hmsc/spatial2", verdict = c(t1_gp_vs_factor, t2_gp_vs_factor), stringsAsFactors = FALSE)
  expect_identical(development_reading(v("b clearly better", "no clear difference"), 0.6)$reading, "shared factors")
  expect_identical(development_reading(v("no clear difference", "a clearly better"), 0.6)$reading, "per-species ranges")
  expect_identical(development_reading(v("b clearly better", "a clearly better"), 0.6)$reading, "trade-off")
  expect_identical(development_reading(v("no clear difference", "no clear difference"), 0.3)$reading, "representation not the bottleneck")
  expect_identical(development_reading(v("no clear difference", "no clear difference"), 0.7)$reading, "no clear difference")
})

test_that("the development reading follows the spec where the first draft of the rule did not", {
  v <- function(t1, t2, extra = NULL) rbind(data.frame(
    cell = c("t1-r0.08-n900", "t2-n900"), a = "occJSDM/gp", b = "Hmsc/spatial2", verdict = c(t1, t2), stringsAsFactors = FALSE), extra)
  # Factors clearly better under both truths: "not clearly worse under truth 2" holds, so the reading is shared factors.
  expect_identical(development_reading(v("b clearly better", "b clearly better"), 0.6)$reading, "shared factors")
  # A clear difference anywhere defeats "representation not the bottleneck", even between two factor arms at 300 sites.
  other <- data.frame(cell = "t2-n300", a = "Hmsc/spatial6", b = "Hmsc/spatial2", verdict = "a clearly better", stringsAsFactors = FALSE)
  expect_false(identical(development_reading(v("no clear difference", "no clear difference", other), 0.3)$reading, "representation not the bottleneck"))
  # The two qualifiers the spec attaches to every reading.
  q <- rbind(v("no clear difference", "no clear difference"),
             data.frame(cell = "t2-n900", a = "Hmsc/spatial6", b = "Hmsc/spatial2", verdict = "a clearly better", stringsAsFactors = FALSE),
             data.frame(cell = "t1-r0.08-n900", a = "Hmsc/spatial2", b = "gllvm/spatial2", verdict = "b clearly better", stringsAsFactors = FALSE))
  r <- development_reading(q, 0.6)
  expect_true(r$factor_dimension); expect_true(r$fitting_method)
  r0 <- development_reading(v("no clear difference", "no clear difference"), 0.6)
  expect_false(r0$factor_dimension); expect_false(r0$fitting_method)
})
cat("Score tests passed.\n")
