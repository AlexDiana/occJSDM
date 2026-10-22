# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-microtiming.R
library(testthat)
source("dev/simstudy/spatial-stage2-comparison/run.R")
# A minimal study directory: the micro-timing summary needs only the config and the dataset table, so no oracle files are hashed.
study <- file.path(tempdir(), "micro-study"); dir.create(file.path(study, "microtiming"), recursive = TRUE)
cfg <- stage2_config(); cfg$sources <- c(x = "0"); saveRDS(cfg, file.path(study, "config.rds"))
ds <- do.call(rbind, lapply(c("t1-r0.08-n300", "t1-r0.08-n900", "t2-n300", "t2-n900"), function(cell) data.frame(
  dataset_id = sprintf("%s-rep%d", cell, 1:10), cell = cell, truth = substr(cell, 1, 2), n = as.integer(sub(".*-n", "", cell)), replicate = 1:10,
  input_file = "x", input_md5 = "x", oracle_file = "x", oracle_md5 = "x", stringsAsFactors = FALSE)))
write.csv(ds, file.path(study, "datasets.csv"), row.names = FALSE)
job_for <- function(key, truth) { j <- job_table(study); j[j$key == key & j$dataset_id == paste0(if (truth == "t1") "t1-r0.08-n900" else "t2-n900", "-rep1"), ] }
fake <- function(key, truth, fit, diag, pred, kind = "bayesian", step = 1L) {
  j <- job_for(key, truth)
  saveRDS(data.frame(job_id = j$job_id, key = key, truth = truth, kind = kind, step = step, signature = timing_signature(key, ladder_settings(cfg, step)),
                     fit_seconds = fit, diagnose_seconds = diag, block_prediction_seconds = pred, gllvm_converged = NA, finished = "now", stringsAsFactors = FALSE), micro_file(study, j, step))
}

test_that("each micro fit lives in its own file and a rerun skips what exists", {
  fake("occJSDM/env", "t1", 20, 2, 0.1); expect_true(file.exists(micro_file(study, job_for("occJSDM/env", "t1"))))
  expect_match(micro_file(study, job_for("occJSDM/gp", "t1"), 8L), "step8")   # measurements at another ladder step never overwrite step 1
  expect_false(file.exists(micro_file(study, job_for("occJSDM/env", "t2"))))
})

test_that("the summary prices measured arms, lists missing ones and marks the projection incomplete", {
  fake("occJSDM/env", "t2", 10, 1, 0.1); fake("occJSDM/gp", "t1", 200, 20, 1); fake("occJSDM/gp", "t2", 100, 10, 0.5); fake("Hmsc/spatial2", "t1", 100, 10, 0.5)
  o <- stage2_microtiming_summary(study)
  expect_setequal(o$missing_arms, c("Hmsc/spatial6", "gllvm/spatial2")); expect_length(o$steps, 8L); expect_false(o$steps[[1]]$complete)
  gp1 <- o$steps[[1]]$per_arm[o$steps[[1]]$per_arm$key == "occJSDM/gp", ]
  # Larger truth is priced: fit = 200 s / 200 iterations x 28,000 iterations x 4.43 = 124,040 s.
  expect_equal(gp1$fit_hours * 3600, 200 / 200 * 28000 * 4.43, tolerance = 1e-9); expect_equal(gp1$jobs, 40L)
  expect_equal(gp1$prediction_hours * 3600, 1 / 150 * (64000 / 10) * 64, tolerance = 1e-9)
  gp2 <- o$steps[[2]]$per_arm[o$steps[[2]]$per_arm$key == "occJSDM/gp", ]; expect_equal(gp2$job_hours / gp1$job_hours, (600 / 900)^2, tolerance = 1e-9)
  h3 <- o$steps[[3]]$per_arm[o$steps[[3]]$per_arm$key == "Hmsc/spatial2", ]; h1 <- o$steps[[1]]$per_arm[o$steps[[1]]$per_arm$key == "Hmsc/spatial2", ]
  expect_equal(h3$job_hours / h1$job_hours, 300 / 900, tolerance = 1e-9)
  expect_true(file.exists(file.path(study, "microtiming.json")))
})

test_that("a measurement at a ladder step's own settings replaces the rescaling for that step", {
  fake("occJSDM/gp", "t1", 10, 1, 0.1, step = 8L); fake("occJSDM/gp", "t2", 8, 1, 0.1, step = 8L); fake("Hmsc/spatial2", "t1", 20, 2, 0.1, step = 8L)
  o <- stage2_microtiming_summary(study)
  gp8 <- o$steps[[8]]$per_arm[o$steps[[8]]$per_arm$key == "occJSDM/gp", ]
  expect_equal(gp8$fit_hours * 3600, 10 / 200 * 28000 * 4.43, tolerance = 1e-9)       # measured at 300 support points, not (300/900)^2 of step 1
  expect_match(o$steps[[8]]$basis[["occJSDM/gp"]], "measured at step 8")
  expect_match(o$steps[[2]]$basis[["occJSDM/gp"]], "rescaled")                           # no measurement at 600 support points
  gp6 <- o$steps[[6]]$per_arm[o$steps[[6]]$per_arm$key == "occJSDM/gp", ]               # step 6 is 300 support points with Hmsc Full: the gp measurement applies
  expect_equal(gp6$fit_hours, gp8$fit_hours, tolerance = 1e-9)
  h8 <- o$steps[[8]]$per_arm[o$steps[[8]]$per_arm$key == "Hmsc/spatial2", ]; expect_equal(h8$fit_hours * 3600, 20 / 200 * 28000 * 4.43, tolerance = 1e-9)
  h7 <- o$steps[[7]]$per_arm[o$steps[[7]]$per_arm$key == "Hmsc/spatial2", ]; expect_equal(h7$fit_hours, h8$fit_hours, tolerance = 1e-9)   # same GPP setting
})

test_that("a gllvm measurement is priced at the iteration cap times six starts", {
  fake("gllvm/spatial2", "t1", 60, 0, 2, kind = "point"); o <- stage2_microtiming_summary(study)
  g <- o$steps[[1]]$per_arm[o$steps[[1]]$per_arm$key == "gllvm/spatial2", ]
  expect_equal(g$fit_hours * 3600, 60 / 200 * 12000 * 6, tolerance = 1e-9); expect_equal(g$prediction_hours * 3600, 2 * 64, tolerance = 1e-9)
  expect_false("gllvm/spatial2" %in% o$missing_arms)
})
cat("Microtiming tests passed.\n")
