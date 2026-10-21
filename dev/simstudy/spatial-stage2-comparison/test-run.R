# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-run.R
library(testthat)
source("dev/simstudy/spatial-stage2-comparison/run.R")
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
study <- file.path(tempdir(), paste0("stage2-run-", as.integer(Sys.time())))
tiny <- list(initial = c(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L), extended = c(nchain = 2L, nburn = 40L, niter = 80L, nthin = 1L))

test_that("the timing ladder steps down in the order the plan describes", {
  cfg <- stage2_config()
  s <- lapply(1:4, function(i) unlist(ladder_settings(cfg, i)[c("support_points_900", "hmsc_method", "hmsc_knots")]))
  expect_equal(unname(s[[1]]), c("900", "Full", NA)); expect_equal(unname(s[[2]]), c("600", "Full", NA))
  expect_equal(unname(s[[3]]), c("900", "GPP", "300")); expect_equal(unname(s[[4]]), c("450", "Full", NA))
  expect_error(ladder_settings(cfg, 9L)); expect_equal(ladder_settings(cfg, 8L)$support_points_900, 300L)
})

test_that("the upper-bound price counts an extended timing fit once, from its initial cost", {
  cfg <- stage2_config()
  # Initial attempt 100 s, extension 340 s: the all-extension price is 4.43 x 100, not 4.43 x 440.
  expect_equal(job_upper_seconds(cfg, "bayesian", c(100, 340), 2L, 10), 100 * 4.43 + 10)
  # An unextended fit's prediction is scaled to the extended draw count (16,000 / 4,000).
  expect_equal(job_upper_seconds(cfg, "bayesian", 100, 1L, 10), 100 * 4.43 + 10 * 4)
  # Point arms: mean start cost times the six-start maximum.
  expect_equal(job_upper_seconds(cfg, "point", c(30, 30, 60), 1L, 5), 40 * 6 + 5)
})

test_that("a timing step measures the initial attempt only and later steps reuse arms whose settings did not change", {
  tstudy <- study; stage2_prepare(tstudy, PILOT, jobs_dir = file.path(tstudy, "jobs-local"))   # one prepare serves every test below
  sl <- 4:10
  p1 <- stage2_timing(tstudy, step = 1L, workers = 2L, mcmc = tiny, slice_sites = 40L, slice_species = sl)
  expect_equal(nrow(p1$timing), 10L); expect_true(all(p1$timing$step_measured == 1L))
  bay <- p1$timing[p1$timing$kind == "bayesian", ]
  expect_true(all(bay$attempts == 1L) && all(bay$status == "initial"))   # the extension is priced, never run, in a timing step
  expect_equal(length(list.files(file.path(tstudy, "jobs-local", "attempts"), pattern = "occJSDM-gp--timing1-2\\.rds$")), 0L)
  expect_true(file.exists(file.path(tstudy, "timing-step1.json"))); expect_true(is.logical(p1$passes))
  p2 <- stage2_timing(tstudy, step = 2L, workers = 2L, mcmc = tiny, slice_sites = 40L, slice_species = sl)
  fresh <- p2$timing[p2$timing$step_measured == 2L, ]; reused <- p2$timing[p2$timing$step_measured == 1L, ]
  expect_setequal(fresh$key, "occJSDM/gp"); expect_equal(nrow(fresh), 2L); expect_equal(nrow(reused), 8L)   # only the gp support points changed
  expect_equal(p2$settings$support_points_900, 600L)
  expect_true(p2$complete)
})

test_that("freeze records settings once and fit refuses to run before it", {
  expect_error(stage2_fit(study, workers = 1L))   # frozen settings are NA
  stage2_freeze(study, list(support_points_900 = 600L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L), timing = list(note = "test"))
  cfg <- readRDS(file.path(study, "config.rds")); expect_equal(cfg$frozen$support_points_900, 600L)
  expect_error(stage2_freeze(study, list(support_points_900 = 900L, hmsc_method = "Full", hmsc_knots = NA, prediction_thin = 1L), timing = list()))
})

test_that("fit on a subset runs every arm, and score writes every table with statuses carried", {
  j <- job_table(study); subset <- j$job_id[j$dataset_id %in% c("t2-n300-rep1", "t2-n300-rep2")]
  r <- stage2_fit(study, workers = 2L, jobs = subset, mcmc = tiny, slice_sites = 40L, slice_species = 4:10)
  expect_equal(length(list.files(file.path(study, "status"), pattern = "^t2-n300-rep[12]--.*\\.json$")), 10L)
  stage2_score(study, species_subset = 4:10)
  sp <- read.csv(file.path(study, "scores", "species.csv")); expect_equal(nrow(sp), 10L * 7L)   # seven species: Hmsc caps factors at the species count, so spatial6 needs six or more
  expect_true(all(sp$status %in% c("pass", "unresolved", "unstable", "unsupported", "no fit")))
  expect_false(any(sp$status %in% c("no fit", "unsupported")))   # every arm, both Hmsc arms included, really fitted and predicted
  agg <- read.csv(file.path(study, "scores", "aggregate.csv")); expect_equal(nrow(agg), 10L)
  pairs <- read.csv(file.path(study, "scores", "pairs.csv")); expect_true(all(pairs$verdict == "inconclusive"))   # two communities only
  expect_true(file.exists(file.path(study, "scores", "reading.json")))
  jobs <- read.csv(file.path(study, "scores", "jobs.csv")); expect_true(all(c("status", "attempts", "selected_attempt", "fit_seconds", "prediction_seconds") %in% names(jobs)))
})

test_that("a missing prediction keeps its row with NA scores and its status instead of dropping the community", {
  j <- job_table(study); id <- j$job_id[j$dataset_id == "t2-n300-rep1" & j$key == "Hmsc/spatial2"]
  expect_true(file.remove(file.path(study, "jobs", "predictions", paste0(id, ".rds"))))
  writeLines(jsonlite::toJSON(list(job_id = id, status = "unsupported"), auto_unbox = TRUE), file.path(study, "status", paste0(id, ".json")))
  stage2_score(study, species_subset = 4:10)
  sp <- read.csv(file.path(study, "scores", "species.csv")); rows <- sp[sp$dataset_id == "t2-n300-rep1" & sp$key == "Hmsc/spatial2", ]
  expect_equal(nrow(rows), 7L); expect_true(all(is.na(rows$headroom_share))); expect_true(all(rows$status == "unsupported"))
  agg <- read.csv(file.path(study, "scores", "aggregate.csv")); expect_equal(nrow(agg), 10L)
  expect_equal(agg$skipped_species[agg$dataset_id == "t2-n300-rep1" & agg$key == "Hmsc/spatial2"], 7L)
})
cat("Run tests passed.\n")
