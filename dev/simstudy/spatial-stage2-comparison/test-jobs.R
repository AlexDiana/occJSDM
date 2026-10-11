# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-jobs.R
library(testthat)
source("dev/simstudy/spatial-stage2-comparison/jobs.R")
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
study <- file.path(tempdir(), paste0("stage2-smoke-", as.integer(Sys.time())))
tiny <- list(initial = c(nchain = 2L, nburn = 20L, niter = 40L, nthin = 1L), extended = c(nchain = 2L, nburn = 40L, niter = 80L, nthin = 1L))

test_that("prepare writes config and the dataset table before any attempt", {
  stage2_prepare(study, PILOT, jobs_dir = file.path(study, "jobs-local"))
  expect_true(file.exists(file.path(study, "config.rds"))); d <- read.csv(file.path(study, "datasets.csv")); expect_equal(nrow(d), 40L)
  j <- job_table(study); expect_equal(nrow(j), 200L); expect_false(anyDuplicated(j$job_id) > 0)
  expect_error(stage2_prepare(study, PILOT, jobs_dir = file.path(study, "jobs-local")))
})

test_that("a Bayesian job runs initial, extends on screen failure, selects, predicts and stores", {
  j <- job_table(study); job <- j[j$dataset_id == "t2-n300-rep1" & j$key == "occJSDM/gp", ]
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L)
  r <- run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  expect_true(r$status %in% c("pass", "unresolved"))   # tiny budgets fail the screen; extension then follows
  expect_equal(length(list.files(file.path(study, "jobs-local", "attempts"), pattern = paste0("^", job$job_id, "-"))), 2L)
  p <- readRDS(file.path(study, "jobs-local", "predictions", paste0(job$job_id, ".rds")))
  expect_equal(dim(p$prediction$p_full), c(1600L, 4L)); expect_equal(p$prediction_thin, 10L)
  expect_identical(run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))$status, r$status)
  expect_equal(length(list.files(file.path(study, "jobs-local", "attempts"), pattern = paste0("^", job$job_id, "-"))), 2L)   # reused, not refitted
  # Selection reads small summary files, never the multi-gigabyte fits; every attempt has one and no partial file is left behind.
  expect_equal(length(list.files(file.path(study, "jobs-local", "summaries"), pattern = paste0("^", job$job_id, "-[0-9]+\\.rds$"))), 2L)
  expect_length(list.files(file.path(study, "jobs-local"), pattern = "\\.tmp$", recursive = TRUE), 0L)
})

test_that("a rerun after a crash between an attempt and its summary rebuilds the summary instead of refitting", {
  j <- job_table(study); job <- j[j$dataset_id == "t2-n300-rep1" & j$key == "occJSDM/gp", ]
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L)
  attempt <- file.path(study, "jobs-local", "attempts", paste0(job$job_id, "-2.rds")); summary <- file.path(study, "jobs-local", "summaries", paste0(job$job_id, "-2.rds"))
  before <- file.info(attempt)$mtime; unlink(summary)   # the process died after the attempt was renamed into place and before its summary was written
  r <- run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  expect_true(file.exists(summary)); expect_identical(file.info(attempt)$mtime, before)
  expect_true(r$status %in% c("pass", "unresolved"))
})

test_that("gp support points follow the frozen setting at 900 sites and the config at 300", {
  j <- job_table(study); cfg <- readRDS(file.path(study, "config.rds"))
  frozen <- list(support_points_900 = 600L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L)
  expect_equal(attempt_settings(j[j$dataset_id == "t2-n900-rep1" & j$key == "occJSDM/gp", ], 1L, cfg, frozen)$support_points, 600L)
  expect_equal(attempt_settings(j[j$dataset_id == "t2-n300-rep1" & j$key == "occJSDM/gp", ], 1L, cfg, frozen)$support_points, 300L)
})

test_that("a gp attempt with fewer support points than sites records the realised count and predicts through the fitted support set", {
  j <- job_table(study); job <- j[j$dataset_id == "t2-n300-rep2" & j$key == "occJSDM/gp", ]; cfg <- readRDS(file.path(study, "config.rds"))
  cfg$support_points_300 <- 20L   # fewer support points than the 40 sliced sites
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L)
  s <- run_attempt(study, job, 1L, cfg, frozen, tiny, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  a <- readRDS(file.path(study, "jobs-local", "attempts", paste0(job$job_id, "-1.rds")))
  expect_equal(a$fit$support_points, 20L); expect_equal(nrow(a$fit$native$infos$list_Xs$X_tilde), 20L)
  lattice_xy <- load_prediction(readRDS(job$input_file))$xy
  bases <- occ_gp_bases(a$fit$native, lattice_xy)   # the kernel is taken against the fitted support set, not the 1,600 lattice points
  expect_equal(ncol(bases[[1]]), 20L); expect_equal(nrow(bases[[1]]), 1600L)
})

test_that("a gllvm job runs three starts and selects or marks unstable", {
  j <- job_table(study); job <- j[j$dataset_id == "t2-n300-rep1" & j$key == "gllvm/spatial2", ]
  frozen <- list(support_points_900 = 300L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 1L)
  r <- run_job(study, job, mcmc = tiny, frozen = frozen, slice_sites = 40L, slice_species = c(4L, 6L, 9L, 12L))
  expect_true(r$status %in% c("pass", "unstable", "no fit"))
  expect_true(length(list.files(file.path(study, "jobs-local", "attempts"), pattern = paste0("^", job$job_id, "-"))) >= 3L)
})
cat("Job tests passed.\n")
