# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-run.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/run.R")
study <- file.path(tempdir(), paste0("pilot-smoke-", as.integer(Sys.time()))); dir.create(study)
tiny <- list(nburn = 10L, niter = 40L, thin = 4L, ext_nburn = 20L, ext_niter = 80L)

test_that("prepare writes config before inputs and inputs have hashes", {
  pilot_prepare(study, communities = 2L, cells = pilot_cells()[pilot_cells()$n == 100L, ],
                jobs_dir = file.path(study, "jobs-local"), workers = 4L)
  cfg <- readRDS(file.path(study, "config.rds"))
  inputs <- list.files(file.path(study, "inputs"), full.names = TRUE)
  expect_true(file.mtime(file.path(study, "config.rds")) < min(file.mtime(inputs)))
  expect_equal(cfg$thresholds, default_thresholds()); expect_equal(cfg$budget, default_budget())
  expect_equal(length(inputs), 4L * 2L)
  man <- read.csv(file.path(study, "input-manifest.csv"))
  expect_equal(nrow(man), 8L); expect_true(all(nchar(man$md5) == 32L))
  expect_false(anyDuplicated(c(man$seed_sites, man$seed_data)) > 0)
  expect_equal(length(list.files(file.path(study, "logs"), pattern = "^generator-conditioning")), 8L)
  expect_error(pilot_prepare(study, communities = 2L, jobs_dir = file.path(study, "jobs-local")))   # never re-prepare
})

test_that("fit runs a tiny budget, reuses existing jobs, and score writes every table", {
  jobs <- job_table(study)
  expect_equal(nrow(jobs), 4L * 2L * 2L)
  r <- pilot_fit(study, truth = "t2", workers = 2L, budget = tiny)
  expect_equal(unname(r[["fitted"]]), 4L)
  done <- list.files(file.path(study, "jobs-local"), pattern = "^t2-")
  expect_equal(length(done), 4L)
  mt <- file.mtime(file.path(study, "jobs-local", done))
  r2 <- pilot_fit(study, truth = "t2", workers = 2L, budget = tiny)
  expect_equal(unname(r2[["existing"]]), 4L)
  expect_identical(file.mtime(file.path(study, "jobs-local", done)), mt)
  r3 <- pilot_fit(study, truth = "t1", workers = 4L, budget = tiny)
  expect_equal(unname(r3[["fitted"]]), 12L)
  expect_length(list.files(file.path(study, "status"), pattern = "\\.error$"), 0L)
  job <- readRDS(file.path(study, "jobs-local", done[1]))
  expect_true(all(c("site_draws", "traces", "prediction", "scores", "screen", "screen_detail", "phases", "factor_records",
                    "conditional_records", "geometry", "budget", "elapsed_seconds") %in% names(job)))
  expect_equal(nrow(job$scores), 12L)
  t1 <- readRDS(file.path(study, "jobs-local", list.files(file.path(study, "jobs-local"), pattern = "^t1-")[1]))
  expect_equal(dim(t1$block_diagnostics$block_proposals), c(14L, 4L))
  expect_equal(dim(t1$block_diagnostics$block_mean_abs_angle), c(14L, 4L))
  expect_null(job$block_diagnostics)   # truth 2 jobs carry none
  pilot_score(study)
  sp <- read.csv(file.path(study, "scores", "species.csv"))
  expect_equal(nrow(sp), 16L * 12L)
  expect_true(all(sp$status %in% c("pass", "unresolved")))
  gate <- read.csv(file.path(study, "scores", "gate.csv"))
  expect_equal(nrow(gate), 3L + 3L)   # three truth 1 cells at n = 100 plus three classes of the truth 2 cell
  expect_true(all(gate$decision %in% c("pass", "fail", "unresolved")))
  for (f in c("aggregate.csv", "map.json", "factors.csv", "conditionals.csv", "geometry.csv", "screen.csv", "jobs.csv"))
    expect_true(file.exists(file.path(study, "scores", f)))
  expect_true(jsonlite::fromJSON(file.path(study, "scores", "map.json"))$branch %in% c("both", "truth1_only", "truth2_only", "stop"))
})

test_that("a prescribed resolution refits only the scoped fits and merges at score time", {
  scope <- list(t1_ranges = 0.08, t1_n = 100L, t2_n = 100L, t2_classes = c("medium", "long"), gate_occupancy = c(0.25, 0.5, 0.75))
  p <- pilot_prescribe_resolution(study, scope = scope, multiplier = 2L)
  expect_true(file.exists(file.path(study, "resolution.json")))
  expect_error(pilot_prescribe_resolution(study, scope = scope))   # never re-prescribe
  expect_equal(p$budget$niter, 2L * default_budget()$ext_niter); expect_equal(p$budget$thin, 2L * 4L * default_budget()$thin)
  # Under the tiny budget every fit is unresolved, so the scope selects all four range 0.08 jobs and all four truth 2 jobs.
  expect_equal(length(p$t1_jobs), 4L); expect_equal(length(p$t2_jobs), 4L)
  sp <- pilot_species(); eligible <- sp$species[sp$class %in% c("medium", "long") & sp$occupancy >= 0.25]
  expect_true(all(vapply(p$t2_refits, function(x) all(x %in% eligible), logical(1))))
  pr <- read_prescription(study)
  expect_identical(pr$t2_refits[[pr$t2_jobs[1]]], p$t2_refits[[pr$t2_jobs[1]]])
  # Run with a tiny budget rather than the prescribed one so that the test stays fast.
  tiny_res <- list(nburn = 10L, niter = 40L, thin = 4L, multiplier = 2L)
  for (j in pr$t1_jobs) expect_identical(run_resolution_job(study, j, NULL, tiny_res), "fitted")
  for (j in pr$t2_jobs) expect_identical(run_resolution_job(study, j, pr$t2_refits[[j]], tiny_res), "fitted")
  expect_identical(run_resolution_job(study, pr$t1_jobs[1], NULL, tiny_res), "existing")
  expect_equal(length(list.files(file.path(study, "jobs-resolution"), pattern = "\\.rds$")), 8L)
  j <- pr$t2_jobs[1]
  r <- readRDS(file.path(study, "jobs-resolution", paste0(j, ".rds"))); o <- readRDS(file.path(study, "jobs-local", paste0(j, ".rds")))
  keep <- !(o$scores$species %in% r$resolution$refit_species)
  expect_true(any(keep) && any(!keep))
  expect_equal(r$scores[keep, ], o$scores[keep, ], ignore_attr = TRUE)
  expect_true(all(r$phases[!keep] == "resolution")); expect_true(all(r$phases[keep] != "resolution"))
  expect_equal(r$site_draws[[which(keep)[1]]], o$site_draws[[which(keep)[1]]])
  expect_false(isTRUE(all.equal(r$site_draws[[which(!keep)[1]]], o$site_draws[[which(!keep)[1]]])))
  expect_equal(nrow(r$screen_detail), 36L); expect_equal(nrow(r$scores), 12L)
  t1 <- readRDS(file.path(study, "jobs-resolution", paste0(pr$t1_jobs[1], ".rds")))
  expect_identical(t1$phases, "resolution"); expect_length(t1$resolution$refit_species, 12L)
  expect_equal(dim(t1$block_diagnostics$block_proposals), c(14L, 4L))
  pilot_score(study)
  jobs <- read.csv(file.path(study, "scores", "jobs.csv"))
  expect_setequal(jobs$job[jobs$resolved], c(pr$t1_jobs, pr$t2_jobs))
  expect_equal(sum(jobs$resolutions[jobs$job %in% pr$t1_jobs]), 4L)   # one joint phase per truth 1 job
  expect_equal(sum(jobs$resolutions[jobs$job %in% pr$t2_jobs]), sum(lengths(pr$t2_refits)))
  expect_true(all(is.na(jobs$resolution_seconds[!jobs$resolved])) && all(jobs$resolution_seconds[jobs$resolved] > 0))
  expect_equal(nrow(read.csv(file.path(study, "scores", "species.csv"))), 16L * 12L)
})
cat("Run tests passed.\n")
