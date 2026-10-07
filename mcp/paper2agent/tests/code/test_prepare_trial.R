# Run from the repository root, with a library containing this source revision:
# Rscript test_prepare_trial.R PINNED_LIBRARY
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
.libPaths(c(normalizePath(args[[1]]), .libPaths()))
library(testthat)
library(occJSDM)

test_that("reference preparation rejects a library missing its pinned package", {
  scratch <- ".superpowers/sdd/2026-10-05-occjsdm-paper2agent-mcp/test-runs"
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  empty_library <- tempfile("empty-library-", tmpdir = scratch)
  dir.create(empty_library)
  output <- tempfile("must-not-be-created-", tmpdir = scratch)
  run <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c("mcp/paper2agent/tests/data/prepare-trial.R", output, empty_library)),
    stdout = TRUE, stderr = TRUE))
  expect_true(!is.null(attr(run, "status")) && attr(run, "status") != 0L)
  expect_match(paste(run, collapse = "\n"), "Pinned occJSDM installation is missing", fixed = TRUE)
  expect_false(dir.exists(output))
})

test_that("trial artifacts agree with native R computations and retain unavailable diagnostics", {
  scratch <- ".superpowers/sdd/2026-10-05-occjsdm-paper2agent-mcp/test-runs"
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  output <- tempfile("occjsdm-trial-", tmpdir = scratch)
  script <- "mcp/paper2agent/tests/data/prepare-trial.R"
  run <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(script, output, args[[1]])), stdout = TRUE, stderr = TRUE))
  status <- attr(run, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(run, collapse = "\n"))
  if (status == 0L) {
    data <- readRDS(file.path(output, "student", "data.rds"))
    expect_identical(nrow(data$info), nrow(data$OTU))
    expect_identical(dim(data$OTU), c(160L, 4L))
    expect_identical(getFromNamespace("inferDataModel", "occJSDM")(data), "two_stage")
    bad <- readRDS(file.path(output, "student", "malformed.rds"))
    expect_identical(nrow(bad$info) - nrow(bad$OTU), 1L)
    no_traits <- readRDS(file.path(output, "student", "no-traits.rds"))
    expect_null(no_traits$traits)
    expect_identical(no_traits$info, data$info)
    expect_identical(no_traits$OTU, data$OTU)
    counts <- readRDS(file.path(output, "student", "unreplicated-counts.rds"))
    expect_error(getFromNamespace("inferDataModel", "occJSDM")(counts),
                 "Counts model not supported yet", fixed = TRUE)

    fit <- readRDS(file.path(output, "student", "fit-short.rds"))
    diagnostics <- read.csv(file.path(output, "grader", "diagnostics.csv"))
    native <- as.data.frame(returnConvergenceDiagnostics(fit))
    expect_identical(names(diagnostics), names(native))
    for (column in c("mean", "sd", "q2.5", "q97.5", "rhat", "ess")) {
      expect_equal(diagnostics[[column]], native[[column]], tolerance = 1e-12)
    }
    expect_identical(diagnostics$param, native$param)
    expect_true(any(is.finite(native$rhat) & native$rhat > 1.01) ||
                any(is.finite(native$ess) & native$ess < 400))
    summary <- read.csv(file.path(output, "grader", "baseline-occupancy.csv"))
    draws <- returnOccupancyRates(fit)
    expect_identical(summary$species, colnames(draws))
    expect_equal(summary$mean, unname(colMeans(draws)), tolerance = 1e-12)
    expect_equal(summary$q2.5, unname(apply(draws, 2, quantile, 0.025)), tolerance = 1e-12)
    expect_equal(summary$q97.5, unname(apply(draws, 2, quantile, 0.975)), tolerance = 1e-12)
    unavailable <- read.csv(file.path(output, "student", "diagnostics-unavailable.csv"))
    expect_true(all(is.na(unavailable$rhat)))
    provenance <- jsonlite::read_json(file.path(output, "grader", "reference.json"))
    expect_identical(provenance$scientific_interpretation, FALSE)
    expect_identical(provenance$fit_settings$MCMCparams$nchain, 2L)
    expect_identical(provenance$package_path, find.package("occJSDM"))
    expect_identical(provenance$input_sha256,
                     digest::digest(file = file.path(output, "student", "data.rds"), algo = "sha256"))
    project <- file.path(output, "student-project")
    expect_true(dir.exists(project), info = "A separate student project must exclude grader materials.")
    if (dir.exists(project)) {
      files <- list.files(project, recursive = TRUE, all.files = TRUE)
      expect_false(any(grepl("^(mcp|tests|dev|grader|\\.git|\\.superpowers)/", files)))
      expect_true(all(c("R/runOccJSDM.R", "src/jsdm.cpp", "trial/data.rds",
                        "teaching/assistant-instructions.md") %in% files))
      manifest <- jsonlite::read_json(file.path(output, "grader", "student-project-manifest.json"))
      expect_setequal(files, names(manifest$files_sha256))
      for (name in files) expect_identical(digest::digest(file = file.path(project, name), algo = "sha256"),
                                           manifest$files_sha256[[name]])
    }
  }
})
