# Check the refreshed data and provenance, optionally against an installation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
repo <- normalizePath(args[[1]], mustWork = TRUE)
evidence <- file.path(repo, "dev/release/2026-10-07-beta")
p <- readRDS(file.path(evidence, "provenance.rds"))
stopifnot(
  p$revision == "a56f5548e68e09d244513d09fcc04ad52f0caa19",
  p$seed == 20261007L,
  identical(p$settings$MCMCparams, list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1)),
  length(p$settings$listPriors) == 0L,
  p$requested_threads$RcppParallel == 1L,
  all(p$requested_threads$environment == "1"),
  p$fit_md5 == unname(tools::md5sum(file.path(repo, "data/sampleresults.rda"))),
  p$input_md5 == unname(tools::md5sum(file.path(repo, "data/sampledata.rda"))),
  identical(p$source_hashes, setNames(tools::md5sum(file.path(repo, names(p$source_hashes))), names(p$source_hashes)))
)
e <- new.env(parent = emptyenv())
load(file.path(repo, "data/sampleresults.rda"), envir = e)
fit <- e$sampleresults
stopifnot(
  length(fit) == 6L, fit$infos$threshold == 1,
  fit$infos$intercept_prior$sd == 1, fit$infos$n == 100L,
  fit$infos$model == "two_stage", fit$infos$n_factors == 2L,
  all(is.finite(fit$results_output$psi_output)),
  identical(dim(fit$results_output$jsdm_output$B_output), c(2L, 10L, 5000L, 2L))
)
if (length(args) == 2L) {
  lib <- normalizePath(args[[2]], mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(occJSDM))
  stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
  stopifnot(identical(occJSDM::sampleresults, fit))
  installed_vignettes <- utils::vignette(package = "occJSDM", lib.loc = lib)$results
  stopifnot(setequal(installed_vignettes[, "Item"], c("occJSDM", "simulateOccJSDMData")))
  diagnostics <- occJSDM::returnConvergenceDiagnostics(fit)
  saved <- read.csv(file.path(evidence, "diagnostics.csv"), check.names = FALSE)
  stopifnot(nrow(diagnostics) == nrow(saved))
  for (column in c("mean", "sd", "q2.5", "q97.5", "rhat", "ess")) {
    stopifnot(isTRUE(all.equal(diagnostics[[column]], saved[[column]], tolerance = 1e-12)))
  }
}
cat("Release fit, input identity, source hashes, seed, defaults and threads verified.\n")
if (length(args) == 2L) cat("Installed data and exactly two beta vignettes verified.\n")
