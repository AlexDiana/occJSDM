# Check the refitted data against its provenance, optionally against an installation.
# Rscript verify-release.R REPO [LIB]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
repo <- normalizePath(args[[1]], mustWork = TRUE)
evidence <- file.path(repo, "dev/release/2026-10-13-parallel-chains-refit")
p <- readRDS(file.path(evidence, "provenance.rds"))
stopifnot(
  p$seed == 20261007L,
  identical(p$MCMCparams, list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1, ncores = 2L)),
  isTRUE(p$identical_apart_from_execution),
  p$execution_two_cores$workers == 2L, p$execution_one_core$workers == 1L,
  identical(p$execution_two_cores$blas_threads_requested,
            c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
              MKL_DOMAIN_NUM_THREADS = "MKL_DOMAIN_ALL=1", VECLIB_MAXIMUM_THREADS = "1")),
  p$fit_md5 == unname(tools::md5sum(file.path(repo, "data/sampleresults.rda"))),
  p$input_md5 == unname(tools::md5sum(file.path(repo, "data/sampledata.rda"))),
  identical(p$source_hashes, stats::setNames(tools::md5sum(file.path(repo, names(p$source_hashes))), names(p$source_hashes))))
e <- new.env(parent = emptyenv()); load(file.path(repo, "data/sampleresults.rda"), envir = e)
fit <- e$sampleresults
stopifnot(length(fit) == 6L, fit$infos$model == "two_stage", fit$infos$n == 100L,
          identical(fit$infos$chain_seeds, p$chain_seeds),
          identical(dim(fit$results_output$jsdm_output$B_output), c(2L, 10L, 5000L, 2L)))
if (length(args) == 2L) {
  lib <- normalizePath(args[[2]], mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(occJSDM))
  stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
  stopifnot(identical(occJSDM::sampleresults, fit))
}
cat("Refit provenance, source hashes, seed, settings and data identity verified.\n")
