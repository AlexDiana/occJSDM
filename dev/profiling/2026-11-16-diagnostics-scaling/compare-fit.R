# Check on a real fit that a build prints the same convergence summary and
# stores the same output as another build.
#
# Usage:
#   Rscript --vanilla compare-fit.R fit <lib> <out.rds>     once per build
#   Rscript --vanilla compare-fit.R compare <a.rds> <b.rds>
#
# "fit" runs the quickstart call on the shipped sampledata (two-stage, two
# occupancy and two collection covariates, a spatial field, two factors) with
# two chains of 100 burn-in and 200 kept iterations from a fixed seed, and
# saves the fitted object's results_output together with everything
# runOccJSDM() printed and warned. "compare" reports whether two such files
# agree exactly.

args <- commandArgs(trailingOnly = TRUE)

if (args[[1L]] == "fit") {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  suppressPackageStartupMessages(library(occJSDM, lib.loc = args[[2L]]))
  msgs <- character()
  warns <- character()
  set.seed(20261116L)
  fit <- withCallingHandlers(
    runOccJSDM(
      data = sampledata,
      listParams = list(n_factors = 2),
      threshold = 1,
      occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
      collCovariates = c("X_theta.1", "X_theta.2"),
      spatCovariates = c("Xs.1", "Xs.2"),
      MCMCparams = list(nchain = 2, nburn = 100, niter = 200, nthin = 1,
                        ncores = 1)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # The summary starts at its first "Parameter Set" line; what comes before
  # it is the sampler's progress output.
  first <- grep("Parameter Set", msgs)[1L]
  saveRDS(list(results_output = fit$results_output,
               summary = msgs[seq(first, length(msgs))],
               warnings = warns),
          args[[3L]])
  cat(msgs[seq(first, length(msgs))], sep = "")
  cat(sprintf("%d warnings\n", length(warns)))
} else {
  a <- readRDS(args[[2L]])
  b <- readRDS(args[[3L]])
  cat("stored results_output identical:",
      identical(a$results_output, b$results_output), "\n")
  cat("printed summary identical:      ", identical(a$summary, b$summary), "\n")
  cat("warnings identical:             ", identical(a$warnings, b$warnings),
      sprintf("(%d)\n", length(a$warnings)))
}
