# Time computeDiagnostics() and record its peak memory on one synthetic block.
#
# Usage:
#   Rscript --vanilla bench-diagnostics.R <label> <n_params> <lib> [niter] [nchain]
#
# <label> is free text for the log ("before", "after"). <lib> is the library
# that holds the occJSDM build to measure, so the same script measures any
# commit: install an export of that commit there (see run-bench.sh).
#
# The block is the only parameter set in a minimal results_output: an
# [n_params, niter, nchain] array stored as jsdm_output$B0_output, which is
# also where computeDiagnostics() reads niter and nchain. The draws are an
# AR(1) series per parameter and chain (coefficient 0.5), built in place so
# that the process holds exactly one copy of the array when the timing starts.
#
# Three memory figures are reported, all in MB:
#   input_mb      the size of the block itself (8 bytes per draw)
#   heap_peak_mb  the most R heap in use at any point of the call, less what
#                 was in use just before it (gc() "max used" after a reset):
#                 what computeDiagnostics() allocated on top of its input
#   rss_base_mb   resident set size of the process just before the call
# run-bench.sh adds the process peak from /usr/bin/time -l, which also counts
# memory that R's heap accounting cannot see (C++ work buffers).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("usage: bench-diagnostics.R <label> <n_params> <lib> [niter] [nchain]")
}
label <- args[[1L]]
n_params <- as.integer(args[[2L]])
lib <- args[[3L]]
niter <- if (length(args) >= 4L) as.integer(args[[4L]]) else 5000L
nchain <- if (length(args) >= 5L) as.integer(args[[5L]]) else 2L

suppressPackageStartupMessages(library(occJSDM, lib.loc = lib))

set.seed(20261116L)
x <- stats::rnorm(as.double(n_params) * niter * nchain)
dim(x) <- c(n_params, niter, nchain)
for (t in seq_len(niter)[-1L]) {
  x[, t, ] <- 0.5 * x[, t - 1L, ] + x[, t, ]
}
results_output <- list(jsdm_output = list(B0_output = x))
rm(x)

heap_mb <- function(g, column) sum(g[, which(colnames(g) == column) + 1L])
rss_mb <- function() {
  as.numeric(system2("ps", c("-o", "rss=", "-p", Sys.getpid()),
                     stdout = TRUE)) / 1024
}

invisible(gc(reset = TRUE))
heap_base <- heap_mb(gc(), "used")
rss_base <- rss_mb()

msgs <- character()
warns <- character()
elapsed <- system.time(
  withCallingHandlers(
    computeDiagnostics(results_output),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
)
heap_peak <- heap_mb(gc(), "max used")

cat(sprintf(paste0(
  "RESULT label=%s occJSDM=%s n_params=%d niter=%d nchain=%d ",
  "elapsed_s=%.2f user_s=%.2f input_mb=%.0f heap_peak_mb=%.0f rss_base_mb=%.0f\n"),
  label, as.character(utils::packageVersion("occJSDM", lib.loc = lib)),
  n_params, niter, nchain, elapsed[["elapsed"]], elapsed[["user.self"]],
  8 * as.double(n_params) * niter * nchain / 2^20,
  heap_peak - heap_base, rss_base))
cat("PRINTED SUMMARY\n", msgs, sep = "")
cat("WARNINGS\n", if (length(warns)) paste0(warns, "\n") else "(none)\n", sep = "")
