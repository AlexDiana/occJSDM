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
# AR(1) series per parameter and chain (coefficient 0.5), built in place, with
# a garbage collection every 100 iterations of the loop that builds them, so
# that when the timing starts the process holds one copy of the array and
# little else.
#
# Memory figures, all in MB:
#   input_mb           the size of the block itself (8 bytes per draw)
#   heap_peak_mb       the most R heap in use at any point of the call, less
#                      what was in use just before it (gc() "max used" after a
#                      reset): what computeDiagnostics() allocated on top of
#                      its input
#   footprint_base_mb  memory footprint of the process just before the call
#                      (footprint(1), what Activity Monitor shows): R, the
#                      loaded packages and the block
# run-bench.sh adds the peak footprint of the process from /usr/bin/time -l,
# which also counts memory that R's heap accounting cannot see (C++ work
# buffers). Peak less base is what the call added to the process.

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
  if (t %% 100L == 0L) gc(FALSE)
}
results_output <- list(jsdm_output = list(B0_output = x))
rm(x)

heap_mb <- function(g, column) sum(g[, which(colnames(g) == column) + 1L])
footprint_mb <- function() {
  out <- suppressWarnings(system2("footprint", as.character(Sys.getpid()),
                                  stdout = TRUE, stderr = FALSE))
  line <- grep("phys_footprint:", out, value = TRUE)[1L]
  if (is.na(line)) return(NA_real_)
  parts <- strsplit(trimws(sub(".*phys_footprint:", "", line)), " +")[[1L]]
  as.numeric(parts[[1L]]) * c(KB = 1 / 1024, MB = 1, GB = 1024)[[parts[[2L]]]]
}

invisible(gc(reset = TRUE))
heap_base <- heap_mb(gc(), "used")
footprint_base <- footprint_mb()

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
  "elapsed_s=%.2f user_s=%.2f input_mb=%.0f heap_peak_mb=%.0f footprint_base_mb=%.0f\n"),
  label, as.character(utils::packageVersion("occJSDM", lib.loc = lib)),
  n_params, niter, nchain, elapsed[["elapsed"]], elapsed[["user.self"]],
  8 * as.double(n_params) * niter * nchain / 2^20,
  heap_peak - heap_base, footprint_base))
cat("PRINTED SUMMARY\n", msgs, sep = "")
cat("WARNINGS\n", if (length(warns)) paste0(warns, "\n") else "(none)\n", sep = "")
