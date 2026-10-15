# Extrapolate the full paper grid's wall time from this pre-study's results file.
# Usage: Rscript estimate.R   (from the run directory)
f <- list.files("results", pattern = "^simstudy-.*[.]rds$", full.names = TRUE)[1]
x <- readRDS(f)
fits <- nrow(unique(x$rows[, c("scenario", "replicate")]))
cat(sprintf("file %s\nelapsed %.2f min | fits %d | failed %d | fits/min %.2f | cores %d | threads/fit %s | sha %s | dirty %s\n",
            basename(f), x$elapsed_min, fits, x$failed_replicates, fits / x$elapsed_min,
            x$provenance$n_cores, x$provenance$rcpp_parallel_threads_in_worker,
            substr(x$provenance$git_sha, 1, 7), x$provenance$git_dirty))
per_fit <- x$elapsed_min / fits
cells <- length(unique(x$rows$scenario))
cat(sprintf("cells %d; wall per fit %.3f min on %d cores\n", cells, per_fit, x$provenance$n_cores))
for (R in c(100, 200, 300, 500))
  cat(sprintf("R = %3d: %5d fits, %5.0f min = %5.1f h\n", R, cells * R, cells * R * per_fit, cells * R * per_fit / 60))
