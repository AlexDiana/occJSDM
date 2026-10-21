# Rscript summarise.R OUTDIR  -> prints the bullet lists for REPORT.md
args <- commandArgs(TRUE)
outdir <- normalizePath(args[[1]], mustWork = TRUE)
r <- read.csv(file.path(outdir, "results.csv"), stringsAsFactors = FALSE)
samples <- function(tag, role) {
  x <- read.csv(file.path(outdir, "rss", paste0(tag, ".csv"))); x <- x[x$role == role, ]
  if (role == "worker") {                          # keep only this fit's workers, from its execution record
    pids <- readRDS(file.path(outdir, "execution", paste0(tag, ".rds")))$worker_pids
    x <- x[x$pid %in% pids, ]
  }
  x
}
worker_pids <- function(tag) readRDS(file.path(outdir, "execution", paste0(tag, ".rds")))$worker_pids
peak_total <- function(tag, role) {
  x <- samples(tag, role)
  if (role == "parent" && !nrow(x)) stop("no parent memory samples for ", tag)
  if (role == "worker") {
    pids <- worker_pids(tag)
    if (is.null(pids)) return(NA_real_)             # in-process fit: no worker existed, reported as such below
    missing <- setdiff(pids, unique(x$pid))
    if (length(missing)) stop("no memory samples for worker(s) ", paste(missing, collapse = ", "), " in ", tag)
  }
  max(tapply(x$rss_kb, x$time, sum)) / 1024       # every sampled process at that second, summed
}
peak_each <- function(tag) {
  if (is.null(worker_pids(tag))) return("none")
  x <- samples(tag, "worker")
  paste(round(tapply(x$rss_kb, x$pid, max) / 1024), collapse = "/")
}
r$parent_peak_mb <- vapply(r$tag, peak_total, numeric(1), role = "parent")
r$workers_peak_mb <- vapply(r$tag, peak_total, numeric(1), role = "worker")
r$worker_peaks_mb <- vapply(r$tag, peak_each, character(1))
base <- mean(r$wall_s[r$label == "serial-ambient"])
for (lab in unique(r$label)) {
  s <- r[r$label == lab, ]
  workers_text <- if (all(is.na(s$workers_peak_mb))) {
    "no worker process (in-process fit)"
  } else {
    sprintf("workers together %.0f MB at most, each worker %s MB in repeat 1 and %s MB in repeat 2",
            max(s$workers_peak_mb), s$worker_peaks_mb[1], s$worker_peaks_mb[2])
  }
  cat(sprintf("- **%s** (ncores %d, %d worker%s, %d thread%s each): wall %.0f s and %.0f s (speed-up %.2fx over serial-ambient), parent peak %.0f MB, %s, identical fit %s, identical RNG %s.\n",
              lab, s$ncores[1], s$workers[1], if (s$workers[1] == 1) "" else "s", s$threads_per_worker[1],
              if (s$threads_per_worker[1] == 1) "" else "s", s$wall_s[1], s$wall_s[2], base / mean(s$wall_s),
              max(s$parent_peak_mb), workers_text,
              paste(s$identical_fit, collapse = "/"), paste(s$identical_rng, collapse = "/")))
}
cat("\nDeclaration basis:", readLines(file.path(outdir, "blas-declaration-basis.txt")), "\n")
