# Summarises the files repro.R wrote for one label.
#
# Usage: Rscript --vanilla summarise.R --label LABEL
#
# A kept sweep is "low" when sigma_b is under 0.02. Per chain: the share of
# kept sweeps low, the longest run of low sweeps, and the number of runs of
# 500 or more. Per fit: the largest measured-coefficient R-hat and the number
# over 1.1, at the first 5,000 kept sweeps and at all of them, the minutes,
# and the acceptance of the sigma_b move's proposals when the fit has the move.
# Last line: the longest low run over all chains of the label.

LOW <- 0.02
LONG <- 500L

args <- commandArgs(trailingOnly = TRUE)
i <- match("--label", args)
if (is.na(i)) stop("missing --label")
label <- args[i + 1L]

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
out_dir <- normalizePath(file.path(dirname(normalizePath(script)), "..", "results",
                                   "sigma-b-move-20261122"))
files <- list.files(out_dir, pattern = sprintf("^%s-d[0-9]+\\.rds$", label), full.names = TRUE)
if (length(files) == 0L) stop("no files for label ", label, " in ", out_dir)
files <- files[order(as.integer(sub(".*-d([0-9]+)\\.rds$", "\\1", files)))]

low_runs <- function(low) {
  r <- rle(low)
  r$lengths[r$values]
}

longest <- 0L
for (f in files) {
  x <- readRDS(f)
  sb <- x$sigma_b
  for (ch in seq_len(ncol(sb))) {
    runs <- low_runs(sb[, ch] < LOW)
    run_max <- if (length(runs)) max(runs) else 0L
    longest <- max(longest, run_max)
    cat(sprintf("%s d%d chain %d: low share %.3f, longest low run %d, runs >= %d: %d, median sigma_b %.4f\n",
                label, x$dataset, ch, mean(sb[, ch] < LOW), run_max, LONG,
                sum(runs >= LONG), stats::median(sb[, ch])))
  }
  rh <- x$rhat
  am <- x$axis_move
  acc <- ""
  if (!is.null(am$sigma_b_acceptance)) {
    a <- am$sigma_b_acceptance
    acc <- sprintf(", sigma_b acceptance adapted %s, large %s, adapted step %s",
                   paste(sprintf("%.2f", a["adapted", ]), collapse = "/"),
                   paste(sprintf("%.2f", a["large", ]), collapse = "/"),
                   paste(sprintf("%.3f", am$sigma_b_step), collapse = "/"))
  }
  cat(sprintf("%s d%d fit: R-hat measured max %.3f (%d over 1.1) at 5,000 kept; %.3f (%d over 1.1) at %d kept; %.1f minutes; low share %.3f%s\n",
              label, x$dataset, rh["first_5000", "max"], as.integer(rh["first_5000", "over_1.1"]),
              rh["all", "max"], as.integer(rh["all", "over_1.1"]), nrow(sb), x$minutes,
              mean(sb < LOW), acc))
}
cat(sprintf("%s: longest low run over all chains %d\n", label, longest))
