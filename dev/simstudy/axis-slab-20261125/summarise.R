# Summarises the files compare.R wrote (evidence step 3 of
# dev/superpowers/specs/2026-11-25-axis-slab-design.md).
#
# Usage, from anywhere:
#   Rscript --vanilla summarise.R [--dir DIR] [--datasets 1:5]
# --dir defaults to dev/simstudy/results/axis-slab-20261125/. It reads
# slab2-d<k>.rds (with the cap) and slabInf-d<k>.rds (without) for each data
# set, says which are missing, and prints what the present files allow.
#
# Per data set: the two numbers that describe the truth (the largest absolute
# true coefficient on an embedding column, and the largest root mean square
# across species of a column's true coefficients), then per arm the RMSE of
# the posterior mean of the embedding's contribution X_emb %*% B_emb, the
# shrinkage slope, the verdict and flagged checks of
# returnGroupedConvergence(), and the minutes.
#
# The four readings, by the design's rules, fixed before the run. Each is
# read only when all the data sets asked for are present in both arms:
#   - Recovery: ratio of the RMSE with the cap to the RMSE without. "no harm"
#     if the ratio is at most 1.10 in at least 4 data sets and at most 1.25 in
#     every one; otherwise "harm".
#   - Shrinkage of the strongest columns: in each data set the 5 embedding
#     columns with the largest true root mean square across species; the
#     slope through the origin of the posterior means on the true values over
#     their species entries, sum(est * truth) / sum(truth^2). "no extra
#     shrinkage" if the slope with the cap is at least 0.95 of the slope
#     without in at least 4 data sets; otherwise "extra shrinkage".
#   - Convergence: a fit fails if returnGroupedConvergence() flags
#     rhat_measured, gradients or disagreement. "no worse" if the number of
#     failing data sets with the cap is at most the number without plus 1;
#     otherwise "worse".
#   - Time: mean minutes per fit. "slower" if the mean with the cap is more
#     than 10% above the mean without; otherwise "not slower".
# The rules assume 5 data sets; with another number they are applied as
# written (at least 4, every one) and say so.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing --", name)
    return(default)
  }
  args[i + 1L]
}
parse_datasets <- function(x) {
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    ends <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(ends[1], ends[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
study_dir <- dirname(normalizePath(script))
dir <- normalizePath(arg("dir", file.path(study_dir, "..", "results", "axis-slab-20261125")),
                     mustWork = TRUE)
datasets <- parse_datasets(arg("datasets", "1:5"))

ARMS <- c(cap = "slab2", none = "slabInf")
N_STRONG <- 5L
FAIL_CHECKS <- c("rhat_measured", "gradients", "disagreement")

files <- outer(ARMS, datasets, function(a, k) file.path(dir, sprintf("%s-d%d.rds", a, k)))
dimnames(files) <- list(names(ARMS), paste0("d", datasets))
present <- matrix(file.exists(files), nrow(files), dimnames = dimnames(files))
cat("Directory:", dir, "\n")
if (any(!present)) {
  cat("Missing:", paste(basename(files[!present]), collapse = ", "), "\n")
} else {
  cat("All", length(files), "files present.\n")
}

col_rms <- function(B) sqrt(rowMeans(B^2))
truth_numbers <- function(B) c(max_abs = max(abs(B)), max_col_rms = max(col_rms(B)))
strong_cols <- function(B) order(col_rms(B), decreasing = TRUE)[seq_len(N_STRONG)]
slope <- function(est, truth) sum(est * truth) / sum(truth^2)

rows <- list()
for (j in seq_along(datasets)) {
  k <- datasets[j]
  x <- lapply(names(ARMS), function(a) if (present[a, j]) readRDS(files[a, j]) else NULL)
  names(x) <- names(ARMS)
  have <- names(ARMS)[present[, j]]
  if (!length(have)) {
    cat(sprintf("\nd%d: no files\n", k))
    next
  }
  truth <- x[[have[1]]]$B_emb_true
  if (length(have) == 2L && !identical(x$cap$B_emb_true, x$none$B_emb_true))
    stop("the two arms of data set ", k, " have different truths")
  tn <- truth_numbers(truth)
  strong <- strong_cols(truth)
  cat(sprintf("\nd%d truth: largest |true coefficient| %.3f; largest column RMS %.3f; strongest %d columns %s (RMS %s)\n",
              k, tn["max_abs"], tn["max_col_rms"], N_STRONG,
              paste(rownames(truth)[strong], collapse = ", "),
              paste(sprintf("%.3f", col_rms(truth)[strong]), collapse = ", ")))
  row <- list(dataset = k, max_abs = tn[["max_abs"]], max_col_rms = tn[["max_col_rms"]])
  for (a in names(ARMS)) {
    r <- x[[a]]
    if (is.null(r)) {
      cat(sprintf("d%d %-7s missing\n", k, ARMS[[a]]))
      next
    }
    est <- r$B_emb_mean[strong, , drop = FALSE]
    s <- slope(est, truth[strong, , drop = FALSE])
    ch <- r$checks
    flagged <- ch$check[ch$flagged]
    fails <- any(ch$flagged[ch$check %in% FAIL_CHECKS])
    val <- function(nm) ch$value[ch$check == nm]
    cat(sprintf("d%d %-7s slab_scale %s: RMSE %.4f; slope %.4f; verdict %s, flagged: %s -> %s; rhat_measured %.3f, gradients %.3f, disagreement %.3f; %.1f minutes; power %s\n",
                k, ARMS[[a]], format(r$slab_scale), r$rmse_contribution, s, r$verdict,
                if (length(flagged)) paste(flagged, collapse = ", ") else "none",
                if (fails) "fails" else "passes",
                val("rhat_measured"), val("gradients"), val("disagreement"),
                r$minutes, r$power_start))
    row[[paste0("rmse_", a)]] <- r$rmse_contribution
    row[[paste0("slope_", a)]] <- s
    row[[paste0("fails_", a)]] <- fails
    row[[paste0("minutes_", a)]] <- r$minutes
  }
  if (length(have) == 2L) {
    cat(sprintf("d%d ratios cap / none: RMSE %.4f; slope %.4f\n", k,
                row$rmse_cap / row$rmse_none, row$slope_cap / row$slope_none))
  }
  rows[[length(rows) + 1L]] <- row
}

cat("\nReadings\n")
if (any(!present)) {
  cat("Not read: the readings need every data set in both arms (", sum(present), " of ",
      length(files), " files present).\n", sep = "")
  quit(save = "no", status = 0)
}
tab <- do.call(rbind, lapply(rows, as.data.frame))
nd <- nrow(tab)
if (nd != 5L) cat("Note: the rules are written for 5 data sets; here there are", nd, "\n")

r_rmse <- tab$rmse_cap / tab$rmse_none
recovery <- if (sum(r_rmse <= 1.10) >= 4L && all(r_rmse <= 1.25)) "no harm" else "harm"
cat(sprintf("Recovery: RMSE ratio cap / none %s; at most 1.10 in %d of %d, largest %.4f -> %s\n",
            paste(sprintf("%.4f", r_rmse), collapse = ", "), sum(r_rmse <= 1.10), nd,
            max(r_rmse), recovery))

r_slope <- tab$slope_cap / tab$slope_none
shrink <- if (sum(r_slope >= 0.95) >= 4L) "no extra shrinkage" else "extra shrinkage"
cat(sprintf("Shrinkage: slopes with the cap %s; without %s; ratio %s; at least 0.95 in %d of %d -> %s\n",
            paste(sprintf("%.4f", tab$slope_cap), collapse = ", "),
            paste(sprintf("%.4f", tab$slope_none), collapse = ", "),
            paste(sprintf("%.4f", r_slope), collapse = ", "), sum(r_slope >= 0.95), nd, shrink))

n_fail_cap <- sum(tab$fails_cap)
n_fail_none <- sum(tab$fails_none)
convergence <- if (n_fail_cap <= n_fail_none + 1L) "no worse" else "worse"
which_sets <- function(f) if (any(f)) paste0("d", tab$dataset[f], collapse = ",") else "none"
cat(sprintf("Convergence: failing fits with the cap %d (%s), without %d (%s) -> %s\n",
            n_fail_cap, which_sets(tab$fails_cap), n_fail_none, which_sets(tab$fails_none),
            convergence))

m_cap <- mean(tab$minutes_cap)
m_none <- mean(tab$minutes_none)
time_reading <- if (m_cap > 1.10 * m_none) "slower" else "not slower"
cat(sprintf("Time: mean minutes with the cap %.1f, without %.1f, ratio %.3f -> %s\n",
            m_cap, m_none, m_cap / m_none, time_reading))
