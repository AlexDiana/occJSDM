#!/usr/bin/env Rscript
# Applies the grouped prior's convergence check, as it stands in the library given, to the pieces
# kept by convergence_check_on_new_fits.R, and lists beside it what says how good each chain is:
# its slope against the truth and its error (root mean square, as a share of the truth's), for the
# embedding's combined effect at the fitted sites. The R-hat of the scales is listed in two parts,
# sigma_emb with sigma_b and the axis scales, because whether the second should decide the verdict
# is the question these fits were run for. (The pieces are kept locally and not tracked, so this
# runs only where the fits were made.)
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/checks/convergence_check_on_new_fits_summary.R <library>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
d <- "dev/simstudy/grouped-mixing-20261118/results/convergence-check"
limits <- occJSDM:::GROUPED_LIMITS
rms <- function(x) sqrt(mean(x^2)); listed <- function(x, f = "%.2f") paste(sprintf(f, x), collapse = " ")
top <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
out <- do.call(rbind, lapply(list.files(d, pattern = "\\.rds$", full.names = TRUE), function(f) {
  x <- readRDS(f); chains <- seq_len(x$nchain)
  E <- vapply(chains, function(k) x$X_emb %*% x$B_emb_mean[, , k], matrix(0, nrow(x$X_emb), ncol(x$true_B_emb)))
  truth <- x$X_emb %*% x$true_B_emb; has_truth <- rms(truth) > 0
  e <- occJSDM:::grouped_effect_checks(E)
  apart <- 0; for (i in chains) for (j in chains) apart <- max(apart, rms(x$B_emb_mean[, , i] - x$B_emb_mean[, , j]))
  dominate <- sum(x$sigma_emb_median >= limits$state)
  rhat_two <- top(c(x$rhat_log_sigma_emb, x$rhat_log_sigma_b)); rhat_meas <- top(x$rhat_measured)
  axis_over <- sum(x$rhat_log_axis > limits$rhat, na.rm = TRUE)
  flag <- function(v) !is.na(v) && v
  others <- dominate > 0 || flag(rhat_meas > limits$rhat) || flag(rhat_two > limits$rhat) ||
    flag(e$agreement < limits$gradients) || flag(e$disagreement > limits$disagreement)
  data.frame(fit = sprintf("%s_s%d", x$truth, x$seed), minutes = round(x$seconds / 60),
    sigma_emb_median = listed(x$sigma_emb_median), leftover_chains = dominate,
    rhat_measured = rhat_meas, measured_over = mean(x$rhat_measured > limits$rhat, na.rm = TRUE),
    rhat_sigma_emb = top(x$rhat_log_sigma_emb), rhat_sigma_b = top(x$rhat_log_sigma_b),
    rhat_axis_max = top(x$rhat_log_axis), axis_over = axis_over,
    n_gradients = e$n_gradients, agreement = e$agreement, disagreement = e$disagreement,
    beyond = if (e$assessed) apart / rms(rowMeans(x$B_emb_mean, dims = 2)) else NA_real_,
    slope_by_chain = if (has_truth) listed(vapply(chains, function(k) sum(E[, , k] * truth) / sum(truth^2), 0)) else NA,
    error_by_chain = if (has_truth) listed(vapply(chains, function(k) rms(E[, , k] - truth) / rms(truth), 0)) else NA,
    error_measured = listed(vapply(chains, function(k) rms(x$B_meas_mean[, , k] - x$true_B_meas) / rms(x$true_B_meas), 0)),
    verdict_axis_in = if (others || axis_over > 0) "review" else "ok",
    verdict_axis_out = if (others) "review" else "ok",
    verdict_when_fitted = x$as_built$verdict, stringsAsFactors = FALSE)
}))
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], round, 3)
options(width = 320); print(out, row.names = FALSE)
utils::write.csv(out, file.path(dirname(d), "convergence_check_on_new_fits.csv"), row.names = FALSE)
