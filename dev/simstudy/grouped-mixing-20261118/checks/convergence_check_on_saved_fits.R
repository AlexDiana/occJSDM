#!/usr/bin/env Rscript
# Runs the two effect checks of the grouped prior's convergence check on the fits of the rotation
# and orientation tests, whose state is already known: which have a chain apart from the others,
# and which have chains where the leftovers carry the signal. The fits themselves were not kept;
# each run's _summary.rds holds every chain's posterior mean of the embedding's combined effect,
# which is all these two checks use. (The summaries are kept locally and not tracked, so this runs
# only where the fits were made.) The limits of the two checks were looked at against these same
# fits. Each chain's slope against the truth and its error (root mean square, as a share of the
# truth's) are listed beside the checks, because they are what says a chain is apart: the checks
# themselves never see the truth.
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/checks/convergence_check_on_saved_fits.R <library>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
d <- "dev/simstudy/grouped-mixing-20261118/results"
files <- c(list.files(file.path(d, "rotation"), pattern = "_summary.rds$", full.names = TRUE),
           list.files(file.path(d, "orientation"), pattern = "_summary.rds$", full.names = TRUE))
limits <- occJSDM:::GROUPED_LIMITS
rms <- function(x) sqrt(mean(x^2)); listed <- function(x) paste(sprintf("%.2f", x), collapse = " ")
out <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f); E <- r$contribution_mean; truth <- r$true_contribution
  x <- occJSDM:::grouped_effect_checks(E); chains <- seq_len(dim(E)[3])
  data.frame(fit = gsub("_b5000_summary.rds", "", basename(f), fixed = TRUE), size = x$size,
             n_gradients = x$n_gradients, agreement = x$agreement, disagreement = x$disagreement,
             flagged_gradients = x$agreement < limits$gradients,
             flagged_disagreement = x$disagreement > limits$disagreement,
             slope_by_chain = listed(vapply(chains, function(k) sum(E[, , k] * truth) / sum(truth^2), 0)),
             error_by_chain = listed(vapply(chains, function(k) rms(E[, , k] - truth) / rms(truth), 0)))
}))
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], round, 3)
options(width = 200); print(out, row.names = FALSE)
utils::write.csv(out, file.path(d, "convergence_check_on_saved_fits.csv"), row.names = FALSE)
