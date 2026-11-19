#!/usr/bin/env Rscript
# Sets the principal-component fits of orientation_test.R beside the fits of rotation_test.R for the
# same data sets. vs_rotated is the error of the embedding contribution over the mean of the rotated fits.
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/orientation_summary.R [--write]
d <- "dev/simstudy/grouped-mixing-20261118/results"
rd <- function(f) { x <- utils::read.csv(f); if (!"active" %in% names(x)) x$active <- 10L; if (!"orient" %in% names(x)) x$orient <- "asis"; x }
files <- c(list.files(file.path(d, "rotation"), pattern = "_b5000[.]csv$", full.names = TRUE), list.files(file.path(d, "orientation"), pattern = "_b5000[.]csv$", full.names = TRUE))
xs <- lapply(files, rd); cols <- Reduce(intersect, lapply(xs, names)); f <- do.call(rbind, lapply(xs, function(x) x[cols]))
f$arm <- ifelse(f$orient == "pca", "principal components", ifelse(f$rotate == 0, "as simulated", paste("rotation", f$rotate)))
f <- f[order(f$active, f$seed, f$arm), ]
rot <- stats::aggregate(err_contribution ~ active + seed, f[grepl("rotation", f$arm), ], mean); names(rot)[3] <- "err_rotated_mean"
f <- merge(f, rot, all.x = TRUE); f$vs_rotated <- f$err_contribution / f$err_rotated_mean
f <- f[order(f$active, f$seed, f$arm), ]; num <- vapply(f, is.numeric, logical(1)); f[num] <- lapply(f[num], signif, 3)
print(f[, c("active","seed","arm","err_contribution","vs_rotated","slope_contribution","err_measured","chain_rms_from_pooled","cells_over","cells_rhat_max","spread_over","meas_over","meas_rhat_max","axis_over","logsig_rhat","sigma_emb_median","chains_in_blend")], row.names = FALSE)
if ("--write" %in% commandArgs(TRUE)) utils::write.csv(f, file.path(d, "orientation", "summary_table.csv"), row.names = FALSE)
