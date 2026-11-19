#!/usr/bin/env Rscript
# Tables of the rotation test from what rotation_test.R wrote under results/rotation: each rotated
# fit against the unrotated fit of the same data set (same truth, same seed).
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/rotation_summary.R [--write]
d <- "dev/simstudy/grouped-mixing-20261118/results/rotation"
read1 <- function(f) { x <- utils::read.csv(f); if (!"active" %in% names(x)) x$active <- 10L; x }
all <- list.files(d, pattern = "_b[0-9]+(_chains)?[.]csv$", full.names = TRUE); chf <- grepl("_chains", all)
bind <- function(files) { xs <- lapply(files, read1); cols <- Reduce(intersect, lapply(xs, names)); do.call(rbind, lapply(xs, function(x) x[cols])) }
f <- bind(all[!chf]); ch <- bind(all[chf])
rng <- stats::aggregate(slope_contribution ~ active + seed + rotate, ch, function(v) sprintf("%.2f to %.2f", min(v), max(v))); names(rng)[4] <- "slope_by_chain"
f <- merge(f, rng)
base <- f[f$rotate == 0, c("active", "seed", "err_contribution", "slope_contribution")]; names(base)[3:4] <- c("err0", "slope0")
f <- merge(f, base, all.x = TRUE); f$err_ratio <- f$err_contribution / f$err0; f$slope_drop <- f$slope0 - f$slope_contribution
f <- f[order(f$active, f$seed, f$rotate), ]
num <- vapply(f, is.numeric, logical(1)); f[num] <- lapply(f[num], signif, 3)
print(f[, c("active","seed","rotate","err_contribution","err_ratio","slope_contribution","slope_drop","slope_by_chain","err_measured","chains_in_blend","sigma_emb_median")], row.names = FALSE)
print(f[, c("active","seed","rotate","chain_rms_from_pooled","cells_over","cells_rhat_max","spread_over","spread_rhat_max","meas_over","meas_rhat_max","axis_over","axis_rhat_max","logsig_rhat")], row.names = FALSE)
rms <- function(x) sqrt(mean(x^2))
for (i in which(f$rotate > 0)) {
  tag0 <- f$tag[f$active == f$active[i] & f$seed == f$seed[i] & f$rotate == 0]
  if (length(tag0) != 1) next
  a <- readRDS(file.path(d, paste0(tag0, "_summary.rds"))); b <- readRDS(file.path(d, paste0(f$tag[i], "_summary.rds")))
  pa <- apply(a$contribution_mean, c(1, 2), mean); pb <- apply(b$contribution_mean, c(1, 2), mean)
  within <- function(x, p) mean(sapply(seq_len(dim(x)[3]), function(k) rms(x[, , k] - p)))
  cat(sprintf("%s against the unrotated fit: rms difference of the posterior mean contributions %.3f; between chains %.3f (unrotated), %.3f (rotated)\n",
              f$tag[i], rms(pa - pb), within(a$contribution_mean, pa), within(b$contribution_mean, pb)))
}
if ("--write" %in% commandArgs(TRUE)) utils::write.csv(f, file.path(d, "summary_table.csv"), row.names = FALSE)
