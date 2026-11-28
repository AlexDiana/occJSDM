# Summarises the files lattrait.R wrote and prints the three readings whose
# rules are fixed in that script's header. Also checks that the 4-trait arm
# repeats the "base" files of repro.R, when those are present.
#
# Usage: Rscript --vanilla lattrait_summarise.R [--dir DIR] [--prefix small-]
# With --prefix small- it reads the files of lattrait.R --truth small. The
# printed readings 1 and 3 are then worded by the rules of the first test; the
# rules for the small truth are in lattrait.R's header and are applied by hand
# from the printed numbers.

LOW <- 0.02
LONG <- 500L

args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
i <- match("--dir", args)
out_dir <- if (is.na(i)) file.path(dirname(normalizePath(script)), "..", "results",
                                   "sigma-b-move-20261122") else args[i + 1L]
out_dir <- normalizePath(out_dir)

read_arm <- function(label) {
  f <- file.path(out_dir, sprintf("%s-d%d.rds", label, 1:5))
  lapply(f, function(x) if (file.exists(x)) readRDS(x) else NULL)
}
low_runs <- function(low) { r <- rle(low); r$lengths[r$values] }
j <- match("--prefix", args)
prefix <- if (is.na(j)) "" else args[j + 1L]
arms <- list(gt2 = read_arm(paste0(prefix, "gt2")), gt4 = read_arm(paste0(prefix, "gt4")))
base <- if (prefix == "") read_arm("base") else vector("list", 5)

chain_stats <- function(x) {
  sb <- x$sigma_b
  t(sapply(seq_len(ncol(sb)), function(ch) {
    runs <- low_runs(sb[, ch] < LOW)
    c(low_share = mean(sb[, ch] < LOW), longest = if (length(runs)) max(runs) else 0,
      long_runs = sum(runs >= LONG), median = stats::median(sb[, ch]))
  }))
}

for (label in names(arms)) for (x in arms[[label]]) {
  if (is.null(x)) next
  cs <- chain_stats(x)
  rh <- x$rhat
  verdict <- if (is.list(x$convergence)) x$convergence$verdict else "not available"
  cat(sprintf("%s d%d: low share by chain %s; longest low stay %s; long stays %d; median sigma_b %s\n",
              label, x$dataset, paste(sprintf("%.3f", cs[, "low_share"]), collapse = " "),
              paste(cs[, "longest"], collapse = " "), as.integer(sum(cs[, "long_runs"])),
              paste(sprintf("%.3f", cs[, "median"]), collapse = " ")))
  cat(sprintf("%s d%d: measured R-hat max %.3f (%d over 1.1) at 5,000 kept; %.3f (%d over 1.1) at all %d; verdict %s; RMSE embedding %.4f, measured %.4f; coverage %.3f, width %.3f; median sigma_emb by chain %s; %.0f minutes\n",
              label, x$dataset, rh["first_5000", "max"], as.integer(rh["first_5000", "over_1.1"]),
              rh["all", "max"], as.integer(rh["all", "over_1.1"]), nrow(x$sigma_b), verdict,
              x$rmse_contribution, x$rmse_measured, x$coverage_measured, x$width_measured,
              paste(sprintf("%.2f", apply(x$sigma_emb, 2, stats::median)), collapse = " "),
              x$minutes))
}

cat("\nDoes the 4-trait arm repeat the base files of repro.R?\n")
for (k in 1:5) {
  if (is.null(arms$gt4[[k]]) || is.null(base[[k]])) { cat("d", k, ": not both present\n", sep = ""); next }
  cat(sprintf("d%d: sigma_b identical %s\n", k, identical(arms$gt4[[k]]$sigma_b, base[[k]]$sigma_b)))
}

done <- function(a) !any(vapply(a, is.null, logical(1)))
if (!done(arms$gt2) || !done(arms$gt4)) {
  cat("\nNot all ten files are present, so no readings.\n")
  quit(save = "no")
}
num <- function(a, f) vapply(a, f, numeric(1))
long2 <- num(arms$gt2, function(x) sum(chain_stats(x)[, "long_runs"]))
long4 <- num(arms$gt4, function(x) sum(chain_stats(x)[, "long_runs"]))
longest2 <- num(arms$gt2, function(x) max(chain_stats(x)[, "longest"]))
longest4 <- num(arms$gt4, function(x) max(chain_stats(x)[, "longest"]))
fail2 <- num(arms$gt2, function(x) x$rhat["first_5000", "over_1.1"]) > 0
fail4 <- num(arms$gt4, function(x) x$rhat["first_5000", "over_1.1"]) > 0
r_emb <- num(arms$gt2, function(x) x$rmse_contribution) / num(arms$gt4, function(x) x$rmse_contribution)
r_meas <- num(arms$gt2, function(x) x$rmse_measured) / num(arms$gt4, function(x) x$rmse_measured)
read_ratio <- function(r) if (sum(r <= 1.10) >= 4) "no harm" else "harm"
f4 <- function(x) paste(sprintf("%.4f", x), collapse = ", ")

cat("\nReadings\n")
cat(sprintf("1. Low stays: long stays with 2 latent traits %d (longest stay by data set %s); with 4, %d (longest %s) -> %s%s\n",
            as.integer(sum(long2)), paste(longest2, collapse = ", "), as.integer(sum(long4)),
            paste(longest4, collapse = ", "),
            if (sum(long2) == 0) "long stays gone" else "long stays remain",
            if (sum(long4) == 0) " (the 4-trait arm has none either, so this does not bear on the suggestion)" else ""))
cat(sprintf("2. Measured coefficients: failing fits with 2 latent traits %d (%s), with 4, %d (%s)\n",
            sum(fail2), paste0("d", which(fail2), collapse = ","), sum(fail4),
            paste0("d", which(fail4), collapse = ",")))
cat(sprintf("3a. Embedding contribution, RMSE ratio 2 over 4: %s -> %s\n", f4(r_emb), read_ratio(r_emb)))
cat(sprintf("3b. Measured coefficients, RMSE ratio 2 over 4: %s -> %s\n", f4(r_meas), read_ratio(r_meas)))
cat(sprintf("    Coverage of the measured coefficients' 95%% intervals: with 2, %s; with 4, %s\n",
            f4(num(arms$gt2, function(x) x$coverage_measured)),
            f4(num(arms$gt4, function(x) x$coverage_measured))))
