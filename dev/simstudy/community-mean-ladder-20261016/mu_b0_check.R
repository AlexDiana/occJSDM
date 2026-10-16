# Decide whether the long-chain step runs (README.md, Mixing). mu_b0 counts as converged in an arm if its
# rank-normalised R-hat is under 1.05 in at least 9 of the 10 refitted replicates and its median bulk ESS is at
# least 100 of the 2,000 draws. Writes draws/mu_b0-check.csv and, if any arm fails, draws/long-arms.txt naming the
# failing arms, comma-separated, then exits with status 1 so that launch.sh runs the long step.
# Usage, from the run directory: Rscript mu_b0_check.R .
run <- commandArgs(trailingOnly = TRUE)[1]
files <- list.files(file.path(run, "draws", "short"), "[.]rds$", full.names = TRUE)
stopifnot(length(files) == 60)
x <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f); a <- r$draws$mu_b0   # niter x nchain
  data.frame(arm = r$arm, replicate = r$replicate, rhat = posterior::rhat(a), ess_bulk = posterior::ess_bulk(a))
}))
s <- do.call(rbind, lapply(split(x, x$arm), function(y) data.frame(
  arm = y$arm[1], replicates = nrow(y), rhat_under_1.05 = sum(y$rhat < 1.05), rhat_max = max(y$rhat),
  ess_bulk_median = median(y$ess_bulk))))
s$converged <- s$rhat_under_1.05 >= 9 & s$ess_bulk_median >= 100
write.csv(s, file.path(run, "draws", "mu_b0-check.csv"), row.names = FALSE)
print(s, digits = 3, row.names = FALSE)
arms <- c("M5_cm", "M10_cm", "M20_cm", "base_rare_cm_m5", "base_rare_cm_m10", "base_rare_cm_m20")
stopifnot(setequal(s$arm, arms))
failing <- arms[arms %in% s$arm[!s$converged]]
if (length(failing)) {
  writeLines(paste(failing, collapse = ","), file.path(run, "draws", "long-arms.txt"))
  message("mu_b0 not converged in: ", paste(failing, collapse = ", "))
  quit(status = 1)
}
message("mu_b0 converged in every arm; no long chains")
