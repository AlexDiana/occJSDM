# Applies README.md's mixing rule to every refit at the study's chain length, after the intervention step. Writes
# mixing-check.csv and draws/long-jobs.txt, the arm and variant labels to refit at 5,000 + 5,000 (empty if none).
# Run from the repository root or from this directory.
here <- if (file.exists("refit_spread.R")) "." else "dev/simstudy/b0-high-replication-20261017"
source(file.path(here, "common.R"))
d <- species_table(read_fits(file.path(here, "draws", c("diag", "interventions"))))
res <- do.call(rbind, lapply(split(d, d$job), function(x) data.frame(
  job = x$job[1], species_fits = nrow(x),
  B0_ess_median = median(x$B0_ess), sd_r_ess_median = median(x$sd_r_ess),
  B0_rhat_over_1.05 = mean(x$B0_rhat > 1.05), sd_r_rhat_over_1.05 = mean(x$sd_r_rhat > 1.05),
  aJ_ess_median = median(x$aJ_ess), m_r_ess_median = median(x$m_r_ess), stringsAsFactors = FALSE)))
res$fails <- res$B0_ess_median < 100 | res$sd_r_ess_median < 100 |
  res$B0_rhat_over_1.05 > 0.10 | res$sd_r_rhat_over_1.05 > 0.10
print(res, row.names = FALSE)
write.csv(res, file.path(here, "mixing-check.csv"), row.names = FALSE)
writeLines(paste(res$job[res$fails], collapse = ","), file.path(here, "draws", "long-jobs.txt"))
message("long chains for: ", if (any(res$fails)) paste(res$job[res$fails], collapse = ", ") else "none")
