# The collapsed shift step at 10 species (README.md).
# Run from the repository root:
#   Rscript dev/simstudy/mu-b0-shift-s10-20261018/analyse.R
# Reads draws/s10 (refit-s100.R with the step on) and the step-off draws of the
# same replicates from the 100-species check's s10 step, whose directory is
# BEFORE_DIR (default: that check's worktree; its draws are not committed).
# Writes:
#   fits.csv          one row per fit and condition: mu_b0's and the mean B0's
#                     R-hat, bulk and tail ESS, posterior mean, sd and Monte
#                     Carlo SE; the community occupancy error; B0 coverage;
#                     the step's acceptance rates and proposal sd; elapsed time
#   before-after.csv  one row per fit: the step-on minus step-off posterior
#                     means of mu_b0 and the mean B0 in combined Monte Carlo
#                     SEs (R3), and the paired mixing statistics
#   rule.csv          R1 (strict and median), R2 and R3 by arm
here <- Sys.getenv("RUN_DIR", "dev/simstudy/mu-b0-shift-s10-20261018")
out <- Sys.getenv("OUT_DIR", here)
before_dir <- Sys.getenv("BEFORE_DIR",
  "/Users/douglasyu/src/occJSDM/.claude/worktrees/nostalgic-mendel-8b7449/dev/simstudy/b0-community-s100-20261018/draws/s10")
rhat <- posterior::rhat; ess_bulk <- posterior::ess_bulk; ess_tail <- posterior::ess_tail
mcse <- posterior::mcse_mean

load_dir <- function(d) {
  f <- sort(list.files(d, "[.]rds$", full.names = TRUE))
  setNames(lapply(f, readRDS), sub("[.]rds$", "", basename(f)))
}
after <- load_dir(file.path(here, "draws", "s10"))
before <- load_dir(before_dir)
stopifnot(length(after) > 0, all(names(after) %in% names(before)))

stats_of <- function(x, condition) {
  d <- x$draws
  mu <- matrix(d$mu_b0, ncol = dim(d$B0)[3])
  mb <- apply(d$B0, c(2, 3), mean)
  b0 <- x$rows[x$rows$block == "B0", ]
  occ <- x$occupancy
  sh <- x$b0_shift
  pick <- function(field, k) if (!is.null(sh[[field]]) && length(sh[[field]]) >= k) sh[[field]][k] else NA_real_
  data.frame(
    fit = sprintf("%s-%02d", x$arm, x$replicate), arm = x$arm, replicate = x$replicate, condition = condition,
    mu_rhat = rhat(mu), mu_ess_bulk = ess_bulk(mu), mu_ess_tail = ess_tail(mu),
    mu_mean = mean(mu), mu_sd = sd(as.vector(mu)), mu_mcse = mcse(mu),
    mb_rhat = rhat(mb), mb_ess_bulk = ess_bulk(mb), mb_mean = mean(mb), mb_sd = sd(as.vector(mb)), mb_mcse = mcse(mb),
    occ_error = mean(occ$occ_est) - mean(occ$occ_true), b0_coverage = mean(b0$covered),
    shift_acceptance_chain1 = pick("acceptance", 1), shift_acceptance_chain2 = pick("acceptance", 2),
    shift_sd_chain1 = pick("sd", 1), shift_sd_chain2 = pick("sd", 2),
    elapsed_minutes = x$elapsed / 60)
}
fits <- rbind(do.call(rbind, lapply(names(after), function(k) stats_of(before[[k]], "step off"))),
              do.call(rbind, lapply(names(after), function(k) stats_of(after[[k]], "step on"))))
rownames(fits) <- NULL
write.csv(fits, file.path(out, "fits.csv"), row.names = FALSE)

off <- fits[fits$condition == "step off", ]; on <- fits[fits$condition == "step on", ]
stopifnot(identical(off$fit, on$fit))
ba <- data.frame(fit = on$fit, arm = on$arm, replicate = on$replicate,
  mu_rhat_off = off$mu_rhat, mu_rhat_on = on$mu_rhat,
  mu_ess_bulk_off = off$mu_ess_bulk, mu_ess_bulk_on = on$mu_ess_bulk,
  mb_ess_bulk_off = off$mb_ess_bulk, mb_ess_bulk_on = on$mb_ess_bulk,
  mu_mean_off = off$mu_mean, mu_mean_on = on$mu_mean, mu_sd = on$mu_sd,
  mu_diff_in_mcse = (on$mu_mean - off$mu_mean) / sqrt(on$mu_mcse^2 + off$mu_mcse^2),
  mu_diff_in_sd = (on$mu_mean - off$mu_mean) / on$mu_sd,
  mb_diff_in_mcse = (on$mb_mean - off$mb_mean) / sqrt(on$mb_mcse^2 + off$mb_mcse^2),
  occ_error_off = off$occ_error, occ_error_on = on$occ_error,
  b0_coverage_off = off$b0_coverage, b0_coverage_on = on$b0_coverage,
  elapsed_off = off$elapsed_minutes, elapsed_on = on$elapsed_minutes)
write.csv(ba, file.path(out, "before-after.csv"), row.names = FALSE)

rule <- do.call(rbind, lapply(unique(fits$arm), function(a) {
  r1 <- function(x) c(rhat_all = all(x$mu_rhat < 1.05), ess_strict = all(x$mu_ess_bulk >= 100),
                      ess_median = median(x$mu_ess_bulk) >= 100)
  o <- off[off$arm == a, ]; n <- on[on$arm == a, ]; b <- ba[ba$arm == a, ]
  data.frame(arm = a,
    R1_strict_off = all(r1(o)[c("rhat_all", "ess_strict")]), R1_median_off = all(r1(o)[c("rhat_all", "ess_median")]),
    R1_strict_on = all(r1(n)[c("rhat_all", "ess_strict")]), R1_median_on = all(r1(n)[c("rhat_all", "ess_median")]),
    median_ess_off = median(o$mu_ess_bulk), median_ess_on = median(n$mu_ess_bulk),
    R2 = median(n$mu_ess_bulk) >= median(o$mu_ess_bulk),
    max_abs_mu_diff_in_mcse = max(abs(b$mu_diff_in_mcse)), max_abs_mb_diff_in_mcse = max(abs(b$mb_diff_in_mcse)),
    R3 = all(abs(b$mu_diff_in_mcse) < 3) && all(abs(b$mb_diff_in_mcse) < 3),
    mean_occ_error_off = mean(o$occ_error), mean_occ_error_on = mean(n$occ_error),
    mean_b0_coverage_off = mean(o$b0_coverage), mean_b0_coverage_on = mean(n$b0_coverage))
}))
write.csv(rule, file.path(out, "rule.csv"), row.names = FALSE)
print(ba[, c("fit", "mu_rhat_off", "mu_rhat_on", "mu_ess_bulk_off", "mu_ess_bulk_on", "mu_diff_in_mcse", "mb_diff_in_mcse")], digits = 3)
print(rule, digits = 3)
