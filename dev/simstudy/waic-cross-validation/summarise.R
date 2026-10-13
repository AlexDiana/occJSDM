# Combine the 22 scores, apply the pre-specified decision rule (DESIGN.md)
# and write compact results. Rscript OUTDIR/scripts/summarise.R OUTDIR RESULTS_DIR
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
outdir <- normalizePath(args[1], mustWork = TRUE)
results_dir <- args[2]
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setup <- readRDS(file.path(outdir, "setup.rds"))
scores <- lapply(setNames(setup$jobs$job, setup$jobs$job), function(job) {
  readRDS(file.path(outdir, "scores", paste0(job, ".rds")))
})
models <- c("one", "two")
sites <- sort(unique(setup$data$info$Site))
n <- length(sites)

# Held-out lppd per site and per chain, from the fold fits.
cv <- lapply(setNames(models, models), function(m) {
  rows <- do.call(rbind, lapply(scores[setup$jobs$job[!is.na(setup$jobs$fold) & setup$jobs$model == m]],
                                function(s) cbind(fold = s$fold, s$summary$sites)))
  rows[match(as.character(sites), rows$site), ]
})
stopifnot(!anyNA(cv$one$lppd), !anyNA(cv$two$lppd))
waic <- lapply(setNames(models, models), function(m) {
  p <- scores[[paste0("full-", m)]]$score$pointwise
  p[match(as.character(sites), p$site), ]
})
sites_table <- data.frame(
  site = sites, fold = cv$one$fold,
  cv_lppd_one = cv$one$lppd, cv_lppd_two = cv$two$lppd,
  waic_elpd_one = waic$one$elpd_waic, waic_elpd_two = waic$two$elpd_waic,
  waic_lppd_one = waic$one$lppd, waic_lppd_two = waic$two$lppd,
  p_waic_one = waic$one$p_waic, p_waic_two = waic$two$p_waic)
stopifnot(identical(cv$one$fold, cv$two$fold))

paired <- function(first, second) {
  d <- first - second
  c(first = sum(first), second = sum(second), difference = sum(d), se = sqrt(n * stats::var(d)))
}
verdict <- function(difference, se) {
  if (abs(difference) >= 4 && abs(difference) >= 2 * se) {
    if (difference > 0) "clear preference: one factor" else "clear preference: two factors"
  } else "no meaningful difference"
}
cv_total <- paired(sites_table$cv_lppd_one, sites_table$cv_lppd_two)
waic_total <- paired(sites_table$waic_elpd_one, sites_table$waic_elpd_two)
fold_mcse <- vapply(scores[setup$jobs$job[!is.na(setup$jobs$fold)]], function(s) s$summary$sum_mcse, numeric(1))
cv_mcse <- sqrt(sum(fold_mcse^2))
chain_columns <- grep("^chain[0-9]+$", names(cv$one), value = TRUE)
per_chain <- do.call(rbind, lapply(chain_columns, function(ch) {
  t <- paired(cv$one[[ch]], cv$two[[ch]])
  data.frame(chain = ch, elpd_one = t[["first"]], elpd_two = t[["second"]],
             difference = t[["difference"]], se = t[["se"]],
             verdict = verdict(t[["difference"]], t[["se"]]))
}))
cv_verdict <- verdict(cv_total[["difference"]], cv_total[["se"]])
waic_verdict <- verdict(waic_total[["difference"]], waic_total[["se"]])
shifted <- vapply(c(-2, 2), function(k) verdict(cv_total[["difference"]] + k * cv_mcse, cv_total[["se"]]), character(1))
robust_to_monte_carlo <- all(shifted == cv_verdict)
chains_agree <- all(per_chain$verdict == cv_verdict)
unclear <- !(robust_to_monte_carlo && chains_agree)
agree <- identical(cv_verdict, waic_verdict)

totals <- data.frame(
  method = c("ten-fold cross-validation", "WAIC (full-data refits)"),
  elpd_one = c(cv_total[["first"]], waic_total[["first"]]),
  elpd_two = c(cv_total[["second"]], waic_total[["second"]]),
  difference_one_minus_two = c(cv_total[["difference"]], waic_total[["difference"]]),
  paired_se = c(cv_total[["se"]], waic_total[["se"]]),
  monte_carlo_se = c(cv_mcse, sqrt(scores[["full-one"]]$waic_mcse^2 + scores[["full-two"]]$waic_mcse^2) / 2),
  verdict = c(cv_verdict, waic_verdict))
totals$deviance_difference <- -2 * totals$difference_one_minus_two

agreement <- data.frame(
  model = models,
  cor_site_elpd = c(stats::cor(sites_table$cv_lppd_one, sites_table$waic_elpd_one),
                    stats::cor(sites_table$cv_lppd_two, sites_table$waic_elpd_two)),
  mean_waic_minus_cv = c(mean(sites_table$waic_elpd_one - sites_table$cv_lppd_one),
                         mean(sites_table$waic_elpd_two - sites_table$cv_lppd_two)),
  flagged_sites_waic = c(sum(sites_table$p_waic_one > 0.4), sum(sites_table$p_waic_two > 0.4)))
paired_site_cor <- stats::cor(sites_table$cv_lppd_one - sites_table$cv_lppd_two,
                              sites_table$waic_elpd_one - sites_table$waic_elpd_two)

fits <- do.call(rbind, lapply(scores, function(s) {
  sc <- s$screen
  by_param <- tapply(sc$flagged, sc$param, sum)
  data.frame(job = s$job, model = s$model, fold = s$fold, seed = s$seed,
             fit_seconds = round(s$fit_seconds, 1), fit_warnings = length(s$fit_warnings),
             max_rhat = max(sc$rhat), min_ess_bulk = min(sc$ess_bulk), min_ess_tail = min(sc$ess_tail),
             flagged = sum(sc$flagged), flagged_excluding_psi = sum(sc$flagged[sc$param != "psi"]),
             flagged_params = paste(names(by_param)[by_param > 0], by_param[by_param > 0], sep = ":", collapse = " "),
             max_quadrature_change = s$max_quadrature_change,
             heldout_sum_mcse = if (is.null(s$summary)) NA_real_ else s$summary$sum_mcse,
             waic_mcse = if (is.null(s$waic_mcse)) NA_real_ else s$waic_mcse)
}))
archive <- do.call(rbind, lapply(models, function(m) {
  a <- scores[[paste0("full-", m)]]$archive_comparison
  if (is.null(a)) NULL else data.frame(model = m, t(a))
}))

decision <- c(
  paste("Cross-validation verdict:", cv_verdict),
  paste("WAIC verdict:", waic_verdict),
  paste("WAIC agrees with cross-validation:", agree),
  paste("Cross-validation verdict unchanged by +/- 2 Monte Carlo SE:", robust_to_monte_carlo,
        paste0("(", paste(shifted, collapse = "; "), ")")),
  paste("Every single-chain recomputation gives the same verdict:", chains_agree),
  paste("Step 1 unclear:", unclear),
  paste("Step 2 needed:", unclear || !agree),
  paste("Largest quadrature change over all fits:", signif(max(fits$max_quadrature_change), 3)))
writeLines(decision, file.path(results_dir, "decision.txt"))
write.csv(sites_table, file.path(results_dir, "sites.csv"), row.names = FALSE)
write.csv(totals, file.path(results_dir, "totals.csv"), row.names = FALSE)
write.csv(per_chain, file.path(results_dir, "per-chain.csv"), row.names = FALSE)
write.csv(agreement, file.path(results_dir, "agreement.csv"), row.names = FALSE)
write.csv(fits, file.path(results_dir, "fits.csv"), row.names = FALSE)
if (!is.null(archive)) write.csv(archive, file.path(results_dir, "archive-comparison.csv"), row.names = FALSE)
saveRDS(list(sites = sites_table, totals = totals, per_chain = per_chain, agreement = agreement,
             paired_site_cor = paired_site_cor, fits = fits, archive = archive, decision = decision,
             cv_mcse = cv_mcse, setup_md5 = unname(tools::md5sum(file.path(outdir, "setup.rds"))),
             package_commit = setup$package_commit, R = setup$R),
        file.path(outdir, "summary.rds"))
print(totals)
print(per_chain)
print(agreement)
cat("Correlation of per-site paired differences, CV vs WAIC:", paired_site_cor, "\n")
writeLines(decision)
