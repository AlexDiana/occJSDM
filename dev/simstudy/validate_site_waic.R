# Run from a source checkout. This uses existing, unchanged teaching fits.
# Rscript dev/simstudy/validate_site_waic.R ARCHIVE_DIRECTORY OUTPUT_DIRECTORY [draws_per_chain]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 2:3)
archive <- normalizePath(args[1])
output_directory <- args[2]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
per_chain <- if (length(args) == 3) as.integer(args[3]) else 100L
source_md5 <- tools::md5sum(c("R/site-waic.R", "src/site_waic.cpp"))
pkgload::load_all(quiet = TRUE)
paths <- c(one_factor = file.path(archive, "prediction-lesson-20260922", "one-factor-fit.rds"),
           two_factors = file.path(archive, "vignette-lesson-20260919", "default-fit.rds"))
results <- list()
for (name in names(paths)) {
  saved <- readRDS(paths[[name]])
  fit <- saved$fit
  dims <- dim(fit$results_output$jsdm_output$B0_output)
  stopifnot(per_chain <= dims[2])
  within_chain <- unique(as.integer(seq(1, dims[2], length.out = per_chain)))
  draws <- unlist(lapply(seq_len(dims[3]), function(chain) within_chain + (chain - 1L) * dims[2]))
  cat(name, length(draws), "draws; input MD5", tools::md5sum(paths[[name]]), "\n")
  warnings <- character()
  started <- proc.time()[[3]]
  score <- withCallingHandlers(
    computeSiteWAIC(fit, threshold = 1, draws = draws),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  elapsed <- proc.time()[[3]] - started
  # A stricter, independently started quadrature sequence on the same draws
  # distinguishes quadrature error from posterior Monte Carlo error.
  stricter <- withCallingHandlers(
    computeSiteWAIC(fit, threshold = 1, draws = draws,
      quadrature = c(31, 51, 81, 121), tolerance = 1e-6),
    warning = function(w) invokeRestart("muffleWarning"))
  stopifnot(max(abs(score$log_lik - stricter$log_lik)) < 1e-4)
  results[[name]] <- list(score = score, warnings = warnings, elapsed = elapsed,
    stricter_waic = stricter$WAIC,
    max_stricter_log_lik_difference = max(abs(score$log_lik - stricter$log_lik)),
    source_fit_md5 = unname(tools::md5sum(paths[[name]])))
  if (requireNamespace("loo", quietly = TRUE)) {
    reference <- suppressWarnings(loo::waic(score$log_lik))
    stopifnot(abs(score$WAIC - reference$estimates["waic", "Estimate"]) < 1e-8,
              abs(score$SE - reference$estimates["waic", "SE"]) < 1e-8,
              abs(score$p_waic - reference$estimates["p_waic", "Estimate"]) < 1e-8)
    cat("Score, penalty and SE agree with loo::waic\n")
  }
  if (requireNamespace("posterior", quietly = TRUE)) {
    # Delta-method influence for WAIC as a function of posterior expectations.
    # Rows stay in iteration-within-chain order. This is approximate Monte
    # Carlo uncertainty, separate from the sampling-across-sites SE.
    log_lik <- score$log_lik
    likelihood <- exp(sweep(log_lik, 2, apply(log_lik, 2, max), "-"))
    centred <- sweep(log_lik, 2, colMeans(log_lik), "-")
    influence <- -2 * rowSums(
      sweep(likelihood, 2, colMeans(likelihood), "/") - 1 -
        nrow(log_lik) / (nrow(log_lik) - 1) *
        sweep(centred^2, 2, colMeans(centred^2), "-"))
    results[[name]]$waic_mcse <- posterior::mcse_mean(
      matrix(influence, length(within_chain), dims[3]))
    cat("Approximate delta-method WAIC MCSE:", results[[name]]$waic_mcse, "\n")
  }
  saveRDS(results[[name]], file.path(output_directory, paste0(name, ".rds")))
  print(c(WAIC = score$WAIC, SE = score$SE, p_waic = score$p_waic,
    flagged_sites = sum(score$pointwise$p_waic > 0.4), elapsed_seconds = elapsed,
    quadrature_difference = results[[name]]$max_stricter_log_lik_difference))
  print(table(score$integration$orders))
}
comparison <- compareSiteWAIC(results$one_factor$score, results$two_factors$score)
print(comparison[c("difference", "SE")])
saveRDS(list(results = results, comparison = comparison, session = sessionInfo(),
  source_md5 = source_md5),
  file.path(output_directory, "validation.rds"))
