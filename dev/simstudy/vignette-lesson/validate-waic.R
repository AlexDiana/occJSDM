# Validate the current three-sample teaching fits. The September validator is unchanged.
# Rscript dev/simstudy/vignette-lesson/validate-waic.R ARCHIVE OUTPUT DRAWS_PER_CHAIN
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
archive <- normalizePath(args[1], mustWork = TRUE)
outdir <- args[2]
per_chain <- as.integer(args[3])
stopifnot(length(per_chain) == 1L, !is.na(per_chain), per_chain >= 2L)
if (dir.exists(outdir)) stop("Use a new output directory; do not overwrite validation.")
library(occJSDM)
source("dev/simstudy/vignette-lesson/helpers.R")
stopifnot(identical(normalizePath(find.package("occJSDM")),
                    normalizePath(file.path(archive, "library/occJSDM"))),
          requireNamespace("loo", quietly = TRUE),
          requireNamespace("posterior", quietly = TRUE))
paths <- c(one_factor = file.path(archive, "prediction", "one-factor-fit.rds"),
           two_factors = file.path(archive, "default-fit.rds"))
stopifnot(all(file.exists(paths)))
dir.create(outdir, recursive = TRUE)
results <- list()
for (name in names(paths)) {
  saved <- readRDS(paths[[name]])
  stopifnot(identical(saved$source_hashes, lesson_source_hashes()))
  fit <- saved$fit
  stopifnot(fit$infos$threshold == 1, fit$infos$n == 100L,
            length(unique(fit$infos$data_info$Sample)) == 300L)
  dims <- dim(fit$results_output$jsdm_output$B0_output)
  stopifnot(per_chain <= dims[2], dims[3] == 4L)
  within_chain <- unique(as.integer(seq(1, dims[2], length.out = per_chain)))
  draws <- unlist(lapply(seq_len(dims[3]), function(chain) within_chain + (chain - 1L) * dims[2]))
  cat(name, length(draws), "draws; input MD5", tools::md5sum(paths[[name]]), "\n")
  warnings <- character()
  started <- proc.time()[[3]]
  score <- withCallingHandlers(computeSiteWAIC(fit, draws = draws), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  stricter <- suppressWarnings(computeSiteWAIC(fit, draws = draws,
    quadrature = c(31, 51, 81, 121), tolerance = 1e-6))
  max_change <- max(abs(score$log_lik - stricter$log_lik))
  stopifnot(max_change < 1e-4)
  reference <- suppressWarnings(loo::waic(score$log_lik))
  stopifnot(abs(score$WAIC - reference$estimates["waic", "Estimate"]) < 1e-8,
            abs(score$SE - reference$estimates["waic", "SE"]) < 1e-8,
            abs(score$p_waic - reference$estimates["p_waic", "Estimate"]) < 1e-8)
  # Delta-method uncertainty in posterior expectations, preserving chain order.
  # This is distinct from the sampling-across-sites SE and quadrature error.
  log_lik <- score$log_lik
  likelihood <- exp(sweep(log_lik, 2, apply(log_lik, 2, max), "-"))
  centred <- sweep(log_lik, 2, colMeans(log_lik), "-")
  influence <- -2 * rowSums(sweep(likelihood, 2, colMeans(likelihood), "/") - 1 -
    nrow(log_lik) / (nrow(log_lik) - 1) *
    sweep(centred^2, 2, colMeans(centred^2), "-"))
  mcse <- posterior::mcse_mean(matrix(influence, length(within_chain), dims[3]))
  results[[name]] <- list(score = score, warnings = warnings,
    elapsed = proc.time()[[3]] - started, stricter_waic = stricter$WAIC,
    max_stricter_log_lik_difference = max_change, waic_mcse = mcse,
    source_fit_md5 = unname(tools::md5sum(paths[[name]])))
  saveRDS(results[[name]], file.path(outdir, paste0(name, ".rds")))
  print(c(WAIC = score$WAIC, site_SE = score$SE, p_waic = score$p_waic,
          flagged_sites = sum(score$pointwise$p_waic > 0.4), mcse = mcse,
          max_loglik_refinement_change = max_change))
  print(table(score$integration$orders))
}
comparison <- compareSiteWAIC(results$one_factor$score, results$two_factors$score)
summary <- do.call(rbind, lapply(names(results), function(name) {
  x <- results[[name]]
  data.frame(model = name, draws = nrow(x$score$log_lik), WAIC = x$score$WAIC,
    site_SE = x$score$SE, p_waic = x$score$p_waic,
    flagged_sites = sum(x$score$pointwise$p_waic > 0.4), waic_mcse = x$waic_mcse,
    stricter_WAIC = x$stricter_waic,
    max_loglik_refinement_change = x$max_stricter_log_lik_difference,
    fit_md5 = x$source_fit_md5)
}))
saveRDS(list(results = results, comparison = comparison, session = sessionInfo(),
  source_md5 = tools::md5sum(c("R/site-waic.R", "src/site_waic.cpp",
    "dev/simstudy/vignette-lesson/validate-waic.R"))), file.path(outdir, "validation.rds"))
write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)
write.csv(comparison$pointwise, file.path(outdir, "paired-sites.csv"), row.names = FALSE)
print(summary)
cat("One minus two WAIC:", comparison$difference, "paired site SE:", comparison$SE,
    "approximate difference MCSE assuming independent fits:", sqrt(sum(summary$waic_mcse^2)), "\n")
