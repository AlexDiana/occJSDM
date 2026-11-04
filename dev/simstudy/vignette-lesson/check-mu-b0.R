# Convergence of the community mean of the occupancy intercepts (mu_b0) and
# of every row of returnConvergenceDiagnostics() in the teaching fits made at
# the community-mean prior (4 November 2026).
# Rscript dev/simstudy/vignette-lesson/check-mu-b0.R ARCHIVE OUT.csv
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
a <- normalizePath(args[1], mustWork = TRUE)
suppressPackageStartupMessages(library(occJSDM))
fits <- c(perfect = "perfect-fit.rds", default = "default-fit.rds",
          alternative = "alternative-fit.rds", alternative_long = "alternative-long-fit.rds",
          prediction_one_factor = "prediction/one-factor-fit.rds",
          unbalanced = "unbalanced/unbalanced-fit.rds")
rows <- list()
for (name in names(fits)) {
  path <- file.path(a, fits[[name]])
  if (!file.exists(path)) { message("missing: ", path); next }
  saved <- readRDS(path)
  fit <- if (!is.null(saved$fit)) saved$fit else saved
  d <- returnConvergenceDiagnostics(fit)
  m <- d[d$param == "mu_b0", ]
  stopifnot(nrow(m) <= 1L)
  # Fixed-prior fits (the alternative arm since 4 November 2026) have no mu_b0.
  if (!nrow(m)) m <- data.frame(mean = NA, q2.5 = NA, q97.5 = NA, rhat = NA, ess = NA)
  mc <- fit$infos$MCMCparams
  rows[[name]] <- data.frame(
    fit = name, file = fits[[name]], model = fit$infos$model,
    b0_prior = if (is.null(fit$results_output$jsdm_output$mu_b0_output)) "fixed" else "community",
    nchain = mc$nchain, nburn = mc$nburn, niter = mc$niter,
    mu_b0_mean = m$mean, mu_b0_q2.5 = m$q2.5, mu_b0_q97.5 = m$q97.5,
    mu_b0_rhat = m$rhat, mu_b0_ess = m$ess,
    max_rhat_all = max(d$rhat, na.rm = TRUE), min_ess_all = min(d$ess, na.rm = TRUE),
    n_rows = nrow(d), n_flagged = sum(d$rhat > 1.01 | d$ess < 400, na.rm = TRUE),
    shift_acceptance = mean(fit$infos$b0_shift$acceptance %||% NA_real_))
  rm(saved, fit); invisible(gc())
}
out <- do.call(rbind, rows)
write.csv(out, args[2], row.names = FALSE)
print(out, row.names = FALSE)
