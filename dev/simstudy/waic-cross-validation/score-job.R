# Score one saved fit: WAIC for a full-data fit, held-out lppd for a fold fit,
# plus the stricter-quadrature check and the convergence screen.
# Rscript OUTDIR/scripts/score-job.R OUTDIR JOB
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
outdir <- normalizePath(args[1], mustWork = TRUE)
job <- args[2]
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "cv-lib.R"))
stopifnot(identical(normalizePath(find.package("occJSDM")),
                    normalizePath(file.path(outdir, "library", "occJSDM"))))
md5 <- function(path) unname(tools::md5sum(path))
setup_path <- file.path(outdir, "setup.rds")
setup <- readRDS(setup_path)
fit_path <- file.path(outdir, "fits", paste0(job, ".rds"))
dir.create(file.path(outdir, "scores"), showWarnings = FALSE)
path <- file.path(outdir, "scores", paste0(job, ".rds"))
if (file.exists(path)) stop("Score already saved: ", path)
saved <- readRDS(fit_path)
fit <- saved$fit
spec <- cv_parse_job(job)
sites <- sort(unique(setup$data$info$Site))
train_sites <- setdiff(sites, saved$heldout_sites)
stopifnot(identical(saved$setup_md5, md5(setup_path)), identical(saved$job, job),
          identical(saved$seed, spec$seed), identical(saved$mcmc, setup$mcmc),
          fit$infos$n_factors == cv_models[[spec$model]],
          setequal(fit$infos$siteNames, train_sites),
          identical(sort(saved$heldout_sites),
                    sort(setup$folds$Site[setup$folds$fold %in% spec$fold])))
# The stored training observations are exactly the declared training sites.
stopifnot(identical(unname(as.matrix(cv_heldout_fit(fit, cv_subset(setup$data, train_sites))$infos$OTU)),
                    unname(as.matrix(fit$infos$OTU))))
dims <- dim(fit$results_output$theta0_output)
iterations <- dims[2]
chains <- dims[3]
stopifnot(iterations == setup$mcmc$niter, chains == setup$mcmc$nchain)
subset <- unlist(lapply(seq_len(chains), function(ch) {
  unique(as.integer(seq(1, iterations, length.out = 250))) + (ch - 1L) * iterations
}))
stricter_quadrature <- c(31L, 51L, 81L, 121L)
started <- proc.time()[[3]]
out <- list(job = job, fold = spec$fold, model = spec$model, seed = spec$seed,
            fit_md5 = md5(fit_path), fit_seconds = as.numeric(difftime(saved$finished, saved$started, units = "secs")),
            fit_warnings = saved$warnings, iterations = iterations, chains = chains)
if (is.na(spec$fold)) {
  waic_warnings <- character()
  score <- withCallingHandlers(occJSDM::computeSiteWAIC(fit), warning = function(w) {
    waic_warnings <<- c(waic_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  stricter <- suppressWarnings(occJSDM::computeSiteWAIC(fit, draws = subset,
    quadrature = stricter_quadrature, tolerance = 1e-6))
  out$max_quadrature_change <- max(abs(stricter$log_lik - score$log_lik[subset, , drop = FALSE]))
  # Delta-method Monte Carlo SE of WAIC, as in vignette-lesson/validate-waic.R.
  ll <- score$log_lik
  likelihood <- exp(sweep(ll, 2, apply(ll, 2, max), "-"))
  centred <- sweep(ll, 2, colMeans(ll), "-")
  influence <- -2 * rowSums(sweep(likelihood, 2, colMeans(likelihood), "/") - 1 -
    nrow(ll) / (nrow(ll) - 1) * sweep(centred^2, 2, colMeans(centred^2), "-"))
  out$waic_mcse <- posterior::mcse_mean(matrix(influence, iterations, chains))
  out$score <- score
  out$waic_warnings <- waic_warnings
  # How the refit compares with the archived fit made with the same seed
  # under R 4.5.0 and the older package (a56f554). Informational only.
  archived <- c(one = "/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-post-pr14-20261007/prediction/one-factor-fit.rds",
                two = "/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-post-pr14-20261007/default-fit.rds")[[spec$model]]
  if (file.exists(archived)) {
    old <- readRDS(archived)
    stopifnot(identical(old$seed, spec$seed))
    compare <- function(a, b) if (identical(dim(a), dim(b))) max(abs(a - b)) else NA_real_
    out$archive_comparison <- c(
      identical_observations = identical(unname(as.matrix(old$fit$infos$OTU)), unname(as.matrix(fit$infos$OTU))),
      B0 = compare(old$fit$results_output$jsdm_output$B0_output, fit$results_output$jsdm_output$B0_output),
      theta0 = compare(old$fit$results_output$theta0_output, fit$results_output$theta0_output),
      p = compare(old$fit$results_output$p_output, fit$results_output$p_output),
      sigmah = compare(old$fit$results_output$jsdm_output$sigmah_output, fit$results_output$jsdm_output$sigmah_output))
    rm(old)
  }
} else {
  heldout <- cv_subset(setup$data, saved$heldout_sites)
  score <- cv_score(fit, heldout)
  stricter <- cv_score(fit, heldout, draws = subset, quadrature = stricter_quadrature, tolerance = 1e-6)
  out$max_quadrature_change <- max(abs(stricter$log_lik - score$log_lik[subset, , drop = FALSE]))
  out$summary <- cv_fold_summary(score, iterations, chains)
  out$log_lik <- score$log_lik
  out$integration_orders <- table(score$integration$orders)
}
out$screen <- cv_screen(fit)
out$score_seconds <- proc.time()[[3]] - started
out$session <- sessionInfo()
partial <- paste0(path, ".partial")
saveRDS(out, partial)
stopifnot(file.rename(partial, path))
cat("Scored", job, "in", round(out$score_seconds, 1), "seconds; max quadrature change",
    signif(out$max_quadrature_change, 3), "; screen flags", sum(out$screen$flagged), "of",
    nrow(out$screen), "\n")
