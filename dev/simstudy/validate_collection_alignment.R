# Focused before/after diagnostic, not a nominal-coverage study.
# From the package root:
# Rscript dev/simstudy/validate_collection_alignment.R [replicates=5] [outdir] [before_ref]
# Both versions use the current compiled backend; this comparison is valid
# only when the backend is unchanged from before_ref.
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 5L
out_dir <- if (length(args) >= 2) args[2] else tempfile("collection-alignment-")
before_ref <- if (length(args) >= 3) args[3] else "6c837bc"
stopifnot(length(n_reps) == 1, !is.na(n_reps), n_reps > 0)
backend_diff <- system2("git", c("diff", "--name-only", shQuote(before_ref),
                                  "--", "src"), stdout = TRUE)
stopifnot(is.null(attr(backend_diff, "status")), length(backend_diff) == 0)
before_source <- system2("git", c("show", shQuote(paste0(before_ref, ":R/runOccJSDM.R"))),
                         stdout = TRUE)
stopifnot(is.null(attr(before_source, "status")))
after_head <- system2("git", "rev-parse HEAD", stdout = TRUE)
after_diff <- system2("git", "diff -- R/runOccJSDM.R", stdout = TRUE)
after_source <- readLines("R/runOccJSDM.R", warn = FALSE)

Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
devtools::load_all(quiet = TRUE)
RcppParallel::setThreadOptions(numThreads = 1)
source("tests/testthat/helper-simstudy.R")
before_env <- new.env(parent = asNamespace("occJSDM"))
eval(parse(text = before_source), envir = before_env)
fitters <- list(before = before_env$runOccJSDM, after = runOccJSDM)
scenarios <- Filter(function(s) s$label %in% c("base", "traits_isolated"),
                     simstudy_scenarios())
mcmc <- list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
slopes <- list()
fits <- list()
started <- Sys.time()

for (scenario in scenarios) {
  for (replicate in seq_len(n_reps)) {
    truth <- draw_truth(scenario, simstudy_seed_for(scenario, replicate))
    sim <- simstudy_simulate(truth)
    fit_seed <- .Random.seed
    info <- sim$data_list$info
    sample_rows <- info[!duplicated(info[c("Site", "Sample")]), ]
    sample_rows <- sample_rows[order(sample_rows$Site, sample_rows$Sample), ]
    stopifnot(scenario$ncov_theta == 1L)
    # The simulator's coefficients apply to raw X_theta; fitting standardises
    # it over samples. Compare slopes on the fitted scale in both versions.
    true_slopes <- sim$true_params$beta_theta_true[2, ] * sd(sample_rows$X_theta)
    true_psi <- plogis(sim$true_params$jsdmParams_true$eta)
    for (version in names(fitters)) {
      assign(".Random.seed", fit_seed, envir = globalenv())
      begin <- proc.time()["elapsed"]
      fit <- suppressMessages(suppressWarnings(fitters[[version]](
        sim$data_list,
        listParams = list(n_factors = scenario$d, n_lattrait = scenario$gt,
                           n_supportpoints = scenario$n_supportpoints),
        occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
        collCovariates = "X_theta",
        spatCovariates = if (scenario$useSpatField) c("Xs.1", "Xs.2") else NULL,
        MCMCparams = mcmc
      )))
      estimated <- apply(fit$results_output$beta_theta_output, c(1, 2), mean)[2, ]
      err <- estimated - true_slopes
      key <- paste(scenario$label, replicate, version, sep = ":")
      slopes[[key]] <- data.frame(
        scenario = scenario$label, replicate, version,
        species = seq_along(true_slopes),
        truth_sign = c("negative", "zero", "positive")[sign(true_slopes) + 2],
        truth = true_slopes, estimate = estimated, bias = err,
        abs_error = abs(err), squared_error = err^2
      )
      psi_error <- fit$results_output$psi_output - true_psi
      fits[[key]] <- data.frame(
        scenario = scenario$label, replicate, version,
        psi_bias = mean(psi_error), psi_mae = mean(abs(psi_error)),
        psi_rmse = sqrt(mean(psi_error^2)),
        B0_bias = mean(apply(fit$results_output$jsdm_output$B0_output, 1, mean) -
                         sim$true_params$jsdmParams_true$B0),
        elapsed_seconds = unname(proc.time()["elapsed"] - begin)
      )
      cat(sprintf("%s replicate %d/%d %s: slope MAE %.3f, %.1f s\n",
                   scenario$label, replicate, n_reps, version, mean(abs(err)),
                   fits[[key]]$elapsed_seconds))
    }
    # Save completed pairs so an interrupted diagnostic retains its evidence.
    saveRDS(list(slopes = do.call(rbind, slopes), fits = do.call(rbind, fits)),
             file.path(out_dir, "paired-results.rds"))
  }
}

slope_rows <- do.call(rbind, slopes)
fit_rows <- do.call(rbind, fits)
# Give each dataset equal weight within a slope-sign group. Its species
# count can differ, so pooling species would define a different estimator.
replicate_summary <- aggregate(cbind(truth, estimate, bias, abs_error, squared_error) ~
                                 scenario + version + truth_sign + replicate,
                                 slope_rows, mean)
summary <- aggregate(cbind(truth, estimate, bias, abs_error, squared_error) ~
                       scenario + version + truth_sign, replicate_summary, mean)
summary$rmse <- sqrt(summary$squared_error)
summary$squared_error <- NULL
# Monte Carlo uncertainty is across datasets, not species from the same fit.
bias_mcse <- aggregate(bias ~ scenario + version + truth_sign, replicate_summary,
                        function(x) sd(x) / sqrt(length(x)))
names(bias_mcse)[names(bias_mcse) == "bias"] <- "bias_mcse"
summary <- merge(summary, bias_mcse)
write.csv(slope_rows, file.path(out_dir, "slopes.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "slope-summary.csv"), row.names = FALSE)
write.csv(fit_rows, file.path(out_dir, "fit-metrics.csv"), row.names = FALSE)
saveRDS(list(before_ref = before_ref,
              after_head = after_head, after_diff = after_diff,
              before_source = before_source, after_source = after_source,
              mcmc = mcmc, replicates = n_reps, threads_per_fit = 1L,
              started = started, finished = Sys.time(), session = sessionInfo()),
         file.path(out_dir, "provenance.rds"))
print(summary, row.names = FALSE)
cat("Results:", normalizePath(out_dir), "\n")
