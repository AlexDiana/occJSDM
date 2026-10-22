# Reuse the existing training community; create independent sites and one fit.
# Rscript prediction-build.R ORIGINAL_ARCHIVE NEW_ARCHIVE prepare|fit
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, args[3] %in% c("prepare", "fit"))
archive <- normalizePath(args[1], mustWork = TRUE)
outdir <- normalizePath(args[2], mustWork = TRUE)
action <- args[3]
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
original <- readRDS(file.path(archive, "input.rds"))
baseline_path <- file.path(archive, "default-fit.rds")
baseline <- readRDS(baseline_path)
md5 <- function(path) unname(tools::md5sum(path))
stopifnot(identical(original, lesson$input),
          identical(original$source_hashes, lesson_source_hashes()),
          identical(md5(baseline_path), lesson$manifests$default$md5),
          identical(normalizePath(find.package("occJSDM")),
                    normalizePath(file.path(archive, "library/occJSDM"))))
validate_lesson_fit_identity(baseline$fit, original)
input_path <- file.path(outdir, "prediction-input.rds")
fit_path <- file.path(outdir, "one-factor-fit.rds")

if (action == "prepare") {
  stopifnot(!file.exists(input_path), !file.exists(fit_path))
  truth <- original$sim$true_params$jsdmParams_true
  scaling <- baseline$fit$infos$list_X_psi_mat
  species <- baseline$fit$infos$speciesNames
  n_new <- 300L
  set.seed(20260925)
  raw_covariates <- matrix(rnorm(n_new * 2L, sd = 10), ncol = 2L,
                          dimnames = list(as.character(101:400), scaling$names_df))
  # Reuse training constants. Never estimate a new scale from test sites.
  standardized <- sweep(sweep(raw_covariates, 2, scaling$mean_df, "-"),
                        2, scaling$sd_df, "/")
  hidden <- matrix(rnorm(n_new * 2L, sd = original$jsdm$sigma_h), ncol = 2L)
  environmental_eta <- sweep(standardized %*% truth$B, 2, truth$B0, "+")
  conditional_probability <- plogis(environmental_eta + hidden %*% truth$L)
  actual_presence <- matrix(rbinom(length(conditional_probability), 1,
                                   as.vector(conditional_probability)), nrow = n_new)
  dimnames(conditional_probability) <- dimnames(actual_presence) <-
    list(rownames(raw_covariates), species)
  training_x <- baseline$fit$infos$data_info |>
    dplyr::distinct(Site, .keep_all = TRUE) |>
    dplyr::select(dplyr::all_of(scaling$names_df))
  outside <- rowSums(vapply(seq_len(ncol(raw_covariates)), function(k) {
    raw_covariates[, k] < min(training_x[[k]]) |
      raw_covariates[, k] > max(training_x[[k]])
  }, logical(n_new))) > 0
  bundle <- list(
    schema = 1L, new_site_seed = 20260925L, fitting_seed = 20260926L,
    public_prediction_seed = 20260927L,
    source_hashes = original$source_hashes,
    original_input_md5 = md5(file.path(archive, "input.rds")),
    original_lesson_md5 = md5("vignettes/teaching-data/nonspatial-lesson.rds"),
    baseline_fit_md5 = md5(baseline_path), scaling = scaling,
    raw_covariates = raw_covariates, standardized = standardized,
    hidden = hidden, environmental_eta = environmental_eta,
    conditional_probability = conditional_probability,
    actual_presence = actual_presence, outside_training_range = outside,
    selected_sites = as.character(101:110),
    design = paste("Compare one and two hidden site factors on identical original PCR training observations.",
                   "Evaluate species-level marginal probabilities at 300 independent new sites.",
                   "Use every saved posterior draw. Do not rank models using current augmented WAIC."),
    generator_md5 = md5("dev/simstudy/vignette-lesson/prediction-build.R")
  )
  saveRDS(bundle, input_path, compress = "xz")
  cat("Saved 300 independent sites before fitting or inspecting predictions.\n")
}

if (action == "fit") {
  stopifnot(file.exists(input_path), !file.exists(fit_path))
  input <- readRDS(input_path)
  stopifnot(identical(input$source_hashes, lesson_source_hashes()),
            identical(input$original_input_md5, md5(file.path(archive, "input.rds"))))
  fit_args <- list(
    data = original$sim$data_list,
    occCovariates = baseline$fit$infos$list_X_psi_mat$names_df,
    collCovariates = "X_theta", spatCovariates = NULL, threshold = 1,
    listParams = list(n_factors = 1L, n_lattrait = 1L),
    MCMCparams = baseline$mcmc, listPriors = baseline$priors,
    summarisedLatentPresences = TRUE
  )
  warnings <- character()
  started <- Sys.time()
  set.seed(input$fitting_seed)
  fit <- withCallingHandlers(do.call(occJSDM::runOccJSDM, fit_args),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
  result <- list(arm = "one_factor", fit = fit, seed = input$fitting_seed,
    mcmc = baseline$mcmc, priors = baseline$priors, warnings = warnings,
    started = started, finished = Sys.time(), session = sessionInfo(),
    source_hashes = lesson_source_hashes(),
    input_md5 = md5(file.path(archive, "input.rds")),
    prediction_input_md5 = md5(input_path))
  saveRDS(result, fit_path, compress = "gzip")
  validate_lesson_fit_identity(fit, original)
  stopifnot(identical(fit$infos$OTU, baseline$fit$infos$OTU),
            identical(fit$X_psi, baseline$fit$X_psi),
            identical(fit$X_theta, baseline$fit$X_theta),
            fit$infos$n_factors == 1L)
  cat("Saved one-factor fit in", as.numeric(difftime(result$finished, started, units = "secs")),
      "seconds; original training observations are identical.\n")
}
