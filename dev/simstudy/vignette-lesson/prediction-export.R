# Rscript prediction-export.R ORIGINAL_ARCHIVE NEW_ARCHIVE
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
archive <- normalizePath(args[1], mustWork = TRUE)
outdir <- normalizePath(args[2], mustWork = TRUE)
library(dplyr)
library(tidyr)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
source("dev/simstudy/vignette-lesson/prediction-math.R")
md5 <- function(p) unname(tools::md5sum(p))
original <- readRDS(file.path(archive, "input.rds"))
input_path <- file.path(outdir, "prediction-input.rds")
input <- readRDS(input_path)
fit_paths <- c(two_factors = file.path(archive, "default-fit.rds"),
               one_factor = file.path(outdir, "one-factor-fit.rds"))
stopifnot(identical(input$source_hashes, lesson_source_hashes()),
          identical(input$original_input_md5, md5(file.path(archive, "input.rds"))),
          identical(input$baseline_fit_md5, md5(fit_paths[["two_factors"]])),
          identical(input$generator_md5, md5("dev/simstudy/vignette-lesson/prediction-build.R")))
truth <- original$sim$true_params$jsdmParams_true
true_spread <- original$jsdm$sigma_h * sqrt(colSums(truth$L^2))
true_marginal <- normal_logistic_matrix(input$environmental_eta,
  rep(0, ncol(input$environmental_eta)), normal_rule(31L))
for (s in seq_along(true_spread)) true_marginal[, s] <- normal_logistic_mean(
  input$environmental_eta[, s], true_spread[s], normal_rule(61L))
truth_table <- expand.grid(Site = rownames(input$raw_covariates),
  species = colnames(input$actual_presence), stringsAsFactors = FALSE) |>
  mutate(truth = c(true_marginal), conditional_truth = c(input$conditional_probability),
         z = c(input$actual_presence),
         outside_training_range = rep(input$outside_training_range, length(true_spread)))
results <- public <- diagnostics <- manifests <- list()
for (label in names(fit_paths)) {
  saved <- readRDS(fit_paths[[label]])
  fit <- saved$fit
  validate_lesson_fit_identity(fit, original)
  stopifnot(identical(saved$source_hashes, original$source_hashes),
            identical(saved$input_md5, input$original_input_md5),
            identical(fit$infos$list_X_psi_mat, input$scaling))
  results[[label]] <- predict_marginal_teaching(fit, input$raw_covariates, label)
  selected <- input$raw_covariates[input$selected_sites, , drop = FALSE]
  set.seed(input$public_prediction_seed)
  native <- occJSDM::predictNewSites(fit, X_psi = as.data.frame(selected),
                                   useSpatial = FALSE, confidence = .95, verbose = FALSE)
  stopifnot(identical(dim(native), c(3L, nrow(selected), fit$infos$S)))
  public[[label]] <- expand.grid(Site = rownames(selected),
    species = fit$infos$speciesNames, stringsAsFactors = FALSE) |>
    mutate(arm = label, lower = c(native[1, , ]), median = c(native[2, , ]),
           upper = c(native[3, , ])) |>
    left_join(truth_table, by = c("Site", "species"), relationship = "one-to-one")
  diagnostics[[label]] <- as.data.frame(occJSDM::returnConvergenceDiagnostics(fit)) |>
    mutate(arm = label, .before = 1)
  manifests[[label]] <- list(file = basename(fit_paths[[label]]),
    md5 = md5(fit_paths[[label]]), seed = saved$seed, mcmc = saved$mcmc,
    priors = saved$priors, warnings = saved$warnings, session = saved$session,
    seconds = as.numeric(difftime(saved$finished, saved$started, units = "secs")),
    waic = occJSDM::extractWAIC(fit), factors = fit$infos$n_factors)
  # Keep the complete per-fit export outside the package for independent checks.
  saveRDS(results[[label]], file.path(outdir, paste0(label, "-predictions.rds")), compress = "xz")
}
cells <- bind_rows(lapply(results, `[[`, "cells")) |>
  left_join(truth_table, by = c("Site", "species"), relationship = "many-to-one") |>
  mutate(error = estimate - truth, absolute_error = abs(error),
    brier = (estimate - z)^2,
    negative_log_score = -if_else(z == 1, log(estimate), log1p(-estimate)))
stopifnot(all(is.finite(cells$negative_log_score)))
scores <- cells |>
  group_by(arm) |>
  summarise(cells = n(), signed_error_pp = 100 * mean(error),
    mean_absolute_error_pp = 100 * mean(absolute_error),
    rmse_pp = 100 * sqrt(mean(error^2)), brier = mean(brier),
    negative_log_score = mean(negative_log_score), .groups = "drop")
site_scores <- cells |>
  group_by(arm, Site) |>
  summarise(brier = mean(brier), negative_log_score = mean(negative_log_score), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = c(brier, negative_log_score))
differences <- bind_rows(lapply(c("brier", "negative_log_score"), function(metric) {
  d <- site_scores[[paste0(metric, "_one_factor")]] - site_scores[[paste0(metric, "_two_factors")]]
  data.frame(metric = metric, sites = length(d), difference = mean(d),
    site_SE = sd(d) / sqrt(length(d)), lower = mean(d) - 1.96 * sd(d) / sqrt(length(d)),
    upper = mean(d) + 1.96 * sd(d) / sqrt(length(d)))
}))
bundle <- list(schema = 1L, input = input, prediction_input_md5 = md5(input_path),
  source_hashes = original$source_hashes, manifests = manifests,
  truth = truth_table, cells = cells, scores = scores, paired_site_differences = differences,
  public = bind_rows(public), diagnostics = bind_rows(diagnostics),
  probability_diagnostics = bind_rows(lapply(results, `[[`, "diagnostics")),
  quadrature = lapply(results, `[[`, "quadrature"),
  source_md5 = tools::md5sum(c("dev/simstudy/vignette-lesson/prediction-build.R",
    "dev/simstudy/vignette-lesson/prediction-math.R",
    "dev/simstudy/vignette-lesson/prediction-export.R")))
saveRDS(bundle, "vignettes/teaching-data/prediction-lesson.rds", compress = "xz")
print(scores)
print(differences)
print(bundle$diagnostics |> filter(is.na(rhat) | is.na(ess) | rhat > 1.01 | ess < 400))
print(bundle$probability_diagnostics)
print(bundle$quadrature)
