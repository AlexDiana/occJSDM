# Export existing fits for the lesson. This script never fits or selects a model.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(2L, 3L))
root <- normalizePath(args[1])
code <- normalizePath(args[2])
# Optional third argument: a CSV of additional fit attempts in the
# all-fit-attempts.csv layout, appended to the runtime table. The revised
# sjSDM export passes its twelve weak-penalty starts here.
extra_attempts <- if (length(args) == 3L) normalizePath(args[3]) else NULL
source(file.path(code, "pilot-math.R"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

training <- readRDS(file.path(root, "inputs/training.rds"))
test_x <- readRDS(file.path(root, "inputs/test-x.rds"))
truth <- readRDS(file.path(root, "truth/truth.rds"))
selection <- readRDS(file.path(root, "results/selection.rds"))$selected
predictions <- readRDS(file.path(root, "results/predictions.rds"))

# Both gradients use training-standard-deviation units, with the other at zero.
grid <- bind_rows(
  tibble(gradient = "environment_1", value = seq(-2, 2, length.out = 51),
         environment_1 = value, environment_2 = 0),
  tibble(gradient = "environment_2", value = seq(-2, 2, length.out = 51),
         environment_1 = 0, environment_2 = value)
)
grid_x <- grid |> select(environment_1, environment_2)

true_parameters <- list(beta = truth$scaled_coefficients, loading = truth$loadings)
true_curves <- point_marginal(true_parameters, grid_x)
curves <- list()
checks <- list()
point_parameters_for_lesson <- list()
for (package in names(selection)) {
  cat("Exporting response curves:", package, "\n")
  saved <- readRDS(file.path(root, "checks", paste0(selection[[package]], "-parameters.rds")))
  if (package %in% c("occJSDM", "Hmsc")) {
    estimate <- matrix(NA_real_, nrow(grid_x), ncol(training$y))
    difference <- 0
    for (rows in split(seq_len(nrow(grid_x)), ceiling(seq_len(nrow(grid_x)) / 17))) {
      coarse <- bayesian_marginal(saved$parameters, grid_x[rows, ], 31L)
      fine <- bayesian_marginal(saved$parameters, grid_x[rows, ], 61L)
      difference <- max(difference, max(abs(coarse - fine)))
      estimate[rows, ] <- apply(fine, 1:2, mean)
    }
  } else {
    point_parameters_for_lesson[[package]] <- saved
    estimate <- point_marginal(saved, grid_x)
    mu <- cbind(1, as.matrix(grid_x)) %*% saved$beta
    ordinary <- logistic_normal(mu, rep(sqrt(colSums(saved$loading^2)), each = nrow(mu)), 121L)
    difference <- max(abs(estimate - ordinary))
    if (package == "sjSDM") estimate <- .999999 * estimate + .0000005
  }
  stopifnot(difference < 1e-4, all(is.finite(estimate)))
  checks[[package]] <- tibble(package, maximum_difference = difference)
  curves[[package]] <- tibble(
    package, gradient = rep(grid$gradient, ncol(training$y)),
    value = rep(grid$value, ncol(training$y)),
    species = rep(colnames(training$y), each = nrow(grid)),
    estimate = as.vector(estimate), truth = as.vector(true_curves)
  )
}

attempts <- read.csv(file.path(code, "pilot-results/all-fit-attempts.csv"))
if (!is.null(extra_attempts)) {
  extra <- read.csv(extra_attempts)
  stopifnot(identical(names(extra), names(attempts)), !any(extra$fit %in% attempts$fit))
  attempts <- rbind(attempts, extra)
}
# When the results root carries a revised selection, keep the original
# selection and the revision's assessment tables in the bundle so the lesson
# can show them. Absent for the original pilot export.
selection_record <- readRDS(file.path(root, "results/selection.rds"))
revision <- if (!is.null(selection_record$original_selected)) {
  read_if <- function(name) {
    path <- file.path(root, "results", name)
    if (file.exists(path)) read.csv(path) else NULL
  }
  list(original_selected = selection_record$original_selected,
       sjSDM_configuration = selection_record$sjSDM_configuration,
       criteria = selection_record$criteria,
       basins = selection_record$basins,
       assessment = selection_record$multistart_assessment,
       comparison = read_if("sjsdm-revision-comparison.csv"),
       prediction_differences = read_if("sjsdm-revision-prediction-differences.csv"),
       species_comparison = read_if("sjsdm-revision-species-comparison.csv"),
       multistart = read.csv(file.path(code, "stability-resolution-results/multistart-summary.csv")),
       curvature = read.csv(file.path(code, "stability-resolution-results/curvature-summary.csv")),
       path_profile = read.csv(file.path(code, "stability-resolution-results/path-profile.csv")),
       species_differences = read.csv(file.path(code, "stability-resolution-results/species-differences.csv")))
} else NULL
diagnostics <- bind_rows(lapply(c("occJSDM", "Hmsc"), function(package) {
  read.csv(file.path(code, "pilot-results", paste0(package, "-diagnostics.csv"))) |>
    mutate(package = package)
}))
input_paths <- c(file.path(root, "inputs/training.rds"),
                 file.path(root, "inputs/test-x.rds"),
                 file.path(root, "truth/truth.rds"),
                 file.path(root, "results/predictions.rds"),
                 file.path(root, "results/selection.rds"),
                 file.path(root, "checks", paste0(unname(selection), "-parameters.rds")),
                 file.path(code, c("export-teaching.R", "pilot-math.R")))

lesson <- list(
  training = training, test_x = test_x, truth = truth,
  predictions = predictions, curves = bind_rows(curves),
  point_parameters = point_parameters_for_lesson,
  integration_checks = bind_rows(checks), attempts = attempts,
  diagnostics = diagnostics, selected = selection, sjsdm_revision = revision,
  provenance = list(seed = truth$seed, source_hashes = tools::md5sum(input_paths),
                    occJSDM_commit = "3a9726760f692ffe0ef21ca125d0c7acd02c7532",
                    sjSDM_commit = "d2ca508853a6e39df493f87c21d9c0136cfe652b",
                    versions = c(occJSDM = "0.1.0", gllvm = "2.0.15", sjSDM = "1.0.7", Hmsc = "3.3-7"))
)
output <- file.path(code, "../../../vignettes/teaching-data/jsdm-comparison.rds")
saveRDS(lesson, output, compress = "xz")
stopifnot(file.copy(
  file.path(root, "results/errors-by-species.csv"),
  file.path(dirname(output), "lesson-N-species-errors.csv"), overwrite = TRUE
))
cat("Saved compact teaching bundle:", normalizePath(output), "\n")
