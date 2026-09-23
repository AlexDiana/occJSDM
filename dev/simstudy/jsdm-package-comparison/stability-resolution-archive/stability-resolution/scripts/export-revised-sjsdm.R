# Export the revised sjSDM predictions and errors into the separately versioned
# "revised/results" directory, alongside the unchanged original results for the
# other three packages, and compare them with the original provisional sjSDM
# result. Runs only after select-sjsdm-revised.R has recorded the selection.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
suppressPackageStartupMessages(library(dplyr))
revised <- file.path(root, "revised")
selection <- readRDS(file.path(revised, "results/selection.rds"))
stopifnot(!file.exists(file.path(revised, "results/predictions.rds")))
input <- readRDS(file.path(root, "inputs/training.rds"))
new_x <- readRDS(file.path(root, "inputs/test-x.rds"))
truth <- readRDS(file.path(root, "truth/truth.rds"))
truth_new <- point_marginal(list(beta = truth$scaled_coefficients, loading = truth$loadings), new_x)
original <- readRDS(file.path(root, "results/predictions.rds"))

sjsdm_rows <- function(p, label, package = "sjSDM") {
  prediction <- point_marginal(p, new_x)
  mu <- cbind(1, as.matrix(new_x)) %*% p$beta
  ordinary <- logistic_normal(mu, rep(sqrt(colSums(p$loading^2)), each = nrow(mu)), 121L)
  marginal_difference <- max(abs(prediction - ordinary))
  a <- joint_integration(p, input$x, input$y, 61L, clamp = TRUE)
  b <- joint_integration(p, input$x, input$y, 121L, clamp = TRUE)
  conditional_difference <- max(abs(a$conditional - b$conditional))
  stopifnot(marginal_difference < 1e-4, conditional_difference < 1e-4)
  prediction <- prediction * .999999 + .0000005
  rows <- list()
  for (target in c("sampled_site_recovery", "new_site_prediction")) {
    x <- if (target == "sampled_site_recovery") input$x else new_x
    estimate <- if (target == "sampled_site_recovery") b$conditional else prediction
    true <- if (target == "sampled_site_recovery") truth$conditional_probability[1:100, ] else truth_new
    y <- if (target == "sampled_site_recovery") input$y else truth$occurrence[101:400, ]
    stopifnot(identical(dim(estimate), dim(true)), all(is.finite(estimate)),
              all(estimate > 0 & estimate < 1))
    rows[[target]] <- data.frame(package = package, fit = label, target = target,
      site = rep(rownames(x), 10), species = rep(colnames(input$y), each = nrow(x)),
      estimate = as.vector(estimate), truth = as.vector(true), observed = as.vector(y))
  }
  cells <- do.call(rbind, rows)
  cells$band <- cut(cells$truth, c(-Inf, .2, .8, Inf), right = FALSE,
                    labels = c("Below 20%", "20% to below 80%", "80% or above"))
  list(cells = cells, checks = data.frame(package = package,
    marginal_max_difference = marginal_difference, conditional_max_difference = conditional_difference))
}

label <- selection$selected[["sjSDM"]]
p <- readRDS(file.path(revised, "checks", paste0(label, "-parameters.rds")))
exported <- sjsdm_rows(p, label)
others <- original[original$package != "sjSDM", ]
cells <- rbind(others, exported$cells)
rownames(cells) <- NULL
stopifnot(nrow(cells) == nrow(original), identical(sort(unique(cells$package)), sort(unique(original$package))))
saveRDS(cells, file.path(revised, "results/predictions.rds"))
write.csv(cells, file.path(revised, "results/predictions.csv"), row.names = FALSE)
checks <- rbind(read.csv(file.path(root, "results/integration-checks.csv")) |> filter(package != "sjSDM"),
                exported$checks)
write.csv(checks, file.path(revised, "results/integration-checks.csv"), row.names = FALSE)

summarise_errors <- function(grouped) {
  grouped |> summarise(cells = n(),
    signed_error_pp = 100 * mean(estimate - truth),
    absolute_error_pp = 100 * mean(abs(estimate - truth)),
    brier_score = mean((estimate - observed)^2),
    negative_log_score = -mean(observed * log(estimate) + (1 - observed) * log1p(-estimate)),
    .groups = "drop")
}
overall <- cells |> group_by(package, target) |> summarise_errors()
by_band <- cells |> group_by(package, target, band, .drop = FALSE) |> summarise_errors()
by_species <- cells |> group_by(package, target, species) |> summarise_errors()
write.csv(overall, file.path(revised, "results/overall-errors.csv"), row.names = FALSE)
write.csv(by_band, file.path(revised, "results/errors-by-band.csv"), row.names = FALSE)
write.csv(by_species, file.path(revised, "results/errors-by-species.csv"), row.names = FALSE)
reloaded <- readRDS(file.path(revised, "results/predictions.rds"))
for (i in seq_len(nrow(overall))) {
  d <- reloaded[reloaded$package == overall$package[i] & reloaded$target == overall$target[i], ]
  stopifnot(abs(100 * sum(abs(d$estimate - d$truth)) / nrow(d) - overall$absolute_error_pp[i]) < 1e-10,
            abs(100 * sum(d$estimate - d$truth) / nrow(d) - overall$signed_error_pp[i]) < 1e-10)
}

# Comparison of three sjSDM results on identical cells: the original provisional
# unpenalised fit, the revised selected fit (basin A) and the best native fit in
# the other local maximum (basin B). Truth is read here only, after selection.
check <- readRDS(file.path(root, "stability-resolution/multistart/multistart-check.rds"))
s <- check$summary
b_start <- s$start[s$basin == "B"][which.max(s$native_penalised_score[s$basin == "B"])]
b_fit <- if (b_start <= 3L) readRDS(file.path(root, "stability/fits",
  sprintf("sjSDM-default-penalty-polish-start-%d.rds", b_start))) else
  readRDS(file.path(root, "stability-resolution/multistart", sprintf("sjSDM-multistart-start-%d.rds", b_start)))
basin_b <- sjsdm_rows(b_fit$parameters, sprintf("sjSDM-multistart-start-%d", b_start), "sjSDM basin B")$cells
variants <- rbind(
  original[original$package == "sjSDM", ] |> mutate(variant = "original unpenalised (provisional)"),
  exported$cells |> mutate(variant = "revised selected, basin A"),
  basin_b |> mutate(package = "sjSDM", variant = "best native fit, basin B"))
comparison <- variants |> group_by(variant, target) |> summarise_errors()
prediction_differences <- variants |>
  select(variant, target, site, species, estimate) |>
  tidyr::pivot_wider(names_from = variant, values_from = estimate) |>
  group_by(target) |>
  summarise(max_abs_difference_original_vs_revised_pp =
              100 * max(abs(`original unpenalised (provisional)` - `revised selected, basin A`)),
            mean_abs_difference_original_vs_revised_pp =
              100 * mean(abs(`original unpenalised (provisional)` - `revised selected, basin A`)),
            max_abs_difference_basin_A_vs_B_pp =
              100 * max(abs(`revised selected, basin A` - `best native fit, basin B`)),
            mean_abs_difference_basin_A_vs_B_pp =
              100 * mean(abs(`revised selected, basin A` - `best native fit, basin B`)), .groups = "drop")
species_comparison <- variants |> group_by(variant, target, species) |> summarise_errors()
write.csv(comparison, file.path(revised, "results/sjsdm-revision-comparison.csv"), row.names = FALSE)
write.csv(prediction_differences, file.path(revised, "results/sjsdm-revision-prediction-differences.csv"), row.names = FALSE)
write.csv(species_comparison, file.path(revised, "results/sjsdm-revision-species-comparison.csv"), row.names = FALSE)
print(overall, width = Inf); print(comparison, width = Inf); print(prediction_differences, width = Inf)
cat("Revised export written under", revised, "\n")
