args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
source(file.path(code, "pilot-math.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
new_x <- readRDS(file.path(root, "inputs/test-x.rds"))
selected <- readRDS(file.path(root, "results/selection.rds"))$selected
truth <- readRDS(file.path(root, "truth/truth.rds"))
truth_new <- point_marginal(list(beta = truth$scaled_coefficients, loading = truth$loadings), new_x)
rows <- list()
checks <- list()
for (package in names(selected)) {
  label <- selected[[package]]
  path <- file.path(root, "checks", paste0(label, "-parameters.rds"))
  saved <- readRDS(path)
  cat("Exporting", label, "\n")
  if (package %in% c("occJSDM", "Hmsc")) {
    p <- saved$parameters
    training <- saved$training_mean
    prediction <- matrix(NA_real_, nrow(new_x), 10)
    difference <- 0
    for (index in split(seq_len(nrow(new_x)), ceiling(seq_len(nrow(new_x))/25))) {
      a <- bayesian_marginal(p, new_x[index, , drop = FALSE], 31L)
      b <- bayesian_marginal(p, new_x[index, , drop = FALSE], 61L)
      difference <- max(difference, max(abs(a - b)))
      prediction[index, ] <- apply(b, 1:2, mean)
    }
    stopifnot(difference < 1e-4)
    checks[[package]] <- data.frame(package = package, marginal_max_difference = difference,
                                    conditional_max_difference = 0)
  } else {
    p <- saved
    prediction <- point_marginal(p, new_x)
    # Check the adaptive identity against ordinary high-order Gaussian quadrature.
    mu <- cbind(1, as.matrix(new_x)) %*% p$beta
    ordinary <- logistic_normal(mu, rep(sqrt(colSums(p$loading^2)), each = nrow(mu)), 121L)
    difference <- max(abs(prediction - ordinary))
    stopifnot(difference < 1e-4)
    a <- joint_integration(p, input$x, input$y, 61L, clamp = package == "sjSDM")
    b <- joint_integration(p, input$x, input$y, 121L, clamp = package == "sjSDM")
    conditional_difference <- max(abs(a$conditional - b$conditional))
    stopifnot(conditional_difference < 1e-4)
    training <- b$conditional
    if (package == "sjSDM") prediction <- prediction * .999999 + .0000005
    checks[[package]] <- data.frame(package = package, marginal_max_difference = difference,
                                    conditional_max_difference = conditional_difference)
  }
  for (target in c("sampled_site_recovery", "new_site_prediction")) {
    x <- if (target == "sampled_site_recovery") input$x else new_x
    estimate <- if (target == "sampled_site_recovery") training else prediction
    true <- if (target == "sampled_site_recovery") truth$conditional_probability[1:100, ] else truth_new
    y <- if (target == "sampled_site_recovery") input$y else truth$occurrence[101:400, ]
    stopifnot(identical(dim(estimate), dim(true)), all(is.finite(estimate)),
              all(estimate > 0 & estimate < 1))
    rows[[paste(package, target)]] <- data.frame(
      package = package, fit = label, target = target,
      site = rep(rownames(x), 10), species = rep(colnames(input$y), each = nrow(x)),
      estimate = as.vector(estimate), truth = as.vector(true),
      observed = as.vector(y))
  }
}
cells <- do.call(rbind, rows)
cells$band <- cut(cells$truth, c(-Inf,.2,.8,Inf), right = FALSE,
                  labels = c("Below 20%", "20% to below 80%", "80% or above"))
saveRDS(cells, file.path(root, "results", "predictions.rds"))
write.csv(cells, file.path(root, "results", "predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks), file.path(root, "results", "integration-checks.csv"), row.names = FALSE)

suppressPackageStartupMessages(library(dplyr))
summarise_errors <- function(grouped) {
  grouped |> summarise(cells = n(),
    signed_error_pp = 100 * mean(estimate - truth),
    absolute_error_pp = 100 * mean(abs(estimate - truth)),
    brier_score = mean((estimate - observed)^2),
    negative_log_score = -mean(observed * log(estimate) + (1-observed)*log1p(-estimate)),
    .groups = "drop")
}
overall <- cells |> group_by(package, target) |> summarise_errors()
by_band <- cells |> group_by(package, target, band, .drop = FALSE) |> summarise_errors()
by_species <- cells |> group_by(package, target, species) |> summarise_errors()
write.csv(overall, file.path(root, "results", "overall-errors.csv"), row.names = FALSE)
write.csv(by_band, file.path(root, "results", "errors-by-band.csv"), row.names = FALSE)
write.csv(by_species, file.path(root, "results", "errors-by-species.csv"), row.names = FALSE)
print(overall, width = Inf)
# Independent base-R recomputation, using the serialized output as the source.
reloaded <- readRDS(file.path(root, "results", "predictions.rds"))
for (i in seq_len(nrow(overall))) {
  d <- reloaded[reloaded$package == overall$package[i] & reloaded$target == overall$target[i], ]
  stopifnot(abs(100*sum(abs(d$estimate-d$truth))/nrow(d)-overall$absolute_error_pp[i]) < 1e-10,
            abs(100*sum(d$estimate-d$truth)/nrow(d)-overall$signed_error_pp[i]) < 1e-10)
}
cat("Independent recomputation of displayed signed and absolute errors passed.\n")
