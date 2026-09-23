# Check teaching exports against the original archive, without fitting models.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
code <- normalizePath(args[2])
repo <- normalizePath(file.path(code, "../../.."))
lesson <- readRDS(file.path(repo, "vignettes/teaching-data/jsdm-comparison.rds"))
source(file.path(code, "pilot-math.R"))
stopifnot(identical(lesson$predictions, readRDS(file.path(root, "results/predictions.rds"))),
          identical(lesson$training, readRDS(file.path(root, "inputs/training.rds"))),
          identical(lesson$truth, readRDS(file.path(root, "truth/truth.rds"))),
          identical(lesson$selected, readRDS(file.path(root, "results/selection.rds"))$selected),
          identical(unname(tools::md5sum(names(lesson$provenance$source_hashes))),
                    unname(lesson$provenance$source_hashes)))
stopifnot(identical(
  read.csv(file.path(repo, "vignettes/teaching-data/lesson-4-species-errors.csv")),
  read.csv(file.path(root, "results/errors-by-species.csv"))
))
keys <- lesson$predictions[c("package", "target", "site", "species")]
stopifnot(nrow(keys) == 16000, !anyDuplicated(keys),
          nrow(lesson$curves) == 4080,
          !anyDuplicated(lesson$curves[c("package", "gradient", "value", "species")]))

# Evaluate the displayed simulator, not a copy of its body in this check.
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
text <- readLines(file.path(repo, "vignettes/occJSDM-lesson-4.Rmd"))
chunk <- function(label) {
  start <- grep(paste0("^```\\{r ", label, "[,}]"), text)
  stopifnot(length(start) == 1)
  end <- which(seq_along(text) > start & text == "```")[1]
  text[(start + 1):(end - 1)]
}
env <- new.env(parent = globalenv())
env$known_truth <- lesson$truth
eval(parse(text = chunk("reproduce-community")), envir = env)
stopifnot(identical(env$presence_absence, lesson$truth$occurrence),
          identical(env$true_probability, lesson$truth$conditional_probability),
          identical(env$training_data$x, lesson$training$x),
          identical(env$training_data$y, lesson$training$y),
          identical(env$test_environment, lesson$test_x))
cat("Displayed simulation exactly reproduces the saved community and inputs.\n")

# Recompute all published error groups using base R, including observed scores.
for (file in c("overall-errors.csv", "errors-by-band.csv", "errors-by-species.csv")) {
  reference <- read.csv(file.path(root, "results", file))
  grouping <- intersect(c("package", "target", "band", "species"), names(reference))
  for (i in seq_len(nrow(reference))) {
    keep <- rep(TRUE, nrow(lesson$predictions))
    for (key in grouping) keep <- keep & as.character(lesson$predictions[[key]]) == reference[[key]][i]
    cells <- lesson$predictions[keep, ]
    stopifnot(nrow(cells) == reference$cells[i])
    values <- c(
      signed_error_pp = 100 * sum(cells$estimate - cells$truth) / nrow(cells),
      absolute_error_pp = 100 * sum(abs(cells$estimate - cells$truth)) / nrow(cells),
      brier_score = sum((cells$estimate - cells$observed)^2) / nrow(cells),
      negative_log_score = -sum(ifelse(cells$observed == 1, log(cells$estimate),
                                      log1p(-cells$estimate))) / nrow(cells)
    )
    stopifnot(max(abs(values - unlist(reference[i, names(values)]))) < 1e-10)
  }
}
cat("All overall, band and species error summaries reproduce.\n")

# Independent direct normal integral for every generating response-curve value.
truth <- lesson$truth
curve_truth <- unique(lesson$curves[c("gradient", "value", "species", "truth")])
for (i in seq_len(nrow(curve_truth))) {
  row <- curve_truth[i, ]
  x <- c(environment_1 = 0, environment_2 = 0)
  x[row$gradient] <- row$value
  raw_x <- x * lesson$training$spread + lesson$training$centre
  s <- match(row$species, truth$parameters$species)
  b <- as.numeric(truth$parameters[s, c("intercept", "environment_1", "environment_2")])
  mu <- sum(c(1, raw_x) * b)
  sd <- sqrt(sum(truth$loadings[, s]^2))
  independently_integrated <- integrate(function(z) plogis(mu + sd*z)*dnorm(z),
                                        -Inf, Inf, rel.tol = 1e-10)$value
  stopifnot(abs(independently_integrated - row$truth) < 1e-8)
}
cat("All 1,020 true response-curve values match raw generating parameters.\n")

# Independently evaluate point-fit curves, and a finer rule for Bayesian curves
# at both endpoints and the midpoint of each gradient (all ten species).
maximum_difference <- 0
for (package in names(lesson$selected)) {
  saved <- readRDS(file.path(root, "checks", paste0(lesson$selected[[package]], "-parameters.rds")))
  rows <- lesson$curves[lesson$curves$package == package, ]
  if (package %in% c("gllvm", "sjSDM")) {
    for (i in seq_len(nrow(rows))) {
      x <- c(environment_1 = 0, environment_2 = 0)
      x[rows$gradient[i]] <- rows$value[i]
      s <- match(rows$species[i], colnames(lesson$training$y))
      mu <- sum(c(1, x) * saved$beta[, s])
      sd <- sqrt(sum(saved$loading[, s]^2))
      value <- integrate(function(z) plogis(mu + sd*z)*dnorm(z), -Inf, Inf,
                         rel.tol = 1e-10)$value
      if (package == "sjSDM") value <- .999999*value + .0000005
      maximum_difference <- max(maximum_difference, abs(value - rows$estimate[i]))
    }
  } else {
    rows <- rows[rows$value %in% c(-2, 0, 2), ]
    p <- saved$parameters
    rule <- normal_rule(121L)
    for (i in seq_len(nrow(rows))) {
      x <- c(environment_1 = 0, environment_2 = 0)
      x[rows$gradient[i]] <- rows$value[i]
      s <- match(rows$species[i], colnames(lesson$training$y))
      beta <- matrix(p$beta[, s, , ], nrow = 3)
      mu <- as.vector(c(1, x) %*% beta)
      sd <- sqrt(colSums(matrix(p$loading[, s, , ], nrow = 2)^2))
      value <- if (package == "Hmsc") mean(pnorm(mu/sqrt(1+sd^2))) else {
        sum(vapply(seq_along(rule$nodes), function(k) {
          rule$weights[k] * mean(plogis(mu + sd*rule$nodes[k]))
        }, numeric(1)))
      }
      maximum_difference <- max(maximum_difference, abs(value - rows$estimate[i]))
    }
  }
}
stopifnot(maximum_difference < 1e-4,
          all(lesson$integration_checks$maximum_difference < 1e-4))
cat("Checked fitted curves; largest independent difference (pp):", 100*maximum_difference, "\n")

env$comparison <- lesson
env$training <- lesson$training
eval(parse(text = chunk("marginal-extraction-example")), envir = env)
expected <- subset(lesson$curves, package == "gllvm" & species == "species_01" &
                     gradient == "environment_1" & value == 0)$estimate
stopifnot(abs(env$fitted_probability - expected) < 1e-8)
cat("Displayed marginal-extraction example matches the fitted curve.\n")
cat("Teaching numerical verification passed.\n")
