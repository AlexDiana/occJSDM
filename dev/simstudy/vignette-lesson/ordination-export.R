# Export exactly the displayed native ordination code from the complete fit.
# No fitting, thinning, replacement fit or package-code modification.
# From repository root with archive/library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/ordination-export.R FULL_FIT_DIRECTORY
# ORDINATION_STAGE optionally redirects these new deliverables for review.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("ORDINATION_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))

library(occJSDM)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
theme_set(ggtern::theme_bw(base_size = 12))
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")

lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
outputs_path <- "vignettes/teaching-data/output-lesson.rds"
snippet_path <- new_path("dev/simstudy/vignette-lesson/ordination-examples.Rmd")
exporter_path <- new_path("dev/simstudy/vignette-lesson/ordination-export.R")
lesson <- readRDS(lesson_path)
outputs <- readRDS(outputs_path)
input <- readRDS(file.path(archive, "input.rds"))
manifest <- lesson$manifests$default
fit_path <- file.path(archive, manifest$file)
stopifnot(
  identical(input, lesson$input),
  identical(input$source_hashes, lesson_source_hashes()),
  identical(outputs$source_hashes, input$source_hashes),
  identical(outputs$lesson_md5, md5(lesson_path)),
  identical(outputs$fit_manifests$default, manifest),
  identical(lesson$input_md5, md5(file.path(archive, "input.rds"))),
  identical(manifest$md5, md5(fit_path)),
  identical(normalizePath(find.package("occJSDM")),
            normalizePath(file.path(archive, "library/occJSDM")))
)
saved <- readRDS(fit_path)
fitmodel <- saved$fit
validate_lesson_fit_identity(fitmodel, input)
stopifnot(
  identical(saved$source_hashes, input$source_hashes),
  identical(saved$input_md5, lesson$input_md5),
  identical(saved$mcmc, manifest$mcmc),
  identical(dim(fitmodel$results_output$jsdm_output$U_output), c(100L, 2L, 6000L, 4L)),
  identical(dim(fitmodel$results_output$jsdm_output$L_output), c(2L, 10L, 6000L, 4L)),
  identical(as.character(fitmodel$infos$siteNames), rownames(input$sim$true_params$z_true)),
  identical(fitmodel$infos$speciesNames, colnames(input$sim$true_params$z_true))
)
known_truth <- input$sim$true_params
true_factors <- known_truth$jsdmParams_true
reconstructed <- sweep(fitmodel$X_psi %*% true_factors$B +
                         true_factors$U %*% true_factors$L,
                       2, true_factors$B0, "+")
stopifnot(max(abs(reconstructed - true_factors$eta)) < 1e-12)

# Also compare the relevant installed R function bodies against this source.
api_functions <- c("returnOrdinationScores", "returnFactorLoadings",
                   "plotOrdinationScores", "plotFactorLoadings", "plotBiplot",
                   "returnFactorScores", "returnFactorLoadings_jsdm",
                   "plotFactorScores", "plotFactorLoadings_jsdm", "reparamFactorModel")
source_environment <- new.env(parent = globalenv())
sys.source("R/jsdmfun.R", envir = source_environment)
sys.source("R/output.R", envir = source_environment)
for (name in api_functions) {
  installed <- getFromNamespace(name, "occJSDM")
  current <- get(name, source_environment)
  stopifnot(identical(formals(installed), formals(current)),
            identical(deparse(body(installed)), deparse(body(current))))
}

# Evaluate exactly the displayed code, including the ordinary native calls.
# purl=FALSE excludes compact-bundle reads, tables and saved-image displays.
code_rmd <- tempfile(fileext = ".Rmd")
code_r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), code_rmd)
knitr::purl(code_rmd, output = code_r, quiet = TRUE, documentation = 0)
student_code_md5 <- md5(code_r)
student_environment <- new.env(parent = globalenv())
student_environment$fitmodel <- fitmodel
student_environment$known_truth <- known_truth
source(code_r, local = student_environment, echo = FALSE, print.eval = FALSE)
unlink(c(code_rmd, code_r))
stopifnot(identical(student_environment$fitmodel, saved$fit), identical(fitmodel, saved$fit))

# Check every draw, not just posterior means or a plotting subset.
original <- fitmodel$results_output$jsdm_output
aligned <- student_environment$fit_for_ordination$results_output$jsdm_output
maximum_error <- c(contribution = 0, species_covariance = 0,
                   site_gram = 0, orthogonality = 0, optimum = 0)
minimum_singular_value <- Inf
for (chain in seq_len(dim(original$L_output)[4])) {
  for (iteration in seq_len(dim(original$L_output)[3])) {
    scores <- original$U_output[, , iteration, chain]
    loadings <- original$L_output[, , iteration, chain]
    scores_aligned <- aligned$U_output[, , iteration, chain]
    loadings_aligned <- aligned$L_output[, , iteration, chain]
    decomposition <- svd(loadings %*% t(true_factors$L))
    rotation <- decomposition$u %*% t(decomposition$v)
    optimum <- sum(loadings^2) + sum(true_factors$L^2) - 2 * sum(decomposition$d)
    errors <- c(
      contribution = max(abs(scores %*% loadings - scores_aligned %*% loadings_aligned)),
      species_covariance = max(abs(crossprod(loadings) - crossprod(loadings_aligned))),
      site_gram = max(abs(tcrossprod(scores) - tcrossprod(scores_aligned))),
      orthogonality = max(abs(crossprod(rotation) - diag(2))),
      optimum = abs(sum((loadings_aligned - true_factors$L)^2) - optimum)
    )
    maximum_error <- pmax(maximum_error, errors)
    minimum_singular_value <- min(minimum_singular_value, decomposition$d)
  }
}
stopifnot(all(maximum_error < 1e-10), minimum_singular_value > 1e-10)

plot_names <- c("sites", "loadings", "biplot")
plots <- setNames(lapply(paste0("native_", plot_names), get,
                        envir = student_environment), plot_names)
figures <- tibble(plot = plot_names,
                  file = paste0("ordination-", plot_names, ".png"),
                  width = c(12, 12, 9), height = c(7, 7, 7))
for (index in seq_len(nrow(figures))) {
  path <- new_path(file.path("vignettes/teaching-data", figures$file[index]))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  set.seed(20260922L)
  ggsave(path, plots[[figures$plot[index]]], width = figures$width[index],
         height = figures$height[index], dpi = 150, bg = "white")
}
figures$md5 <- map_chr(figures$file, ~ md5(new_path(file.path("vignettes/teaching-data", .x))))
plot_records <- lapply(plots, function(plot) {
  built <- ggplot_build(plot)
  list(data = if (is.data.frame(plot$data)) as.data.frame(plot$data) else NULL,
       layer_data = lapply(plot$layers, function(layer) {
         if (is.data.frame(layer$data)) as.data.frame(layer$data) else NULL
       }),
       layers = lapply(built$data, as.data.frame),
       layout = built$layout$layout)
})

ordination_examples <- list(
  schema = 1L,
  truth = list(sites = student_environment$site_truth,
               loadings = student_environment$loading_truth),
  selected_sites = student_environment$shown_sites,
  score_comparison = student_environment$score_comparison,
  quantiles = list(sites = student_environment$site_quantiles,
                   loadings = student_environment$loading_quantiles),
  ordinary_quantiles = list(sites = student_environment$ordinary_site_quantiles,
                            loadings = student_environment$ordinary_loading_quantiles),
  arrow_multiplier = student_environment$arrow_multiplier,
  plots = plot_records, figures = figures,
  checks = list(draws = 24000L, maximum_error = maximum_error,
                minimum_alignment_singular_value = minimum_singular_value),
  provenance = list(
    source_hashes = lesson_source_hashes(), lesson_md5 = md5(lesson_path),
    outputs_md5 = md5(outputs_path), input_md5 = lesson$input_md5,
    fit_manifest = manifest, snippet_md5 = md5(snippet_path),
    exporter_md5 = md5(exporter_path), student_code_md5 = student_code_md5,
    package_library = normalizePath(find.package("occJSDM")),
    package_files_md5 = tools::md5sum(list.files(find.package("occJSDM"),
                                               recursive = TRUE, full.names = TRUE)),
    seed = 20260922L, session = sessionInfo(),
    note = paste("Complete default fit: 6000 iterations x 4 chains.",
                 "Each draw rotated/reflected against known generating species loadings",
                 "using orthogonal Procrustes, with scores transformed jointly.",
                 "Original fitted object unchanged; no scaling, refitting or thinning.",
                 "Truth-assisted orientation is a simulation diagnostic only.",
                 "Circle areas encode products of marginal 95% widths, not joint coverage.")
  )
)
saveRDS(ordination_examples, new_path("vignettes/teaching-data/ordination-examples.rds"),
        compress = "xz")
print(maximum_error)
cat("Exported three native ordination plots from 24000 aligned draws.\n")
