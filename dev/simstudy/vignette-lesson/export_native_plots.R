# Export native package plots from the existing full fit; never run MCMC.
# Run from the repository root with the archive's library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/export_native_plots.R FULL_FIT_DIRECTORY
# Optional NATIVE_PLOT_STAGE redirects only this deliverable's paths for review.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("NATIVE_PLOT_STAGE", ".")
new_path <- function(path) file.path(stage, path)

library(occJSDM)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")

md5 <- function(path) unname(tools::md5sum(path))
lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
output_path <- "vignettes/teaching-data/output-lesson.rds"
snippet_path <- new_path("dev/simstudy/vignette-lesson/native-plot-examples.Rmd")
exporter_path <- new_path("dev/simstudy/vignette-lesson/export_native_plots.R")
lesson <- readRDS(lesson_path)
outputs <- readRDS(output_path)
input <- readRDS(file.path(archive, "input.rds"))
manifest <- lesson$manifests$default
fit_path <- file.path(archive, manifest$file)
stopifnot(
  identical(lesson$input, input),
  identical(input$source_hashes, lesson_source_hashes()),
  identical(outputs$source_hashes, input$source_hashes),
  identical(outputs$lesson_md5, md5(lesson_path)),
  identical(manifest, outputs$fit_manifests$default),
  identical(manifest$md5, md5(fit_path)),
  identical(lesson$input_md5, md5(file.path(archive, "input.rds"))),
  identical(normalizePath(find.package("occJSDM")),
            normalizePath(file.path(archive, "library/occJSDM")))
)
saved <- readRDS(fit_path)
fitmodel <- saved$fit
validate_lesson_fit_identity(fitmodel, input)
stopifnot(identical(saved$source_hashes, input$source_hashes),
          identical(saved$input_md5, lesson$input_md5),
          identical(saved$mcmc, manifest$mcmc),
          identical(dim(fitmodel$results_output$p_output)[3:4], c(6000L, 4L)))

species <- fitmodel$infos$speciesNames
primers <- as.character(fitmodel$infos$primerNames)
truth <- input$sim$true_params
beta_collection <- truth$beta_theta_true
scale_info <- fitmodel$infos$list_X_theta_mat
beta_collection[1, ] <- beta_collection[1, ] +
  scale_info$mean_df * beta_collection[2, ]
beta_collection[2, ] <- scale_info$sd_df * beta_collection[2, ]
theta <- plogis(beta_collection[1, ])
p <- effective_detection_rate(input$params$p, input$params$mu1, input$params$sigma1)
q <- effective_detection_rate(input$params$q, input$params$mu0, input$params$sigma0)
stopifnot(identical(species, colnames(input$sim$data_list$OTU)),
          identical(primers, as.character(seq_len(input$settings$P))))

# Exact Poisson-binomial distribution: add one independent species indicator
# at a time. mass[k + 1] is the chance of detecting exactly k species.
count_distribution <- function(probabilities) {
  mass <- 1
  for (probability in probabilities) {
    mass <- c(mass * (1 - probability), 0) + c(0, mass * probability)
  }
  mass
}
survey_counts <- expand_grid(M = 1:4, K = 1:6) |>
  mutate(
    species_probability = map2(M, K, function(M, K) {
      1 - (1 - theta * (1 - apply((1 - p)^K, 2, prod)))^M
    }),
    mass = map(species_probability, count_distribution),
    lower = map_dbl(mass, ~ which(cumsum(.x) >= .025)[1] - 1),
    median = map_dbl(mass, ~ which(cumsum(.x) >= .5)[1] - 1),
    upper = map_dbl(mass, ~ which(cumsum(.x) >= .975)[1] - 1),
    expectation = map_dbl(species_probability, sum)
  )

native_truth <- list(
  environment = outputs$coefficients |>
    filter(arm == "default", block == "Environment") |>
    select(covariate, species = term, truth),
  collection = outputs$collection |>
    select(covariate, species = term, truth),
  occupancy_rates = tibble(species, truth = plogis(truth$jsdmParams_true$B0)),
  collection_rates = tibble(species, truth = theta),
  laboratory = expand_grid(Primer = primers, species = species, rate = c("p", "q")) |>
    mutate(truth = pmap_dbl(list(Primer, species, rate), function(Primer, species, rate) {
      values <- if (rate == "p") p else q
      values[match(Primer, primers), match(species, fitmodel$infos$speciesNames)]
    })),
  correlations = outputs$correlations |>
    filter(arm == "default", match(species1, species) > match(species2, species)) |>
    mutate(truth_label = if_else(is.na(truth), "NA", sprintf("%.2f", truth))) |>
    select(species1, species2, truth, truth_label),
  survey_counts = survey_counts
)

# Evaluate the exact displayed plotting bodies. Change eval only for the
# explicitly exportable chunks. All image-reading chunks have purl=FALSE.
native_examples <- list(truth = native_truth)
snippet <- readLines(snippet_path)
snippet <- sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE", snippet, fixed = TRUE)
temporary_rmd <- tempfile(fileext = ".Rmd")
temporary_r <- tempfile(fileext = ".R")
writeLines(snippet, temporary_rmd)
knitr::purl(temporary_rmd, output = temporary_r, quiet = TRUE, documentation = 0)
student_code_md5 <- md5(temporary_r)
environment <- new.env(parent = globalenv())
environment$fitmodel <- fitmodel
environment$native_examples <- native_examples
source(temporary_r, local = environment, echo = FALSE, print.eval = FALSE)
unlink(c(temporary_rmd, temporary_r))

plot_names <- c("environment", "collection", "occupancy_rates", "collection_rates",
                "primer_1", "primer_2", "correlations", "effort_k", "effort_m")
plots <- setNames(lapply(paste0("native_", plot_names), get, envir = environment), plot_names)
figures <- tibble(
  plot = plot_names,
  file = paste0("native-plot-", gsub("_", "-", plot_names), ".png"),
  width = 9, height = c(4.8, 4.8, 5.2, 5.2, 5.5, 5.5, 7, 7.5, 9)
)
for (i in seq_len(nrow(figures))) {
  ggsave(new_path(file.path("vignettes/teaching-data", figures$file[i])),
         plots[[figures$plot[i]]], width = figures$width[i], height = figures$height[i],
         dpi = 150, bg = "white")
}
figures$md5 <- map_chr(figures$file, ~ md5(new_path(file.path("vignettes/teaching-data", .x))))

# Keep only plotted tables and built coordinates, not plot expression
# environments (which capture the fit) or any posterior arrays.
plot_records <- lapply(plots, function(plot) {
  built <- ggplot_build(plot)
  list(data = as.data.frame(plot$data), layers = lapply(built$data, as.data.frame),
       x_limits = lapply(built$layout$panel_scales_x, function(scale) scale$get_limits()),
       y_limits = lapply(built$layout$panel_scales_y, function(scale) scale$get_limits()),
       layout = built$layout$layout)
})
native_examples <- list(
  schema = 1L, truth = native_truth, plots = plot_records, figures = figures,
  provenance = list(
    source_hashes = lesson_source_hashes(),
    lesson_md5 = md5(lesson_path), output_md5 = md5(output_path),
    input_md5 = lesson$input_md5, fit_manifest = manifest,
    snippet_md5 = md5(snippet_path), exporter_md5 = md5(exporter_path),
    student_code_md5 = student_code_md5,
    package_library = normalizePath(find.package("occJSDM")),
    package_files_md5 = tools::md5sum(list.files(find.package("occJSDM"),
                                               recursive = TRUE, full.names = TRUE)),
    seed = 20260922L, session = sessionInfo(),
    note = paste("Native calls used the unchanged complete 6000 x 4 fit.",
                 "The cumulative routine internally selects 500 draws and simulates surveys.",
                 "Truth quantiles use the exact Poisson-binomial distribution.",
                 "Ordination awaits a declared, validated joint factor alignment.")
  )
)
saveRDS(native_examples, new_path("vignettes/teaching-data/native-plots.rds"), compress = "xz")
cat("Exported", nrow(figures), "native plots and compact truth/provenance tables.\n")
