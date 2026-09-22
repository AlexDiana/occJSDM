# Export the exact displayed native calls without rerunning MCMC.
# Run from the repository root with the archive library first in R_LIBS.
# Rscript dev/simstudy/vignette-lesson/remaining-plots-export.R FULL_FIT_DIRECTORY
# REMAINING_PLOT_STAGE optionally redirects this deliverable's input/output.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("REMAINING_PLOT_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))
library(occJSDM)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")

lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
output_path <- "vignettes/teaching-data/output-lesson.rds"
snippet_path <- new_path("dev/simstudy/vignette-lesson/remaining-plots-examples.Rmd")
exporter_path <- new_path("dev/simstudy/vignette-lesson/remaining-plots-export.R")
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

# Verify that the actual loaded plotting code is the declared source.
source_environment <- new.env(parent = globalenv())
sys.source("R/output.R", source_environment)
sys.source("R/jsdmfun.R", source_environment)
api_names <- c("returnOccupancyGradient", "plotOccupancyGradient", "plotSpeciesRates",
               "plotStage1FPRates", "plotStage2FPRates", "plotDetectionRates",
               "plotCovariateEffect", "returnCovariateEffect_base", "plotCovariateEffect_base",
               "create_covariates_matrix")
stopifnot(all(vapply(api_names, function(name) {
  identical(body(get(name, asNamespace("occJSDM"))), body(get(name, source_environment))) &&
    identical(formals(get(name, asNamespace("occJSDM"))), formals(get(name, source_environment)))
}, logical(1))))

species <- fitmodel$infos$speciesNames
primers <- as.character(fitmodel$infos$primerNames)
b <- input$sim$true_params$jsdmParams_true
covariates <- colnames(fitmodel$X_psi)
scaling <- tibble(
  covariate = covariates,
  mean = unname(fitmodel$infos$list_X_psi_mat$mean_df),
  sd = unname(fitmodel$infos$list_X_psi_mat$sd_df)
)
gradients <- map_dfr(seq_along(covariates), function(k) {
  x <- seq(quantile(fitmodel$X_psi[, k], .02),
           quantile(fitmodel$X_psi[, k], .98), length.out = 40)
  design <- matrix(apply(fitmodel$X_psi, 2, median), nrow = 40,
                   ncol = length(covariates), byrow = TRUE)
  design[, k] <- x
  probability <- plogis(sweep(design %*% b$B, 2, b$B0, "+"))
  expand_grid(species = species, x = x) |>
    mutate(covariate = covariates[k], truth = as.vector(probability),
           raw_x = scaling$mean[k] + scaling$sd[k] * x,
           species = factor(species, levels = fitmodel$infos$speciesNames))
})
rate_truth <- expand_grid(Species = species, Primer = primers) |>
  mutate(
    p = map2_dbl(Species, Primer, ~ input$params$p[match(.y, primers), match(.x, species)] *
                   pnorm(log(1.5), input$params$mu1, input$params$sigma1, lower.tail = FALSE)),
    q = map2_dbl(Species, Primer, ~ input$params$q[match(.y, primers), match(.x, species)] *
                   pnorm(log(1.5), input$params$mu0, input$params$sigma0, lower.tail = FALSE)),
    Primer = factor(Primer, levels = primers)
  )
remaining_examples <- list(truth = list(
  gradients = gradients,
  stage1_fp = tibble(Species = species, truth = input$params$theta0),
  stage2_fp = rate_truth |> select(Species, Primer, truth = q),
  detection = rate_truth |> select(Species, Primer, truth = p)
))

# Only explicitly exportable chunks are evaluated. The fit is never thinned.
rmd <- tempfile(fileext = ".Rmd")
r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), rmd)
knitr::purl(rmd, output = r, quiet = TRUE, documentation = 0)
student_code_md5 <- md5(r)
environment <- new.env(parent = globalenv())
environment$fitmodel <- fitmodel
environment$remaining_examples <- remaining_examples
source(r, local = environment, echo = FALSE, print.eval = FALSE)
unlink(c(rmd, r))
plot_names <- c("gradient_1", "gradient_2", "stage1_fp", "stage2_fp", "detection")
plots <- setNames(lapply(paste0("remaining_", plot_names), get, envir = environment), plot_names)
figures <- tibble(plot = plot_names,
                  file = paste0("remaining-plots-", gsub("_", "-", plot_names), ".png"),
                  width = 10, height = c(8, 8, 5.5, 6, 6))
for (i in seq_len(nrow(figures))) {
  ggsave(new_path(file.path("vignettes/teaching-data", figures$file[i])),
         plots[[i]], width = figures$width[i], height = figures$height[i], dpi = 150, bg = "white")
}
figures$md5 <- map_chr(figures$file, ~ md5(new_path(file.path("vignettes/teaching-data", .x))))
plot_records <- lapply(plots, function(plot) {
  built <- ggplot_build(plot)
  native_data <- if (is.data.frame(plot$data)) plot$data else plot$layers[[1]]$data
  list(data = as.data.frame(native_data), layers = lapply(built$data, as.data.frame),
       x_limits = lapply(built$layout$panel_scales_x, function(s) s$get_limits()),
       y_limits = lapply(built$layout$panel_scales_y, function(s) s$get_limits()),
       layout = built$layout$layout,
       primer_colours = if (is.null(built$plot$scales$get_scales("colour"))) NULL else
         setNames(built$plot$scales$get_scales("colour")$map(primers), primers))
})
remaining_examples$schema <- 1L
remaining_examples$scaling <- scaling
remaining_examples$plots <- plot_records
remaining_examples$figures <- figures
remaining_examples$provenance <- list(
  source_hashes = lesson_source_hashes(), lesson_md5 = md5(lesson_path),
  output_md5 = md5(output_path), input_md5 = lesson$input_md5,
  fit_manifest = manifest, snippet_md5 = md5(snippet_path), exporter_md5 = md5(exporter_path),
  student_code_md5 = student_code_md5, api_names = api_names,
  package_library = normalizePath(find.package("occJSDM")),
  package_files_md5 = tools::md5sum(list.files(find.package("occJSDM"), recursive = TRUE, full.names = TRUE)),
  session = sessionInfo(),
  note = "Native plots use the unchanged complete 6000 x 4 fit. Site factors equal zero in gradient curves."
)
saveRDS(remaining_examples, new_path("vignettes/teaching-data/remaining-plots-data.rds"), compress = "xz")
cat("Exported five native plots, matching truth, scaling, and exact-code provenance.\n")
