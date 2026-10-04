# Export the displayed plotCovariateEffect() call without rerunning MCMC.
# Run from the repository root with the archive library first in R_LIBS.
# Rscript dev/simstudy/vignette-lesson/covariate-effect-export.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
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
remaining_path <- "vignettes/teaching-data/remaining-plots-data.rds"
snippet_path <- "dev/simstudy/vignette-lesson/covariate-effect-examples.Rmd"
exporter_path <- "dev/simstudy/vignette-lesson/covariate-effect-export.R"
bundle_path <- "vignettes/teaching-data/covariate-effect-data.rds"
lesson <- readRDS(lesson_path)
outputs <- readRDS(output_path)
remaining <- readRDS(remaining_path)
input <- readRDS(file.path(archive, "input.rds"))
manifest <- lesson$manifests$default
fit_path <- file.path(archive, manifest$file)
stopifnot(
  identical(lesson$input, input),
  identical(input$source_hashes, lesson_source_hashes()),
  identical(outputs$source_hashes, input$source_hashes),
  identical(outputs$lesson_md5, md5(lesson_path)),
  identical(manifest, outputs$fit_manifests$default),
  identical(remaining$provenance$fit_manifest, manifest),
  identical(remaining$provenance$source_hashes, input$source_hashes),
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
api_names <- c("plotCovariateEffect", "returnCovariateEffect", "returnCovariateEffect_base",
               "covariate_response_grid", "plot_covariate_response", "create_covariates_matrix",
               "validate_covariate_response", "stopIfNoRawCovariates")
stopifnot(all(vapply(api_names, function(name) {
  identical(body(get(name, asNamespace("occJSDM"))), body(get(name, source_environment))) &&
    identical(formals(get(name, asNamespace("occJSDM"))), formals(get(name, source_environment)))
}, logical(1))))

# The generating curve on the function's own grid. The simulator applies B to
# standardized covariates, and the fit's X0_psi is already standardized.
covariate <- "X_psi.EnvCov.1"
other <- "X_psi.EnvCov.2"
scaling <- fitmodel$infos$list_X_psi_mat
mean_1 <- unname(scaling$mean_df[[covariate]])
sd_1 <- unname(scaling$sd_df[[covariate]])
stopifnot(identical(scaling$names_df, c(covariate, other)),
          all(unlist(scaling$is_numeric)),
          isTRUE(all.equal(remaining$scaling$mean, unname(scaling$mean_df), tolerance = 1e-12)),
          isTRUE(all.equal(remaining$scaling$sd, unname(scaling$sd_df), tolerance = 1e-12)))
X0 <- fitmodel$infos$X0_psi
standardized_grid <- seq(min(X0[[covariate]]), max(X0[[covariate]]), length.out = 200L)
x <- standardized_grid * sd_1 + mean_1
other_median <- median(X0[[other]])
b <- input$sim$true_params$jsdmParams_true
species <- fitmodel$infos$speciesNames
truth <- map_dfr(seq_along(species), function(s) {
  tibble(Species = species[s], x = x,
         truth = plogis(b$B0[s] + b$B[1, s] * (x - mean_1) / sd_1 + b$B[2, s] * other_median))
})
covariate_effect_examples <- list(truth = truth)

# Only explicitly exportable chunks are evaluated. The fit is never thinned.
rmd <- tempfile(fileext = ".Rmd")
r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), rmd)
knitr::purl(rmd, output = r, quiet = TRUE, documentation = 0)
student_code_md5 <- md5(r)
environment <- new.env(parent = globalenv())
environment$fitmodel <- fitmodel
environment$covariate_effect_examples <- covariate_effect_examples
source(r, local = environment, echo = FALSE, print.eval = FALSE)
unlink(c(rmd, r))
plot <- get("covariate_effect_1", envir = environment)
stopifnot(inherits(plot, "ggplot"))
figures <- tibble(file = "covariate-effect-1.png", width = 10, height = 8)
ggsave(file.path("vignettes/teaching-data", figures$file), plot,
       width = figures$width, height = figures$height, dpi = 150, bg = "white")
figures$md5 <- md5(file.path("vignettes/teaching-data", figures$file))

# How well each species' generating curve is covered by the plotted band.
plotted <- as_tibble(plot$data)
stopifnot(nrow(plotted) == nrow(truth), setequal(plotted$Species, species))
joined <- plotted |>
  mutate(point = row_number(), .by = Species) |>
  inner_join(truth |> mutate(point = row_number(), .by = Species), by = c("Species", "point"))
stopifnot(nrow(joined) == nrow(truth), max(abs(joined$x.x - joined$x.y)) < 1e-10)
checks <- joined |>
  summarise(inside_share = mean(truth >= lower & truth <= upper),
            max_gap = max(abs(median - truth)), .by = Species) |>
  arrange(match(Species, species))

covariate_effect_examples$schema <- 1L
covariate_effect_examples$figures <- figures
covariate_effect_examples$checks <- checks
covariate_effect_examples$provenance <- list(
  source_hashes = lesson_source_hashes(), lesson_md5 = md5(lesson_path),
  output_md5 = md5(output_path), remaining_md5 = md5(remaining_path),
  input_md5 = lesson$input_md5, fit_manifest = manifest,
  snippet_md5 = md5(snippet_path), exporter_md5 = md5(exporter_path),
  student_code_md5 = student_code_md5, api_names = api_names,
  package_library = normalizePath(find.package("occJSDM")),
  package_files_md5 = tools::md5sum(list.files(find.package("occJSDM"), recursive = TRUE, full.names = TRUE)),
  session = sessionInfo(),
  note = "The native plot uses the unchanged complete 6000 x 4 fit. Gradient 2 is held at its standardized median; site factors equal zero."
)
saveRDS(covariate_effect_examples, bundle_path, compress = "xz")
cat(sprintf("Exported plotCovariateEffect() for %s: one PNG, %d truth points, coverage checks and exact-code provenance.\n",
            covariate, nrow(truth)))
