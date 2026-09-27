# Native trait plots from the unchanged complete default fit. No MCMC.
# Run from repository root with archive/library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/native-traits-export.R FULL_FIT_DIRECTORY
# NATIVE_TRAITS_STAGE optionally redirects these new deliverables for review.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("NATIVE_TRAITS_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))
library(occJSDM)
library(dplyr)
library(ggplot2)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
outputs_path <- "vignettes/teaching-data/output-lesson.rds"
snippet_path <- new_path("dev/simstudy/vignette-lesson/native-traits-examples.Rmd")
exporter_path <- new_path("dev/simstudy/vignette-lesson/native-traits-export.R")
lesson <- readRDS(lesson_path)
outputs <- readRDS(outputs_path)
input <- readRDS(file.path(archive, "input.rds"))
manifest <- lesson$manifests$default
fit_path <- file.path(archive, manifest$file)
stopifnot(
  identical(lesson$input, input),
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
stopifnot(identical(saved$source_hashes, input$source_hashes),
          identical(saved$input_md5, lesson$input_md5),
          identical(saved$mcmc, manifest$mcmc),
          identical(dim(fitmodel$results_output$jsdm_output$G_output), c(2L, 2L, 6000L, 4L)),
          identical(colnames(fitmodel$Tr), colnames(input$sim$data_list$traits)),
          identical(colnames(fitmodel$X_psi),
                    grep("^X_psi", names(input$sim$data_list$info), value = TRUE)),
          max(abs(fitmodel$Tr - scale(input$sim$data_list$traits))) < 1e-12)

# The simulator standardizes environment before generating B. Reconstruct eta
# to check that no extra environmental rescaling belongs in the trait truth.
factors <- input$sim$true_params$jsdmParams_true
eta <- sweep(fitmodel$X_psi %*% factors$B + factors$U %*% factors$L,
             2, factors$B0, "+")
stopifnot(max(abs(eta - factors$eta)) < 1e-12)
trait_sd <- apply(input$sim$data_list$traits, 2, sd)
trait_truth <- sweep(factors$G, 1, trait_sd, "*")
expected <- outputs$coefficients |>
  filter(arm == "default", block == "Trait")
for (index in seq_len(nrow(expected))) {
  row <- expected[index, ]
  trait <- match(row$term, colnames(fitmodel$Tr))
  environment <- match(row$covariate, colnames(fitmodel$X_psi))
  stopifnot(abs(row$truth - trait_truth[trait, environment]) < 1e-12)
}

# Require installed functions to match the checked current source as well.
source_environment <- new.env(parent = globalenv())
sys.source("R/jsdmfun.R", envir = source_environment)
sys.source("R/output.R", envir = source_environment)
for (name in c("plotTraitsCoefficients", "returnTraitsCoeff", "plotCoefficient")) {
  installed <- getFromNamespace(name, "occJSDM")
  current <- get(name, source_environment)
  stopifnot(identical(formals(installed), formals(current)),
            identical(deparse(body(installed)), deparse(body(current))))
}

code_rmd <- tempfile(fileext = ".Rmd")
code_r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), code_rmd)
knitr::purl(code_rmd, output = code_r, quiet = TRUE, documentation = 0)
student_code_md5 <- md5(code_r)
student_environment <- new.env(parent = globalenv())
student_environment$fitmodel <- fitmodel
student_environment$outputs <- outputs
source(code_r, local = student_environment, echo = FALSE, print.eval = FALSE)
unlink(c(code_rmd, code_r))
stopifnot(identical(student_environment$fitmodel, saved$fit))
plots <- list(gradient_1 = student_environment$native_traits_1,
              gradient_2 = student_environment$native_traits_2)
figures <- data.frame(plot = names(plots),
                      file = c("native-traits-gradient-1.png", "native-traits-gradient-2.png"),
                      width = 8, height = 4.8)
for (index in seq_len(nrow(figures))) {
  path <- new_path(file.path("vignettes/teaching-data", figures$file[index]))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plots[[figures$plot[index]]], width = figures$width[index],
         height = figures$height[index], dpi = 150, bg = "white")
}
figures$md5 <- vapply(figures$file, function(file) {
  md5(new_path(file.path("vignettes/teaching-data", file)))
}, character(1))
plot_records <- lapply(plots, function(plot) {
  built <- ggplot_build(plot)
  list(data = as.data.frame(plot$data),
       layers = lapply(built$data, as.data.frame),
       x_labels = built$layout$panel_scales_x[[1]]$get_limits())
})
result <- list(
  schema = 1L, truth = student_environment$native_trait_truth,
  scaling = data.frame(trait = names(trait_sd), sd = unname(trait_sd)),
  plots = plot_records, figures = figures,
  provenance = list(
    source_hashes = lesson_source_hashes(), lesson_md5 = md5(lesson_path),
    outputs_md5 = md5(outputs_path), input_md5 = lesson$input_md5,
    fit_manifest = manifest, snippet_md5 = md5(snippet_path),
    exporter_md5 = md5(exporter_path), student_code_md5 = student_code_md5,
    package_library = normalizePath(find.package("occJSDM")),
    package_files_md5 = tools::md5sum(list.files(find.package("occJSDM"),
                                               recursive = TRUE, full.names = TRUE)),
    session = sessionInfo(),
    note = "Native trait intervals use all 6000 x 4 draws from the unchanged default fit."
  )
)
saveRDS(result, new_path("vignettes/teaching-data/native-traits-data.rds"), compress = "xz")
cat("Exported both native trait plots from all 24000 retained draws.\n")
