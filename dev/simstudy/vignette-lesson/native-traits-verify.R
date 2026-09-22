# Check native plotted intervals directly from G draws and original truth.
# Run from repository root with archive/library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/native-traits-verify.R FULL_FIT_DIRECTORY
# NATIVE_TRAITS_STAGE optionally points to staged deliverables.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("NATIVE_TRAITS_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))
close <- function(left, right) {
  stopifnot(length(left) == length(right), all(is.finite(left)), all(is.finite(right)),
            max(abs(as.numeric(left) - as.numeric(right))) < 1e-10)
}
library(occJSDM)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
examples <- readRDS(new_path("vignettes/teaching-data/native-traits-data.rds"))
provenance <- examples$provenance
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
input <- readRDS(file.path(archive, "input.rds"))
outputs <- readRDS("vignettes/teaching-data/output-lesson.rds")
snippet_path <- new_path("dev/simstudy/vignette-lesson/native-traits-examples.Rmd")
stopifnot(
  identical(lesson$input, input),
  identical(provenance$source_hashes, lesson_source_hashes()),
  identical(provenance$source_hashes, input$source_hashes),
  identical(provenance$source_hashes, outputs$source_hashes),
  identical(provenance$lesson_md5, md5("vignettes/teaching-data/nonspatial-lesson.rds")),
  identical(provenance$outputs_md5, md5("vignettes/teaching-data/output-lesson.rds")),
  identical(provenance$input_md5, md5(file.path(archive, "input.rds"))),
  identical(provenance$input_md5, lesson$input_md5),
  identical(provenance$fit_manifest, lesson$manifests$default),
  identical(provenance$fit_manifest, outputs$fit_manifests$default),
  identical(provenance$snippet_md5, md5(snippet_path)),
  identical(provenance$exporter_md5,
            md5(new_path("dev/simstudy/vignette-lesson/native-traits-export.R"))),
  identical(provenance$package_files_md5, tools::md5sum(names(provenance$package_files_md5))),
  identical(normalizePath(find.package("occJSDM")), provenance$package_library)
)
fit_path <- file.path(archive, provenance$fit_manifest$file)
stopifnot(identical(md5(fit_path), provenance$fit_manifest$md5))
saved <- readRDS(fit_path)
fit <- saved$fit
validate_lesson_fit_identity(fit, input)
stopifnot(identical(saved$source_hashes, provenance$source_hashes),
          identical(saved$input_md5, provenance$input_md5),
          identical(saved$mcmc, provenance$fit_manifest$mcmc),
          identical(dim(fit$results_output$jsdm_output$G_output), c(2L, 2L, 6000L, 4L)))
traits <- colnames(input$sim$data_list$traits)
environments <- grep("^X_psi", names(input$sim$data_list$info), value = TRUE)
stopifnot(identical(colnames(fit$Tr), traits),
          identical(colnames(fit$X_psi), environments))
raw_traits <- input$sim$data_list$traits
trait_means <- colMeans(raw_traits)
# Direct sample SD and transformed design independently of exporter sweep.
trait_sd <- sqrt(colSums((raw_traits - rep(trait_means, each = nrow(raw_traits)))^2) /
                   (nrow(raw_traits) - 1))
close(fit$Tr, (raw_traits - rep(trait_means, each = nrow(raw_traits))) /
        rep(trait_sd, each = nrow(raw_traits)))
close(examples$scaling$sd, trait_sd)
stopifnot(identical(examples$scaling$trait, traits))
factors <- input$sim$true_params$jsdmParams_true
close(fit$X_psi %*% factors$B + factors$U %*% factors$L +
        rep(factors$B0, each = nrow(fit$X_psi)), factors$eta)

for (environment in seq_along(environments)) {
  record <- examples$plots[[environment]]
  bars <- record$layers[[1]]
  crosses <- record$layers[[3]]
  stopifnot(nrow(bars) == length(traits), nrow(crosses) == length(traits),
            all(crosses$shape == 4), all(record$layers[[2]]$yintercept == 0),
            identical(sort(as.character(record$x_labels)), sort(traits)))
  for (row in seq_len(nrow(bars))) {
    trait <- match(record$x_labels[as.integer(bars$x[row])], traits)
    stopifnot(!is.na(trait))
    draws <- fit$results_output$jsdm_output$G_output[trait, environment, , ]
    close(c(bars$ymin[row], bars$ymax[row]), quantile(draws, c(.025, .975)))
  }
  for (row in seq_len(nrow(crosses))) {
    trait_name <- record$x_labels[as.integer(crosses$x[row])]
    trait <- match(trait_name, traits)
    value <- factors$G[trait, environment] * trait_sd[trait]
    close(crosses$y[row], value)
    table_row <- examples$truth[examples$truth$covariate == environments[environment] &
                                 examples$truth$trait == trait_name, ]
    stopifnot(nrow(table_row) == 1L)
    close(table_row$truth, value)
  }
}

# Exported bodies must match both canonical snippet and integrated Lesson 3.
code_rmd <- tempfile(fileext = ".Rmd")
code_r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), code_rmd)
knitr::purl(code_rmd, output = code_r, quiet = TRUE, documentation = 0)
stopifnot(identical(md5(code_r), provenance$student_code_md5))
unlink(c(code_rmd, code_r))
exportable_chunks <- function(lines) {
  starts <- grep("^```\\{r native-traits-.*purl=TRUE", lines)
  setNames(lapply(starts, function(start) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }), sub("^```\\{r ([^,]+),.*", "\\1", lines[starts]))
}
expected <- exportable_chunks(readLines(snippet_path))
actual <- exportable_chunks(readLines("vignettes/occJSDM-lesson-3.Rmd"))
stopifnot(length(expected) == 3L)
if (any(names(expected) %in% names(actual))) {
  stopifnot(all(names(expected) %in% names(actual)),
            identical(expected, actual[names(expected)]))
  cat("Lesson 3 displays the exact exported native trait code.\n")
}
for (index in seq_len(nrow(examples$figures))) {
  figure <- examples$figures[index, ]
  path <- new_path(file.path("vignettes/teaching-data", figure$file))
  stopifnot(identical(md5(path), figure$md5))
  stopifnot(identical(dim(png::readPNG(path, native = TRUE)),
                      as.integer(c(figure$height, figure$width) * 150)))
}
cat("PASS: native intervals from all 24000 draws, trait/environment identities,\n")
cat("generating trait-scale conversion, plotted truth coordinates, code and PNG evidence.\n")
