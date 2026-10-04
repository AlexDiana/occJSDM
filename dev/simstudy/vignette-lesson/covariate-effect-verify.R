# Independent numerical and graphical checks of the plotCovariateEffect() export,
# using the unchanged full fit. Run from the repository root with the archive
# library first in R_LIBS.
# Rscript dev/simstudy/vignette-lesson/covariate-effect-verify.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
md5 <- function(path) unname(tools::md5sum(path))
close <- function(a, b, tolerance = 1e-10) {
  stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)),
            max(abs(as.numeric(a) - as.numeric(b))) < tolerance)
}
library(occJSDM)
library(ggplot2)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
bundle_path <- "vignettes/teaching-data/covariate-effect-data.rds"
x <- readRDS(bundle_path)
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
input <- lesson$input
p <- x$provenance
snippet_path <- "dev/simstudy/vignette-lesson/covariate-effect-examples.Rmd"
stopifnot(
  identical(x$schema, 1L),
  identical(p$source_hashes, lesson_source_hashes()),
  identical(p$lesson_md5, md5("vignettes/teaching-data/nonspatial-lesson.rds")),
  identical(p$output_md5, md5("vignettes/teaching-data/output-lesson.rds")),
  identical(p$remaining_md5, md5("vignettes/teaching-data/remaining-plots-data.rds")),
  identical(p$snippet_md5, md5(snippet_path)),
  identical(p$exporter_md5, md5("dev/simstudy/vignette-lesson/covariate-effect-export.R")),
  identical(p$input_md5, md5(file.path(archive, "input.rds"))),
  identical(p$fit_manifest, lesson$manifests$default),
  identical(p$package_files_md5, tools::md5sum(names(p$package_files_md5))),
  identical(normalizePath(find.package("occJSDM")), p$package_library)
)
fit_path <- file.path(archive, p$fit_manifest$file)
stopifnot(identical(md5(fit_path), p$fit_manifest$md5))
saved <- readRDS(fit_path)
f <- saved$fit
validate_lesson_fit_identity(f, input)
stopifnot(identical(saved$source_hashes, p$source_hashes),
          identical(saved$input_md5, p$input_md5),
          identical(saved$mcmc, p$fit_manifest$mcmc),
          identical(dim(f$results_output$p_output)[3:4], c(6000L, 4L)))
source_environment <- new.env(parent = globalenv())
sys.source("R/output.R", source_environment)
sys.source("R/jsdmfun.R", source_environment)
stopifnot(all(c("plotCovariateEffect", "returnCovariateEffect", "returnCovariateEffect_base",
                "covariate_response_grid", "plot_covariate_response",
                "create_covariates_matrix") %in% p$api_names))
for (name in p$api_names) {
  stopifnot(identical(body(get(name, asNamespace("occJSDM"))), body(get(name, source_environment))),
            identical(formals(get(name, asNamespace("occJSDM"))), formals(get(name, source_environment))))
}

# Prove that the displayed chunks are those evaluated by the exporter, and
# rebuild the exported plot from them for the checks below.
rmd <- tempfile(fileext = ".Rmd")
r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), rmd)
knitr::purl(rmd, output = r, documentation = 0, quiet = TRUE)
stopifnot(identical(md5(r), p$student_code_md5))
environment <- new.env(parent = globalenv())
environment$fitmodel <- f
environment$covariate_effect_examples <- x
source(r, local = environment, echo = FALSE, print.eval = FALSE)
unlink(c(rmd, r))
plot <- get("covariate_effect_1", envir = environment)
exportable_chunks <- function(lines) {
  starts <- grep("^```\\{r covariate-effect-.*purl=TRUE", lines)
  setNames(lapply(starts, function(start) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }), sub("^```\\{r ([^,]+),.*", "\\1", lines[starts]))
}
expected_chunks <- exportable_chunks(readLines(snippet_path))
lesson_chunks <- exportable_chunks(readLines("vignettes/occJSDM-lesson-3.Rmd"))
stopifnot(length(expected_chunks) == 2L)
if (length(lesson_chunks)) {
  stopifnot(identical(sort(names(lesson_chunks)), sort(names(expected_chunks))),
            identical(expected_chunks, lesson_chunks[names(expected_chunks)]))
  cat("Lesson 3 displays the exact exported covariate-effect plotting body.\n")
}

# The scaling, the simulator's scale, and the grid, from the raw survey table.
covariate <- "X_psi.EnvCov.1"
covariates <- colnames(f$X_psi)
stopifnot(identical(covariates, c(covariate, "X_psi.EnvCov.2")))
sp <- f$infos$speciesNames
b <- input$sim$true_params$jsdmParams_true
info <- input$sim$data_list$info
raw <- as.matrix(info[!duplicated(info$Site), covariates])
mean_1 <- f$infos$list_X_psi_mat$mean_df[[covariate]]
sd_1 <- f$infos$list_X_psi_mat$sd_df[[covariate]]
close(c(mean_1, sd_1), c(mean(raw[, 1]), sd(raw[, 1])))
close(as.matrix(f$infos$X0_psi), scale(raw))
close(f$X_psi, scale(raw))
close(f$X_psi %*% b$B + b$U %*% b$L + rep(b$B0, each = nrow(raw)), b$eta)
stopifnot(identical(sp, colnames(input$sim$data_list$OTU)))
other_median <- median(scale(raw)[, 2])

d <- as.data.frame(plot$data)
stopifnot(nrow(d) == 200L * length(sp), setequal(d$Species, sp))
j <- f$results_output$jsdm_output
built <- ggplot_build(plot)
for (s in seq_along(sp)) {
  rows <- which(d$Species == sp[s])
  grid <- d$x[rows]
  close(range(grid), range(raw[, 1]))
  close(diff(grid), rep((max(raw[, 1]) - min(raw[, 1])) / 199, 199), 1e-9)
  # Independently form each linear predictor from all 24,000 matched draws.
  z <- (grid - mean(raw[, 1])) / sd(raw[, 1])
  eta <- outer(z, as.vector(j$B_output[1, s, , ])) +
    matrix(as.vector(j$B0_output[s, , ]) + other_median * as.vector(j$B_output[2, s, , ]),
           nrow = length(z), ncol = length(j$B0_output[s, , ]), byrow = TRUE)
  expected <- t(apply(plogis(eta), 1, quantile, probs = c(.5, .025, .975), names = FALSE))
  close(d$median[rows], expected[, 1])
  close(d$lower[rows], expected[, 2])
  close(d$upper[rows], expected[, 3])
  trows <- which(x$truth$Species == sp[s])
  close(x$truth$x[trows], grid)
  close(x$truth$truth[trows], plogis(b$B0[s] + b$B[1, s] * z + b$B[2, s] * other_median))
}
stopifnot(nrow(x$truth) == nrow(d))
cat("All 2000 medians and 95% limits recomputed from the draws; truth recomputed from raw units.\n")

# Check the ribbon, the median line and the dashed truth against each panel's label.
key <- function(species, value) paste(species, sprintf("%.10f", value), sep = "/")
data_key <- key(d$Species, d$x)
truth_key <- key(x$truth$Species, x$truth$x)
layout <- built$layout$layout
stopifnot(length(built$data) == 3L, setequal(layout$Species, sp))
for (layer in 1:3) {
  points <- built$data[[layer]]
  keys <- key(layout$Species[match(points$PANEL, layout$PANEL)], points$x)
  stopifnot(nrow(points) == nrow(d), !anyDuplicated(keys), setequal(keys, data_key))
  if (layer == 1L) {
    close(points$ymin, d$lower[match(keys, data_key)])
    close(points$ymax, d$upper[match(keys, data_key)])
  } else if (layer == 2L) {
    close(points$y, d$median[match(keys, data_key)])
  } else {
    close(points$y, x$truth$truth[match(keys, truth_key)])
    stopifnot(all(points$linetype == "dashed"), all(points$colour == "black"))
  }
}

# The recorded coverage checks, recomputed.
inside <- sapply(sp, function(s) {
  rows <- d$Species == s
  t <- x$truth$truth[match(data_key[rows], truth_key)]
  c(mean(t >= d$lower[rows] & t <= d$upper[rows]), max(abs(d$median[rows] - t)))
})
stopifnot(identical(x$checks$Species, sp))
close(x$checks$inside_share, inside[1, ])
close(x$checks$max_gap, inside[2, ])
cat(sprintf("Truth inside the band at every grid point for %d of %d species; largest median-truth gap %.3f.\n",
            sum(inside[1, ] == 1), length(sp), max(inside[2, ])))

file <- file.path("vignettes/teaching-data", x$figures$file)
stopifnot(nrow(x$figures) == 1L, identical(md5(file), x$figures$md5), file.info(file)$size > 5000,
          identical(as.integer(readBin(file, "raw", n = 8)), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)),
          file.info(bundle_path)$size < 1e6)
cat("Covariate-effect figure, truth, independent medians and exact displayed code verified.\n")
