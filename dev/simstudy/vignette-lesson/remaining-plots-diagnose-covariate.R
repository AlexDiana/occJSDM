# Cheap reproduction of the existing plotCovariateEffect numeric response bug.
# No MCMC and no changes to source, saved fit, package, or its namespace.
# Run from repository root with the archive library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/remaining-plots-diagnose-covariate.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
library(occJSDM)
library(dplyr)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
path <- file.path(archive, lesson$manifests$default$file)
stopifnot(identical(unname(tools::md5sum(path)), lesson$manifests$default$md5),
          identical(lesson$input$source_hashes, lesson_source_hashes()),
          identical(normalizePath(find.package("occJSDM")), normalizePath(file.path(archive, "library/occJSDM"))))
saved <- readRDS(path)
f <- saved$fit
validate_lesson_fit_identity(f, lesson$input)
stopifnot(identical(saved$source_hashes, lesson$input$source_hashes),
          identical(dim(f$results_output$p_output)[3:4], c(6000L, 4L)))
source_environment <- new.env(parent = globalenv())
sys.source("R/jsdmfun.R", source_environment)
sys.source("R/output.R", source_environment)
for (name in c("plotCovariateEffect", "plotCovariateEffect_base", "returnCovariateEffect_base", "create_covariates_matrix")) {
  stopifnot(identical(body(get(name, asNamespace("occJSDM"))), body(get(name, source_environment))))
}
covariate <- "X_psi.EnvCov.1"
old <- plotCovariateEffect(f, covNames = covariate, idx_species = c(1, 10))[[1]]$data
print(old |>
        group_by(Species) |>
        summarise(minimum_median = min(mean), maximum_median = max(mean),
                  medians_outside_probability_scale = sum(mean < 0 | mean > 1),
                  grid_points = n(), .groups = "drop"))
stopifnot(all(old$mean[old$Species == "OTU_1"] > 1),
          all(old$mean[old$Species == "OTU_10"] < 0))

# Numeric helper uses B0 + plogis(X_sub * B), mixing incompatible scales.
j <- f$results_output$jsdm_output
scaling <- f$infos$list_X_psi_mat
first <- old[old$Species == "OTU_1", ][1, ]
standardized_again <- (first$x - scaling$mean_df[[covariate]]) / scaling$sd_df[[covariate]]
wrong_draws <- as.vector(j$B0_output[1, , ]) + plogis(standardized_again * as.vector(j$B_output[1, 1, , ]))
stopifnot(isTRUE(all.equal(unname(quantile(wrong_draws, c(.025, .5, .975))),
                          unname(unlist(first[c("lower", "mean", "upper")])), tolerance = 1e-12)))

# The purported raw X0 grid is already standardized in this saved fit.
info <- f$infos$data_info
raw <- info[!duplicated(info$Site), covariate]
stored <- f$infos$X0_psi[[covariate]]
stopifnot(isTRUE(all.equal(stored, f$X_psi[, covariate], check.attributes = FALSE)),
          !isTRUE(all.equal(stored, raw, check.attributes = FALSE)))
print(data.frame(scale = c("original predictor", "stored X0 predictor", "effect helper x"),
                 min = c(min(raw), min(stored), min(old$x)),
                 max = c(max(raw), max(stored), max(old$x))))
cat("Confirmed: inverse-logit applied before adding B0; X0 holds standardized rather than raw values.\n")
cat("The wrapper also omits confidence when calling its base helper, and the numeric helper does not add other predictors.\n")
cat("Use plotOccupancyGradient for the checked zero-site-factor response curves.\n")
