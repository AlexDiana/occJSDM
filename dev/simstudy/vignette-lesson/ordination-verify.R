# Independent verification of identity, orientation, quantiles and native plots.
# Run from repository root with archive/library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/ordination-verify.R FULL_FIT_DIRECTORY
# ORDINATION_STAGE optionally points to the staged new deliverables.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("ORDINATION_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))
close <- function(left, right, tolerance = 1e-9) {
  stopifnot(length(left) == length(right), all(is.finite(left)), all(is.finite(right)),
            max(abs(as.numeric(left) - as.numeric(right))) < tolerance)
}
library(occJSDM)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
examples <- readRDS(new_path("vignettes/teaching-data/ordination-examples.rds"))
provenance <- examples$provenance
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
outputs <- readRDS("vignettes/teaching-data/output-lesson.rds")
input <- readRDS(file.path(archive, "input.rds"))
snippet_path <- new_path("dev/simstudy/vignette-lesson/ordination-examples.Rmd")
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
            md5(new_path("dev/simstudy/vignette-lesson/ordination-export.R"))),
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
          identical(saved$mcmc, provenance$fit_manifest$mcmc))
original <- fit$results_output$jsdm_output
stopifnot(identical(dim(original$U_output), c(100L, 2L, 6000L, 4L)),
          identical(dim(original$L_output), c(2L, 10L, 6000L, 4L)))
truth <- input$sim$true_params
factors <- truth$jsdmParams_true
sites <- as.character(fit$infos$siteNames)
species <- fit$infos$speciesNames
stopifnot(identical(sites, rownames(truth$z_true)),
          identical(species, colnames(truth$z_true)),
          identical(examples$truth$sites$site, sites),
          identical(examples$truth$loadings$species, species),
          identical(examples$selected_sites, sites[seq(1, 100, by = 10)]))
close(examples$truth$sites$x, factors$U[, 1])
close(examples$truth$sites$y, factors$U[, 2])
close(examples$truth$loadings$x, factors$L[1, ])
close(examples$truth$loadings$y, factors$L[2, ])
close(sweep(fit$X_psi %*% factors$B + factors$U %*% factors$L,
            2, factors$B0, "+"), factors$eta)

# Check the exact exported student code and any installed Lesson 3 copy.
code_rmd <- tempfile(fileext = ".Rmd")
code_r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), code_rmd)
knitr::purl(code_rmd, output = code_r, quiet = TRUE, documentation = 0)
stopifnot(identical(md5(code_r), provenance$student_code_md5))
unlink(c(code_rmd, code_r))
exportable_chunks <- function(lines) {
  starts <- grep("^```\\{r ordination-.*purl=TRUE", lines)
  setNames(lapply(starts, function(start) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }), sub("^```\\{r ([^,]+),.*", "\\1", lines[starts]))
}
expected_chunks <- exportable_chunks(readLines(snippet_path))
lesson_chunks <- exportable_chunks(readLines("vignettes/occJSDM-lesson-3.Rmd"))
stopifnot(length(expected_chunks) == 5L)
if (any(names(expected_chunks) %in% names(lesson_chunks))) {
  stopifnot(all(names(expected_chunks) %in% names(lesson_chunks)),
            identical(expected_chunks, lesson_chunks[names(expected_chunks)]))
  cat("Lesson 3 displays the exact exported ordination code.\n")
}

# Solve the two-dimensional problem analytically, independently of SVD.
# Compare the best rotation (determinant +1) and reflection (determinant -1).
# Both maximize trace(Q' C), C = L L_true', within their determinant class.
# The larger maximum gives the minimum loading distance without rescaling.
score_draws <- original$U_output
loading_draws <- original$L_output
maximum_error <- c(contribution = 0, species_covariance = 0,
                   site_gram = 0, orthogonality = 0)
for (chain in seq_len(dim(loading_draws)[4])) {
  for (iteration in seq_len(dim(loading_draws)[3])) {
    scores <- original$U_output[, , iteration, chain]
    loadings <- original$L_output[, , iteration, chain]
    cross_matrix <- loadings %*% t(factors$L)
    rotation_pair <- c(cross_matrix[1, 1] + cross_matrix[2, 2],
                       cross_matrix[2, 1] - cross_matrix[1, 2])
    reflection_pair <- c(cross_matrix[1, 1] - cross_matrix[2, 2],
                         cross_matrix[2, 1] + cross_matrix[1, 2])
    rotation_score <- sqrt(sum(rotation_pair^2))
    reflection_score <- sqrt(sum(reflection_pair^2))
    if (rotation_score >= reflection_score) {
      direction <- rotation_pair / rotation_score
      rotation <- matrix(c(direction[1], direction[2], -direction[2], direction[1]), 2)
    } else {
      direction <- reflection_pair / reflection_score
      rotation <- matrix(c(direction[1], direction[2], direction[2], -direction[1]), 2)
    }
    stopifnot(all(is.finite(rotation)))
    new_scores <- scores %*% rotation
    new_loadings <- crossprod(rotation, loadings)
    score_draws[, , iteration, chain] <- new_scores
    loading_draws[, , iteration, chain] <- new_loadings
    maximum_error <- pmax(maximum_error, c(
      contribution = max(abs(scores %*% loadings - new_scores %*% new_loadings)),
      species_covariance = max(abs(crossprod(loadings) - crossprod(new_loadings))),
      site_gram = max(abs(tcrossprod(scores) - tcrossprod(new_scores))),
      orthogonality = max(abs(crossprod(rotation) - diag(2)))
    ))
  }
}
stopifnot(all(maximum_error < 1e-9), examples$checks$draws == 24000L,
          all(examples$checks$maximum_error < 1e-10))
site_quantiles <- apply(score_draws, c(1, 2), quantile, probs = c(.025, .5, .975))
loading_quantiles <- apply(loading_draws, c(1, 2), quantile, probs = c(.025, .5, .975))
close(examples$quantiles$sites, site_quantiles)
close(examples$quantiles$loadings, loading_quantiles)
close(examples$ordinary_quantiles$sites,
      apply(original$U_output, c(1, 2), quantile, probs = c(.025, .5, .975)))
close(examples$ordinary_quantiles$loadings,
      apply(original$L_output, c(1, 2), quantile, probs = c(.025, .5, .975)))
stopifnot(identical(fit, saved$fit))
for (index in seq_len(nrow(examples$score_comparison))) {
  row <- examples$score_comparison[index, ]
  site <- match(row$site, sites)
  close(row$truth, factors$U[site, 1])
  close(unlist(row[c("lower", "estimate", "upper")]), site_quantiles[, site, 1])
}

# The native median labels, circle radii, truth crosses and connectors all
# get checked in actual built plot coordinates, not just source tables.
for (type in c("sites", "loadings")) {
  record <- examples$plots[[type]]
  labels <- as.character(record$layers[[2]]$label)
  indices <- match(labels, if (type == "sites") sites else species)
  stopifnot(!anyNA(indices))
  if (type == "sites") {
    stopifnot(identical(labels, examples$selected_sites))
    q <- site_quantiles[, indices, , drop = FALSE]
    centers <- q[2, , ]
    widths <- q[3, , ] - q[1, , ]
    true_xy <- factors$U[indices, , drop = FALSE]
  } else {
    stopifnot(identical(labels, species))
    q <- loading_quantiles[, , indices, drop = FALSE]
    centers <- t(q[2, , ])
    widths <- t(q[3, , ] - q[1, , ])
    true_xy <- t(factors$L[, indices, drop = FALSE])
  }
  close(record$layers[[2]]$x, centers[, 1])
  close(record$layers[[2]]$y, centers[, 2])
  close(record$data$r, sqrt(widths[, 1] * widths[, 2] / pi))
  close(record$layers[[3]]$x, centers[, 1])
  close(record$layers[[3]]$y, centers[, 2])
  close(record$layers[[3]]$xend, true_xy[, 1])
  close(record$layers[[3]]$yend, true_xy[, 2])
  close(record$layers[[4]]$x, true_xy[, 1])
  close(record$layers[[4]]$y, true_xy[, 2])
  stopifnot(all(record$layers[[4]]$shape == 4))
  circles <- record$layers[[1]]
  # Each native circle and truth overlay must share its labelled panel.
  label_panels <- as.character(record$layers[[2]]$PANEL)
  stopifnot(length(unique(label_panels)) == nrow(centers),
            identical(as.character(record$layers[[3]]$PANEL), label_panels),
            identical(as.character(record$layers[[4]]$PANEL), label_panels))
  for (index in seq_len(nrow(centers))) {
    circle <- circles[as.character(circles$PANEL) == label_panels[index], ]
    stopifnot(nrow(circle) > 300)
    close((circle$x - centers[index, 1])^2 + (circle$y - centers[index, 2])^2,
          rep(record$data$r[index]^2, nrow(circle)))
  }
}

biplot <- examples$plots$biplot
site_medians <- site_quantiles[2, , ]
loading_medians <- loading_quantiles[2, , ]
multiplier <- .8 * max(sqrt(rowSums(site_medians^2))) /
  max(sqrt(colSums(loading_medians^2)))
close(examples$arrow_multiplier, multiplier)
close(biplot$layers[[1]]$x, site_medians[, 1])
close(biplot$layers[[1]]$y, site_medians[, 2])
stopifnot(identical(as.character(biplot$layers[[3]]$label), species))
close(biplot$layers[[2]]$xend, loading_medians[1, ] * multiplier)
close(biplot$layers[[2]]$yend, loading_medians[2, ] * multiplier)
close(biplot$layers[[4]]$xend, factors$L[1, ] * multiplier)
close(biplot$layers[[4]]$yend, factors$L[2, ] * multiplier)
close(biplot$layers[[4]]$x, rep(0, length(species)))
close(biplot$layers[[4]]$y, rep(0, length(species)))
for (index in seq_len(nrow(examples$figures))) {
  figure <- examples$figures[index, ]
  path <- new_path(file.path("vignettes/teaching-data", figure$file))
  stopifnot(identical(md5(path), figure$md5))
  dimensions <- dim(png::readPNG(path, native = TRUE))
  stopifnot(identical(dimensions, as.integer(c(figure$height, figure$width) * 150)))
}
cat("PASS: identities, 24000 joint alignments, preserved contributions and Gram matrices,\n")
cat("marginal quantiles, exact teaching code, native plotted quantities and truth overlays.\n")
print(maximum_error)
