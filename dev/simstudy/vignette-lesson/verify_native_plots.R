# Independent numerical and graphical checks; no fitting or package changes.
# Run from the repository root with the archive's library first in R_LIBS:
# Rscript dev/simstudy/vignette-lesson/verify_native_plots.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("NATIVE_PLOT_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))
close <- function(a, b, tolerance = 1e-10) {
  stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)),
            max(abs(as.numeric(a) - as.numeric(b))) < tolerance)
}

library(occJSDM)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
x <- readRDS(new_path("vignettes/teaching-data/native-plots.rds"))
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
input <- lesson$input
truth <- input$sim$true_params
provenance <- x$provenance
snippet_path <- new_path("dev/simstudy/vignette-lesson/native-plot-examples.Rmd")
stopifnot(
  identical(provenance$source_hashes, lesson_source_hashes()),
  identical(provenance$lesson_md5, md5("vignettes/teaching-data/nonspatial-lesson.rds")),
  identical(provenance$output_md5, md5("vignettes/teaching-data/output-lesson.rds")),
  identical(provenance$snippet_md5, md5(snippet_path)),
  identical(provenance$exporter_md5, md5(new_path("dev/simstudy/vignette-lesson/export_native_plots.R"))),
  identical(provenance$input_md5, md5(file.path(archive, "input.rds"))),
  identical(provenance$fit_manifest, lesson$manifests$default),
  identical(provenance$package_files_md5, tools::md5sum(names(provenance$package_files_md5))),
  identical(normalizePath(find.package("occJSDM")), provenance$package_library)
)
fit_path <- file.path(archive, provenance$fit_manifest$file)
stopifnot(identical(md5(fit_path), provenance$fit_manifest$md5))
saved <- readRDS(fit_path)
f <- saved$fit
validate_lesson_fit_identity(f, input)
stopifnot(identical(saved$source_hashes, provenance$source_hashes),
          identical(saved$input_md5, provenance$input_md5),
          identical(dim(f$results_output$p_output)[3:4], c(6000L, 4L)))

# Re-extract the student bodies to detect a stale figure/code pairing.
rmd <- tempfile(fileext = ".Rmd")
r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), rmd)
knitr::purl(rmd, output = r, documentation = 0, quiet = TRUE)
stopifnot(identical(md5(r), provenance$student_code_md5))
unlink(c(rmd, r))

# When this source appendix is copied into Lesson 3, verify that its displayed
# plotting bodies still match those actually evaluated for the saved figures.
exportable_chunks <- function(lines) {
  starts <- grep("^```\\{r native-.*purl=TRUE", lines)
  setNames(lapply(starts, function(start) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }), sub("^```\\{r ([^,]+),.*", "\\1", lines[starts]))
}
expected_chunks <- exportable_chunks(readLines(snippet_path))
lesson_chunks <- exportable_chunks(readLines("vignettes/occJSDM-lesson-3.Rmd"))
stopifnot(length(expected_chunks) == 10L)
if (any(names(expected_chunks) %in% names(lesson_chunks))) {
  stopifnot(all(names(expected_chunks) %in% names(lesson_chunks)),
            identical(expected_chunks, lesson_chunks[names(expected_chunks)]))
  cat("Lesson 3 displays the exact exported native plotting bodies.\n")
}

sp <- f$infos$speciesNames
j <- f$results_output$jsdm_output
b <- truth$jsdmParams_true
close(f$X_psi %*% b$B + b$U %*% b$L + rep(b$B0, each = nrow(f$X_psi)), b$eta)
beta <- truth$beta_theta_true
beta[1, ] <- beta[1, ] + f$infos$list_X_theta_mat$mean_df * beta[2, ]
beta[2, ] <- beta[2, ] * f$infos$list_X_theta_mat$sd_df
samples <- f$infos$data_info[!duplicated(f$infos$data_info$Sample), ]
close(f$X_theta %*% beta, cbind(1, samples$X_theta) %*% truth$beta_theta_true)

# Test the displayed coordinates as well as values. This catches missing or
# reordered truth crosses on the native factor axes after ggplot builds them.
check_crosses <- function(record, table, layer) {
  built <- record$layers[[layer]]
  stopifnot(nrow(built) == nrow(table), all(is.finite(built$x)),
            all(is.finite(built$y)), all(as.character(built$PANEL) == "1"))
  labels <- record$x_limits[[1]][as.integer(built$x)]
  stopifnot(identical(as.character(labels), as.character(table$species)))
  close(built$y, table$truth)
}

for (name in c("environment", "collection")) {
  rec <- x$plots[[name]]
  expected_cov <- if (name == "environment") "X_psi.EnvCov.1" else "X_theta"
  table <- x$truth[[name]][x$truth[[name]]$covariate == expected_cov, ]
  check_crosses(rec, table, 3)
  for (i in seq_len(nrow(rec$data))) {
    s <- match(rec$data$Output[i], sp)
    draws <- if (name == "environment") j$B_output[1, s, , ] else
      f$results_output$beta_theta_output[2, s, , ]
    close(unlist(rec$data[i, c("2.5%", "97.5%")]), quantile(draws, c(.025, .975)))
    close(table$truth[table$species == sp[s]],
          if (name == "environment") b$B[1, s] else beta[2, s])
  }
}
for (name in c("occupancy_rates", "collection_rates")) {
  rec <- x$plots[[name]]
  table <- x$truth[[name]]
  check_crosses(rec, table, 2)
  for (i in seq_len(nrow(rec$data))) {
    s <- match(rec$data$Species[i], sp)
    draws <- if (name == "occupancy_rates") j$B0_output[s, , ] else
      f$results_output$beta_theta_output[1, s, , ]
    close(unlist(rec$data[i, c("Min", "Max")]), quantile(plogis(draws), c(.025, .975)))
    close(table$truth[table$species == sp[s]],
          plogis(if (name == "occupancy_rates") b$B0[s] else beta[1, s]))
  }
}
for (primer in 1:2) {
  rec <- x$plots[[paste0("primer_", primer)]]
  table <- x$truth$laboratory[x$truth$laboratory$Primer == as.character(primer), ]
  check_crosses(rec, table, 3)
  for (i in seq_len(nrow(rec$data))) {
    s <- match(rec$data$Species[i], sp)
    close(unlist(rec$data[i, c("p1", "p2")]),
          quantile(f$results_output$p_output[primer, s, , ], c(.025, .975)))
    close(unlist(rec$data[i, c("q1", "q2")]),
          quantile(f$results_output$q_output[primer, s, , ], c(.025, .975)))
    for (rate in c("p", "q")) {
      meanlog <- if (rate == "p") input$params$mu1 else input$params$mu0
      sdlog <- if (rate == "p") input$params$sigma1 else input$params$sigma0
      expected <- input$params[[rate]][primer, s] *
        plnorm(1.5, meanlog, sdlog, lower.tail = FALSE)
      close(table$truth[table$species == sp[s] & table$rate == rate], expected)
    }
  }
}
cat("All coefficient/rate intervals, standardized truth, species and primer axes verified.\n")

# Independently reconstruct every residual-correlation interval from loading
# draws and map BOTH native x markers and truth labels back to named pairs.
rec <- x$plots$correlations
tiles <- rec$layers[[1]]
markers <- rec$layers[[2]]
labels <- rec$layers[[3]]
stopifnot(nrow(tiles) == 45L, nrow(labels) == 45L, nrow(markers) == 45L)
pair_at <- function(x_coord, y_coord) {
  paste(rec$x_limits[[1]][as.integer(round(x_coord))],
        rec$y_limits[[1]][as.integer(round(y_coord))], sep = "/")
}
tile_pairs <- pair_at(tiles$x, tiles$y)
stopifnot(identical(tile_pairs, paste(rec$data$Var1, rec$data$Var2, sep = "/")))
marker_pairs <- pair_at(markers$x, markers$y)
label_pairs <- pair_at(labels$x, labels$y - .27)
truth_pairs <- paste(x$truth$correlations$species1, x$truth$correlations$species2, sep = "/")
stopifnot(identical(label_pairs, truth_pairs), !anyDuplicated(label_pairs),
          setequal(tile_pairs, label_pairs), all(labels$label == x$truth$correlations$truth_label))
for (i in seq_len(nrow(rec$data))) {
  s1 <- match(rec$data$Var1[i], sp)
  s2 <- match(rec$data$Var2[i], sp)
  l1 <- matrix(j$L_output[, s1, , ], nrow = 2)
  l2 <- matrix(j$L_output[, s2, , ], nrow = 2)
  draws <- colSums(l1 * l2) / sqrt(colSums(l1^2) * colSums(l2^2))
  interval <- quantile(draws, c(.025, .5, .975))
  close(rec$data$value[i], round(interval[2], 2))
  crosses <- interval[1] < 0 && interval[3] > 0
  stopifnot((tile_pairs[i] %in% marker_pairs) == crosses)
  a <- x$truth$correlations[match(tile_pairs[i], truth_pairs), ]
  norm <- sqrt(sum(b$L[, s1]^2) * sum(b$L[, s2]^2))
  if (norm == 0) stopifnot(is.na(a$truth), a$truth_label == "NA") else
    close(a$truth, sum(b$L[, s1] * b$L[, s2]) / norm)
}
cat("All 45 native correlation medians, uncertainty markers and truth labels match named pairs.\n")

# Independent truth calculation: average over the Binomial number of collected
# field samples, then enumerate all 2^10 possible species detection patterns.
patterns <- as.matrix(expand.grid(rep(list(0:1), length(sp))))
counts <- rowSums(patterns)
survey <- x$truth$survey_counts
p_truth <- input$params$p * plnorm(1.5, input$params$mu1, input$params$sigma1,
                                  lower.tail = FALSE)
for (i in seq_len(nrow(survey))) {
  a <- survey[i, ]
  probability <- vapply(seq_along(sp), function(s) {
    collected <- 0:a$M
    sum(dbinom(collected, a$M, plogis(beta[1, s])) *
          (1 - prod(1 - p_truth[, s])^(a$K * collected)))
  }, numeric(1))
  pattern_probability <- apply(patterns, 1, function(observed) {
    prod(ifelse(observed == 1, probability, 1 - probability))
  })
  mass <- vapply(0:length(sp), function(count) sum(pattern_probability[counts == count]), numeric(1))
  close(a$species_probability[[1]], probability)
  close(a$mass[[1]], mass)
  close(sum(mass), 1)
  close(a$expectation, sum((0:length(sp)) * mass))
  q <- vapply(c(.025, .5, .975), function(level) min(which(cumsum(mass) >= level)) - 1, numeric(1))
  close(unlist(a[c("lower", "median", "upper")]), q)
}
set.seed(provenance$seed)
recomputed <- plotCumulativeSpeciesDetections(f, M = 4, K = 6, primer = 0, byK = TRUE)
stopifnot(identical(as.data.frame(recomputed$data), x$plots$effort_k$data),
          identical(x$plots$effort_k$data, x$plots$effort_m$data))
for (direction in c("k", "m")) {
  rec <- x$plots[[paste0("effort_", direction)]]
  interval <- rec$layers[[2]]
  median <- rec$layers[[3]]
  axis <- if (direction == "k") "K" else "M"
  facet <- if (direction == "k") "M" else "K"
  stopifnot(nrow(interval) == 24L, nrow(median) == 24L,
            all(is.finite(interval$x)), all(is.finite(interval$ymin)),
            all(is.finite(interval$ymax)), all(is.finite(median$y)))
  close(interval$x, survey[[axis]] + .12)
  close(median$x, interval$x)
  close(interval$ymin, survey$lower)
  close(interval$ymax, survey$upper)
  close(median$y, survey$median)
  panels <- match(as.character(interval$PANEL), as.character(rec$layout$PANEL))
  close(rec$layout[[facet]][panels], survey[[facet]])
}
stopifnot(all(survey$K <= 6L), all(survey$M <= 4L))
cat("Exact survey-count truth verified by exhaustive enumeration; both native directions and all facets verified.\n")

for (i in seq_len(nrow(x$figures))) {
  file <- new_path(file.path("vignettes/teaching-data", x$figures$file[i]))
  stopifnot(identical(md5(file), x$figures$md5[i]), file.info(file)$size > 5000)
  signature <- readBin(file, what = "raw", n = 8)
  stopifnot(identical(as.integer(signature), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
}
stopifnot(file.info(new_path("vignettes/teaching-data/native-plots.rds"))$size < 1e6)
cat("Nine PNGs, compact export, full-fit hash, source hashes and exact teaching-code provenance verified.\n")
