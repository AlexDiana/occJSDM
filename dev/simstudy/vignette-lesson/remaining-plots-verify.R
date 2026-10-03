# Independent numerical and graphical checks, using the unchanged full fit.
# Run from repository root with the archive library first in R_LIBS.
# Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
stage <- Sys.getenv("REMAINING_PLOT_STAGE", ".")
new_path <- function(path) file.path(stage, path)
md5 <- function(path) unname(tools::md5sum(path))
close <- function(a, b, tolerance = 1e-10) {
  stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)),
            max(abs(as.numeric(a) - as.numeric(b))) < tolerance)
}
library(occJSDM)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
x <- readRDS(new_path("vignettes/teaching-data/remaining-plots-data.rds"))
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
input <- lesson$input
p <- x$provenance
snippet_path <- new_path("dev/simstudy/vignette-lesson/remaining-plots-examples.Rmd")
stopifnot(
  identical(p$source_hashes, lesson_source_hashes()),
  identical(p$lesson_md5, md5("vignettes/teaching-data/nonspatial-lesson.rds")),
  identical(p$output_md5, md5("vignettes/teaching-data/output-lesson.rds")),
  identical(p$snippet_md5, md5(snippet_path)),
  identical(p$exporter_md5, md5(new_path("dev/simstudy/vignette-lesson/remaining-plots-export.R"))),
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
for (name in p$api_names) {
  stopifnot(identical(body(get(name, asNamespace("occJSDM"))), body(get(name, source_environment))),
            identical(formals(get(name, asNamespace("occJSDM"))), formals(get(name, source_environment))))
}

# Prove that the displayed chunks are those evaluated by the exporter.
rmd <- tempfile(fileext = ".Rmd")
r <- tempfile(fileext = ".R")
writeLines(sub("eval=FALSE, purl=TRUE", "eval=TRUE, purl=TRUE",
               readLines(snippet_path), fixed = TRUE), rmd)
knitr::purl(rmd, output = r, documentation = 0, quiet = TRUE)
stopifnot(identical(md5(r), p$student_code_md5))
unlink(c(rmd, r))
exportable_chunks <- function(lines) {
  starts <- grep("^```\\{r remaining-plots-.*purl=TRUE", lines)
  setNames(lapply(starts, function(start) {
    end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
    lines[(start + 1L):(end - 1L)]
  }), sub("^```\\{r ([^,]+),.*", "\\1", lines[starts]))
}
expected_chunks <- exportable_chunks(readLines(snippet_path))
lesson_chunks <- exportable_chunks(readLines("vignettes/occJSDM-lesson-3.Rmd"))
stopifnot(length(expected_chunks) == 6L)
if (any(names(expected_chunks) %in% names(lesson_chunks))) {
  stopifnot(all(names(expected_chunks) %in% names(lesson_chunks)),
            identical(expected_chunks, lesson_chunks[names(expected_chunks)]))
  cat("Lesson 3 displays the exact exported remaining plotting bodies.\n")
}

sp <- f$infos$speciesNames
primers <- as.character(f$infos$primerNames)
covariates <- colnames(f$X_psi)
b <- input$sim$true_params$jsdmParams_true
j <- f$results_output$jsdm_output
info <- input$sim$data_list$info
raw <- as.matrix(info[!duplicated(info$Site), covariates])
close(x$scaling$mean, colMeans(raw))
close(x$scaling$sd, apply(raw, 2, sd))
stopifnot(identical(x$scaling$covariate, covariates))
close(f$X_psi, scale(raw))
close(f$X_psi %*% b$B + b$U %*% b$L + rep(b$B0, each = nrow(raw)), b$eta)
stopifnot(identical(sp, colnames(input$sim$data_list$OTU)), identical(primers, c("1", "2")))

curve_key <- function(species, value) paste(species, sprintf("%.12f", value), sep = "/")
for (k in seq_along(covariates)) {
  rec <- x$plots[[paste0("gradient_", k)]]
  truth <- x$truth$gradients[x$truth$gradients$covariate == covariates[k], ]
  grid <- seq(quantile(f$X_psi[, k], .02), quantile(f$X_psi[, k], .98), length.out = 40)
  close(unique(rec$data$x), grid)
  close(unique(truth$x), grid)
  close(truth$raw_x, x$scaling$mean[k] + x$scaling$sd[k] * truth$x)
  stopifnot(nrow(rec$data) == 400L, nrow(truth) == 400L,
            setequal(as.character(rec$data$species), sp))
  for (s in seq_along(sp)) {
    rows <- which(as.character(rec$data$species) == sp[s])
    trows <- which(truth$species == sp[s])
    # Independently form each linear predictor from all 24,000 matched draws.
    for (g in seq_along(grid)) {
      design <- apply(f$X_psi, 2, median)
      design[k] <- grid[g]
      eta <- as.vector(j$B0_output[s, , ])
      for (c in seq_along(covariates)) eta <- eta + design[c] * as.vector(j$B_output[c, s, , ])
      expected <- quantile(plogis(eta), c(.025, .5, .975))
      close(unlist(rec$data[rows[g], c("low", "med", "high")]), expected)
      close(truth$truth[trows[g]], plogis(b$B0[s] + sum(design * b$B[, s])))
    }
  }
  # Check every native ribbon/median and truth line against the actual panel label.
  data_key <- curve_key(rec$data$species, rec$data$x)
  truth_key <- curve_key(truth$species, truth$x)
  for (layer in c(1L, 2L, 4L)) {
    points <- rec$layers[[layer]]
    panel_species <- rec$layout$species[match(points$PANEL, rec$layout$PANEL)]
    keys <- curve_key(panel_species, points$x)
    stopifnot(nrow(points) == 400L, !anyDuplicated(keys),
              setequal(keys, data_key), all(is.finite(points$x)))
    if (layer == 1L) {
      close(points$ymin, rec$data$low[match(keys, data_key)])
      close(points$ymax, rec$data$high[match(keys, data_key)])
    } else if (layer == 2L) {
      close(points$y, rec$data$med[match(keys, data_key)])
    } else {
      close(points$y, truth$truth[match(keys, truth_key)])
      stopifnot(all(points$linetype == "dashed"), all(points$colour == "black"))
    }
  }
  rug <- rec$layers[[3]]
  stopifnot(nrow(rug) == nrow(raw) * length(sp))
  for (panel in unique(rug$PANEL)) close(sort(rug$x[rug$PANEL == panel]), sort(f$X_psi[, k]))
}
cat("Both gradients: 800 intervals, 800 truth points, all facets, rugs, and raw/standardized scaling verified.\n")

for (name in c("stage1_fp", "stage2_fp", "detection")) {
  rec <- x$plots[[name]]
  bars <- rec$layers[[1]]
  crosses <- rec$layers[[2]]
  truth <- x$truth[[name]]
  stopifnot(all(is.finite(bars$x)), all(is.finite(bars$ymin)), all(is.finite(bars$ymax)),
            all(is.finite(crosses$x)), all(is.finite(crosses$y)), all(crosses$shape == 4))
  bar_species <- as.character(rec$x_limits[[1]][round(bars$x)])
  cross_species <- as.character(rec$x_limits[[1]][round(crosses$x)])
  if (name == "stage1_fp") {
    stopifnot(nrow(bars) == 10L, nrow(crosses) == 10L, setequal(bar_species, sp))
    close(crosses$y, input$params$theta0[match(cross_species, sp)])
    close(truth$truth, input$params$theta0[match(truth$Species, sp)])
    for (i in seq_len(nrow(bars))) {
      s <- match(bar_species[i], sp)
      close(c(bars$ymin[i], bars$ymax[i]), quantile(f$results_output$theta0_output[s, , ], c(.025, .975)))
    }
  } else {
    stopifnot(nrow(bars) == 20L, nrow(crosses) == 20L)
    bar_primer <- names(rec$primer_colours)[match(bars$colour, rec$primer_colours)]
    # For two primers, native dodge puts primer 1 left and primer 2 right.
    cross_primer <- primers[ifelse(crosses$x < round(crosses$x), 1, 2)]
    bar_key <- paste(bar_species, bar_primer)
    cross_key <- paste(cross_species, cross_primer)
    stopifnot(!anyNA(bar_primer), !anyDuplicated(bar_key), !anyDuplicated(cross_key),
              setequal(bar_key, cross_key))
    close(crosses$x, bars$x[match(cross_key, bar_key)])
    rate <- if (name == "stage2_fp") "q" else "p"
    mu <- if (rate == "q") input$params$mu0 else input$params$mu1
    sigma <- if (rate == "q") input$params$sigma0 else input$params$sigma1
    expected_truth <- input$params[[rate]] * plnorm(1.5, mu, sigma, lower.tail = FALSE)
    for (i in seq_len(nrow(bars))) {
      s <- match(bar_species[i], sp)
      primer <- match(bar_primer[i], primers)
      draws <- f$results_output[[paste0(rate, "_output")]][primer, s, , ]
      close(c(bars$ymin[i], bars$ymax[i]), quantile(draws, c(.025, .975)))
      close(crosses$y[match(bar_key[i], cross_key)], expected_truth[primer, s])
      trow <- truth$Species == sp[s] & as.character(truth$Primer) == primers[primer]
      close(truth$truth[trow], expected_truth[primer, s])
    }
  }
}
cat("All 50 rate intervals and truth crosses match species, primers, dodge coordinates, and threshold targets.\n")
for (i in seq_len(nrow(x$figures))) {
  file <- new_path(file.path("vignettes/teaching-data", x$figures$file[i]))
  stopifnot(identical(md5(file), x$figures$md5[i]), file.info(file)$size > 5000,
            identical(as.integer(readBin(file, "raw", n = 8)), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
}
stopifnot(file.info(new_path("vignettes/teaching-data/remaining-plots-data.rds"))$size < 1e6)
cat("Five PNGs, compact bundle, full-fit/library/source hashes, and exact displayed code verified.\n")
