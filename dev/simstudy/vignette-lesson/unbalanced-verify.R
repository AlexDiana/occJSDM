# Independent checks of the reduced observations and saved fit, without MCMC.
# Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R FULL_FIT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
archive <- normalizePath(args[1], mustWork = TRUE)
library(dplyr)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
lesson <- readRDS("vignettes/teaching-data/nonspatial-lesson.rds")
compact <- readRDS("vignettes/teaching-data/unbalanced-lesson.rds")
input <- readRDS(file.path(archive, "input.rds"))
saved <- readRDS(file.path(archive, compact$fit_manifest$file))
fit <- saved$fit
original <- lesson$input$sim$data_list
dat <- input$sim$data_list
stopifnot(
  identical(input$sim$true_params, lesson$input$sim$true_params),
  identical(input$source_hashes, lesson_source_hashes()),
  identical(saved$source_hashes, lesson_source_hashes()),
  identical(compact$source_hashes, lesson_source_hashes()),
  identical(compact$original_lesson_md5,
            unname(tools::md5sum("vignettes/teaching-data/nonspatial-lesson.rds"))),
  identical(compact$exporter_hashes, tools::md5sum(names(compact$exporter_hashes))),
  identical(compact$input_md5, unname(tools::md5sum(file.path(archive, "input.rds")))),
  identical(saved$input_md5, compact$input_md5),
  identical(compact$fit_manifest$md5,
            unname(tools::md5sum(file.path(archive, compact$fit_manifest$file))))
)

# An independent base-R composite-key selection must give exactly the same rows.
key <- function(x) paste(x$Site, x$Sample, sep = "/")
removed <- compact$removal$removed_samples
keep <- !(key(original$info) %in% key(removed))
stopifnot(identical(which(keep), compact$removal$retained_rows),
          identical(dat$info, original$info[keep, , drop = FALSE]),
          identical(dat$OTU, original$OTU[keep, , drop = FALSE]),
          identical(dat$traits, original$traits),
          sum(!keep) == 36L, sum(keep) == 2364L,
          identical(as.numeric(removed$Site), c(12, 31, 52)),
          identical(as.numeric(removed$Sample), c(24, 61, 103)),
          nrow(unique(dat$info[c("Site", "Sample")])) == 197L,
          identical(sort(unique(dat$info$Site)), sort(unique(original$info$Site))),
          identical(colnames(dat$OTU), colnames(original$OTU)),
          all(table(interaction(dat$info$Site, dat$info$Sample,
                                dat$info$Primer, drop = TRUE)) == 6L),
          !anyNA(dat$OTU))
validate_lesson_fit_identity(fit, input)
stopifnot(sum(fit$infos$M == 1L) == 3L, sum(fit$infos$M == 2L) == 97L,
          length(fit$infos$K) == 394L, all(fit$infos$K == 6L),
          identical(as.character(fit$infos$siteNames),
                    rownames(input$sim$true_params$z_true)),
          fit$infos$ps == 0L, fit$infos$n_factors == 2L)

# Match all 1,000 probabilities to all unchanged site/species truths by identity.
cells <- compact$cells
truth <- input$sim$true_params
truth_index <- cbind(match(cells$Site, rownames(truth$z_true)),
                     match(cells$species, colnames(truth$z_true)))
stopifnot(!anyNA(truth_index), nrow(cells) == 1000L,
          identical(cells$truth, plogis(truth$jsdmParams_true$eta[truth_index])),
          identical(cells$z, truth$z_true[truth_index]),
          max(abs(cells$estimate - c(fit$results_output$psi_output))) < 1e-10,
          identical(cells$conditional, c(fit$results_output$z_output)))

# Reconstruct every retained draw for OTU_1 independently of the shared scorer.
j <- fit$results_output$jsdm_output
species <- match("OTU_1", fit$infos$speciesNames)
iterations <- dim(j$B0_output)[2]
chains <- dim(j$B0_output)[3]
probabilities <- array(NA_real_, c(100L, iterations, chains))
for (chain in seq_len(chains)) {
  for (iteration in seq_len(iterations)) {
    eta <- rep(j$B0_output[species, iteration, chain], 100L)
    for (covariate in seq_len(ncol(fit$X_psi))) {
      eta <- eta + fit$X_psi[, covariate] *
        j$B_output[covariate, species, iteration, chain]
    }
    for (factor in seq_len(dim(j$U_output)[2])) {
      eta <- eta + j$U_output[, factor, iteration, chain] *
        j$L_output[factor, species, iteration, chain]
    }
    probabilities[, iteration, chain] <- plogis(eta)
  }
}
cell_index <- match(paste(fit$infos$siteNames, "OTU_1"),
                    paste(cells$Site, cells$species))
stopifnot(max(abs(apply(probabilities, 1, mean) - cells$estimate[cell_index])) < 1e-12,
          max(abs(apply(probabilities, 1, quantile, probs = 0.025) -
                    cells$lower[cell_index])) < 1e-12,
          max(abs(apply(probabilities, 1, quantile, probs = 0.975) -
                    cells$upper[cell_index])) < 1e-12)
for (band in compact$groups$band) {
  selected <- switch(band, All = rep(TRUE, nrow(cells)),
                     Low = cells$truth < 0.2, High = cells$truth > 0.8,
                     Middle = cells$truth >= 0.2 & cells$truth <= 0.8)
  expected <- compact$groups[compact$groups$band == band, ]
  errors <- cells$estimate[selected] - cells$truth[selected]
  stopifnot(abs(mean(errors) - expected$signed_error) < 1e-12,
            abs(mean(abs(errors)) - expected$mae) < 1e-12)
}
stopifnot(isTRUE(all.equal(compact$diagnostics,
                          as.data.frame(occJSDM::returnConvergenceDiagnostics(fit)))),
          identical(compact$warnings, saved$warnings),
          nrow(compact$samples) == 1970L,
          !any(key(compact$samples) %in% key(removed)))
samples <- compact$samples
sample_truth_index <- cbind(match(as.character(samples$Sample), rownames(truth$w_true)),
                            match(samples$species, colnames(truth$w_true)))
site_truth_index <- cbind(match(as.character(samples$Site), rownames(truth$z_true)),
                          match(samples$species, colnames(truth$z_true)))
fitted_samples <- unique(fit$infos$data_info[c("Site", "Sample")])
sample_fit_index <- cbind(match(key(samples), key(fitted_samples)),
                          match(samples$species, fit$infos$speciesNames))
stopifnot(!anyNA(c(sample_truth_index, site_truth_index, sample_fit_index)),
          identical(samples$w, truth$w_true[sample_truth_index]),
          identical(samples$z, truth$z_true[site_truth_index]),
          identical(samples$sample_probability,
                    fit$results_output$w_output[sample_fit_index]))

# Execute exactly the displayed preparation chunks, preserving all identities.
chunk <- function(path, label) {
  lines <- readLines(path, warn = FALSE)
  start <- which(startsWith(lines, paste0("```{r ", label)))
  stopifnot(length(start) == 1L)
  end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
  parse(text = lines[(start + 1L):(end - 1L)])
}
displayed <- new.env(parent = globalenv())
displayed$survey_data <- original
displayed$known_truth <- truth
for (label in c("unbalanced-select-samples", "unbalanced-paired-rows")) {
  code <- chunk("vignettes/occJSDM-lesson-0.Rmd", label)
  stopifnot(identical(code, chunk("dev/simstudy/vignette-lesson/unbalanced-lesson-0.Rmd", label)))
  eval(code, displayed)
}
stopifnot(identical(displayed$unbalanced_data, dat),
          identical(displayed$removed_samples, removed),
          identical(displayed$known_truth, truth))

# Replace only the model call with a capturing function. Never run a second fit.
displayed$capture_fit_call <- function(...) list(...)
load_code <- chunk("vignettes/occJSDM-lesson-1.Rmd", "unbalanced-load")
stopifnot(identical(load_code, chunk("dev/simstudy/vignette-lesson/unbalanced-lesson-1.Rmd",
                                    "unbalanced-load")))
previous_directory <- setwd("vignettes")
tryCatch(eval(load_code, displayed), finally = setwd(previous_directory))
stopifnot(identical(displayed$unbalanced_data, dat))
optional <- chunk("vignettes/occJSDM-lesson-1.Rmd", "unbalanced-fit-optional")
stopifnot(identical(optional, chunk("dev/simstudy/vignette-lesson/unbalanced-lesson-1.Rmd",
                                   "unbalanced-fit-optional")),
          optional[[1]][[2]] == saved$seed)
stopifnot(identical(optional[[2]][[3]][[1]], quote(occJSDM::runOccJSDM)))
optional[[2]][[3]][[1]] <- quote(capture_fit_call)
eval(optional, displayed)
called <- displayed$unbalanced_fit
stopifnot(identical(called$data, dat),
          identical(called$MCMCparams, saved$mcmc),
          identical(called$listPriors, saved$priors),
          identical(called$listParams, list(n_factors = 2L, n_lattrait = 1L)),
          identical(called$occCovariates, colnames(fit$X_psi)),
          identical(called$collCovariates, "X_theta"),
          is.null(called$spatCovariates), called$threshold == 1,
          isTRUE(called$summarisedLatentPresences), saved$seed == 20260924L)
cat("Verified paired deletion, all unchanged truths, fitted identities, all saved psi means,\n",
    "OTU_1 draw reconstruction and intervals, errors, diagnostics, and displayed code.\n")
