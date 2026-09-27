# Run from the repository root with the matching archived library in R_LIBS.
# Rscript dev/simstudy/vignette-lesson/unbalanced-build.R OUTDIR prepare|fit|summarise
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, args[2] %in% c("prepare", "fit", "summarise"))
outdir <- normalizePath(args[1], mustWork = TRUE)
action <- args[2]
library(dplyr)
source("dev/simstudy/vignette-lesson/helpers.R")
source("dev/simstudy/vignette-lesson/score_lesson.R")
lesson_path <- "vignettes/teaching-data/nonspatial-lesson.rds"
lesson <- readRDS(lesson_path)
original <- lesson$input
stopifnot(identical(original$source_hashes, lesson_source_hashes()))
input_path <- file.path(outdir, "input.rds")
fit_path <- file.path(outdir, "unbalanced-fit.rds")

if (action == "prepare") {
  stopifnot(!file.exists(input_path), !file.exists(fit_path))
  survey_data <- original$sim$data_list
  sample_keys <- survey_data$info |>
    distinct(Site, Sample) |>
    arrange(Site, Sample)
  # Declare this random deletion before inspecting any fitted results.
  set.seed(3947)
  removed_sites <- sample(sort(unique(sample_keys$Site)), size = 3)
  removed_samples <- sample_keys |>
    filter(Site %in% removed_sites) |>
    group_by(Site) |>
    slice_sample(n = 1) |>
    ungroup() |>
    arrange(Site)
  # One explicit key removes all 12 PCR rows belonging to a field sample.
  retained_rows <- survey_data$info |>
    mutate(original_row = row_number()) |>
    anti_join(removed_samples, by = c("Site", "Sample")) |>
    pull(original_row)
  unbalanced_data <- survey_data
  unbalanced_data$info <- survey_data$info[retained_rows, , drop = FALSE]
  unbalanced_data$OTU <- survey_data$OTU[retained_rows, , drop = FALSE]
  samples_per_site <- unbalanced_data$info |>
    distinct(Site, Sample) |>
    count(Site, name = "samples")
  pcrs_per_primer <- unbalanced_data$info |>
    count(Site, Sample, Primer, name = "PCRs")
  stopifnot(nrow(removed_samples) == 3L,
            n_distinct(removed_samples$Site) == 3L,
            length(retained_rows) == 2364L,
            nrow(survey_data$info) - length(retained_rows) == 36L,
            n_distinct(unbalanced_data$info$Sample) == 197L,
            nrow(samples_per_site) == 100L,
            sum(samples_per_site$samples == 1L) == 3L,
            sum(samples_per_site$samples == 2L) == 97L,
            nrow(pcrs_per_primer) == 394L,
            all(pcrs_per_primer$PCRs == 6L),
            identical(sort(unique(unbalanced_data$info$Primer)), c(1L, 2L)),
            identical(colnames(unbalanced_data$OTU), colnames(survey_data$OTU)),
            ncol(unbalanced_data$OTU) == 10L,
            identical(unbalanced_data$traits, survey_data$traits),
            !anyNA(unbalanced_data$OTU))
  # Keep the original simulation settings and ALL original latent truths.
  # Only the observed tables change. No site truth is regenerated or discarded.
  input <- original
  input$sim$data_list <- unbalanced_data
  input$removal <- list(
    seed = 3947L, fit_seed = 20260924L, rng_kind = RNGkind(),
    removed_samples = removed_samples, retained_rows = retained_rows,
    samples_per_site = samples_per_site, pcrs_per_primer = pcrs_per_primer,
    original_lesson_md5 = unname(tools::md5sum(lesson_path)),
    selection = "Three uniformly sampled distinct sites; one uniformly sampled whole field sample at each, in numeric site order"
  )
  stopifnot(identical(input$sim$true_params, original$sim$true_params))
  saveRDS(input, input_path, compress = "xz")
  write.csv(removed_samples, file.path(outdir, "removed-samples-before-fitting.csv"),
            row.names = FALSE)
  print(removed_samples)
  print(count(samples_per_site, samples, name = "sites"))
  cat("Retained 197 samples and 2364 PCR rows at all 100 sites, for all 10 species.\n")
}

if (action == "fit") {
  stopifnot(!file.exists(fit_path))
  input <- readRDS(input_path)
  stopifnot(identical(input$source_hashes, lesson_source_hashes()),
            identical(input$sim$true_params, original$sim$true_params))
  dat <- input$sim$data_list
  mcmc <- list(nchain = 4L, nburn = 3000L, niter = 6000L, nthin = 1L)
  priors <- list()
  stopifnot(identical(mcmc, lesson$manifests$default$mcmc),
            identical(priors, lesson$manifests$default$priors))
  fit_args <- list(
    data = dat, listParams = list(n_factors = 2L, n_lattrait = 1L),
    threshold = 1, occCovariates = grep("^X_psi", names(dat$info), value = TRUE),
    collCovariates = "X_theta", spatCovariates = NULL,
    MCMCparams = mcmc, listPriors = priors, summarisedLatentPresences = TRUE
  )
  warnings <- character()
  set.seed(20260924)
  started <- Sys.time()
  fit <- withCallingHandlers(
    do.call(occJSDM::runOccJSDM, fit_args),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))
  )
  result <- list(
    arm = "unbalanced", fit = fit, seed = 20260924L, mcmc = mcmc,
    priors = priors, warnings = warnings, started = started,
    finished = Sys.time(), session = sessionInfo(),
    source_hashes = lesson_source_hashes(),
    input_md5 = unname(tools::md5sum(input_path))
  )
  # Save before validation so a diagnostic issue never loses the expensive fit.
  saveRDS(result, fit_path, compress = "gzip")
  validate_lesson_fit_identity(fit, input)
  cat("Saved the single approved fit in",
      as.numeric(difftime(result$finished, started, units = "secs")), "seconds.\n")
}

if (action == "summarise") {
  input <- readRDS(input_path)
  result <- readRDS(fit_path)
  stopifnot(identical(result$source_hashes, lesson_source_hashes()),
            identical(result$input_md5, unname(tools::md5sum(input_path))),
            identical(input$sim$true_params, original$sim$true_params))
  score <- score_lesson_fit(result, input)
  saveRDS(score, file.path(outdir, "unbalanced-summary.rds"), compress = "xz")
  # The scorer matches original global Sample IDs to w_true. Unequal sample
  # counts need no adaptation of the shared helper or truth reindexing.
  compact <- c(list(
    schema = 1L, removal = input$removal,
    original_lesson_md5 = unname(tools::md5sum(lesson_path)),
    input_md5 = result$input_md5, source_hashes = result$source_hashes,
    source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    fit_manifest = list(
      file = basename(fit_path), md5 = unname(tools::md5sum(fit_path)),
      seed = result$seed, mcmc = result$mcmc, priors = result$priors,
      warnings = result$warnings, session = result$session,
      seconds = as.numeric(difftime(result$finished, result$started, units = "secs"))
    ),
    exporter_hashes = tools::md5sum(c(
      "dev/simstudy/vignette-lesson/unbalanced-build.R",
      "dev/simstudy/vignette-lesson/helpers.R",
      "dev/simstudy/vignette-lesson/score_lesson.R"
    ))
  ), score)
  saveRDS(compact, "vignettes/teaching-data/unbalanced-lesson.rds", compress = "xz")
  print(score$groups)
  print(score$diagnostics |>
          filter(is.na(rhat) | is.na(ess) | rhat > 1.01 | ess < 400))
  cat("Fit warnings:", paste(result$warnings, collapse = " | "), "\n")
  cat("Maximum reconstructed psi difference:", score$reconstruction_difference, "\n")
}
