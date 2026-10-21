# Job driver for stage 2: attempts, selection, prediction and storage. One job is one dataset and arm.
D2 <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("upstream/adapter-contract.R", "upstream/basis.R", "upstream/adapter-occjsdm.R", "upstream/adapter-hmsc.R",
            "upstream/adapter-gllvm.R", "upstream/diagnostics.R", "upstream/select.R", "config.R", "data.R", "score.R")) source(file.path(D2, f))

stage2_prepare <- function(study, pilot_study, jobs_dir, config = stage2_config()) {
  dir.create(study, showWarnings = FALSE, recursive = TRUE)
  stopifnot(!file.exists(file.path(study, "config.rds")))
  dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)
  for (d in c("attempts", "predictions", "summaries")) dir.create(file.path(jobs_dir, d), showWarnings = FALSE)
  if (!file.exists(file.path(study, "jobs"))) file.symlink(normalizePath(jobs_dir), file.path(study, "jobs"))
  for (d in c("status", "scores", "logs")) dir.create(file.path(study, d), showWarnings = FALSE)
  datasets <- stage2_datasets(pilot_study, config)
  config$pilot_study <- pilot_study; config$created <- Sys.time(); config$git <- system("git rev-parse HEAD", intern = TRUE)
  config$r_version <- R.version.string
  config$packages <- sapply(c("occJSDM", "Hmsc", "gllvm", "posterior"), function(p) as.character(packageVersion(p)))
  config$sources <- tools::md5sum(list.files(D2, pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
  saveRDS(config, file.path(study, "config.rds"))
  write.csv(datasets, file.path(study, "datasets.csv"), row.names = FALSE)
  invisible(config)
}

job_table <- function(study) {
  cfg <- readRDS(file.path(study, "config.rds")); d <- read.csv(file.path(study, "datasets.csv"), stringsAsFactors = FALSE)
  jobs <- merge(d, cfg$arms, by = NULL)
  jobs$key <- arm_key(jobs$package, jobs$arm)
  jobs$job_id <- sprintf("%s--%s-%s", jobs$dataset_id, jobs$package, jobs$arm)
  jobs$kind <- ifelse(jobs$package == "gllvm", "point", "bayesian")
  jobs[order(jobs$job_id), ]
}

stage2_seed <- function(job_id, attempt_id, kind) {
  b <- as.numeric(utf8ToInt(paste0("stage2:", kind, ":", job_id)))
  as.integer((sum(b * seq_along(b)) * 7919 + sum(b) + 97 * attempt_id) %% 2147483646 + 1)
}

attempt_settings <- function(job, attempt_id, config, frozen, mcmc = config$mcmc) {
  sched <- if (job$kind == "bayesian" && attempt_id == 2L) mcmc$extended else mcmc$initial
  list(seed = stage2_seed(job$job_id, attempt_id, "fit"), prediction_seed = stage2_seed(job$job_id, 1L, "prediction"),
       job_id = job$job_id, attempt_id = attempt_id, source_hashes = config$sources, mcmc = sched, gllvm = config$gllvm,
       support_points = if (job$n == 900L) frozen$support_points_900 else config$support_points_300,
       hmsc_method = frozen$hmsc_method, hmsc_knots = frozen$hmsc_knots)
}

# Test-only slicing: `sites` is the number of leading sites to keep and `species` the species columns.
slice_training <- function(tr, sites, species) {
  if (is.null(sites)) return(tr)
  sites <- seq_len(sites)
  tr$site_id <- tr$site_id[sites]; tr$species_id <- tr$species_id[species]; tr$y <- tr$y[sites, species, drop = FALSE]
  tr$environment <- tr$environment[sites, , drop = FALSE]; tr$xy <- tr$xy[sites, , drop = FALSE]; tr
}

# Attempts are written atomically, so a completed attempt file is never partial; its small summary lives apart in
# jobs/summaries so that selection, scoring and timing never read a multi-gigabyte fit. A missing summary is rebuilt.
save_attempt <- function(study, job, attempt_id, record) {
  path <- file.path(study, "jobs", "attempts", sprintf("%s-%d.rds", job$job_id, attempt_id))
  tmp <- paste0(path, ".tmp"); saveRDS(record, tmp); file.rename(tmp, path)
  tmp <- paste0(summary_path(study, job, attempt_id), ".tmp"); saveRDS(record$summary, tmp); file.rename(tmp, summary_path(study, job, attempt_id))
  record$summary
}
summary_path <- function(study, job, attempt_id) file.path(study, "jobs", "summaries", sprintf("%s-%d.rds", job$job_id, attempt_id))

run_attempt <- function(study, job, attempt_id, config, frozen, mcmc, slice_sites = NULL, slice_species = NULL) {
  path <- file.path(study, "jobs", "attempts", sprintf("%s-%d.rds", job$job_id, attempt_id))
  if (file.exists(summary_path(study, job, attempt_id))) return(readRDS(summary_path(study, job, attempt_id)))
  if (file.exists(path)) { s <- readRDS(path)$summary; saveRDS(s, summary_path(study, job, attempt_id)); return(s) }
  cm <- readRDS(job$input_file); tr <- slice_training(load_training(cm, job$dataset_id, job$input_md5), slice_sites, slice_species)
  basis <- fit_basis(tr, job$arm); settings <- attempt_settings(job, attempt_id, config, frozen, mcmc)
  # A test slice has fewer sites than the production support-point settings; production never slices.
  if (!is.null(slice_sites)) settings$support_points <- min(settings$support_points, nrow(tr$y))
  t0 <- Sys.time()
  fit <- tryCatch(switch(job$package, occJSDM = fit_occjsdm(tr, basis, settings), Hmsc = fit_hmsc(tr, basis, settings), gllvm = fit_gllvm(tr, basis, settings)),
                  error = function(e) structure(list(message = conditionMessage(e)), class = "fit_error"))
  if (inherits(fit, "fit_error")) {
    summary <- list(attempt_id = attempt_id, status = "failed", acceptable = FALSE, error = fit$message, elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"))
    return(save_attempt(study, job, attempt_id, list(summary = summary, fit = NULL, diagnostics = NULL)))
  }
  diag <- diagnose_fit(fit, tr, basis, study_config())
  training_probability <- if (job$kind == "point") { p <- predict_gllvm(fit, c(tr[c("schema_version", "site_id", "species_id", "environment", "xy")]), basis); p$p_full } else NULL
  summary <- list(attempt_id = attempt_id, status = "complete", acceptable = !any(diag$failed), objective = if (job$kind == "point") fit$native$objective else NA_real_,
                  training_probability = training_probability, elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), warnings = fit$warnings)
  save_attempt(study, job, attempt_id, list(summary = summary, fit = fit, diagnostics = diag, basis = basis, training = tr))
}

select_job <- function(study, job, config) {
  files <- list.files(file.path(study, "jobs", "summaries"), pattern = paste0("^", job$job_id, "-[0-9]+\\.rds$"), full.names = TRUE)
  summaries <- lapply(files, readRDS)
  select_attempts(summaries, study_config(), kind = job$kind)
}

predict_job <- function(study, job, selected_attempt, frozen) {
  out <- file.path(study, "jobs", "predictions", paste0(job$job_id, ".rds"))
  if (file.exists(out)) return(readRDS(out))
  a <- readRDS(file.path(study, "jobs", "attempts", sprintf("%s-%d.rds", job$job_id, selected_attempt)))
  cm <- readRDS(job$input_file); pr <- load_prediction(cm); pr$species_id <- a$training$species_id
  thin <- if (is.na(frozen$prediction_thin)) 1L else frozen$prediction_thin
  fit <- a$fit
  if (job$kind == "bayesian" && thin > 1L) fit <- thin_fit(fit, thin)
  t0 <- Sys.time()
  prediction <- switch(job$package, occJSDM = predict_occjsdm(fit, pr, a$basis), Hmsc = predict_hmsc(fit, pr, a$basis), gllvm = predict_gllvm(fit, pr, a$basis))
  record <- list(job_id = job$job_id, selected_attempt = selected_attempt, prediction = prediction, prediction_thin = thin,
                 elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), finished = Sys.time())
  tmp <- paste0(out, ".tmp"); saveRDS(record, tmp); file.rename(tmp, out); record
}

# Keeps every thin-th retained draw of a Bayesian fit for prediction only; selection and diagnostics used all draws.
thin_fit <- function(fit, thin) {
  if (fit$package == "occJSDM") {
    js <- fit$native$results_output$jsdm_output; keep <- seq(1L, dim(js$B0_output)[2], by = thin)
    for (nm in c("B0_output", "B_output", "Bs_output", "sigmab_output", "sigmabs_output", "idx_ls_output")) if (!is.null(js[[nm]])) {
      d <- dim(js[[nm]]); it <- length(d) - 1L
      js[[nm]] <- if (it == 1L) js[[nm]][keep, , drop = FALSE] else if (it == 2L) js[[nm]][, keep, , drop = FALSE] else js[[nm]][, , keep, , drop = FALSE]
    }
    fit$native$results_output$jsdm_output <- js
  } else if (fit$package == "Hmsc") {
    fit$native$postList <- lapply(fit$native$postList, function(ch) ch[seq(1L, length(ch), by = thin)])
  }
  fit$prediction_thin <- thin; fit
}

run_job <- function(study, job, config = readRDS(file.path(study, "config.rds")), frozen = config$frozen, mcmc = config$mcmc,
                    slice_sites = NULL, slice_species = NULL) {
  repeat {
    sel <- select_job(study, job, config)
    if (sel$status == "pending") { for (id in sel$next_attempts) run_attempt(study, job, id, config, frozen, mcmc, slice_sites, slice_species); next }
    break
  }
  status <- switch(sel$status, selected = "pass", unresolved_diagnostics = "unresolved", unstable = "unstable", failed = "no fit", budget_amendment_required = "no fit")
  prediction <- if (!is.na(sel$selected_attempt)) tryCatch(predict_job(study, job, sel$selected_attempt, frozen),
                                                         error = function(e) { status <<- "unsupported"; list(error = conditionMessage(e)) }) else NULL
  result <- list(job_id = job$job_id, status = status, selection = sel, prediction_error = prediction$error, finished = Sys.time())
  writeLines(jsonlite::toJSON(result[c("job_id", "status")], auto_unbox = TRUE), file.path(study, "status", paste0(job$job_id, ".json")))
  result
}
