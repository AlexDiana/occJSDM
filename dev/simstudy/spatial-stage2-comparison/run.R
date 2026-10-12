# Study driver for stage 2. Always through run-r:
#   run-r -f run.R --args --mode=prepare     --study=DIR --pilot=PILOT_DIR --jobs=VOLUME_DIR
#   run-r -f run.R --args --mode=microtiming --study=DIR             (minutes; micro fits are discarded)
#   run-r -f run.R --args --mode=timing      --study=DIR --step=1    (ladder step 1, 2, ...)
#   run-r -f run.R --args --mode=freeze      --study=DIR --step=1    (freeze the settings of a timing step that passed)
#   run-r -f run.R --args --mode=fit         --study=DIR --workers=8
#   run-r -f run.R --args --mode=score       --study=DIR
source("dev/simstudy/spatial-stage2-comparison/jobs.R")

# Ladder steps are ordered by total reduction, then by the Hmsc method: step 1 is 900 support points with Hmsc Full,
# step 2 is 600 with Full, step 3 is 900 with GPP at 300 knots, and so on down.
ladder_settings <- function(config, step) {
  sp <- config$ladder$support_points_900; hm <- config$ladder$hmsc
  grid <- expand.grid(sp = seq_along(sp), hm = seq_along(hm))
  grid <- grid[order(grid$sp + grid$hm, grid$hm), ]
  stopifnot(step >= 1L, step <= nrow(grid))
  g <- grid[step, ]
  list(support_points_900 = sp[g$sp], hmsc_method = hm[[g$hm]][1], hmsc_knots = as.integer(hm[[g$hm]][2]), prediction_thin = 10L)
}

# Upper-bound price of one job from its measured attempts: Bayesian arms at the initial attempt's cost times the
# extension factor (initial plus extension), point arms at the mean start cost times the maximum number of starts;
# prediction at the extended draw count (an unextended fit's prediction is scaled up by the draw ratio).
job_upper_seconds <- function(cfg, kind, attempt_seconds, selected_attempt, prediction_seconds) {
  draw_ratio <- unname(cfg$mcmc$extended[["niter"]] / cfg$mcmc$initial[["niter"]])
  if (kind == "bayesian") {
    fit <- attempt_seconds[1] * cfg$extension_factor
    pred <- prediction_seconds * (if (identical(as.integer(selected_attempt), 1L)) draw_ratio else 1)
  } else {
    fit <- mean(attempt_seconds) * cfg$diagnostics$max_starts; pred <- prediction_seconds
  }
  fit + pred
}

read_summaries <- function(study, job_id) {
  files <- list.files(file.path(study, "jobs", "summaries"), pattern = paste0("^", job_id, "-[0-9]+\\.rds$"), full.names = TRUE)
  files <- files[order(as.integer(sub(".*-([0-9]+)\\.rds$", "\\1", files)))]
  lapply(files, readRDS)
}

# Settings that a timing measurement depends on: a step only re-measures the arms whose setting changed.
timing_signature <- function(key, frozen) {
  if (key == "occJSDM/gp") paste0("support", frozen$support_points_900)
  else if (startsWith(key, "Hmsc/")) paste0("hmsc-", frozen$hmsc_method, "-", frozen$hmsc_knots)
  else "fixed"
}
timing_cols <- c("job_id", "key", "truth", "kind", "signature", "step_measured", "status", "attempts", "fit_seconds", "initial_seconds", "selected_attempt", "prediction_seconds", "error")

# One ladder step. A Bayesian timing fit is the initial attempt plus prediction on it: the spec prices the all-extension
# projection from the initial cost times 4.43, so the extension is priced and never run. Point arms run their starts.
# Arms whose setting is unchanged since an earlier step reuse that step's measurement.
stage2_timing <- function(study, step = 1L, workers = 5L, mcmc = NULL, slice_sites = NULL, slice_species = NULL) {
  cfg <- readRDS(file.path(study, "config.rds")); frozen <- ladder_settings(cfg, step); if (is.null(mcmc)) mcmc <- cfg$mcmc
  jobs <- job_table(study); pick <- jobs[jobs$n == 900L & jobs$replicate == 1L, ]   # 2 truths x 5 arms = 10 jobs
  pick$signature <- vapply(pick$key, timing_signature, "", frozen = frozen)
  usable_status <- c("initial", "pass", "unresolved", "unstable")
  earlier <- do.call(rbind, lapply(seq_len(step - 1L), function(k) { f <- file.path(study, sprintf("timing-step%d.json", k))
    if (file.exists(f)) { t <- jsonlite::fromJSON(f)$timing; t[t$status %in% usable_status & is.finite(t$prediction_seconds), ] } else NULL }))
  reused <- function(i) { if (is.null(earlier)) return(NULL); m <- earlier[earlier$key == pick$key[i] & earlier$truth == pick$truth[i] & earlier$signature == pick$signature[i], ]
    if (nrow(m)) m[nrow(m), timing_cols] else NULL }
  pick$job_id <- paste0(pick$job_id, "--timing", step)
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(pick)), function(i) tryCatch({
    old <- reused(i); if (!is.null(old)) return(old)
    job <- pick[i, ]
    if (job$kind == "bayesian") {
      s <- run_attempt(study, job, 1L, cfg, frozen, mcmc, slice_sites, slice_species)
      rec <- if (identical(s$status, "complete")) predict_job(study, job, 1L, frozen) else NULL
      status <- if (is.null(rec)) "failed" else "initial"; secs <- s$elapsed_seconds; sel <- 1L; pred <- if (is.null(rec)) NA_real_ else rec$elapsed_seconds; n_att <- 1L
    } else {
      r <- run_job(study, job, config = cfg, frozen = frozen, mcmc = mcmc, slice_sites = slice_sites, slice_species = slice_species)
      secs <- read_attempt_seconds(study, job$job_id); n_att <- length(secs); status <- r$status
      sel <- if (is.na(r$selection$selected_attempt)) NA_integer_ else as.integer(r$selection$selected_attempt)
      p <- file.path(study, "jobs", "predictions", paste0(job$job_id, ".rds")); pred <- if (file.exists(p)) readRDS(p)$elapsed_seconds else NA_real_
    }
    data.frame(job_id = job$job_id, key = job$key, truth = job$truth, kind = job$kind, signature = job$signature, step_measured = step, status = status,
               attempts = n_att, fit_seconds = sum(secs), initial_seconds = secs[1], selected_attempt = sel, prediction_seconds = pred, error = NA_character_, stringsAsFactors = FALSE)
  }, error = function(e) data.frame(job_id = pick$job_id[i], key = pick$key[i], truth = pick$truth[i], kind = pick$kind[i], signature = pick$signature[i], step_measured = step,
               status = "error", attempts = NA_integer_, fit_seconds = NA_real_, initial_seconds = NA_real_, selected_attempt = NA_integer_, prediction_seconds = NA_real_,
               error = conditionMessage(e), stringsAsFactors = FALSE)), mc.cores = workers, mc.preschedule = FALSE)
  timing <- do.call(rbind, res)
  usable <- timing[timing$status %in% usable_status & is.finite(timing$prediction_seconds), ]
  complete <- all(arm_key(cfg$arms$package, cfg$arms$arm) %in% usable$key)
  per_job <- if (nrow(usable)) data.frame(usable, job_upper = vapply(seq_len(nrow(usable)), function(i) job_upper_seconds(cfg, usable$kind[i],
    read_attempt_seconds(study, usable$job_id[i]), usable$selected_attempt[i], usable$prediction_seconds[i]), numeric(1))) else NULL
  per_arm <- if (complete) aggregate(job_upper ~ key, per_job, max) else NULL
  n_jobs <- table(jobs$key)
  cpu_hours <- if (complete) sum(per_arm$job_upper * as.numeric(n_jobs[per_arm$key])) / 3600 else NA_real_
  projection <- list(step = step, settings = frozen, timing = timing, per_arm = per_arm, complete = complete,
                     projected_cpu_hours_upper = cpu_hours, projected_wall_hours_upper = cpu_hours / cfg$workers, limit_wall_hours = cfg$wall_hours,
                     passes = complete && cpu_hours / cfg$workers < cfg$wall_hours, wall_seconds = as.numeric(Sys.time() - t0, units = "secs"),
                     note = "Upper bound: every job priced at its arm's 900-site cost, the larger of the two truths; Bayesian arms at 4.43 times the measured initial attempt; point arms at the maximum starts; prediction at the extended draw count. Arms whose setting is unchanged since an earlier step reuse its measurement (step_measured).")
  writeLines(jsonlite::toJSON(projection, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"), file.path(study, sprintf("timing-step%d.json", step)))
  projection
}

read_attempt_seconds <- function(study, job_id) vapply(read_summaries(study, job_id), function(x) x$elapsed_seconds, numeric(1))

# Micro-timing: one chain, 50 warm-up and 150 retained draws (gllvm: one start, maxit 200) at ladder step 1, on the two 900-site
# replicate-1 datasets. Nothing here is a study attempt. Each micro fit saves its own result under microtiming/ so that one stuck or
# killed fit never loses the others, and a rerun skips what exists. Extrapolation assumptions are written into the JSON.
micro_schedule <- c(nchain = 1L, nburn = 50L, niter = 150L, nthin = 1L)
micro_file <- function(study, job, step = 1L) file.path(study, "microtiming", paste0(job$job_id, if (step == 1L) "" else paste0("--step", step), ".rds"))

stage2_microtime_job <- function(study, job, cfg, frozen, step = 1L) {
  cm <- readRDS(job$input_file); tr <- load_training(cm, job$dataset_id, job$input_md5); pr <- load_prediction(cm)
  basis <- fit_basis(tr, job$arm); settings <- attempt_settings(job, 1L, cfg, frozen, list(initial = micro_schedule, extended = micro_schedule))
  if (job$kind == "point") settings$gllvm$maxit <- 200L
  t_fit <- system.time(fit <- switch(job$package, occJSDM = fit_occjsdm(tr, basis, settings), Hmsc = fit_hmsc(tr, basis, settings), gllvm = fit_gllvm(tr, basis, settings)))[["elapsed"]]
  t_diag <- system.time(suppressWarnings(diagnose_fit(fit, tr, basis, study_config())))[["elapsed"]]
  block <- subset_prediction(pr, seq_len(cfg$prediction_block))
  t_pred <- system.time(switch(job$package, occJSDM = predict_occjsdm(fit, block, basis), Hmsc = predict_hmsc(fit, block, basis), gllvm = predict_gllvm(fit, block, basis)))[["elapsed"]]
  data.frame(job_id = job$job_id, key = job$key, truth = job$truth, kind = job$kind, step = step, signature = timing_signature(job$key, frozen), fit_seconds = t_fit, diagnose_seconds = t_diag, block_prediction_seconds = t_pred,
             gllvm_converged = if (job$kind == "point") isTRUE(fit$native$convergence) else NA, finished = format(Sys.time()), stringsAsFactors = FALSE)
}

stage2_microtiming <- function(study, workers = 5L, keys = NULL, step = 1L) {
  cfg <- readRDS(file.path(study, "config.rds")); frozen <- ladder_settings(cfg, step)
  dir.create(file.path(study, "microtiming"), showWarnings = FALSE)
  jobs <- job_table(study); pick <- jobs[jobs$n == 900L & jobs$replicate == 1L, ]
  if (!is.null(keys)) pick <- pick[pick$key %in% keys, ]
  pick <- pick[!file.exists(vapply(seq_len(nrow(pick)), function(i) micro_file(study, pick[i, ], step), "")), ]
  parallel::mclapply(seq_len(nrow(pick)), function(i) tryCatch({
    r <- stage2_microtime_job(study, pick[i, ], cfg, frozen, step); tmp <- paste0(micro_file(study, pick[i, ], step), ".tmp"); saveRDS(r, tmp); file.rename(tmp, micro_file(study, pick[i, ], step)); "done"
  }, error = function(e) { writeLines(conditionMessage(e), paste0(micro_file(study, pick[i, ], step), ".error")); "error" }), mc.cores = workers, mc.preschedule = FALSE)
  stage2_microtiming_summary(study)
}

# Prices every arm that has a micro measurement. An arm with no measurement is listed as missing and makes the projection incomplete.
stage2_microtiming_summary <- function(study) {
  cfg <- readRDS(file.path(study, "config.rds")); frozen <- ladder_settings(cfg, 1L)
  files <- list.files(file.path(study, "microtiming"), pattern = "\\.rds$", full.names = TRUE)
  normalise <- function(f) { r <- readRDS(f); if (!"step" %in% names(r)) r$step <- 1L; if (!"signature" %in% names(r)) r$signature <- timing_signature(r$key, frozen); r }
  micro_tab <- if (length(files)) do.call(rbind, lapply(files, normalise)) else NULL
  errors <- list.files(file.path(study, "microtiming"), pattern = "\\.error$", full.names = TRUE)
  jobs <- job_table(study); n_jobs <- table(jobs$key); micro <- micro_schedule; draws_micro <- micro[["niter"]] * micro[["nchain"]]
  blocks <- ceiling(1600 / cfg$prediction_block)
  draws_prod_init <- cfg$mcmc$initial[["nchain"]] * cfg$mcmc$initial[["niter"]]; draws_prod_ext <- cfg$mcmc$extended[["nchain"]] * cfg$mcmc$extended[["niter"]]
  iters_init <- cfg$mcmc$initial[["nchain"]] * (cfg$mcmc$initial[["nburn"]] + cfg$mcmc$initial[["niter"]])
  per_job_seconds <- function(r, frozen) {
    if (r$kind == "bayesian") {
      per_iter <- r$fit_seconds / (micro[["nburn"]] + micro[["niter"]]); diag_per_draw <- r$diagnose_seconds / draws_micro
      c(fit = per_iter * iters_init * cfg$extension_factor, diagnose = diag_per_draw * (draws_prod_init + draws_prod_ext),
        prediction = r$block_prediction_seconds / draws_micro * (draws_prod_ext / frozen$prediction_thin) * blocks)
    } else c(fit = r$fit_seconds / 200 * cfg$gllvm$maxit * cfg$diagnostics$max_starts, diagnose = 0, prediction = r$block_prediction_seconds * blocks)
  }
  all_keys <- arm_key(cfg$arms$package, cfg$arms$arm); measured <- if (is.null(micro_tab)) character() else unique(micro_tab$key); missing <- setdiff(all_keys, measured)
  # For each arm and step: the measurement taken at that step's own setting if there is one, else the step 1 measurement rescaled.
  price <- function(st, rescale) {
    basis <- character(); arms <- lapply(measured, function(k) {
      own <- micro_tab[micro_tab$key == k & micro_tab$signature == timing_signature(k, st), ]
      if (nrow(own)) { rows <- own; f <- 1; basis[[k]] <<- paste("measured at step", max(own$step)) }
      else { rows <- micro_tab[micro_tab$key == k & micro_tab$step == 1L, ]; f <- rescale(k); basis[[k]] <<- "rescaled from step 1" }
      if (!nrow(rows)) return(NULL)
      rows <- rows[order(rows$step), ]; rows <- rows[!duplicated(rows$truth, fromLast = TRUE), ]
      secs <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) per_job_seconds(rows[i, ], frozen))); s <- secs[which.max(rowSums(secs)), ]
      data.frame(key = k, fit_hours = s[["fit"]] * f / 3600, diagnose_hours = s[["diagnose"]] * f / 3600, prediction_hours = s[["prediction"]] * f / 3600,
                 job_hours = sum(s) * f / 3600, jobs = as.integer(n_jobs[k]), worker_hours = sum(s) * f / 3600 * as.integer(n_jobs[k]), stringsAsFactors = FALSE)
    })
    t <- do.call(rbind, arms)
    list(per_arm = t, basis = as.list(basis), worker_hours_measured_arms = if (is.null(t)) 0 else sum(t$worker_hours),
         wall_hours_at_workers = if (is.null(t)) 0 else sum(t$worker_hours) / cfg$workers, complete = length(missing) == 0L)
  }
  steps <- lapply(1:8, function(s) {
    st <- ladder_settings(cfg, s)
    rescale <- function(k) if (k == "occJSDM/gp") (st$support_points_900 / frozen$support_points_900)^2 else
      if (grepl("^Hmsc/", k) && identical(st$hmsc_method, "GPP")) st$hmsc_knots / 900 else 1
    c(list(step = s, settings = st), price(st, rescale))
  })
  out <- list(micro_fits = micro_tab, missing_arms = missing, error_files = basename(errors), assumptions = list(
    micro_schedule = "one chain, 50 warm-up, 150 retained; gllvm one start with maxit 200",
    bayesian_fit = "per-iteration cost times 4 chains times 7,000 iterations times the extension factor 4.43",
    diagnostics = "per-retained-draw cost times (4 x 4,000 + 4 x 16,000) draws",
    prediction = "per-draw, per-25-site-block cost times 4 x 16,000 / prediction thinning draws times 64 blocks",
    gllvm = "cost per iteration at the cap times maxit 12,000 times 6 starts: a worst case that assumes no start converges early",
    jobs_priced = "every job at the 900-site cost of its arm, the larger of the two truths (upper bound for the 300-site jobs)",
    missing = "an arm with no measurement is excluded from the totals and the projection is marked incomplete"), steps = steps)
  writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"), file.path(study, "microtiming.json"))
  invisible(out)
}

stage2_freeze <- function(study, settings, timing) {
  cfg <- readRDS(file.path(study, "config.rds"))
  stopifnot(all(is.na(unlist(cfg$frozen))), setequal(names(settings), names(cfg$frozen)))
  cfg$frozen <- settings; cfg$frozen_at <- Sys.time(); cfg$frozen_from <- timing
  saveRDS(cfg, file.path(study, "config.rds")); invisible(cfg)
}

stage2_fit <- function(study, workers = 8L, jobs = NULL, mcmc = NULL, slice_sites = NULL, slice_species = NULL) {
  cfg <- readRDS(file.path(study, "config.rds"))
  if (any(is.na(unlist(cfg$frozen)[c("support_points_900", "hmsc_method", "prediction_thin")]))) stop("Settings are not frozen; run timing and freeze first")
  all_jobs <- job_table(study); if (!is.null(jobs)) all_jobs <- all_jobs[all_jobs$job_id %in% jobs, ]
  all_jobs <- all_jobs[order(-all_jobs$n, all_jobs$job_id), ]   # the 900-site jobs dominate the cost, so start them first
  if (is.null(mcmc)) mcmc <- cfg$mcmc
  res <- parallel::mclapply(seq_len(nrow(all_jobs)), function(i) tryCatch(run_job(study, all_jobs[i, ], cfg, cfg$frozen, mcmc, slice_sites, slice_species)$status,
    error = function(e) { writeLines(conditionMessage(e), file.path(study, "status", paste0(all_jobs$job_id[i], ".error"))); "error" }),
    mc.cores = workers, mc.preschedule = FALSE)
  table(unlist(res))
}

stage2_score <- function(study, species_subset = NULL) {
  cfg <- readRDS(file.path(study, "config.rds")); jobs <- job_table(study)
  done <- jobs[file.exists(file.path(study, "status", paste0(jobs$job_id, ".json"))), ]
  status_of <- function(id) jsonlite::fromJSON(file.path(study, "status", paste0(id, ".json")))$status
  rows <- list(); job_rows <- list()
  for (ds in unique(done$dataset_id)) {
    cm <- readRDS(done$input_file[done$dataset_id == ds][1]); truth <- load_truth(cm)
    if (!is.null(species_subset)) for (nm in c("psi_lattice", "field_lattice", "eta_lattice", "psi_training", "field_training")) truth[[nm]] <- truth[[nm]][, species_subset, drop = FALSE]
    if (!is.null(species_subset)) { truth$occupancy <- truth$occupancy[species_subset]; truth$species <- truth$species[species_subset]; truth$spatial_sd <- truth$spatial_sd[species_subset] }
    oracle <- load_oracle_reference(done$oracle_file[done$dataset_id == ds][1])
    if (!is.null(species_subset)) oracle$probability_mean <- oracle$probability_mean[, species_subset, drop = FALSE]
    load_pred <- function(id) { p <- file.path(study, "jobs", "predictions", paste0(id, ".rds")); if (file.exists(p)) readRDS(p) else NULL }
    env_id <- done$job_id[done$dataset_id == ds & done$key == "occJSDM/env"]; env <- if (length(env_id)) load_pred(env_id) else NULL
    for (i in which(done$dataset_id == ds)) {
      j <- done[i, ]; st <- status_of(j$job_id); p <- load_pred(j$job_id)
      rows[[length(rows) + 1L]] <- score_arm(if (is.null(p)) NULL else p$prediction, if (is.null(env)) NULL else env$prediction, oracle, truth, ds, j$package, j$arm,
                                             if (is.null(env) && st %in% c("pass", "unresolved", "unstable")) "unsupported" else st)
      a <- list.files(file.path(study, "jobs", "attempts"), pattern = paste0("^", j$job_id, "-[0-9]+\\.rds$"), full.names = TRUE)
      sel <- select_job(study, j, cfg)
      job_rows[[length(job_rows) + 1L]] <- data.frame(job_id = j$job_id, dataset_id = ds, cell = j$cell, key = j$key, status = st, attempts = length(a),
        selected_attempt = sel$selected_attempt, fit_seconds = sum(read_attempt_seconds(study, j$job_id)),
        prediction_seconds = if (is.null(p)) NA_real_ else p$elapsed_seconds, oracle_status = oracle$status, stringsAsFactors = FALSE)
    }
  }
  species <- do.call(rbind, rows); write.csv(species, file.path(study, "scores", "species.csv"), row.names = FALSE)
  agg <- aggregate_arm(species, cfg$gate_occupancy); write.csv(agg, file.path(study, "scores", "aggregate.csv"), row.names = FALSE)
  pairs <- paired_comparisons(agg, cfg); write.csv(pairs, file.path(study, "scores", "pairs.csv"), row.names = FALSE)
  finite <- agg[is.finite(agg$headroom_share), ]
  max_share <- if (nrow(finite)) max(aggregate(headroom_share ~ cell + key, finite, median)$headroom_share) else NA_real_
  writeLines(jsonlite::toJSON(development_reading(pairs, max_share), auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA), file.path(study, "scores", "reading.json"))
  write.csv(do.call(rbind, job_rows), file.path(study, "scores", "jobs.csv"), row.names = FALSE)
  invisible(pairs)
}

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  get <- function(key, default = NULL) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), args, value = TRUE)); if (length(v)) v else default }
  mode <- get("mode"); study <- get("study")
  switch(mode,
    prepare = stage2_prepare(study, get("pilot"), jobs_dir = get("jobs")),
    microtiming = { k <- get("keys"); print(lapply(stage2_microtiming(study, keys = if (is.null(k)) NULL else strsplit(k, ",")[[1]], step = as.integer(get("step", "1")))$steps, function(s) s[c("step", "worker_hours_measured_arms", "wall_hours_at_workers", "complete")])) },
    timing = print(stage2_timing(study, step = as.integer(get("step", "1")))[c("settings", "complete", "projected_cpu_hours_upper", "projected_wall_hours_upper", "passes")]),
    freeze = { s <- as.integer(get("step")); t <- jsonlite::fromJSON(file.path(study, sprintf("timing-step%d.json", s))); stopifnot(isTRUE(t$passes))
               stage2_freeze(study, ladder_settings(readRDS(file.path(study, "config.rds")), s), timing = t[c("step", "projected_wall_hours_upper")]) },
    fit = print(stage2_fit(study, workers = as.integer(get("workers", "8")))),
    score = print(stage2_score(study)),
    stop("Unknown mode ", mode))
  cat("Completed mode", mode, "at", format(Sys.time()), "\n")
}
