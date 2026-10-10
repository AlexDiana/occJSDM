# Driver for the spatial recoverability pilot. Always run through run-r:
#   run-r -f run.R --args --mode=prepare --study=DIR --jobs=VOLUME_DIR [--workers=8]
#   run-r -f run.R --args --mode=timing  --study=DIR
#   run-r -f run.R --args --mode=fit     --study=DIR --truth=t2 --workers=8
#   run-r -f run.R --args --mode=score   --study=DIR
PILOT_DIR <- "dev/simstudy/spatial-recoverability-pilot"
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R", "predict.R", "score.R")) source(file.path(PILOT_DIR, f))

default_budget <- function() list(nburn = 3000L, niter = 4000L, thin = 10L, ext_nburn = 8000L, ext_niter = 16000L)

pilot_prepare <- function(study, communities = 10L, cells = pilot_cells(), jobs_dir, workers = 1L) {
  dir.create(study, showWarnings = FALSE, recursive = TRUE)
  stopifnot(!file.exists(file.path(study, "config.rds")))
  dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(file.path(study, "jobs"))) file.symlink(normalizePath(jobs_dir), file.path(study, "jobs"))
  for (d in c("inputs", "status", "scores", "logs")) dir.create(file.path(study, d), showWarnings = FALSE)
  config <- list(created = Sys.time(), communities = communities, cells = cells, thresholds = default_thresholds(),
                 budget = default_budget(), r_version = R.version.string, r_home = R.home(),
                 packages = sapply(c("Rcpp", "RcppArmadillo", "posterior"), function(p) as.character(packageVersion(p))),
                 sources = tools::md5sum(c(list.files(PILOT_DIR, pattern = "\\.(R|cpp)$", full.names = TRUE),
                                           "dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp")),
                 git = system("git rev-parse HEAD", intern = TRUE))
  saveRDS(config, file.path(study, "config.rds"))
  Sys.sleep(1)   # config strictly predates every input
  grid <- merge(cells, data.frame(replicate = seq_len(communities)), by = NULL)
  rows <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    cm <- make_community(g$truth, g$n, g$replicate, range = g$range)
    path <- file.path(study, "inputs", sprintf("%s-rep%d.rds", g$cell, g$replicate))
    saveRDS(cm, path)
    write.csv(cbind(cell = g$cell, replicate = g$replicate, cm$conditioning),
              file.path(study, "logs", sprintf("generator-conditioning-%s-rep%d.csv", g$cell, g$replicate)), row.names = FALSE)
    data.frame(cell = g$cell, truth = g$truth, n = g$n, range = g$range, replicate = g$replicate, file = basename(path),
               md5 = unname(tools::md5sum(path)), seed_sites = cm$seeds[["sites"]], seed_data = cm$seeds[["data"]],
               rejections = cm$arrangements$rejections, stringsAsFactors = FALSE)
  }, mc.cores = workers, mc.preschedule = FALSE)
  bad <- vapply(rows, function(r) inherits(r, "try-error") || !is.data.frame(r), logical(1))
  if (any(bad)) stop("Input generation failed for ", sum(bad), " communities: ", paste(unlist(rows[bad]), collapse = "; "))
  manifest <- do.call(rbind, rows)
  write.csv(manifest, file.path(study, "input-manifest.csv"), row.names = FALSE)
  invisible(config)
}

job_table <- function(study) {
  man <- read.csv(file.path(study, "input-manifest.csv"), stringsAsFactors = FALSE)
  jobs <- merge(man, data.frame(arrangement = c("uniform", "clustered"), stringsAsFactors = FALSE), by = NULL)
  jobs$job <- sprintf("%s-rep%d-%s", jobs$cell, jobs$replicate, jobs$arrangement)
  jobs[order(jobs$job), ]
}

fit_with_extension <- function(budget, fit_fun, seed_label) {
  fit <- fit_fun(nburn = budget$nburn, niter = budget$niter, thin = budget$thin, seed = pilot_seed(seed_label))
  sc <- screen_traces(fit$traces); phase <- "initial"
  if (!all(sc$pass)) {
    fit <- fit_fun(nburn = budget$ext_nburn, niter = budget$ext_niter, thin = budget$thin * 4L,
                   seed = pilot_seed(paste0(seed_label, "-extension")))
    sc <- screen_traces(fit$traces); phase <- "extension"
  }
  list(fit = fit, screen = sc, phase = phase, status = if (all(sc$pass)) "pass" else "unresolved")
}

run_job <- function(study, job, budget) {
  out <- file.path(study, "jobs", paste0(job$job, ".rds"))
  if (file.exists(out)) return(invisible("existing"))
  cm <- readRDS(file.path(study, "inputs", job$file)); arrangement <- job$arrangement
  t0 <- Sys.time(); S <- nrow(cm$species)
  if (cm$truth == "t1") {
    r <- fit_with_extension(budget, function(nburn, niter, thin, seed) run_joint_oracle(cm, arrangement, nburn, niter, thin, seed),
                            pilot_label("oracle", job$job))
    screen <- data.frame(species = seq_len(S),
                         status = vapply(seq_len(S), function(s) if (all(r$screen$pass[r$screen$species == s])) "pass" else "unresolved", character(1)),
                         stringsAsFactors = FALSE)
    oracle <- list(joint = r$fit); phases <- r$phase; factor_records <- r$fit$factor_record
    draws <- r$fit$site_draws; proposals <- r$fit$mean_proposals; screen_detail <- r$screen; traces <- r$fit$traces
  } else {
    fits <- lapply(seq_len(S), function(s) fit_with_extension(budget,
      function(nburn, niter, thin, seed) run_species_oracle(cm, arrangement, s, nburn, niter, thin, seed),
      pilot_label("oracle", job$job, cm$species$species[s])))
    screen <- data.frame(species = seq_len(S), status = vapply(fits, `[[`, character(1), "status"), stringsAsFactors = FALSE)
    oracle <- list(species = seq_len(S), oracles = lapply(fits, `[[`, "fit")); phases <- vapply(fits, `[[`, character(1), "phase")
    factor_records <- do.call(rbind, lapply(fits, function(f) f$fit$factor_record))
    draws <- lapply(fits, function(f) f$fit$site_draws[[1]]); proposals <- vapply(fits, function(f) mean(f$fit$mean_proposals), numeric(1))
    screen_detail <- do.call(rbind, lapply(seq_len(S), function(s) { d <- fits[[s]]$screen; d$species <- s; d }))
    traces <- lapply(fits, function(f) f$fit$traces)
  }
  prediction <- predict_lattice(cm, arrangement, oracle, seed = pilot_seed(pilot_label("lattice", job$job)))
  scores <- score_species(cm, arrangement, prediction, screen)
  geometry <- layout_geometry(cm, arrangement)
  result <- list(job = job$job, cell = cm$cell, replicate = cm$replicate, arrangement = arrangement, truth = cm$truth, n = cm$n,
                 site_draws = draws, traces = traces,
                 prediction = prediction[c("probability_mean", "field_mean", "site_field_mean", "site_probability_mean")],
                 scores = scores, screen = screen, screen_detail = screen_detail, phases = phases, mean_proposals = proposals,
                 factor_records = factor_records, conditional_records = prediction$conditional_records, geometry = geometry,
                 budget = budget, elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), finished = Sys.time())
  tmp <- paste0(out, ".tmp"); saveRDS(result, tmp); file.rename(tmp, out)
  writeLines(jsonlite::toJSON(list(job = job$job, status = if (all(screen$status == "pass")) "pass" else "unresolved",
                                   elapsed_seconds = round(result$elapsed_seconds, 1), phases = paste(phases, collapse = ",")), auto_unbox = TRUE),
             file.path(study, "status", paste0(job$job, ".json")))
  invisible("fitted")
}

pilot_fit <- function(study, truth, workers = 8L, budget = readRDS(file.path(study, "config.rds"))$budget, jobs = NULL) {
  compile_pilot_samplers(".")
  all_jobs <- job_table(study); all_jobs <- all_jobs[all_jobs$truth == truth, ]
  if (!is.null(jobs)) all_jobs <- all_jobs[all_jobs$job %in% jobs, ]
  res <- parallel::mclapply(seq_len(nrow(all_jobs)), function(i) tryCatch(run_job(study, all_jobs[i, ], budget),
    error = function(e) { writeLines(conditionMessage(e), file.path(study, "status", paste0(all_jobs$job[i], ".error"))); "error" }),
    mc.cores = workers, mc.preschedule = FALSE)
  table(unlist(res))
}

pilot_timing <- function(study) {
  compile_pilot_samplers(".")
  jobs <- job_table(study); budget <- readRDS(file.path(study, "config.rds"))$budget
  pick <- c(jobs$job[jobs$truth == "t2" & jobs$n == 900L & jobs$arrangement == "uniform"][1],
            jobs$job[jobs$truth == "t1" & jobs$n == 900L & !is.na(jobs$range) & jobs$range == 0.15 & jobs$arrangement == "uniform"][1])
  times <- lapply(pick, function(j) { t0 <- Sys.time(); r <- run_job(study, jobs[jobs$job == j, ], budget)
    list(result = r, seconds = as.numeric(Sys.time() - t0, units = "secs")) })
  names(times) <- pick
  n_jobs <- table(jobs$truth)
  secs <- vapply(times, `[[`, numeric(1), "seconds")
  projection <- list(timing = times, jobs = as.list(n_jobs),
                     projected_cpu_hours_upper = unname((secs[1] * n_jobs[["t2"]] + secs[2] * n_jobs[["t1"]]) / 3600),
                     note = "Upper bound: every job priced at the 900-site cost; extensions add up to 3x for the jobs that need them.")
  writeLines(jsonlite::toJSON(projection, auto_unbox = TRUE, pretty = TRUE, digits = 4), file.path(study, "timing.json"))
  projection
}

pilot_score <- function(study) {
  files <- list.files(file.path(study, "jobs"), pattern = "\\.rds$", full.names = TRUE)
  results <- lapply(files, function(f) { r <- readRDS(f); r$site_draws <- NULL; r$traces <- NULL; r })
  species <- do.call(rbind, lapply(results, `[[`, "scores"))
  write.csv(species, file.path(study, "scores", "species.csv"), row.names = FALSE)
  agg <- do.call(rbind, lapply(split(species, species$truth), aggregate_scores)); rownames(agg) <- NULL
  write.csv(agg, file.path(study, "scores", "aggregate.csv"), row.names = FALSE)
  th <- readRDS(file.path(study, "config.rds"))$thresholds
  uni <- agg[agg$arrangement == "uniform", ]
  uni$key <- paste(uni$cell, uni$group)
  gate <- do.call(rbind, lapply(split(uni, uni$key), function(rows) cbind(rows[1, c("truth", "cell", "n", "range", "group")], gate_cells(rows, th))))
  rownames(gate) <- NULL
  gate$class <- ifelse(gate$truth == "t2", gate$group, NA_character_)
  write.csv(gate, file.path(study, "scores", "gate.csv"), row.names = FALSE)
  writeLines(jsonlite::toJSON(read_map(gate), auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(study, "scores", "map.json"))
  write.csv(do.call(rbind, lapply(results, function(r) cbind(job = r$job, r$factor_records))), file.path(study, "scores", "factors.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, function(r) cbind(job = r$job, r$conditional_records))), file.path(study, "scores", "conditionals.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, function(r) cbind(job = r$job, cell = r$cell, replicate = r$replicate, r$geometry))), file.path(study, "scores", "geometry.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(results, function(r) cbind(job = r$job, r$screen_detail))), file.path(study, "scores", "screen.csv"), row.names = FALSE)
  status <- do.call(rbind, lapply(results, function(r) data.frame(job = r$job, truth = r$truth, cell = r$cell, replicate = r$replicate,
    arrangement = r$arrangement, status = if (all(r$screen$status == "pass")) "pass" else "unresolved",
    extensions = sum(r$phases == "extension"), elapsed_seconds = r$elapsed_seconds, mean_proposals = mean(r$mean_proposals), stringsAsFactors = FALSE)))
  write.csv(status, file.path(study, "scores", "jobs.csv"), row.names = FALSE)
  invisible(gate)
}

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  get <- function(key, default = NULL) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), args, value = TRUE)); if (length(v)) v else default }
  mode <- get("mode"); study <- get("study")
  switch(mode,
    prepare = pilot_prepare(study, jobs_dir = get("jobs"), workers = as.integer(get("workers", "1"))),
    timing = print(pilot_timing(study)),
    fit = print(pilot_fit(study, truth = get("truth"), workers = as.integer(get("workers", "8")))),
    score = print(pilot_score(study)),
    stop("Unknown mode ", mode))
  cat("Completed mode", mode, "at", format(Sys.time()), "\n")
}
