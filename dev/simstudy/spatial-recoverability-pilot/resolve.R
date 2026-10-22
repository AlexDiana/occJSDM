# Prescribed resolution extension for the recoverability pilot (spec revision record,
# 11 October 2026). Refits, at a fixed multiple of the extension budget, the fits whose
# convergence screen failed in the cells that decide the map: truth 1 at range 0.08 (the
# whole job is refit jointly) and truth 2 at n = 900 in the medium and long classes (each
# unresolved gate species is refit alone; resolved species are carried over unchanged).
# Results go to jobs-resolution/ and replace the corresponding fits at score time. The
# original job files are never modified. Sourced by run.R after its other sources.

resolution_budget <- function(budget, multiplier = 4L) {
  list(nburn = as.integer(budget$ext_nburn * multiplier), niter = as.integer(budget$ext_niter * multiplier),
       thin = as.integer(budget$thin * 4L * multiplier), multiplier = as.integer(multiplier))
}

default_resolution_scope <- function() list(t1_ranges = 0.08, t1_n = c(100L, 300L, 900L), t2_n = 900L,
                                            t2_classes = c("medium", "long"), gate_occupancy = c(0.25, 0.50, 0.75))

resolution_dir <- function(study) {
  link <- file.path(study, "jobs-resolution")
  if (!file.exists(link)) {
    target <- Sys.readlink(file.path(study, "jobs"))
    dir <- if (nzchar(target)) file.path(dirname(target), "jobs-resolution") else file.path(study, "jobs-resolution-local")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    if (normalizePath(dir) != normalizePath(link, mustWork = FALSE)) file.symlink(normalizePath(dir), link)
  }
  link
}

# Writes resolution.json from the current species.csv. Refuses to run twice.
pilot_prescribe_resolution <- function(study, scope = default_resolution_scope(), multiplier = 4L) {
  path <- file.path(study, "resolution.json")
  stopifnot(!file.exists(path))
  species <- read.csv(file.path(study, "scores", "species.csv"), stringsAsFactors = FALSE)
  species$job <- sprintf("%s-rep%d-%s", species$cell, species$replicate, species$arrangement)
  unresolved_gate <- species$occupancy %in% scope$gate_occupancy & species$status == "unresolved"
  t1 <- species[species$truth == "t1" & species$range %in% scope$t1_ranges & species$n %in% scope$t1_n & unresolved_gate, ]
  t2 <- species[species$truth == "t2" & species$n %in% scope$t2_n & species$class %in% scope$t2_classes & unresolved_gate, ]
  t2_jobs <- sort(unique(t2$job))
  t2_refits <- lapply(t2_jobs, function(j) sort(t2$species[t2$job == j])); names(t2_refits) <- t2_jobs
  budget <- readRDS(file.path(study, "config.rds"))$budget
  prescription <- list(created = format(Sys.time()), git = system("git rev-parse HEAD", intern = TRUE),
                       scope = scope, multiplier = as.integer(multiplier), budget = resolution_budget(budget, multiplier),
                       rule = paste("Truth 1: every job in the scoped cells with any unresolved gate species is refit jointly at the resolution budget.",
                                    "Truth 2: every unresolved gate species of the scoped classes in the scoped cells is refit alone at the resolution budget;",
                                    "resolved species are carried over unchanged. Thresholds, screen and seeds-by-label scheme unchanged; seed labels carry the suffix 'resolution'."),
                       t1_jobs = sort(unique(t1$job)), t2_jobs = t2_jobs, t2_refits = t2_refits,
                       counts = list(t1_jobs = length(unique(t1$job)), t2_jobs = length(t2_jobs), t2_species = nrow(t2)))
  writeLines(jsonlite::toJSON(prescription, auto_unbox = TRUE, pretty = TRUE), path)
  resolution_dir(study)
  invisible(prescription)
}

read_prescription <- function(study) {
  p <- jsonlite::fromJSON(file.path(study, "resolution.json"), simplifyVector = TRUE)
  p$t1_jobs <- as.character(unlist(p$t1_jobs)); p$t2_jobs <- as.character(unlist(p$t2_jobs))
  p$t2_refits <- lapply(p$t2_refits, as.character)
  p
}

run_resolution_job <- function(study, job_name, refit_species, budget) {
  out <- file.path(resolution_dir(study), paste0(job_name, ".rds"))
  if (file.exists(out)) return(invisible("existing"))
  original_path <- file.path(study, "jobs", paste0(job_name, ".rds"))
  orig <- readRDS(original_path)
  jobs <- job_table(study); job <- jobs[jobs$job == job_name, ]
  cm <- readRDS(file.path(study, "inputs", job$file)); arrangement <- job$arrangement; S <- nrow(cm$species)
  t0 <- Sys.time()
  result <- orig
  if (cm$truth == "t1") {
    fit <- run_joint_oracle(cm, arrangement, budget$nburn, budget$niter, budget$thin,
                            seed = pilot_seed(pilot_label("oracle", job_name, "resolution")))
    sc <- screen_traces(fit$traces)
    prediction <- predict_lattice(cm, arrangement, list(joint = fit), seed = pilot_seed(pilot_label("lattice", job_name, "resolution")))
    result$site_draws <- fit$site_draws; result$traces <- fit$traces
    result$screen <- data.frame(species = seq_len(S), status = vapply(seq_len(S), function(s) if (all(sc$pass[sc$species == s])) "pass" else "unresolved", character(1)),
                                stringsAsFactors = FALSE)
    result$screen_detail <- sc; result$phases <- "resolution"; result$mean_proposals <- fit$mean_proposals
    result$block_diagnostics <- fit[c("block_proposals", "block_mean_abs_angle")]
    result$factor_records <- fit$factor_record; result$conditional_records <- prediction$conditional_records
    result$prediction <- prediction[c("probability_mean", "field_mean", "site_field_mean", "site_probability_mean")]
    refit <- seq_len(S)
  } else {
    refit <- match(refit_species, cm$species$species)
    stopifnot(!anyNA(refit), length(refit) >= 1L)
    fits <- lapply(refit, function(s) run_species_oracle(cm, arrangement, s, budget$nburn, budget$niter, budget$thin,
                                                         seed = pilot_seed(pilot_label("oracle", job_name, cm$species$species[s], "resolution"))))
    prediction <- predict_lattice(cm, arrangement, list(species = refit, oracles = fits),
                                  seed = pilot_seed(pilot_label("lattice", job_name, "resolution")))
    for (k in seq_along(refit)) {
      s <- refit[k]
      result$site_draws[[s]] <- fits[[k]]$site_draws[[1]]; result$traces[[s]] <- fits[[k]]$traces
      sc <- screen_traces(fits[[k]]$traces); sc$species <- s
      result$screen$status[result$screen$species == s] <- if (all(sc$pass)) "pass" else "unresolved"
      result$screen_detail <- rbind(result$screen_detail[result$screen_detail$species != s, ], sc[names(result$screen_detail)])
      result$phases[s] <- "resolution"; result$mean_proposals[s] <- mean(fits[[k]]$mean_proposals)
      for (nm in names(result$prediction)) result$prediction[[nm]][, s] <- prediction[[nm]][, k]
    }
    result$factor_records <- rbind(orig$factor_records, do.call(rbind, lapply(fits, function(f) { r <- f$factor_record; r$what <- paste0(r$what, "-resolution"); r })))
    result$conditional_records <- rbind(orig$conditional_records, within(prediction$conditional_records, what <- paste0(what, "-resolution")))
  }
  result$scores <- score_species(cm, arrangement, result$prediction, result$screen)
  result$resolution <- list(refit_species = cm$species$species[refit], budget = budget, original_job_md5 = unname(tools::md5sum(original_path)),
                            elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), finished = Sys.time())
  tmp <- paste0(out, ".tmp"); saveRDS(result, tmp); file.rename(tmp, out)
  writeLines(jsonlite::toJSON(list(job = job_name, status = if (all(result$screen$status == "pass")) "pass" else "unresolved",
                                   refit = cm$species$species[refit], elapsed_seconds = round(result$resolution$elapsed_seconds, 1)), auto_unbox = TRUE),
             file.path(study, "status", paste0(job_name, ".resolution.json")))
  invisible("fitted")
}

pilot_resolve <- function(study, workers = 8L) {
  compile_pilot_samplers(".")
  p <- read_prescription(study)
  jobs <- c(p$t1_jobs, p$t2_jobs)
  res <- parallel::mclapply(seq_along(jobs), function(i) tryCatch(
    run_resolution_job(study, jobs[i], if (jobs[i] %in% p$t2_jobs) p$t2_refits[[jobs[i]]] else NULL, p$budget),
    error = function(e) { writeLines(conditionMessage(e), file.path(study, "status", paste0(jobs[i], ".resolution.error"))); "error" }),
    mc.cores = workers, mc.preschedule = FALSE)
  table(unlist(res))
}
