#!/usr/bin/env Rscript
# Lesson 2 sweep driver. Modes: prepare, oracle, freeze, pilot, initial, long.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit) > 1L) stop("Repeated option: ", name)
  if (!length(hit)) { if (is.null(default)) stop("Missing --", name); return(default) }
  substring(hit, nchar(name) + 4L)
}
known <- c("repo", "study", "mode", "workers", "keys")
stopifnot(all(sub("^--([^=]+)=.*$", "\\1", args) %in% known))
repo <- normalizePath(option("repo", ".")); study <- option("study")
mode <- option("mode"); workers <- as.integer(option("workers", "2"))
stopifnot(mode %in% c("prepare", "oracle", "freeze", "pilot", "initial", "long"), workers %in% 1:8)
dir.create(study, showWarnings = FALSE, recursive = TRUE); study <- normalizePath(study)
scripts <- file.path(repo, "dev/simstudy/spatial-design-sweep")
source(file.path(scripts, "generator.R"))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
atomic <- function(object, path) {
  tmp <- paste0(path, ".tmp"); saveRDS(object, tmp); stopifnot(file.rename(tmp, path))
}
say <- function(...) { cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console() }
# Write a csv once; if it already exists, require it to equal the regenerated frame
# after a read.csv round trip (numeric formatting differs, so compare with all.equal).
publish_csv <- function(frame, path, drop_paths = character()) {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp))
  if (!file.exists(path)) { write.csv(frame, path, row.names = FALSE); return(invisible(FALSE)) }
  write.csv(frame, tmp, row.names = FALSE)
  old <- read.csv(path, stringsAsFactors = FALSE); new <- read.csv(tmp, stringsAsFactors = FALSE)
  for (column in drop_paths) { old[[column]] <- basename(old[[column]]); new[[column]] <- basename(new[[column]]) }
  stopifnot(isTRUE(all.equal(old, new)))
  invisible(TRUE)
}
arrangements <- c("spread", "pairs", "clustered", "grid")
generator_hash <- tools::md5sum(file.path(scripts, "generator.R"))

# ---- prepare --------------------------------------------------------------
dir.create(file.path(study, "inputs"), showWarnings = FALSE)
inputs <- character()
for (r in 1:3) {
  community <- sprintf("rep%02d", r)
  path <- file.path(study, "inputs", paste0(community, ".rds"))
  input <- make_sweep_input(r)
  if (file.exists(path)) stopifnot(identical(readRDS(path), input)) else atomic(input, path)
  inputs[community] <- path
}
statistics <- do.call(rbind, lapply(inputs, function(p) readRDS(p)$statistics))
rownames(statistics) <- NULL
stat_path <- file.path(study, "design-statistics.csv")
publish_csv(statistics, stat_path)
if (any(!statistics$in_grid)) say("WARNING: a design's standardised range lies outside the fitter's grid; see design-statistics.csv")
manifest <- expand.grid(replicate = 1:3, arrangement = arrangements, arm = c("binary", "two_stage"),
                        stringsAsFactors = FALSE)
manifest$community <- sprintf("rep%02d", manifest$replicate)
manifest$key <- sprintf("%s-%s-%s", manifest$community, manifest$arrangement, manifest$arm)
manifest$input_md5 <- unname(tools::md5sum(inputs[manifest$community]))
manifest$fit_seed <- mapply(function(a, r, arr) sweep_seed(paste0("spatial-design-fit-", a, "-", arr), r),
                            manifest$arm, manifest$replicate, manifest$arrangement)
manifest <- manifest[order(manifest$replicate, match(manifest$arrangement, arrangements), manifest$arm), ]
rownames(manifest) <- NULL
manifest_path <- file.path(study, "input-manifest.rds")
if (file.exists(manifest_path)) stopifnot(identical(readRDS(manifest_path), manifest)) else atomic(manifest, manifest_path)
publish_csv(manifest, file.path(study, "input-manifest.csv"))
if (mode == "prepare") { say("Prepared 3 communities, 4 arrangements, 24 fit jobs; generator md5", generator_hash); quit(status = 0) }

# ---- oracle ---------------------------------------------------------------
if (mode == "oracle") {
  source(file.path(scripts, "oracle.R"))
  out <- file.path(study, "oracle"); dir.create(out, showWarnings = FALSE)
  hashes <- compile_oracle(repo, cache = file.path(out, "cpp-cache"))
  jobs <- expand.grid(community = names(inputs), arrangement = arrangements, species = 1:8,
                      stringsAsFactors = FALSE)
  work <- function(i) {
    job <- jobs[i, ]; input <- readRDS(inputs[job$community]); tr <- input$surveys[[job$arrangement]]$truth
    key <- sprintf("%s-%s-species%02d", job$community, job$arrangement, job$species)
    seed <- sweep_seed(paste0("spatial-design-oracle-", job$arrangement, "-", job$species), input$replicate)
    lattice_seed <- sweep_seed(paste0("spatial-design-oracle-lattice-", job$arrangement, "-", job$species), input$replicate)
    perform <- function(phase) {
      dest <- file.path(out, paste0(key, "-", phase, ".rds"))
      nburn <- if (phase == "initial") 1000L else 2000L
      niter <- if (phase == "initial") 2000L else 4000L
      if (file.exists(dest)) {
        res <- readRDS(dest)
        stopifnot(identical(unname(res$hashes), unname(hashes)), res$seed == seed, res$phase == phase,
                  res$nburn == nburn, res$niter == niter, res$key == key,
                  identical(unname(res$generator_hash), unname(generator_hash)),
                  identical(res$input_md5, unname(tools::md5sum(inputs[job$community]))))
        return(res)
      }
      started <- Sys.time()
      fit <- run_oracle(input, job$arrangement, job$species, seed, nburn, niter)
      scored <- score_oracle(fit$draws, tr$field[, job$species], tr$psi[, job$species], fit$offset)
      lattice_index <- input$landscape$index$lattice
      offset_lattice <- input$landscape$B0[job$species] + input$landscape$B[job$species] * input$landscape$environment[lattice_index]
      lat <- oracle_lattice(fit$draws, tr$xy, input$landscape$points[lattice_index, ], input$landscape$range,
                            offset_lattice, thin = 4L, seed = lattice_seed)
      res <- c(list(job = job, key = key, phase = phase, seed = seed, nburn = nburn, niter = niter, nchain = 4L,
                    hashes = hashes, generator_hash = generator_hash, input_md5 = unname(tools::md5sum(inputs[job$community])),
                    elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")),
                    draws = fit$draws, mean_proposals = fit$mean_proposals, lattice = lat), scored)
      atomic(res, dest); say(key, phase, "seconds", round(res$elapsed), "flags", res$metrics$flag_count); res
    }
    a <- perform("initial"); initial_flags <- a$metrics$flag_count
    if (initial_flags > 0L) a <- perform("long")
    data.frame(job, key = key, phase = a$phase, initial_flag_count = initial_flags, a$metrics,
               file = file.path(out, paste0(key, "-", a$phase, ".rds")))
  }
  results <- if (workers == 1L) lapply(seq_len(nrow(jobs)), work) else
    parallel::mclapply(seq_len(nrow(jobs)), work, mc.cores = min(workers, 2L), mc.preschedule = FALSE, mc.set.seed = FALSE)
  stopifnot(!any(vapply(results, inherits, logical(1), "try-error")))
  table <- do.call(rbind, results); table$result_md5 <- unname(tools::md5sum(table$file))
  stopifnot(nrow(table) == 96L)
  publish_csv(table, file.path(out, "oracle-selected.csv"), drop_paths = "file")
  say("oracle finished; 96 posteriors;", sum(table$flag_count > 0), "still flagged after the doubled rerun")
  quit(status = 0)
}
# ---- freeze ---------------------------------------------------------------
if (mode == "freeze") {
  stopifnot(!dir.exists(file.path(study, "library")))
  dirty <- system2("git", c("-C", shQuote(repo), "status", "--porcelain", "--", "R", "src", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)
  stopifnot(length(dirty) == 0L)
  rev <- system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)
  src <- file.path(study, "source-main"); dir.create(src)
  for (d in c("R", "src", "man", "data", "inst")) if (dir.exists(file.path(repo, d))) file.copy(file.path(repo, d), src, recursive = TRUE)
  file.copy(file.path(repo, c("DESCRIPTION", "NAMESPACE")), src)
  unlink(list.files(file.path(src, "src"), pattern = "\\.(o|so|dll)$", full.names = TRUE))
  dir.create(file.path(study, "library"))
  status <- system2("R", c("CMD", "INSTALL", "--preclean", paste0("--library=", shQuote(file.path(study, "library"))), shQuote(src)),
                    stdout = file.path(study, "install.log"), stderr = file.path(study, "install.log"))
  stopifnot(status == 0L)
  writeLines(rev, file.path(study, "source-revision.txt"))
  cat(sprintf("- **%s, amendment 2.** Production code frozen at main revision `%s` and installed into the study library before the first full fit.\n",
              format(Sys.Date(), "%d %B %Y"), rev), file = file.path(scripts, "PLAN.md"), append = TRUE)
  say("Frozen revision", rev, "installed into", file.path(study, "library"))
  quit(status = 0)
}

# ---- fits -----------------------------------------------------------------
.libPaths(c(file.path(study, "library"), .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == normalizePath(file.path(study, "library/occJSDM")))
RcppParallel::setThreadOptions(numThreads = 1)
source(file.path(scripts, "score.R"))
production <- c(file.path(study, "source-main", c("DESCRIPTION", "NAMESPACE")),
                list.files(file.path(study, "source-main/R"), pattern = "\\.R$", full.names = TRUE),
                list.files(file.path(study, "source-main/src"), pattern = "\\.(cpp|h)$|^Makevars", full.names = TRUE),
                file.path(find.package("occJSDM"), "libs/occJSDM.so"))
fit_hashes <- tools::md5sum(c(production, file.path(scripts, c("generator.R", "run.R"))))
score_hash <- score_source_hashes(repo)
jobs <- split(manifest, manifest$key)[manifest$key]
keys <- option("keys", "")
if (nzchar(keys)) { keys <- strsplit(keys, ",", fixed = TRUE)[[1]]; stopifnot(all(keys %in% names(jobs))); jobs <- jobs[keys] } else
if (mode == "pilot") jobs <- jobs["rep01-clustered-two_stage"] else
if (mode == "long") stop("Long runs require explicit --keys from summarise.R --mode=select")
mcmc <- switch(mode, pilot = list(nchain = 2L, nburn = 40L, niter = 60L, nthin = 1L),
               initial = list(nchain = 2L, nburn = 3000L, niter = 5000L, nthin = 1L),
               long = list(nchain = 4L, nburn = 6000L, niter = 12000L, nthin = 1L))
out <- file.path(study, mode); dir.create(out, showWarnings = FALSE)
settings <- list(revision = readLines(file.path(study, "source-revision.txt")), fit_hashes = fit_hashes,
                 score_hash = score_hash, mcmc = mcmc, mode = mode, workers = workers, threads_per_fit = 1L,
                 priors = "Unchanged defaults", session = sessionInfo())
settings_file <- file.path(out, "settings.rds")
if (file.exists(settings_file)) { old <- readRDS(settings_file)
  stopifnot(identical(old$fit_hashes, fit_hashes), identical(old$mcmc, mcmc), identical(old$score_hash, score_hash)) } else atomic(settings, settings_file)
run_job <- function(job) {
  input_file <- inputs[job$community]
  stopifnot(identical(unname(tools::md5sum(input_file)), job$input_md5))
  dest <- file.path(out, paste0(job$key, "-result.rds")); fitfile <- file.path(out, paste0(job$key, "-fit.rds"))
  if (file.exists(dest)) { old <- readRDS(dest); stopifnot(identical(old$job, job), identical(old$mcmc, mcmc)); return(job$key) }
  input <- readRDS(input_file); data <- input$surveys[[job$arrangement]][[job$arm]]
  say(job$key, "started")
  if (file.exists(fitfile)) { saved <- readRDS(fitfile); stopifnot(identical(saved$job, job), identical(saved$mcmc, mcmc))
    fit <- saved$fit; warnings <- saved$warnings; started <- saved$started; finished <- saved$finished } else {
    set.seed(job$fit_seed); warnings <- character(); started <- Sys.time()
    fit <- withCallingHandlers(suppressMessages(occJSDM::runOccJSDM(data,
      listParams = list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = 100L), listPriors = list(), threshold = 1,
      occCovariates = "environment", collCovariates = if (job$arm == "binary") NULL else "collection",
      spatCovariates = c("longitude", "latitude"), MCMCparams = mcmc)),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    finished <- Sys.time()
    atomic(list(fit = fit, job = job, mcmc = mcmc, fit_hashes = fit_hashes, warnings = warnings, started = started, finished = finished), fitfile)
  }
  scores <- score_sweep_fit(fit, input, job$arrangement, job$arm)
  atomic(c(list(job = job, mcmc = mcmc, fit_hashes = fit_hashes, score_hash = score_hash, warnings = warnings,
                started = started, finished = finished), scores), dest)
  say(job$key, "complete; fit seconds", round(as.numeric(difftime(finished, started, units = "secs"))),
      "; warnings", length(warnings), "; max group Rhat", round(max(scores$groups$rhat, na.rm = TRUE), 3))
  job$key
}
work <- function(job) tryCatch(run_job(job), error = function(e) { say(job$key, "ERROR:", conditionMessage(e)); list(key = job$key, error = conditionMessage(e)) })
status <- if (workers == 1L) lapply(jobs, work) else parallel::mclapply(jobs, work, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
atomic(status, file.path(out, paste0("status-", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")))
stopifnot(all(vapply(status, is.character, logical(1))))
say("Completed", length(status), mode, "fits.")
