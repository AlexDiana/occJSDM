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
stop("Mode ", mode, " is added in Task 3")
