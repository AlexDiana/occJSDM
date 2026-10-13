# Time and measure the quickstart fit at several core budgets, twice each, and
# check identity against the in-process fit. Run detached through launch.sh.
# Rscript benchmark.R OUTDIR [LABEL ...]
# With no labels all configurations run; otherwise only the named ones do, and
# an unknown label stops the run before any fit.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L)
outdir <- normalizePath(args[[1]], mustWork = TRUE)
lib <- file.path(outdir, "library")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"),
          nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")), nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS_BASIS")))
writeLines(Sys.getenv("OCCJSDM_REFERENCE_BLAS_BASIS"), file.path(outdir, "blas-declaration-basis.txt"))
sampler <- file.path(outdir, "scripts", "sample-rss.sh")
reference_fit <- file.path(outdir, "reference-fit.rds")       # kept on disk, never resident during a measured fit
reference_rng <- file.path(outdir, "reference-rng.rds")
# A run directory is used once: a rerun must not compare a new build against an
# earlier run's reference or inherit its sampler stop markers.
stopifnot(!file.exists(reference_fit), !file.exists(reference_rng), !file.exists(file.path(outdir, "results.csv")),
          length(list.files(file.path(outdir, "rss"))) == 0L)
configs <- list(
  list(label = "serial-ambient", ncores = 1L, pin_one_thread = FALSE),
  list(label = "serial-1thread", ncores = 1L, pin_one_thread = TRUE),
  list(label = "ncores-2", ncores = 2L, pin_one_thread = FALSE),
  list(label = "ncores-4", ncores = 4L, pin_one_thread = FALSE),
  list(label = "ncores-8", ncores = 8L, pin_one_thread = FALSE))
requested <- args[-1L]
if (length(requested)) {
  known <- vapply(configs, function(cfg) cfg$label, "")
  unknown <- setdiff(requested, known)
  if (length(unknown)) stop("unknown configuration label(s): ", paste(unknown, collapse = ", "), "; known labels: ", paste(known, collapse = ", "))
  configs <- Filter(function(cfg) cfg$label %in% requested, configs)
}
fit_once <- function(ncores) {
  set.seed(20261007)
  fit <- runOccJSDM(data = sampledata, listParams = list(n_factors = 2), threshold = 1,
                    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
                    collCovariates = c("X_theta.1", "X_theta.2"),
                    spatCovariates = c("Xs.1", "Xs.2"),
                    MCMCparams = list(nchain = 4, nburn = 5000, niter = 5000, nthin = 1, ncores = ncores),
                    summarisedLatentPresences = TRUE)
  list(fit = fit, seed_state = get(".Random.seed", envir = globalenv()), draws = list(rnorm(3), runif(3), sample(100, 3)))
}
rows <- list()
for (cfg in configs) for (rep in 1:2) {
  tag <- sprintf("%s-rep%d", cfg$label, rep)
  if (cfg$pin_one_thread) Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")   # the variable is all RcppParallel consults
  rss <- file.path(outdir, "rss", paste0(tag, ".csv"))
  stopifnot(!file.exists(rss), !file.exists(paste0(rss, ".stop")))
  system2("sh", c(sampler, Sys.getpid(), shQuote(rss)), wait = FALSE)
  t0 <- Sys.time()
  res <- suppressMessages(suppressWarnings(fit_once(cfg$ncores)))
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  file.create(paste0(rss, ".stop"))               # ends this configuration's sampler before anything else is loaded
  Sys.sleep(2)
  if (cfg$pin_one_thread) Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
  if (!file.exists(reference_fit)) {
    saveRDS(occJSDM:::drop_execution(res$fit), reference_fit)
    saveRDS(list(seed_state = res$seed_state, draws = res$draws), reference_rng)
  }
  reference <- readRDS(reference_fit); rng_ref <- readRDS(reference_rng)   # outside the measured interval
  identical_fit <- identical(occJSDM:::drop_execution(res$fit), reference)
  identical_rng <- identical(res$seed_state, rng_ref$seed_state) && identical(res$draws, rng_ref$draws)
  rm(reference, rng_ref)
  ex <- res$fit$infos$execution
  # Memory results are reported only when the samples exist: the parent must
  # have been sampled, and every worker the fit records must appear.
  sampled <- read.csv(rss)
  if (!any(sampled$role == "parent")) stop("no parent memory samples for ", tag)
  missing_workers <- setdiff(ex$worker_pids, sampled$pid[sampled$role == "worker"])
  if (length(missing_workers)) stop("no memory samples for worker(s) ", paste(missing_workers, collapse = ", "), " in ", tag)
  rows[[tag]] <- data.frame(tag = tag, label = cfg$label, rep = rep, ncores = cfg$ncores, workers = ex$workers,
                            threads_per_worker = ex$threads_per_worker, wall_s = wall,
                            chain_s = paste(round(ex$elapsed, 1), collapse = ";"),
                            identical_fit = identical_fit, identical_rng = identical_rng,
                            blas = ex$blas, stringsAsFactors = FALSE)
  saveRDS(res$fit$infos$execution, file.path(outdir, "execution", paste0(tag, ".rds")))
  write.csv(do.call(rbind, rows), file.path(outdir, "results.csv"), row.names = FALSE)
  cat("DONE", tag, round(wall, 1), "s identical_fit", identical_fit, "\n")
  rm(res); invisible(gc())
}
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
cat("ALL DONE\n")
