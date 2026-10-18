# Refit fits of the 100-species mixing check (b0-community-s100-20261018)
# with the same chains as that run, storing every nthin-th draw together with
# the latent terms of the occupancy linear predictor, so that the partner of
# the mu_b0 shift can be found (README.md).
# Usage: Rscript refit.R <package root> <output directory> <cores> <jobs> <nburn> <niter> <nthin>
#   jobs: comma-separated arm:replicate pairs, for example
#         base_rare_cm_s100:3,base_rare_cm_s100:4,base_cm_s100:2
# Each job saves <output directory>/<arm>-<replicate>.rds as it finishes, so
# a killed run keeps its finished jobs; a rerun skips those. Workers are set
# up as dev/simstudy/theta-community-mean-20261016/refit.R sets up its own
# (clusterSetRNGStream(cl, 20260727L), one RcppParallel thread, then
# devtools::load_all()), and the fit is simstudy_fit()'s runOccJSDM() call
# with summarisedLatentPresences = FALSE added, so that z and psi are stored
# per draw. A chain runs nburn + niter * nthin iterations and stores every
# nthin-th draw after burn-in, so nburn = 5000, niter = 500, nthin = 10 is the
# main run's chain (5,000 + 5,000) with every 10th draw kept. Storage draws no
# random numbers (check-storage.R), so the stored draws equal the main run's
# at every 10th iteration; analyse.R checks that against the main run's draws.
args <- commandArgs(trailingOnly = TRUE)
pkg_root <- normalizePath(args[1], mustWork = TRUE)
out_dir <- args[2]
cores <- as.integer(args[3])
job_spec <- strsplit(strsplit(args[4], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
mcmc <- list(nchain = 2, nburn = as.integer(args[5]), niter = as.integer(args[6]),
             nthin = as.integer(args[7]))
run_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
library(parallel)

jobs <- data.frame(arm = vapply(job_spec, `[`, "", 1),
                   replicate = as.integer(vapply(job_spec, `[`, "", 2)),
                   stringsAsFactors = FALSE)
jobs$file <- file.path(out_dir, sprintf("%s-%02d.rds", jobs$arm, jobs$replicate))
jobs <- jobs[!file.exists(jobs$file), ]
message(sprintf("%s: %d jobs to run at %d + %d x %d", format(Sys.time()), nrow(jobs),
                mcmc$nburn, mcmc$niter, mcmc$nthin))

cl <- makeCluster(cores, outfile = "")
clusterSetRNGStream(cl, 20260727L)
clusterExport(cl, c("pkg_root", "run_dir"))
invisible(clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  RcppParallel::setThreadOptions(numThreads = 1)
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))
  source(file.path(run_dir, "terms.R"))
  NULL
}))

# simstudy_fit() with summarisedLatentPresences = FALSE added: the same call
# otherwise, so the chain seeds drawn inside runOccJSDM() are the main run's.
fit_keeping_latents <- function(sim, scenario, MCMCparams) {
  a <- simstudy_fit_args(sim, scenario)
  suppressMessages(suppressWarnings(
    runOccJSDM(a$data, listParams = a$listParams, listPriors = a$listPriors,
               occCovariates = a$occCovariates, collCovariates = a$collCovariates,
               spatCovariates = a$spatCovariates, MCMCparams = MCMCparams,
               summarisedLatentPresences = FALSE)))
}

one <- function(k) {
  job <- jobs[k, ]
  sc <- simstudy_scenarios()
  s <- sc[[which(vapply(sc, `[[`, "", "label") == job$arm)]]
  started <- proc.time()[["elapsed"]]
  truth <- draw_truth(s, simstudy_seed_for(s, job$replicate))
  sim <- simstudy_simulate(truth)
  fit <- fit_keeping_latents(sim, s, mcmc)
  # simstudy_occupancy_rows() expects the n x S means that
  # summarisedLatentPresences = TRUE returns.
  fit_means <- fit
  fit_means$results_output$psi_output <- apply(fit$results_output$psi_output, c(1, 2), mean)
  fit_means$results_output$z_output <- apply(fit$results_output$z_output, c(1, 2), mean)
  res <- c(
    list(arm = job$arm, replicate = job$replicate, mcmc = mcmc,
         rows = simstudy_rows(fit, sim, truth, s, job$replicate),
         occupancy = simstudy_occupancy_rows(fit_means, sim, s, job$replicate)),
    keep_fit(fit, sim),
    list(elapsed = proc.time()[["elapsed"]] - started))
  saveRDS(res, job$file)
  msg <- sprintf("%s %s replicate %d in %.0f s", format(Sys.time(), "%H:%M:%S"),
                 job$arm, job$replicate, res$elapsed)
  message(msg)
  msg
}
clusterExport(cl, c("jobs", "one", "mcmc", "fit_keeping_latents"))
done <- parLapplyLB(cl, seq_len(nrow(jobs)), function(k) {
  tryCatch(one(k), error = function(e) paste("FAILED", jobs$arm[k], jobs$replicate[k],
                                             conditionMessage(e)))
})
message(paste(unlist(done), collapse = "\n"))
message(sprintf("%s: finished, %d files in %s", format(Sys.time()),
                length(list.files(out_dir, "[.]rds$")), out_dir))
stopCluster(cl)
