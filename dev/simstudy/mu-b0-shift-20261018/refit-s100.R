# dev/simstudy/theta-community-mean-20261016/refit.R as it stands, with one
# change: the fit's infos$b0_shift (the shift step's proposal sd and
# acceptance rates per chain) is saved with the draws. It refits replicates
# of named arms, keeping the draws that the runner discards, for the rerun of
# the 100-species mixing check with the shift step on (README.md).
# Usage: Rscript refit.R <package root> <output directory> <cores> <arms> <replicates> <nburn> <niter>
#   arms: comma-separated labels; replicates: first-last, for example 1-10.
# Each job saves <output directory>/<arm>-<replicate>.rds as it finishes, so a
# killed run keeps its finished jobs; a rerun skips those. Workers are set up
# as dev/simstudy/run_study.R sets up its own: clusterSetRNGStream(cl,
# 20260727L), one RcppParallel thread, then devtools::load_all(). Under that
# L'Ecuyer generator the set.seed() in draw_truth() gives the runner's datasets
# and chain seeds, so at the study's chain length each refit should reproduce
# its arm's stored rows exactly.
args <- commandArgs(trailingOnly = TRUE)
pkg_root <- normalizePath(args[1], mustWork = TRUE)
out_dir <- args[2]
cores <- as.integer(args[3])
arms <- strsplit(args[4], ",", fixed = TRUE)[[1]]
reps <- as.integer(strsplit(args[5], "-", fixed = TRUE)[[1]])
mcmc <- list(nchain = 2, nburn = as.integer(args[6]), niter = as.integer(args[7]), nthin = 1)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
library(parallel)

jobs <- expand.grid(replicate = seq(reps[1], reps[2]), arm = arms, stringsAsFactors = FALSE)
jobs$file <- file.path(out_dir, sprintf("%s-%02d.rds", jobs$arm, jobs$replicate))
jobs <- jobs[!file.exists(jobs$file), ]
message(sprintf("%s: %d jobs to run at %d + %d", format(Sys.time()), nrow(jobs), mcmc$nburn, mcmc$niter))

cl <- makeCluster(cores, outfile = "")
clusterSetRNGStream(cl, 20260727L)
clusterExport(cl, c("pkg_root"))
invisible(clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  RcppParallel::setThreadOptions(numThreads = 1)
  suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
  source(file.path(pkg_root, "tests", "testthat", "helper-simstudy.R"))
  NULL
}))

one <- function(k) {
  job <- jobs[k, ]
  sc <- simstudy_scenarios()
  s <- sc[[which(vapply(sc, `[[`, "", "label") == job$arm)]]
  started <- proc.time()[["elapsed"]]
  truth <- draw_truth(s, simstudy_seed_for(s, job$replicate))
  sim <- simstudy_simulate(truth)
  fit <- simstudy_fit(sim, s, MCMCparams = mcmc)
  ro <- fit$results_output
  jo <- ro$jsdm_output
  res <- list(
    arm = job$arm, replicate = job$replicate, mcmc = mcmc,
    rows = simstudy_rows(fit, sim, truth, s, job$replicate),
    occupancy = simstudy_occupancy_rows(fit, sim, s, job$replicate),
    draws = list(mu_b0 = jo$mu_b0_output, mu_theta = ro$mu_theta_output, B0 = jo$B0_output,
                 beta_theta0 = if (is.null(ro$beta_theta_output)) NULL else
                   ro$beta_theta_output[1, , , , drop = TRUE],
                 theta0 = ro$theta0_output, p = ro$p_output),
    intercept_prior = fit$infos$intercept_prior,
    collection_intercept_prior = fit$infos$collection_intercept_prior,
    b0_shift = fit$infos$b0_shift,
    elapsed = proc.time()[["elapsed"]] - started)
  saveRDS(res, job$file)
  msg <- sprintf("%s %s replicate %d in %.0f s", format(Sys.time(), "%H:%M:%S"),
                 job$arm, job$replicate, res$elapsed)
  message(msg)
  msg
}
clusterExport(cl, c("jobs", "one", "mcmc"))
done <- parLapplyLB(cl, seq_len(nrow(jobs)), function(k) {
  tryCatch(one(k), error = function(e) paste("FAILED", jobs$arm[k], jobs$replicate[k],
                                             conditionMessage(e)))
})
message(paste(unlist(done), collapse = "\n"))
message(sprintf("%s: finished, %d files in %s", format(Sys.time()),
                length(list.files(out_dir, "[.]rds$")), out_dir))
stopCluster(cl)
