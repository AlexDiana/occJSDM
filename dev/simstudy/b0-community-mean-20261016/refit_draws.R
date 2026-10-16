# Refit replicates 1 to 10 of the six community-mean arms and their six
# default-prior counterparts at the study's settings, keeping the draws that
# the runner discards, for the mixing check (README.md, C12 and C13).
# Usage: Rscript refit_draws.R <package root> <output directory> <cores>
# Each job saves <output directory>/<arm>-<replicate>.rds as it finishes, so a
# killed run keeps its finished jobs; a rerun skips those. Workers are set up
# as dev/simstudy/run_study.R sets up its own: clusterSetRNGStream(cl,
# 20260727L), one RcppParallel thread, then devtools::load_all(). Under that
# L'Ecuyer generator the set.seed() in draw_truth() gives the runner's datasets
# and chain seeds, so each refit should reproduce its arm's rows exactly.
args <- commandArgs(trailingOnly = TRUE)
pkg_root <- normalizePath(args[1], mustWork = TRUE)
out_dir <- args[2]
cores <- as.integer(args[3])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
library(parallel)

arms <- c("base_cm", "occupancy_cm", "binary_cm",
          "base_rare_cm", "occupancy_rare_cm", "binary_rare_cm",
          "base", "occupancy", "binary", "base_rare", "occupancy_rare", "binary_rare")
jobs <- expand.grid(replicate = 1:10, arm = arms, stringsAsFactors = FALSE)
jobs$file <- file.path(out_dir, sprintf("%s-%02d.rds", jobs$arm, jobs$replicate))
jobs <- jobs[!file.exists(jobs$file), ]
message(sprintf("%s: %d jobs to run", format(Sys.time()), nrow(jobs)))

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
  fit <- simstudy_fit(sim, s)
  ro <- fit$results_output
  jo <- ro$jsdm_output
  res <- list(
    arm = job$arm, replicate = job$replicate,
    rows = simstudy_rows(fit, sim, truth, s, job$replicate),
    occupancy = simstudy_occupancy_rows(fit, sim, s, job$replicate),
    draws = list(mu_b0 = jo$mu_b0_output, B0 = jo$B0_output,
                 beta_theta0 = if (is.null(ro$beta_theta_output)) NULL else
                   ro$beta_theta_output[1, , , , drop = TRUE],
                 theta0 = ro$theta0_output, p = ro$p_output),
    intercept_prior = fit$infos$intercept_prior,
    elapsed = proc.time()[["elapsed"]] - started)
  saveRDS(res, job$file)
  message(sprintf("%s %s replicate %d in %.0f s", format(Sys.time(), "%H:%M:%S"),
                  job$arm, job$replicate, res$elapsed))
  sprintf("%s %s replicate %d in %.0f s", format(Sys.time(), "%H:%M:%S"),
          job$arm, job$replicate, res$elapsed)
}
clusterExport(cl, c("jobs", "one"))
done <- parLapplyLB(cl, seq_len(nrow(jobs)), function(k) {
  tryCatch(one(k), error = function(e) paste("FAILED", jobs$arm[k], jobs$replicate[k],
                                             conditionMessage(e)))
})
message(paste(unlist(done), collapse = "\n"))
message(sprintf("%s: finished, %d files in %s", format(Sys.time()),
                length(list.files(out_dir, "[.]rds$")), out_dir))
stopCluster(cl)
