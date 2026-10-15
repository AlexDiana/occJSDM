# Refit base and binary, replicates 1 to 20, at the study's code (de4389c plus
# the calibration record) with the experimental half-Cauchy spatial amplitude
# prior. Labels are unchanged, so the seeds, data and truth are those of the
# production cells. Usage: Rscript halfcauchy.R <tree> <out.rds>
a <- commandArgs(TRUE)
library(parallel)
suppressMessages(devtools::load_all(a[1], quiet = TRUE))
cl <- makeCluster(5)
clusterExport(cl, "a")
# As run_study.R does: L'Ecuyer streams on the workers, so that set.seed() in
# draw_truth() reproduces the production datasets.
clusterSetRNGStream(cl, 20260727L)
clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  suppressMessages(devtools::load_all(a[1], quiet = TRUE))
  RcppParallel::setThreadOptions(numThreads = 1)
  source(file.path(a[1], "tests", "testthat", "helper-simstudy.R"))
  NULL
})
jobs <- expand.grid(cell = c("base", "binary"), r = 1:20, stringsAsFactors = FALSE)
res <- parLapplyLB(cl, seq_len(nrow(jobs)), function(i, jobs) {
  sc <- Filter(function(s) s$label == jobs$cell[i], simstudy_scenarios())[[1]]
  sc$listPriors <- list(sigma_bs_prior = "half_cauchy")
  simstudy_replicate(sc, jobs$r[i], list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1))$params
}, jobs = jobs)
stopCluster(cl)
saveRDS(list(rows = do.call(rbind, res)), a[2])
