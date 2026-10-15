# Fitted spatial range and amplitude for replicates 1 to 10 of the two
# production parents and the four range-inside-grid cells. The runner does not
# save them, so this refits those replicates exactly as run_study.R does and
# keeps the range and amplitude draws.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/range-inside-grid-20261015/range-diag.R <out.rds>
#
# As run_study.R does: L'Ecuyer streams on the workers, so set.seed() in
# draw_truth() reproduces the production datasets, and the MCMC stream that
# follows is the runner's too. Each refit is therefore the same fit as the
# matching replicate in the run (single thread, same code and BLAS); the
# summary checks this by comparing B and B0 posterior means with the run's.
a <- commandArgs(TRUE)
stopifnot(length(a) == 1, !file.exists(a[1]))
library(parallel)
pkg <- normalizePath(".")
cl <- makeCluster(5)
clusterExport(cl, "pkg")
clusterSetRNGStream(cl, 20260727L)
clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  RcppParallel::setThreadOptions(numThreads = 1)
  suppressMessages(devtools::load_all(pkg, quiet = TRUE))
  source(file.path(pkg, "tests", "testthat", "helper-simstudy.R"))
  NULL
})
cells <- c("base", "binary", "base_ls05", "binary_ls05", "base_ls075", "binary_ls075")
jobs <- expand.grid(cell = cells, r = 1:10, stringsAsFactors = FALSE)
res <- parLapplyLB(cl, seq_len(nrow(jobs)), function(i, jobs) {
  sc <- Filter(function(s) s$label == jobs$cell[i], simstudy_scenarios())[[1]]
  r <- jobs$r[i]
  truth <- draw_truth(sc, simstudy_seed_for(sc, r))
  sim <- simstudy_simulate(truth)
  fit <- simstudy_fit(sim, sc, list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1))
  jo <- fit$results_output$jsdm_output
  grid <- fit$infos$l_s_grid
  info <- sim$data_list$info
  if ("Site" %in% names(info)) info <- info[!duplicated(info$Site), ]
  sds <- apply(as.matrix(info[c("Xs.1", "Xs.2")]), 2, sd)
  list(scenario = sc$label, replicate = r, l_s_raw = truth$jsdmParams$l_s,
       l_s_std = truth$jsdmParams$l_s / sds, l_s_grid = grid,
       l_s_draws = grid[as.vector(jo$idx_ls_output)],
       sigma_bs_draws = as.vector(jo$sigmabs_output),
       rows = simstudy_rows(fit, sim, truth, sc, r))
}, jobs = jobs)
stopCluster(cl)
saveRDS(res, a[1])
