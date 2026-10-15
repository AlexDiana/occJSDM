# Fitted spatial amplitude, range and spatial share for replicates 1 to 10 of
# base and binary, at the production range and at 0.075 raw, each under the
# default InvGamma(10, 1) amplitude prior and under the half-Cauchy prior
# (8 cells, 80 fits). The runner saves none of these, so this refits those
# replicates exactly as run_study.R does and keeps them.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/amplitude-prior-20261015/amp-diag.R <out.rds>
#
# As run_study.R does: L'Ecuyer streams on the workers, so set.seed() in
# draw_truth() reproduces the production datasets and the MCMC stream that
# follows is the runner's. Each refit is therefore the same fit as the matching
# replicate in the corresponding run; amp-diag-summary.R checks this.
a <- commandArgs(TRUE)
stopifnot(length(a) == 1, !file.exists(a[1]))
library(parallel)
pkg <- normalizePath(".")
cl <- makeCluster(5)
clusterExport(cl, "pkg")
clusterSetRNGStream(cl, 20260727L)
invisible(clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  RcppParallel::setThreadOptions(numThreads = 1)
  suppressMessages(devtools::load_all(pkg, quiet = TRUE))
  source(file.path(pkg, "tests", "testthat", "helper-simstudy.R"))
  NULL
}))
cells <- c("base", "base_hc", "base_ls075", "base_ls075_hc",
           "binary", "binary_hc", "binary_ls075", "binary_ls075_hc")
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
  jp <- sim$true_params$jsdmParams_true
  # varPart_output is S x 4 x niter x nchain; column 2 is the spatial share.
  vp <- jo$varPart_output
  list(scenario = sc$label, replicate = r,
       prior = sc$listPriors$sigma_bs_prior %||% "inverse_gamma",
       l_s_raw = truth$jsdmParams$l_s, l_s_std = truth$jsdmParams$l_s / sds,
       l_s_grid = grid,
       l_s_draws = grid[as.vector(jo$idx_ls_output)],
       sigma_bs_draws = as.vector(jo$sigmabs_output),
       spatial_share_fit = apply(vp[, 2, , , drop = FALSE], 1, mean),
       spatial_share_true = jp$varPart[, "Spatial"],
       field_sd_true = sd(as.vector(jp$spatField)),
       rows = simstudy_rows(fit, sim, truth, sc, r))
}, jobs = jobs)
stopCluster(cl)
saveRDS(res, a[1])
