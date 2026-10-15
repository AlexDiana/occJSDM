# Spatial range and amplitude posteriors for binary replicates 1 to 5 at one
# commit. Usage: Rscript spatial-diag.R <tree> <out.rds>
a <- commandArgs(TRUE)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
suppressMessages(devtools::load_all(a[1], quiet = TRUE))
RcppParallel::setThreadOptions(numThreads = 1)
source(file.path(a[1], "tests", "testthat", "helper-simstudy.R"))
sc <- Filter(function(s) s$label == "binary", simstudy_scenarios())[[1]]
out <- lapply(1:5, function(r) {
  truth <- draw_truth(sc, simstudy_seed_for(sc, r))
  sim <- simstudy_simulate(truth)
  fit <- simstudy_fit(sim, sc)
  jo <- fit$results_output$jsdm_output
  grid <- fit$infos$l_s_grid
  ls <- grid[as.vector(jo$idx_ls_output)]
  b0 <- jo$B0_output
  list(replicate = r, l_s = ls, l_s_grid = grid,
       sigma_bs = as.vector(jo$sigmabs_output),
       B0_sd = apply(matrix(b0, dim(b0)[1]), 1, sd),
       Bs_sd = sd(as.vector(jo$Bs_output)))
})
saveRDS(out, a[2])

# This script runs in a plain Rscript process, so set.seed() in draw_truth()
# uses the default Mersenne-Twister generator rather than the L'Ecuyer streams
# that run_study.R sets on its workers. Its datasets therefore follow the
# binary cell's design but are not the production replicates; the comparison
# between commits is still paired, on the same five datasets.
#
# spatial-diagnostic-binary.csv was summarised from the two output files with:
#
# out <- do.call(rbind, lapply(c("54b5cca", "f1d2e12"), function(t) {
#   d <- readRDS(paste0("diag-", t, ".rds"))
#   do.call(rbind, lapply(d, function(x) data.frame(
#     commit = t, scenario = "binary", replicate = x$replicate,
#     l_s_mean = mean(x$l_s), l_s_sd = sd(x$l_s),
#     l_s_share_at_grid_max = mean(x$l_s == max(x$l_s_grid)),
#     sigma_bs_mean = mean(x$sigma_bs), B0_posterior_sd_mean = mean(x$B0_sd))))
# }))
