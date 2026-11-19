# A small grouped binary fit with 5,000 burn-in sweeps, where the hold applies. usage: <lib> <worktree> <out.rds>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
for (h in c("helper-fixtures.R", "helper-grouped-prior.R")) source(file.path(a[2], "tests/testthat", h))
sim <- simulate_grouped_fixture(model = "binary", p_emb = 6L, S = 8L); dat <- sim$data_list
set.seed(5L)
fit <- suppressMessages(suppressWarnings(runOccJSDM(dat, listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
  occCovariates = grouped_occ_covariates(dat), spatCovariates = NULL, MCMCparams = list(nchain = 2, nburn = 5000, niter = 50, nthin = 1))))
saveRDS(list(jsdm = fit$results_output$jsdm_output, hold = fit$infos$axis_move$hold_sweeps), a[3])
