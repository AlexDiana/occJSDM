# Fits that must not change between the baseline and the prototype: every fixed-prior fit, and
# grouped fits where the axis scale move is off. Plus one grouped binary fit as a positive control
# (it must differ). usage: Rscript identity.R <lib> <worktree> <out.rds>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
for (h in c("helper-fixtures.R", "helper-grouped-prior.R")) source(file.path(a[2], "tests/testthat", h))
mc <- list(nchain = 2, nburn = 40, niter = 40, nthin = 1)
grouped <- function(model, embedding = TRUE, priors = list(), seed = 5L, ...) {
  sim <- simulate_grouped_fixture(model = model, p_emb = 6L, S = 8L, ...); dat <- sim$data_list
  set.seed(seed)
  fit <- suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = c(list(n_factors = 1L), if (embedding) list(embeddingVars = dat$embeddingVars)),
    listPriors = priors, occCovariates = grouped_occ_covariates(dat), spatCovariates = NULL, MCMCparams = mc)))
  list(jsdm = fit$results_output$jsdm_output, coef_prior = fit$infos$coef_prior$type, axis = !is.null(fit$infos$axis_move))
}
out <- list(
  fixed_binary = grouped("binary", embedding = FALSE),
  fixed_occupancy = grouped("occupancy", embedding = FALSE),
  fixed_two_stage = grouped("two_stage", embedding = FALSE),
  fixed_continuous = grouped("continuous", embedding = FALSE),
  grouped_occupancy = grouped("occupancy"),
  grouped_two_stage = grouped("two_stage"),
  grouped_continuous = grouped("continuous"),
  grouped_binary_no_embedding = grouped("binary", embedding = FALSE, priors = list(b_prior = "grouped")),
  grouped_binary = grouped("binary"))
saveRDS(out, a[3])
