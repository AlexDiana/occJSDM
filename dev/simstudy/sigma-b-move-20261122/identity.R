# Checks that fits without the sigma_b move draw exactly the random numbers
# they drew before it. Fits three small simulated data sets with the package
# loaded from --lib and saves their B_output arrays:
#   1. a binary fit under the fixed prior (spatial, as the tier-1 fixture),
#   2. a grouped fit of model type "occupancy" (no scale moves),
#   3. a grouped binary fit (the moves, and with this branch the sigma_b move).
# Run against a library built from a2a64a1 and one built from this branch,
# 1 and 2 must be identical and 3 must differ.
#
# Usage: Rscript --vanilla identity.R --lib PATH --out FILE
# The simulator calls are copied from tests/testthat/helper-fixtures.R and
# tests/testthat/helper-grouped-prior.R.

Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) stop("missing --", name)
  args[i + 1L]
}
lib <- normalizePath(arg("lib"), mustWork = TRUE)
out <- arg("out")
library(occJSDM, lib.loc = lib)
if (normalizePath(dirname(system.file(package = "occJSDM"))) != lib)
  stop("occJSDM was not loaded from ", lib)
RcppParallel::setThreadOptions(numThreads = 1)

MCMC <- list(nchain = 2, nburn = 20, niter = 20, nthin = 1)

# From helper-fixtures.R: simulate_fixture() and fit_fixture().
simulate_fixture <- function(model = "two_stage", n = 40L, S = 4L, M = 2L, P = 1L, K = 2L,
                             g = 2L, gt = 2L, d = 2L, useSpatField = TRUE, seed = 1L) {
  set.seed(seed)
  M_vec <- if (model %in% c("occupancy", "two_stage")) rep(M, n) else rep(1L, n)
  N     <- sum(M_vec)
  K_vec <- if (model == "two_stage") rep(K, N * P) else rep(1L, N * P)
  datasettings <- list(n = n, S = S, g = g, M = M_vec, P = P, K = K_vec,
                       ncov_psi = 2L, ncov_theta = 1L)
  jsdmParams <- list(gt = gt, d = d, ds = 0,
                     sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                     sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
                     tau = rep(1, S), useSpatField = useSpatField)
  params <- list(p              = matrix(stats::runif(P * S, 0.3, 0.6), P, S),
                 q              = matrix(stats::runif(P * S, 0.01, 0.05), P, S),
                 theta0         = stats::runif(S, 0.02, 0.1),
                 theta_baseline = stats::runif(S, 0.2, 0.5))
  simulateOccJSDMData(datasettings, params, jsdmParams, model = model)
}
fit_fixture <- function(sim) {
  dat <- sim$data_list
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = list(n_factors = 2L, n_supportpoints = 8L),
    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
    spatCovariates = c("Xs.1", "Xs.2"), MCMCparams = MCMC)))
}

# From helper-grouped-prior.R: simulate_grouped_fixture() and
# grouped_occ_covariates().
simulate_grouped_fixture <- function(model = "binary", n = 60L, S = 12L, p_meas = 2L,
                                     p_emb = 8L, gt = 2L, seed = 1L, embedding = list(),
                                     M = 2L, K = 2L) {
  set.seed(seed)
  replicated <- model %in% c("occupancy", "two_stage")
  M_vec <- if (replicated) rep(M, n) else rep(1L, n)
  K_vec <- if (model == "two_stage") rep(K, sum(M_vec)) else rep(1L, sum(M_vec))
  datasettings <- list(n = n, S = S, g = 2L, M = M_vec, P = 1L,
                       K = K_vec, ncov_psi = p_meas, ncov_theta = 1L)
  jsdmParams <- list(gt = gt, d = 1L, ds = 0,
                     sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                     sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
                     tau = rep(1, S), useSpatField = FALSE,
                     embedding = modifyList(list(n_emb = p_emb), embedding))
  params <- list()
  if (replicated)
    params <- list(p              = matrix(stats::runif(S, 0.3, 0.6), 1L, S),
                   q              = matrix(stats::runif(S, 0.01, 0.05), 1L, S),
                   theta0         = stats::runif(S, 0.02, 0.1),
                   theta_baseline = stats::runif(S, 0.2, 0.5))
  simulateOccJSDMData(datasettings, params, jsdmParams, model = model)
}
grouped_occ_covariates <- function(dat) {
  c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars),
    dat$embeddingVars)
}
fit_grouped <- function(sim) {
  dat <- sim$data_list
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
    occCovariates = grouped_occ_covariates(dat), MCMCparams = MCMC)))
}

seeded <- function(sim, fit) {
  set.seed(11L); occJSDM:::setOccJSDMSeed(11L)
  fit(sim)$results_output$jsdm_output$B_output
}

res <- list(
  fixed_binary = seeded(simulate_fixture(model = "binary"), fit_fixture),
  grouped_occupancy = seeded(simulate_grouped_fixture(model = "occupancy", p_emb = 6L, S = 8L), fit_grouped),
  grouped_binary = seeded(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L), fit_grouped))
saveRDS(res, out)
cat("saved", out, "\n")
