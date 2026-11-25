# Checks that a grouped fit with slab_scale = Inf, and a fixed-prior fit, draw
# exactly the random numbers they drew at the baseline commit a03a1b8. Fits four
# small simulated data sets (2 chains of 20 + 20) with the package loaded from
# --lib and saves, for each, B_output, sigmab_output, sigmaemb_output and
# lambda_output:
#   1. a binary fit under the fixed prior (spatial, as the tier-1 fixture),
#   2. a grouped fit of model type "binary",
#   3. a grouped fit of model type "occupancy",
#   4. a grouped fit of model type "two_stage".
# Run against a library built from a03a1b8 (no --slab) and one built from this
# branch with --slab Inf; all four fits must be identical. With --slab Inf the
# grouped fits get listPriors = list(slab_scale = Inf); the fixed-prior fit
# never gets slab_scale. The script stops if a grouped fit is not grouped.
#
# Usage: Rscript --vanilla identity.R --lib PATH --out FILE [--slab Inf]
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
slab <- if (is.na(match("--slab", args))) NULL else as.numeric(arg("slab"))
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
  priors <- if (is.null(slab)) list() else list(slab_scale = slab)
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
    listPriors = priors,
    occCovariates = grouped_occ_covariates(dat), MCMCparams = MCMC)))
}

seeded <- function(sim, fit, grouped) {
  set.seed(11L); occJSDM:::setOccJSDMSeed(11L)
  jd <- fit(sim)$results_output$jsdm_output
  if (grouped && (length(jd$sigmaemb_output) == 0L || length(jd$lambda_output) == 0L ||
                  anyNA(jd$sigmaemb_output)))
    stop("a fit meant to be grouped has no embedding scales")
  if (!grouped && !is.null(jd$sigmaemb_output)) stop("the fixed-prior fit is grouped")
  jd[c("B_output", "sigmab_output", "sigmaemb_output", "lambda_output")]
}

res <- list(
  fixed_binary = seeded(simulate_fixture(model = "binary"), fit_fixture, FALSE),
  grouped_binary = seeded(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L),
                          fit_grouped, TRUE),
  grouped_occupancy = seeded(simulate_grouped_fixture(model = "occupancy", p_emb = 6L, S = 8L),
                             fit_grouped, TRUE),
  grouped_two_stage = seeded(simulate_grouped_fixture(model = "two_stage", p_emb = 6L, S = 8L),
                             fit_grouped, TRUE))
for (n in names(res)) cat(n, ": sigmaemb_output length", length(res[[n]]$sigmaemb_output),
                          ", lambda_output length", length(res[[n]]$lambda_output),
                          ", B_output length", length(res[[n]]$B_output), "\n")
saveRDS(res, out)
cat("saved", out, "\n")
