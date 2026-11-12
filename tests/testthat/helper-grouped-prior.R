# Helpers for the grouped covariate prior tests.
seed_both <- function(s) { set.seed(s); occJSDM:::setOccJSDMSeed(as.integer(s)) }

fixed_prior_reference_fit <- function(model) {
  i <- match(model, c("binary", "continuous", "occupancy", "two_stage"))
  spatial <- model %in% c("binary", "two_stage")   # covers both range functions
  sim <- simulate_fixture(model = model, useSpatField = spatial, seed = 100L + i)
  set.seed(200L + i)
  fit <- fit_fixture(sim, spatial = spatial)
  fit$results_output$jsdm_output[c("B0_output", "B_output", "G_output", "A_output",
                                   "C_output", "L_output", "sigmab_output")]
}

# Small simulated dataset with an embedding block beside the measured covariates.
# gt is the simulator's number of latent traits, p_meas its ncov_psi, d one factor.
simulate_grouped_fixture <- function(model = "binary", n = 60L, S = 12L, p_meas = 2L,
                                     p_emb = 8L, gt = 2L, seed = 1L, embedding = list()) {
  set.seed(seed)
  datasettings <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L,
                       K = rep(1L, n), ncov_psi = p_meas, ncov_theta = 1L)
  jsdmParams <- list(gt = gt, d = 1L, ds = 0,
                     sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                     sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
                     tau = rep(1, S), useSpatField = FALSE,
                     embedding = modifyList(list(n_emb = p_emb), embedding))
  simulateOccJSDMData(datasettings, list(), jsdmParams, model = model)
}

grouped_occ_covariates <- function(dat) {
  c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars),
    dat$embeddingVars)
}
