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
