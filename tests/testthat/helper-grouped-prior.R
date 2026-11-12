# Helpers for the grouped covariate prior tests.
seed_both <- function(s) { set.seed(s); occJSDM:::setOccJSDMSeed(as.integer(s)) }

# One of the four fits that the first test of test-grouped-prior.R compares
# with fixtures/fixed-prior-reference.rds. The saved fits were made by this
# function with the package loaded from a checkout of commit 4f2ee6c, the last
# one before the grouped prior. They must only ever be regenerated that way,
# never at the branch head, where the file would record the code under test;
# fixtures/make-fixed-prior-reference.R gives the steps.
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
# The replicated model types get the replicate counts and detection parameters
# that simulate_fixture() gives them: M samples per site for "occupancy" and
# "two_stage", K PCR replicates per sample for "two_stage", one primer. The
# detection parameters are drawn only for those two model types, so a "binary"
# or "continuous" call keeps the random stream it had before they were added.
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

# Inputs of the spatial range step for a simulated fixture. The spatial side
# is built as test-spatial-range.R builds it for its own references (and as
# runOccJSDM() does): the location covariates through create_covariates_matrix(),
# FIXTURE_KNOTS support points and the basis at three candidate ranges. Rows
# of one site share a location label and a basis row. The prior means, factor
# scores and Polya-Gamma weights are arbitrary draws, and kappa is the
# detection matrix minus one half.
range_reference_setup <- function(sim) {
  info <- sim$data_list$info
  n <- nrow(info); S <- ncol(sim$data_list$OTU); d <- 1L
  X <- unname(apply(as.matrix(info[, fixture_occ_covariates()]), 2,
                    function(v) (v - mean(v)) / stats::sd(v)))
  set.seed(31)
  Xs <- occJSDM:::create_covariates_matrix(info[, fixture_spat_covariates()],
                                           spline_vars = FALSE, remove_intercept = TRUE)$X
  list_Xs <- occJSDM:::computeSpatialSummaries(Xs, FIXTURE_KNOTS)
  ps <- list_Xs$ps
  list(X = X, U = matrix(rnorm(n * d), n, d),
       M_B = matrix(rnorm(ncol(X) * S, mean = .4), ncol(X), S),
       M_Bs = matrix(rnorm(ps * S, mean = -.3), ps, S), sigma_bs = 0.7,
       kappa = unname((sim$data_list$OTU > 0) - .5),
       Omega = matrix(rexp(n * S) / 4, n, S), Xs_centers = list_Xs$Xs_centers,
       list_SoRSummaries = suppressMessages(
         occJSDM:::precomputeSORmatrices(c(.05, .15, .25), list_Xs)),
       location = list_Xs$Xs_index)
}

# Dense reference for the range log-weights with one prior standard deviation
# per covariate: for each species 0.5 * h' Q^-1 h - 0.5 * log det Q, with
# Q = Z' diag(omega) Z + P and h = Z' kappa + P m, summed over species, plus
# the Gamma(a_l_s, b_l_s) log density of the range at that grid point.
range_reference <- function(setup, sd_b, a_l_s = 1, b_l_s = 1) {
  with(setup, {
    n <- nrow(X); d <- ncol(U); ps <- nrow(M_Bs); S <- ncol(Omega)
    grid <- list_SoRSummaries$l_s_grid
    P <- diag(c(1, 1 / sd_b^2, rep(1, d), rep(1 / sigma_bs^2, ps)))
    m <- rbind(0, M_B, matrix(0, d, S), M_Bs)
    vapply(seq_along(grid), function(j) {
      Ks <- matrix(list_SoRSummaries$Ks_all[, , j], nrow = n)
      H <- matrix(0, n, ps)
      for (i in seq_len(n)) H[i, Xs_centers[i, ]] <- Ks[i, ]
      Z <- cbind(1, X, U, H)
      total <- 0
      for (s in seq_len(S)) {
        Q <- t(Z) %*% (Omega[, s] * Z) + P
        h <- t(Z) %*% kappa[, s] + P %*% m[, s]
        total <- total + 0.5 * sum(h * solve(Q, h)) -
          0.5 * as.numeric(determinant(Q, logarithm = TRUE)$modulus)
      }
      total + stats::dgamma(grid[j], a_l_s, b_l_s, log = TRUE)
    }, numeric(1))
  })
}
