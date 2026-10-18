# The terms of the occupancy linear predictor, rebuilt from stored draws the
# way computePsiCoef() (R/jsdmfun.R) and KsBproduct() (src/jsdm.cpp) build
# them: eta = B0 + X B + Ks Bs + U L. Sourced by check-storage.R, which checks
# that logistic(eta) equals the stored psi, and by analyse.R.
spatial_term <- function(Ks, Bs, Xs_centers) {
  # KsBproduct(): SE[i, s] = sum_j Ks[i, j] * Bs[Xs_centers[i, j], s]
  n <- nrow(Ks); m <- ncol(Ks)
  vapply(seq_len(ncol(Bs)), function(s) rowSums(Ks * matrix(Bs[Xs_centers, s], n, m)),
         numeric(n))
}
eta_terms <- function(draws, design, t, chain) {
  n <- nrow(design$X_psi)
  B0 <- draws$B0[, t, chain]
  S <- length(B0)
  XB <- design$X_psi %*% matrix(draws$B[, , t, chain], ncol = S)
  Ks <- matrix(design$Ks_all[, , draws$idx_ls[t, chain]], nrow = n)
  SE <- spatial_term(Ks, matrix(draws$Bs[, , t, chain], ncol = S), design$Xs_centers)
  UL <- matrix(draws$U[, , t, chain], nrow = n) %*% matrix(draws$L[, , t, chain], ncol = S)
  list(B0 = B0, XB = XB, SE = SE, UL = UL,
       eta = matrix(B0, n, S, byrow = TRUE) + XB + SE + UL)
}
# What analyse.R and check-storage.R keep of a fit, from its runOccJSDM()
# result and the simulation that produced it.
keep_fit <- function(fit, sim) {
  ro <- fit$results_output
  jo <- ro$jsdm_output
  z <- ro$z_output
  storage.mode(z) <- "integer"
  # The spatial basis at every range of the grid, as the sampler builds it
  # from this fit's support points; idx_ls picks the range per draw.
  Ks_all <- suppressMessages(precomputeSORmatrices(fit$infos$l_s_grid, fit$infos$list_Xs))$Ks_all
  list(
    draws = list(mu_b0 = jo$mu_b0_output, B0 = jo$B0_output, B = jo$B_output,
                 Bs = jo$Bs_output, idx_ls = jo$idx_ls_output,
                 U = jo$U_output, L = jo$L_output,
                 sigmab = jo$sigmab_output, sigmabs = jo$sigmabs_output,
                 sigmah = jo$sigmah_output,
                 beta_theta0 = ro$beta_theta_output[1, , , , drop = TRUE],
                 theta0 = ro$theta0_output, p = ro$p_output,
                 z = z, psi = ro$psi_output),
    design = list(X_psi = fit$X_psi, Tr = fit$Tr,
                  Xs_centers = fit$infos$list_Xs$Xs_centers,
                  Ks_all = Ks_all, l_s_grid = fit$infos$l_s_grid),
    truth = sim$true_params,
    OTU = sim$data_list$OTU, info = sim$data_list$info,
    chain_seeds = fit$infos$chain_seeds,
    intercept_prior = fit$infos$intercept_prior)
}
