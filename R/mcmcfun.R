
logistic <- function(x){
  1 / (1 + exp(-x))
}

computeTheta <- function(X_theta, beta_theta){
  logistic(X_theta %*% beta_theta)
}

# Conjugate draw of the community mean of the collection intercepts, for the
# experimental theta_intercept_prior = "community": beta_theta[1, s] ~
# Normal(mu_theta, sd^2) for every species and mu_theta ~ Normal(0,
# mu_theta_sd^2).
sample_mu_theta <- function(beta_theta0, sd, mu_theta_sd) {
  precision <- 1/mu_theta_sd^2 + length(beta_theta0)/sd^2
  rnorm(1, sum(beta_theta0)/sd^2/precision, sqrt(1/precision))
}

# sample_betatheta_cpp_parallel() leaves a species' collection coefficients
# unchanged when it is absent from every site, although with no collection
# data their conditional is their prior. Under theta_intercept_prior =
# "community" a stale intercept would feed mu_theta's draw, so this draws
# such species from the prior the sampler used. B_betatheta is diagonal, as
# runOccJSDM() builds it. Default fits do not call this.
sample_betatheta_unoccupied <- function(beta_theta, z, idx_z, b_betatheta,
                                        B_betatheta) {
  unoccupied <- which(colSums(z[idx_z, , drop = FALSE] == 1) == 0)
  for (s in unoccupied) {
    beta_theta[, s] <- rnorm(nrow(beta_theta), b_betatheta, sqrt(diag(B_betatheta)))
  }
  beta_theta
}

# OK
sample_theta0 <- function(z, w, idx_z, a_theta0, b_theta0){

  S <- ncol(z)

  z_all <- z[idx_z, , drop = FALSE]

  theta0 <- rep(NA, S)

  for (s in 1:S) {

    z0_1 <- sum(z_all[,s] == 0 & w[,s] == 1)
    z0_0 <- sum(z_all[,s] == 0 & w[,s] == 0)

    theta0[s] <- rbeta(1, a_theta0 + z0_1, b_theta0 + z0_0)

  }

  theta0
}



# WAIC CALCULATIONS -----

update_waic_summary <- function(logliks, list_waic, iter){

  M2 <- list_waic$M2
  mean_log <- list_waic$mean_log
  mean_lik <- list_waic$mean_lik

  delta_log <- logliks - mean_log
  mean_log <- mean_log + delta_log / iter

  delta2_log <- logliks - mean_log
  M2 <- M2 + delta_log * delta2_log

  delta_lik <- exp(logliks) - mean_lik
  mean_lik <- mean_lik + delta_lik / iter

  list("M2" = M2,
       "mean_log" = mean_log,
       "mean_lik" = mean_lik)

}

compute_waic <- function(list_waic, numIters){

  mean_lik <- list_waic$mean_lik
  var_loglik <- list_waic$M2 / (numIters - 1)

  lppd <- sum(log(mean_lik), na.rm = T)
  p_waic <- sum(var_loglik , na.rm = T)

  WAIC <- - 2 * (lppd - p_waic)

  WAIC
}
