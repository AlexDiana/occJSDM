
logistic <- function(x){
  1 / (1 + exp(-x))
}

computeTheta <- function(X_theta, beta_theta){
  logistic(X_theta %*% beta_theta)
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
