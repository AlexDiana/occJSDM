createDistance <- function(points) {

  # Convert to matrix if data.frame
  if (is.data.frame(points)) {
    points <- as.matrix(points)
  }

  # Use R's built-in dist function
  dist_matrix <- as.matrix(dist(points, method = "euclidian"))

  return(dist_matrix)
}

simulateData <- function(n, S, M, L, K, d,
                         ncov_psi, ncov_theta, ncov_ord,
                         a_sigma, b_sigma,
                         edith_params){


  M <- rep(M, n)
  N <- sum(M)
  L <- rep(L, N)
  N2 <- sum(L)
  K <- rep(K, N2)
  N3 <- sum(K)

  sumM <- c(0, cumsum(M)[-n])
  sumL <- c(0, cumsum(L)[-N])
  sumK <- c(0, cumsum(K)[-N2])

  idx_z <- rep(1:n, M)
  L_rep <- rep(1:N, L) # right only in same cases
  LK_rep <- rep(L_rep, K)

  X_psi <- matrix(rnorm(n * ncov_psi), n, ncov_psi)
  X_theta <- cbind(1, matrix(rnorm(N * (ncov_theta - 1)), N, ncov_theta - 1))
  X_ord <- matrix(rnorm(n * ncov_ord), n, ncov_ord)

  beta_psi_true <- matrix(sample(c(-1,1,0), ncov_psi * S, replace = T), ncov_psi, S)
  beta_theta_true <- matrix(1, ncov_theta, S)
  (beta_ord_true <- matrix(sample(c(-1,1,0), ncov_ord * d, replace = T), ncov_ord, d))

  E_true <- matrix(rnorm(n * d, sd = .1), n, d)
  U_true <- X_ord %*% beta_ord_true + E_true
  L_true <- matrix(rnorm(d * S, sd = 1), d, S)
  if(d == 1){
    L_true[1,1] <- 1
  } else {
    for (i in 1:d) {
      L_true[i,i] <- 1
    }
  }
  L_true[seq_len(d), seq_len(d)][lower.tri(L_true[seq_len(d), seq_len(d)])] <- 0
  sigma <- rgamma(S, a_sigma, b_sigma)

  Xbeta_logb <- X_psi %*% beta_psi_true + U_true %*% L_true
  logb_true <- t(sapply(1:n, function(i){
    sapply(1:S, function(j){
      0#rnorm(1, mean = Xbeta_logb[i,j], sd = sigma[j])
    })
  }))

  # river params
  tau <- exp(seq(log(0.25), log(20), length.out = n_rivers)) # smaller tau, less transportation

  conc_true <- compute_C_hat_cpp(edith_params$A_si,
                                 edith_params$L_ij,
                                 edith_params$nu_ij,
                                 logb_true,
                                 tau,
                                 edith_params$idx_region,
                                 edith_params$idx_ancestors,
                                 edith_params$Q_j)

  beta_conc <- rep(1, S)
  conc_true_samples <- conc_true[idx_z,]
  theta_true <-
    logistic(X_theta %*% beta_theta_true + conc_true_samples * matrix(beta_conc, N, S, byrow = T))
  # theta0_true <- rep(0.05, S)
  p_true <- matrix(.9, max(L), S)
  # q_true <- matrix(sample(c(0.05,0.05), size = max(L) * S, replace = T, prob = c(0.95,0.05)),
  #                  max(L), S)

  w <- matrix(0, N, S)
  y <- matrix(NA, N3, S)

  w <- sapply(1:S, function(s){
    sapply(1:N, function(i){
      rbinom(1, 1, theta_true[i,s])
    })
  })

  w_rep <- w[LK_rep,]

  primer_rep <- rep(rep(1:L[1], eac = K[1]), times = N)

  y <- sapply(1:S, function(s){
    sapply(1:N3, function(i){
      if(w_rep[i,s] == 1){
        rbinom(1, 1, p_true[primer_rep[i],s])
      } else {
        0
      }
    })
  })

  true_params <- list(
    "beta_psi_true" = beta_psi_true,
    "beta_theta_true" = beta_theta_true,
    "beta_ord_true" = beta_ord_true,
    "logb_true" = logb_true,
    "w_true" = w,
    "sigma_true" = sigma,
    "tau_true" = tau,
    "theta_true" = theta_true,
    "p_true" = p_true,
    "L_true" = L_true,
    "U_true" = U_true
  )

  data_list <- list(
    "X_psi" = X_psi,
    "X_theta" = X_theta,
    "X_ord" = X_ord,
    "y" = y
  )


  list(true_params = true_params,
       data_list = data_list)
}
