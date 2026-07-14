

logit <- function(x) log(x / (1-x))

logistic <- function(x) 1 / (1 + exp(-x))

simulateOccPlusData <- function(list_datasettings,
                         list_params,
                         list_jsdmParams){
  # read data settings
  {
    n <- list_datasettings$n
    S <- list_datasettings$S
    g <- list_datasettings$g
    M <- list_datasettings$M
    L <- list_datasettings$L
    K <- list_datasettings$K
    ncov_psi <- list_datasettings$ncov_psi
    ncov_theta <- list_datasettings$ncov_theta
  }

  # read jsdm params
  {
    gt <- list_jsdmParams$gt
    d <- list_jsdmParams$d
    ds <- list_jsdmParams$ds
    sigma_b <- list_jsdmParams$sigma_b
    sigma_bs <- list_jsdmParams$sigma_bs
    sigma_ts <- list_jsdmParams$sigma_ts
    sigma_h <- list_jsdmParams$sigma_h
    l_s <- list_jsdmParams$l_s
  }

  # read param settings
  {
    p <- list_params$p
    q <- list_params$q
    theta_baseline <- list_params$theta_baseline
  }

  N <- sum(M)
  N2 <- sum(L)
  N3 <- sum(K)
  sumM <- c(0, cumsum(M)[-n])
  sumL <- c(0, cumsum(L)[-N])
  sumK <- c(0, cumsum(K)[-N2])

  idx_z <- rep(1:n, M)
  L_rep <- rep(1:N, L) # right only in same cases
  LK_rep <- rep(L_rep, K)

  list_simjSDMData <- simulateData(
    n, S, ncov_psi,
    g, gt, d, tau = NULL, ds,
    sigma_b, sigma_bs, sigma_ts, sigma_h, l_s,
    useSpatField = F, usingSplines = F, model = "binary")

  # data
  {
    z <- list_simjSDMData$data$z
    X_psi <- list_simjSDMData$data$X
    Tr <- list_simjSDMData$data$Tr
    Xs <- list_simjSDMData$data$Xs
  }

  X_theta <- cbind(1, matrix(rnorm(N * (ncov_theta - 1)), N, ncov_theta - 1))

  # beta_psi_true <- matrix(sample(c(-1,1,0), ncov_psi * S, replace = T), ncov_psi, S)
  beta_theta_true <- matrix(1, ncov_theta, S)
  beta_theta_true[1,] <- logit(theta_baseline)

  theta_true <- logistic(X_theta %*% beta_theta_true)

  theta0_true <- rep(0.05, S)
  p_true <- matrix(.9, max(L), S)
  q_true <- matrix(sample(c(0.05,0.05), size = max(L) * S, replace = T, prob = c(0.95,0.05)),
                   max(L), S)

  z_rep <- z[rep(1:n, M),]
  delta <- matrix(0, N, S)
  y <- matrix(NA, N3, S)

  w <- sapply(1:S, function(s){
    sapply(1:N, function(i){
      if(z_rep[i,s] == 1){
        rbinom(1, 1, theta_true[i,s])
      } else {
        rbinom(1, 1, theta0_true[s])
      }
    })
  })

  L_rep <- rep(1:N, L)
  LK_rep <- rep(L_rep, K)
  w_rep <- w[LK_rep,]

  idx_z <- rep(1:n, M)
  # idx_k <- as.numeric(as.factor(data_info_sample))

  primer_rep <- rep(rep(1:L[1], each = K[1]), times = N) # only right in some cases

  cimk_true <- sapply(1:S, function(s){
    sapply(1:N3, function(i){
      if(w_rep[i,s] == 1){
        rbinom(1, 1, p_true[primer_rep[i],s])
      } else {
        2 * rbinom(1, 1, q_true[primer_rep[i],s])
      }
    })
  })

  y <- sapply(1:S, function(s){
    sapply(1:N3, function(i){
      if(cimk_true[i,s] == 1){

        1

      } else if (cimk_true[i,s] == 2) {

        1

      } else {

        0

      }
    })
  })

  speciesNames <- paste0("OTU_", 1:S)

  colnames(y) <- speciesNames
  rownames(Tr) <- speciesNames

  true_params <- list(
    "beta_theta_true" = beta_theta_true,
    "jsdmParams_true" = list_simjSDMData$trueParams,
    "z_true" = z,
    "w_true" = w,
    "theta_true" = theta_true,
    "p_true" = p_true,
    "q_true" = q_true
  )

  data_list <- list(
    "X_psi" = X_psi,
    "X_theta" = X_theta,
    "Xs" = Xs,
    "Tr" = Tr,
    "y" = y
  )

  list(true_params = true_params,
       data_list = data_list)
}

