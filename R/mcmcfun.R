
createDescendants <- function(idx_ancestors){

  n <- length(idx_ancestors)
  idx_descendants <- vector("list", n)

  # Invert the relationship
  for (i in seq_along(idx_ancestors)) {
    for (a in idx_ancestors[[i]]) {
      idx_descendants[[a]] <- c(idx_descendants[[a]], i)
    }
  }

  idx_descendants
}

logistic <- function(x){
  1 / (1 + exp(-x))
}

computeU <- function(X_ord, beta_ord, E){
  X_ord %*% beta_ord + E
}

computePsi <- function(X_psi, beta_psi, U, LL){

  logit_psi = X_psi %*% beta_psi + U %*% LL

  logistic(logit_psi)

}

computeTheta <- function(X_theta, beta_theta){
  logistic(X_theta %*% beta_theta)
}

# OK
# update_tau <- function(tau, w, sd_tau){
#
#   tau <- rep(0, n_regions)
#
#   for (i in seq_len(n_regions)) {
#
#     tau_star <- rnorm(1, tau, sd = sd_tau)
#
#     conc_old <- compu
#     conc_star <-
#
#   }
#
# }

update_logb <- function(logb, w, X_psi, beta_psi, U, LL,
                        X_theta, beta_theta, sigma,
                        tau, M, sumM, edith_params,
                        sigma_proposals = .02){

  prior_logb <- X_psi %*% beta_psi + U %*% LL
  Xtheta_betatheta <- X_theta %*% beta_theta

  A_si <- edith_params$A_si
  Q_j <- edith_params$Q_j
  L_ij <- edith_params$L_ij
  nu_ij <- edith_params$nu_ij
  idx_ancestors <- edith_params$idx_ancestors
  idx_descendants <- edith_params$idx_descendants
  idx_region <- edith_params$idx_region

  logb <- update_logz_cpp(
    logb,
    w,
    prior_logb,
    sigma,
    Xtheta_betatheta,
    tau,
    M,
    sumM,
    sigma_proposals = 0.02,
    A_si,
    L_ij,
    nu_ij,
    idx_region,
    idx_ancestors,
    idx_descendants,
    Q_j
  )

  logb

}

update_tau <- function(tau, logb, w, X_theta, betatheta, M_site, sumM,
                      edith_params, sd_proposal = .05, a_tau, b_tau){

  Xtheta_betatheta <- X_theta %*% beta_theta

  A_si <- edith_params$A_si
  Q_j <- edith_params$Q_j
  L_ij <- edith_params$L_ij
  nu_ij <- edith_params$nu_ij
  idx_ancestors <- edith_params$idx_ancestors
  idx_descendants <- edith_params$idx_descendants
  idx_region <- edith_params$idx_region

  tau <- update_tau_cpp(tau, logb, w, Xtheta_betatheta, M_site, sumM,
                        A_si, L_ij, nu_ij, idx_region, idx_ancestors, idx_descendants,
                    Q_j, sd_proposal, a_tau, b_tau)

  tau

}

sample_sigma <- function(logb, X_psi, beta_psi, U, LL, a_sigma, b_sigma) {

  logb <- logb - X_psi %*% beta_psi - U %*% LL

  sigma <- update_sigma_cpp(sigma, logb, a_sigma, b_sigma)

  sigma

}

compute_edith_contribution_element_r <- function(
  A_s_i, L_i, v_i, logb_i, tau, Q_j, idx_list
){

  sum_term <- 0

  for (k in seq_along(idx_list)) {

    sum_term <- sum_term + compute_edith_contribution(
      A_s_i[idx_list[k]],
      L_i[idx_list[k]],
      v_i[idx_list[k]],
      logb_i[idx_list[k]],
      tau,
      Q_j
    )

  }

  sum_term

}

update_logz_r <- function(logb,
                          w,
                          prior_logb,
                          sigma,
                          Xtheta_betatheta,
                          tau,
                          M_site,
                          sigma_proposals,
                          A_s_i,
                          L_ij,
                          nu_ij,
                          idx_region,
                          idx_ancestors,
                          idx_descendants,
                          Q_j){

    n <- length(M)
    S <- ncol(log_b)

    for (i in 1:n) {

      idx_anc <- idx_ancestors[[i]]
      idx_desc <- idx_descendants[[i]]

      L_i <- L_ij[i,]
      nu_i <- nu_ij[i,]
      logb_i <- logb[,i]

      for (j in 1:S) {

        logb_current <- logb[i,j]

        logb_new <- rnorm(1, logb_current, sigma_proposals)

        loglik_new <- 0
        loglik_old <- 0

        for (k in seq_along(idx_desc)) { # loop through the sites affected

          idx_descendant <- idx_desc[k]

          conc_old =
            compute_edith_contribution_element_r(
              A_s_i,
              L_i,
              v_i,
              logb_i,
              tau[idx_region[i]],
              Q_j[i],
              idx_anc
            )

          conc_element_old =
            compute_edith_contribution(
              A_s_i[i],
              L_ij(i,idx_descendant),
              v_ij(i,idx_descendant),
              logb_current,
              tau[idx_region[i]],
              Q_j[i]
            )

          conc_element_new =
            compute_edith_contribution(
              A_s_i[i],
              L_ij(i,idx_descendant),
              v_ij(i,idx_descendant),
              logb_new,
              tau[idx_region[i] - 1],
              Q_j[i]
            )

          conc_new <-  conc_old -
            conc_element_old +
            conc_element_new;

          # for(m in 1:M_site[idx_descendant]){
          #
          #   p_current = 1 / (1 + exp(- Xtheta_betatheta(sumM[idx_descendant] + m, j) - conc_old));
          #   p_star = 1 / (1 + exp(- Xtheta_betatheta(sumM[idx_descendant] + m, j) - conc_new));
          #
          #   loglik_new += R::dbinom(w(sum_m + m, j), 1, p_star, 1);
          #   loglik_old += R::dbinom(w(sum_m + m, j), 1, p_current, 1);
          #
          # }

        }

        # double logprior = R::dnorm(logb_new, prior_logb(i,j), sigma[j], 1) -
        #   R::dnorm(logb_current, prior_logb(i,j), sigma[j], 1);
        #
        # double logposterior = logprior + loglik_new - loglik_old;
        #
        # if(R::runif(0,1) < exp(logposterior)){
        #
        #   logb(i, j) = logb_new;
        #
        # }

      }


    }

    logb

}

# OK
sample_z <- function(w, psi, theta, theta0){

  S <- ncol(psi)
  n <- nrow(psi)

  z <- matrix(NA, n, S)

  for (s in 1:S) {

    for (i in 1:n) {

      p_zsequal1 <-
        sum(
          sapply(1:M[i], function(m){
            dbinom(w[sumM[i] + m,s], 1, theta[sumM[i] + m,s], log = T)
          })
        ) + dbinom(1, 1, psi[i,s], log = T)

      p_zsequal0 <-
        sum(
          sapply(1:M[i], function(m){
            dbinom(w[sumM[i] + m,s], 1, theta0[s], log = T)
          })
        ) + dbinom(0, 1, psi[i,s], log = T)

      p_1 <- exp(p_zsequal1) / (exp(p_zsequal1) + exp(p_zsequal0))

      z[i,s] <- rbinom(1, 1, p_1)

    }

  }

  z
}

# OK
sample_w <- function(y, theta, theta0, p, q,
                     M, K, sumL, sumM, sumK, maxL){

  S <- ncol(theta)
  N <- nrow(theta)
  w <- matrix(NA, nrow = N, ncol = S)

  for (s in 1:S) {

    for (i in 1:n) {

      for (m in 1:M[i]) {

        p_wsequal1 <-
          sum(
            sapply(1:maxL, function(l){
              sapply(1:K[sumL[sumM[i] + m] + l], function(k){
                dbinom(y[sumK[sumL[sumM[i] + m] + l] + k,s], 1,
                       p[l,s], log = T)
              })
            })
          )

        p_wsequal0 <-
          sum(
            sapply(1:maxL, function(l){
              sapply(1:K[sumL[sumM[i] + m] + l], function(k){
                dbinom(y[sumK[sumL[sumM[i] + m] + l] + k,s], 1,
                       q[l,s], log = T)
              })
            })
          )

        if(z[i, s] == 1){

          p_wsequal1 <- p_wsequal1 + dbinom(1, 1, theta[sumM[i] + m,s], log = T)
          p_wsequal0 <- p_wsequal0 + dbinom(0, 1, theta[sumM[i] + m,s], log = T)

        } else {

          p_wsequal1 <- p_wsequal1 + dbinom(1, 1, theta0[s], log = T)
          p_wsequal0 <- p_wsequal0 + dbinom(0, 1, theta0[s], log = T)

        }

        p_ws1 <- exp(p_wsequal1) / (exp(p_wsequal1) + exp(p_wsequal0))

        w[sumM[i] + m,s] <- rbinom(1, 1, p_ws1)

      }

    }

  }

  w
}

# OK
sample_pq <- function(c_imk, w, primerIdx, idx_k, maxL, a_p,
                      b_p, a_q, b_q){

  p <- matrix(NA, maxL, S)
  q <- matrix(NA, maxL, S)

  w_all <- w[idx_k,,drop=F]

  for (s in 1:S) {

    for (l in 1:maxL) {

      w1_primerl_cases_1 <- sum(primerIdx == l & w_all[,s] == 1 & c_imk[,s] == 1)
      w1_primerl_cases_0 <- sum(primerIdx == l & w_all[,s] == 1 & c_imk[,s] == 0)

      w0_primerl_cases_1 <- sum(primerIdx == l & w_all[,s] == 0 & c_imk[,s] == 2)
      w0_primerl_cases_0 <- sum(primerIdx == l & w_all[,s] == 0 & c_imk[,s] == 0)

      p[l,s] <- rbeta(1, a_p + w1_primerl_cases_1, b_p + w1_primerl_cases_0)
      q[l,s] <- rbeta(1, a_q + w0_primerl_cases_1, b_q + w0_primerl_cases_0)


    }

  }

  list("p" = p,
       "q" = q)

}

# OK
sample_betatheta <- function(w, z, beta_theta, idx_z, X_theta,
                             b_betatheta, B_betatheta){

  ncov_theta <- nrow(beta_theta)
  S <- ncol(beta_theta)

  z_all <- z[idx_z,]

  for (s in 1:S) {

    k <- as.vector(t(w[z_all[,s]==1,s])) - .5

    n <- rep(1, length(k))

    X_thetasubset <- X_theta[z_all[,s]==1,,drop=F]

    beta_theta[,s] <- sample_beta_nocov_cpp(beta_theta[,s], X_thetasubset,
                                            b_betatheta, B_betatheta, n, k)

  }

  beta_theta
}

# OK
sample_theta0 <- function(z, w, idx_z, a_theta0, b_theta0){

  S <- ncol(z)

  z_all <- z[idx_z,]

  theta0 <- rep(NA, S)

  for (s in 1:S) {

    z0_1 <- sum(z_all[,s] == 0 & w[,s] == 1)
    z0_0 <- sum(z_all[,s] == 0 & w[,s] == 0)

    theta0[s] <- rbeta(1, a_theta0 + z0_1, b_theta0 + z0_0)

  }

  theta0
}

# OK
sample_betapsiLL <- function(beta_psi, LL, Omega, X_psi, k, U,
                             prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  d <- nrow(LL)
  S <- ncol(LL)

  beta_all <- matrix(NA, ncov_psi + d, S)

  for (s in 1:S) {

    elemZero <- max(d - s, 0)
    elem1 <- ifelse(s <= d, 1, 0)
    # elemZero <- elem1 <- 0
    elemNonZero <- d - elemZero - elem1

    if(elem1 > 0){

      k_current <- k[,s] - Omega[,s] * U[,s]

    } else {

      k_current <- k[,s]

    }

    X_all <- cbind(X_psi, U[,seq_len(elemNonZero)])

    b_current <- rep(prior_beta_psi, ncov_psi + elemNonZero)
    B_current <- diag(prior_beta_psi_sd, nrow = ncov_psi + elemNonZero)

    # betapntsiLL_current <- c(beta_psi[,s], LL[,s])
    if(ncol(X_all) > 0){
      betapsiLL <- sample_beta_cpp(X_all, B_current, b_current, Omega[,s], k_current)
      beta_all[seq_len(ncov_psi),s] <- betapsiLL[seq_len(ncov_psi)]
    }

    if(elemNonZero > 0){

      beta_all[ncov_psi + seq_len(elemNonZero),s] <-
        betapsiLL[ncov_psi + seq_len(elemNonZero)]

    }


  }

  if(T){

    # set lower diagonal to 0
    beta_all[ncov_psi + seq_len(d),seq_len(d)][!upper.tri(beta_all[ncov_psi + seq_len(d),seq_len(d)])] <- 0

    # set diagonal to 1
    if(d > 1){

      diag(beta_all[ncov_psi + seq_len(d),seq_len(d)]) <- 1
    } else {
      beta_all[ncov_psi + 1,1] <- 1
    }

  }


  beta_psi <- beta_all[seq_len(ncov_psi),,drop=F]
  LL <- beta_all[ncov_psi + seq_len(d),,drop=F]

  list("beta_psi" = beta_psi,
       "LL" = LL)

}

# OK
sample_LL <- function(beta_psi, LL, X_psi, logb, U,
                             prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  d <- nrow(LL)
  S <- ncol(LL)
  n <- nrow(U)

  LL <- matrix(NA, d, S)

  logb <- logb - X_psi %*% beta_psi

  for (s in 1:S) {

    elemZero <- max(d - s, 0)
    elem1 <- ifelse(s <= d, 1, 0)
    # elemZero <- elem1 <- 0
    elemNonZero <- d - elemZero - elem1

    if(elem1 > 0){

      logb_current <- logb[,s] -  U[,s]

    } else {

      logb_current <- logb[,s]

    }



    # betapntsiLL_current <- c(beta_psi[,s], LL[,s])

    if(elemNonZero > 0){

      X_all <- matrix(U[,seq_len(elemNonZero)], n, elemNonZero)

      b_current <- rep(prior_beta_psi, elemNonZero)
      B_current <- diag(prior_beta_psi_sd, nrow = elemNonZero)

      LL_s <- sample_beta_cpp(X_all, B_current, b_current, rep(1, nrow(X_all)), logb_current)
      LL[seq_len(elemNonZero),s] <- LL_s

    }

  }

  if(T){

    # set lower diagonal to 0
    LL[seq_len(d),seq_len(d)][!upper.tri(LL[seq_len(d),seq_len(d)])] <- 0

    # set diagonal to 1
    if(d > 1){
      for (i in 1:d) {
        LL[i,i] <- 1
      }
      # diag(LL[ncov_psi + seq_len(d),seq_len(d)]) <- 1
    } else {
      LL[1,1] <- 1
    }

  }

  LL

}

sample_betapsi <- function(beta_psi, LL, X_psi, logb, U,
                             prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  d <- nrow(LL)
  S <- ncol(LL)

  beta_psi <- matrix(NA, ncov_psi, S)

  for (s in 1:S) {

    logb_current <- logb[,s] - U %*% LL[,s]

    b_current <- rep(prior_beta_psi, ncov_psi)
    B_current <- diag(prior_beta_psi_sd, nrow = ncov_psi)


    beta_psi[,s] <- sample_beta_cpp(X_psi, B_current, b_current,
                                    rep(1, nrow(X_psi)), as.vector(logb_current))

  }

  beta_psi

}

# OK
sample_betaord <- function(logb, LL, X_ord, X_psi, beta_psi, E){

  ncov_ord <- ncol(X_ord)
  d <- nrow(LL)
  S <- ncol(LL)

  Omega_vec <- as.vector(Omega)

  logb_new <- logb - (X_psi %*% beta_psi + E %*% LL)

  Xtilde <- kronecker(t(LL), X_ord)

  beta_ord <- sample_beta_cpp(Xtilde,
                              diag(2, nrow = ncol(Xtilde)),
                              rep(0, ncol(Xtilde)),
                              rep(1, nrow(Xtilde)),
                              logb_new)

  matrix(beta_ord, ncov_ord, d, byrow = F)
}

sample_E <- function(logb, X_psi, beta_psi, X_ord, beta_ord, LL){

  d <- nrow(LL)
  n <- nrow(logb)

  logb_new <- logb - (X_psi %*% beta_psi + X_ord %*% beta_ord %*% LL)

  E <- matrix(NA, n, d)

  B_current <- diag(2, nrow = d)
  b_current <- rep(0, d)

  for (i in 1:n) {

    E[i,] <- sample_beta_cpp(t(LL), B_current, b_current,
                             rep(1, length(logb_new[i,])),
                             logb_new[i,])

  }

  E

}

sample_betaordE <- function(k, Omega, X_psi, beta_psi, X_ord, LL){

  d <- nrow(LL)
  n <- nrow(k)

  k_new <- k - X_psi %*% beta_psi

  B_current <- diag(1, nrow = (ncov_ord + n) * d)
  b_current <- rep(0, (ncov_ord + n) * d)

  Xtilde <- cbind(X_ord, diag(1, nrow = n))
  X_current <- kronecker(t(LL), Xtilde)

  beta_ordE <- sample_beta_cpp(X_current, B_current, b_current,
                               as.vector(Omega), as.vector(k_new))
  beta_ordE_mat <- matrix(beta_ordE, ncov_ord + n, d)

  beta_ord <- beta_ordE_mat[1:ncov_ord,,drop=F]
  E <- beta_ordE_mat[ncov_ord + 1:n,,drop=F]

  list("beta_ord" = beta_ord,
       "E" = E)

}

sample_psivars <- function(logb, X_psi, beta_psi, X_ord, beta_ord, E, LL,
                           prior_beta_psi, prior_beta_psi_sd){

  ncov_psi <- ncol(X_psi)
  ncov_ord <- ncol(X_ord)

  if(ncov_psi > 0){
    beta_psi <- sample_betapsi(beta_psi, LL, X_psi, logb, U,
                               prior_beta_psi, prior_beta_psi_sd)
  }

  LL <- sample_LL(beta_psi, LL, X_psi, logb, U,
                  prior_beta_psi, prior_beta_psi_sd)

  E <- sample_E(logb, X_psi, beta_psi, X_ord, beta_ord, LL)

  if(ncov_ord > 0){
    # list_betaordE <- sample_betaordE(k, Omega, X_psi, beta_psi, X_ord, LL)
    # beta_ord <- list_betaordE$beta_ord
    # E <- list_betaordE$E
    beta_ord <- sample_betaord(logb, LL, X_ord, X_psi, beta_psi, E)
  }

  list("beta_psi" = beta_psi,
       "beta_ord" = beta_ord,
       "LL" = LL,
       "E" = E)

}

sample_musigma <- function(sigma1,
                            logy1, c_imk,
                            mu_mu1, sd_mu1,
                            a_sigma1, b_sigma1, tpfp){

  idx_tpfp <- ifelse(tpfp, 1, 2)

  logy1_cimk1 <- logy1[c_imk == idx_tpfp]
  n_samples <- length(logy1_cimk1)

  if(n_samples > 0){

    sample_mean <- mean(logy1_cimk1)

    ww <- (1/sd_mu1^2) / ((1/sd_mu1^2) + (n_samples / sigma1^2))
    mu_n <- ww * mu_mu1 + (1-ww) * sample_mean
    sigma2_n <- 1 / ((1 / sd_mu1^2) + (n_samples / sigma1^2))


  } else {

    mu_n <- mu_mu1
    sigma2_n <- sd_mu1

  }

  mu1 <- rnorm(1, mean = mu_n, sd = sqrt(sigma2_n))

  if(n_samples > 0){

    alpha_n <- a_sigma1 + n_samples / 2
    beta_n <- b_sigma1 + 0.5 * sum((logy1_cimk1 - mu1)^2)

  } else {

    alpha_n <- a_sigma1
    beta_n <- b_sigma1

  }

  sigma1 <- sqrt(1 / rgamma(1, shape = alpha_n, rate = beta_n))

  list("mu1" = mu1,
       "sigma1" = sigma1)
}

sample_sigma0 <- function(logy1, c_imk,
                          a_sigma1, b_sigma1){

  logy1_cimk1 <- logy1[c_imk == 2]
  n_samples <- length(logy1_cimk1)

  if(n_samples > 0){

    alpha_n <- a_sigma1 + n_samples / 2
    beta_n <- b_sigma1 + 0.5 * sum((logy1_cimk1)^2)

  } else {

    alpha_n <- a_sigma1
    beta_n <- b_sigma1

  }

  sigma0 <- sqrt(1 / rgamma(1, shape = alpha_n, rate = beta_n))

  sigma0
}

sample_mu0sigma0 <- function(y, c_imk,
                            a_pi0, b_pi0,
                            a_sigma0, b_sigma0){

  logy1_cimk0 <- y[c_imk == 0]
  logy1_cimk0_pos <- logy1_cimk0[logy1_cimk0 > 0]

  n_samples <- length(logy1_cimk0)

  num0 <- sum(logy1_cimk0 == 0)
  num1 <- sum(logy1_cimk0 > 0)

  pi0 <- rbeta(1, a_pi0 + num0, b_pi0 + (n_samples - num0))

  alpha_n <- a_sigma0 + num1 / 2
  beta_n <- b_sigma0 + 0.5 * sum(logy1_cimk0_pos^2)

  sigma0 <- sqrt(1 / rgamma(1, shape = alpha_n, rate = beta_n))

  list("pi0" = pi0,
       "sigma0" = sigma0)
}

sample_cimk <- function(logy1, mu1, sigma1, pi0, sigma0,
                        p, q, idx_k, primerIdx){

  N3 <- length(idx_k)
  S <- ncol(p)

  c_imk <- matrix(NA, N3, S)

  for (i in 1:N3) {
    for (s in 1:S) {

      term1_loglik <- dnorm(logy1[i,s], mu1, sigma1, log = T)
      term2_loglik <- ifelse(logy1[i,s] == 0, log(pi0),
                             dnorm(y[i,s], 0, sigma0, log = T) - log(.5))
      # term2_loglik <- ifelse(logy1[i,s] == 0, log(pi0),
      #                        dlaplace(logy1[i,s], 0, sigma0, log = T) - log(.5))
                             # dnorm(logy1[i,s], 0, sigma0, log = T) - log(.5))

        # log(pi0 * dnorm(logy1[i,s], 0, sigma0, log = T) +
        #                     (1 - pi0))

      if(w[idx_k[i],s] == 1){
        term1_prior <- dbinom(1, 1, p[primerIdx[i],s], log = T)
        term2_prior <- dbinom(0, 1, p[primerIdx[i],s], log = T)
      } else {
        term1_prior <- dbinom(1, 1, q[primerIdx[i],s], log = T)
        term2_prior <- dbinom(0, 1, q[primerIdx[i],s], log = T)
      }

      term12_diff <- (term1_loglik + term1_prior) - (term2_loglik + term2_prior)

      # p_cimk1 <- exp(term12_diff) / (exp(term12_diff) + 1)
      p_cimk1 <- 1 / (exp(-term12_diff) + 1)


      c_imk[i,s] <- rbinom(1, 1, p_cimk1)
    }
  }

  c_imk
}

loglik_sigma1 <- function(w, logy1){

  sum(
    sapply(1:S, function(s){
      sapply(1:N3, function(i){

        p[primerIdx[]]

      })
    })
  )

}
