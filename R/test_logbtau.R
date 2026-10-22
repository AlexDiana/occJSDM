{

  edith_params$idx_descendants <- createDescendants(edith_params$idx_ancestors)

  data_info <- as.data.frame(data$info)
  OTU <- data$OTU

  data_info <- data_info %>%
    dplyr::arrange(Site, Sample, Primer)

  M_df <- data_info %>%
    dplyr::group_by(Site, Sample) %>%
    slice(1) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(M = n())

  M <- M_df$M
  names(M) <- M_df$Site

  n <- length(M)

  sumM <- c(0, cumsum(M)[-n])

}

logb_true <- list_simulatedData$true_params$logb_true
beta_psi_true <- list_simulatedData$true_params$beta_psi_true
beta_ord_true <- list_simulatedData$true_params$beta_ord_true
beta_theta_true <- list_simulatedData$true_params$beta_theta_true
L_true <- list_simulatedData$true_params$L_true
U_true <- list_simulatedData$true_params$U_true
w_true <- list_simulatedData$true_params$w_true
sigma_true <- list_simulatedData$true_params$sigma_true
tau_true <- list_simulatedData$true_params$tau_true

X_psi <- list_simulatedData$data_list$X_psi
X_ord <- list_simulatedData$data_list$X_ord
X_theta <- list_simulatedData$data_list$X_theta


logb <- logb_true
U <- U_true
LL <- L_true
beta_psi <- beta_psi_true
beta_theta <- beta_theta_true
w <- w_true
sigma <- sigma_true
tau <- tau_true

niter <- 5000
tau_output <- array(NA, dim = c(n_rivers, niter))
logb_output <- array(NA, dim = c(n, S, niter))

prior_beta_psi <- 0
prior_beta_psi_sd <- 1

for (iter in 1:niter) {

  print(iter)

  tau <- update_tau(tau, logb, w, X_theta, beta_theta,
                    M, sumM, edith_params, sd_proposal = .05,
                    a_tau = 1,
                    b_tau = 1)
  logb <- update_logb(logb, w, X_psi, beta_psi, U, LL,
                      X_theta, beta_theta, sigma,
                      tau, M, sumM, edith_params,
                      sigma_proposals = .05)

  tau_output[,iter] <- tau
  logb_output[,,iter] <- logb

}

# tau
{
  tau_true_vec <- as.vector(tau_true)

  tau_CI <- apply(tau_output, 1, function(x){ quantile(x, probs = c(0.025, 0.975))})

  ggplot(data = NULL,
         aes(x = 1:ncol(tau_CI),
             y = tau_true_vec,
             ymin = tau_CI[1,],
             ymax = tau_CI[2,])) +
    geom_point(size = 1, color = "red") +
    geom_errorbar()

  qplot(1:niter, tau_output[2,])

}

# log b
{
  logb_true_vec <- as.vector(logb_true)

  logb_output_vec <- apply(logb_output, 3, as.vector)
  logb_CI <- apply(logb_output_vec, 1, function(x){ quantile(x, probs = c(0.025, 0.975))})

  ggplot(data = NULL,
         aes(x = 1:ncol(logb_CI),
             y = logb_true_vec,
             ymin = logb_CI[1,],
             ymax = logb_CI[2,])) +
    geom_point(size = 1, color = "red") +
    geom_errorbar()
}
