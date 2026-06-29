computeMinESS <- function(results_output){

  beta_ord_output <- results_output$beta_ord_output
  beta_psi_output <- results_output$beta_psi_output
  beta_theta_output <- results_output$beta_theta_output
  LL_output <- results_output$LL_output

  ESS_betaord <- matrix(NA, dim(beta_ord_output)[1],dim(beta_ord_output)[2])
  for (x in 1:nrow(ESS_betaord)) {
    for (y in 1:ncol(ESS_betaord)) {
      chains <- lapply(1:dim(beta_ord_output)[4], function(i) coda::mcmc(beta_ord_output[x, y, , i]))
      mcmc_list <- coda::mcmc.list(chains)
      ESS_betaord[x,y] <- coda::effectiveSize(mcmc_list)
    }
  }

  ESS_betapsi <- matrix(NA, dim(beta_psi_output)[1],dim(beta_psi_output)[2])
  for (x in 1:nrow(ESS_betapsi)) {
    for (y in 1:ncol(ESS_betapsi)) {
      chains <- lapply(1:dim(beta_psi_output)[4], function(i) coda::mcmc(beta_psi_output[x, y, , i]))
      mcmc_list <- coda::mcmc.list(chains)
      ESS_betapsi[x,y] <- coda::effectiveSize(mcmc_list)
    }
  }

  ESS_betatheta <- matrix(NA, dim(beta_theta_output)[1],dim(beta_theta_output)[2])
  for (x in 1:nrow(ESS_betatheta)) {
    for (y in 1:ncol(ESS_betatheta)) {
      chains <- lapply(1:dim(beta_theta_output)[4], function(i) coda::mcmc(beta_theta_output[x, y, , i]))
      mcmc_list <- coda::mcmc.list(chains)
      ESS_betatheta[x,y] <- coda::effectiveSize(mcmc_list)
    }
  }

  ESS_LL <- matrix(NA, dim(LL_output)[1],dim(LL_output)[2])
  for (x in 1:nrow(ESS_LL)) {
    for (y in 1:ncol(ESS_LL)) {
      chains <- lapply(1:dim(LL_output)[4], function(i) coda::mcmc(LL_output[x, y, , i]))
      mcmc_list <- coda::mcmc.list(chains)
      ESS_LL[x,y] <- coda::effectiveSize(mcmc_list)
    }
  }

  return(min(ESS_betaord, ESS_betatheta, ESS_betapsi, ESS_LL[ESS_LL != 0]))

}

