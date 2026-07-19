computeESSparams <- function(param_output){

  dim1 <- dim(param_output)[1]
  dim2 <- dim(param_output)[2]

  ESS_param <- matrix(NA, dim1, dim2)
  if(dim1 > 0 & dim2 > 0){
    for (x in 1:nrow(ESS_param)) {
      for (y in 1:ncol(ESS_param)) {
        chains <- lapply(1:dim(param_output)[4], function(i) coda::mcmc(param_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_param[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  }


  ESS_param

}

computeMinESS <- function(results_output){

  beta0_psi_output <- results_output$jsdm_output$B0_output
  beta_psi_output <- results_output$jsdm_output$B_output
  beta_theta_output <- results_output$beta_theta_output
  L_output <- results_output$jsdm_output$L_output

  ESS_beta0psi <- rep(NA, dim(beta0_psi_output)[1])
  for (x in 1:length(ESS_beta0psi)) {
    chains <- lapply(1:dim(beta0_psi_output)[3],
                     function(i) coda::mcmc(beta0_psi_output[x, , i]))
    mcmc_list <- coda::mcmc.list(chains)
    ESS_beta0psi[x] <- coda::effectiveSize(mcmc_list)
  }

  ESS_betapsi <- matrix(NA, dim(beta_psi_output)[1],dim(beta_psi_output)[2])
  if(dim(beta_psi_output)[1] > 1){
    for (x in 1:nrow(ESS_betapsi)) {
      for (y in 1:ncol(ESS_betapsi)) {
        chains <- lapply(1:dim(beta_psi_output)[4], function(i) coda::mcmc(beta_psi_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_betapsi[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  }

  if(!is.null(beta_theta_output)){
    ESS_betatheta <- matrix(NA, dim(beta_theta_output)[1],dim(beta_theta_output)[2])
    for (x in 1:nrow(ESS_betatheta)) {
      for (y in 1:ncol(ESS_betatheta)) {
        chains <- lapply(1:dim(beta_theta_output)[4], function(i) coda::mcmc(beta_theta_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_betatheta[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  } else {
    ESS_betatheta <- NULL
  }

  ESS_L <- matrix(NA, dim(L_output)[1],dim(L_output)[2])
  if(dim(L_output)[1] > 0){
    for (x in 1:nrow(ESS_L)) {
      for (y in 1:ncol(ESS_L)) {
        chains <- lapply(1:dim(L_output)[4], function(i) coda::mcmc(L_output[x, y, , i]))
        mcmc_list <- coda::mcmc.list(chains)
        ESS_L[x,y] <- coda::effectiveSize(mcmc_list)
      }
    }
  }


  return(min(ESS_beta0psi,
             ESS_betatheta,
             ESS_betapsi[!is.na(ESS_betapsi)],
             ESS_L[ESS_L != 0]))

}

