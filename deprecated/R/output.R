
logit <- function(x){
  log(x / (1 - x))
}

logistic <- function(x) 1 / (1 + exp(-x))

#' returnOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Returns the occupancy covariates coefficients posterior sample
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnOccupancyCovariates <- function(fitmodel){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_psi <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  occCovNames <- colnames(fitmodel$X_psi)
  # idxcov <- which(occCovNames == covName)

  # samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  # samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_psi]

  beta_psi_output0 <-
    matrix_of_draws[,grepl("beta_psi\\[", colnames(matrix_of_draws))]

  niter <- nrow(beta_psi_output0)

  beta_psi_output <- array(NA, dim = c(niter, ncov_psi, S))
  for(iter in 1:niter){
    beta_psi_output[iter,,] <- matrix(beta_psi_output0[iter,], ncov_psi, S, byrow = F)
  }

  dimnames(beta_psi_output)[[2]] <- occCovNames
  dimnames(beta_psi_output)[[3]] <- speciesNames

  beta_psi_output

}

#' plotOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyCovariates <- function(fitmodel,
                                    covName = NULL,
                                    idx_species = NULL
){

  beta_psi_output <- returnOccupancyCovariates(fitmodel)

  if(is.null(covName)){
    stop("No name provided")
  }

  # matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_psi <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  occCovNames <- colnames(fitmodel$X_psi)
  idxcov <- which(occCovNames == covName)

  if(length(idxcov) == 0){
    stop("Covariate name not found. If you are using a categorical covariates,
         the name might have changed to code the level. Use
         colnames(fitmodel$X_occ) to find the new names")
  }

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  # samples_subset <- beta_psi_output[,,]

  samples_subset <- matrix(beta_psi_output[,idxcov, idx_species],
                           dim(beta_psi_output)[1], length(idx_species))

  # samples_subset <- beta_ord_output[,idxcov,]
  # samples_subset <- samples_subset[,idx_species,drop=F]

  # param <- "beta_psi"
  #
  # samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  # samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_psi]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames[idx_species]) %>%
    mutate(speciesOrder = order(`2.5%`)) #%>%
    # filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  plot_occcovs <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + geom_hline(aes(yintercept = 0), color = "red")

  plot_occcovs

}


#' returnOrdinationCovariates
#'
#' Ordination covariate coefficients.
#'
#' @details
#' Returns the 95% credible interval of the ordination covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnOrdinationCovariatesOutput <- function(fitmodel){

  # if(is.null(covName)){
  #   stop("No name provided")
  # }

  matrix_of_draws <- fitmodel$matrix_of_draws

  d <- fitmodel$infos$d
  ncov_ord <- fitmodel$infos$ncov_ord
  speciesNames <- fitmodel$infos$speciesNames
  ordCovNames <- colnames(fitmodel$X_ord)
  nsites <- length(fitmodel$infos$siteNames)
  nspecies <- length(fitmodel$infos$speciesNames)

  # if(length(idxcov) == 0){
  #   stop("Covariate name not found. If you are using a categorical covariates,
  #        the name might have changed to code the level. Use
  #        colnames(fitmodel$X_ord) to find the new names")
  # }
  #
  # if(is.null(idx_factors)){
  #   idx_factors <- 1:d
  # }

  {
    U_output0 <-
      matrix_of_draws[,grepl("U\\[", colnames(matrix_of_draws))]
    L_output0 <-
      matrix_of_draws[,grepl("L\\[", colnames(matrix_of_draws))]
    beta_ord_output0 <-
      matrix_of_draws[,grepl("beta_ord\\[", colnames(matrix_of_draws))]

    niter <- nrow(U_output0)

    U_output <- array(NA, dim = c(niter, nsites, d))
    for(iter in 1:niter){
      U_output[iter,,] <- matrix(U_output0[iter,], nsites, d, byrow = F)
    }

    L_output <- array(NA, dim = c(niter, d, nspecies))
    for(iter in 1:niter){
      L_output[iter,,] <- matrix(L_output0[iter,], d, nspecies, byrow = F)
    }

    beta_ord_output <- array(NA, dim = c(niter, ncov_ord, d))
    for(iter in 1:niter){
      beta_ord_output[iter,,] <- matrix(beta_ord_output0[iter,], ncov_ord, d, byrow = F)
    }

    niter <- nrow(L_output)

    L_output_reparam <- L_output
    U_output_reparam <- U_output
    # E_output_reparam <- E_output
    beta_ord_output_reparam <- beta_ord_output

    d <- dim(L_output)[2]

    for (iter in 1:niter) {
      # print(iter)

      if(d == 1){

        L1 <- L_output[iter,1,1]
        L_output_reparam[iter,1,] <- L_output[iter,1,] / L1
        U_output_reparam[iter,,1] <- U_output[iter,,1] * L1
        beta_ord_output_reparam[iter,,1] <- beta_ord_output[iter,,1] * L1

      } else {

        L_current <- L_output[iter,,]
        # E_current <- E_output[iter,,]
        U_current <- U_output[iter,,]
        beta_ord_current <- beta_ord_output[iter,,]

        qr_decomp <- qr(L_current)
        Q_current <- qr.Q(qr_decomp)
        R_current <- qr.R(qr_decomp)

        Q2 <- Q_current %*% diag(diag(R_current), nrow = d)
        invQ2 <- diag(1 / diag(R_current), nrow = d) %*% t(Q_current)

        betapsiord_new <- beta_ord_current %*% Q2
        # E_new <- E_current %*% Q2
        L_new <- invQ2 %*% L_current
        U_new <- U_current %*% Q2

        L_output_reparam[iter,,] <- L_new
        # E_output_reparam[iter,] <- E_new
        U_output_reparam[iter,,] <- U_new
        beta_ord_output_reparam[iter,,] <- betapsiord_new
      }

    }
  }

  beta_ord_output_reparam

  # samples_subset <- beta_ord_output_reparam
  # samples_subset <- samples_subset[,idxcov + 0:(d - 1)*ncov_ord,drop=F]

  # samples_subset

}

#' plotOrdinationCovariates
#'
#' Ordination covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the ordination covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOrdinationCovariates <- function(fitmodel,
                                     covName = NULL,
                                     idx_factors = NULL
){

  # if(is.null(covName)){
  #   stop("No name provided")
  # }
  #
  # matrix_of_draws <- fitmodel$matrix_of_draws
  #
  # d <- fitmodel$infos$d
  # ncov_ord <- fitmodel$infos$ncov_ord
  # speciesNames <- fitmodel$infos$speciesNames
  # ordCovNames <- colnames(fitmodel$X_ord)
  # idxcov <- which(ordCovNames == covName)
  # nsites <- nrow(fitmodel$X_psi)
  #
  # if(length(idxcov) == 0){
  #   stop("Covariate name not found. If you are using a categorical covariates,
  #        the name might have changed to code the level. Use
  #        colnames(fitmodel$X_ord) to find the new names")
  # }
  #
  # if(is.null(idx_factors)){
  #   idx_factors <- 1:d
  # }
  #
  # param <- "beta_ord"
  #
  # {
  #   U_output0 <-
  #     matrix_of_draws[,grepl("U\\[", colnames(matrix_of_draws))]
  #   L_output0 <-
  #     matrix_of_draws[,grepl("L\\[", colnames(matrix_of_draws))]
  #   beta_ord_output0 <-
  #     matrix_of_draws[,grepl("beta_ord\\[", colnames(matrix_of_draws))]
  #
  #   niter <- nrow(U_output0)
  #
  #   U_output <- array(NA, dim = c(niter, nsites, d))
  #   for(iter in 1:niter){
  #     U_output[iter,,] <- matrix(U_output0[iter,], nsites, d, byrow = F)
  #   }
  #
  #   L_output <- array(NA, dim = c(niter, d, nspecies))
  #   for(iter in 1:niter){
  #     L_output[iter,,] <- matrix(L_output0[iter,], d, nspecies, byrow = F)
  #   }
  #
  #   beta_ord_output <- array(NA, dim = c(niter, ncov_ord, d))
  #   for(iter in 1:niter){
  #     beta_ord_output[iter,,] <- matrix(beta_ord_output0[iter,], ncov_ord, d, byrow = F)
  #   }
  #
  #   L_output_reparam <- L_output
  #   U_output_reparam <- U_output
  #   # E_output_reparam <- E_output
  #   beta_ord_output_reparam <- beta_ord_output
  #
  #   d <- dim(L_output)[2]
  #
  #   for (iter in 1:niter) {
  #     print(iter)
  #
  #     if(d == 1){
  #
  #       L1 <- L_output[iter,1,1]
  #       L_output_reparam[iter,1,] <- L_output[iter,1,] / L1
  #       U_output_reparam[iter,,1] <- U_output[iter,,1] * L1
  #       beta_ord_output_reparam[iter,,1] <- beta_ord_output[iter,,1] * L1
  #
  #     } else {
  #
  #       L_current <- L_output[iter,,]
  #       # E_current <- E_output[iter,,]
  #       U_current <- U_output[iter,,]
  #       beta_ord_current <- beta_ord_output[iter,,]
  #
  #       qr_decomp <- qr(L_current)
  #       Q_current <- qr.Q(qr_decomp)
  #       R_current <- qr.R(qr_decomp)
  #
  #       Q2 <- Q_current %*% diag(diag(R_current), nrow = d)
  #       invQ2 <- diag(1 / diag(R_current), nrow = d) %*% t(Q_current)
  #
  #       betapsiord_new <- beta_ord_current %*% Q2
  #       # E_new <- E_current %*% Q2
  #       L_new <- invQ2 %*% L_current
  #       U_new <- U_current %*% Q2
  #
  #       L_output_reparam[iter,,] <- L_new
  #       # E_output_reparam[iter,] <- E_new
  #       U_output_reparam[iter,,] <- U_new
  #       beta_ord_output_reparam[iter,,] <- betapsiord_new
  #     }
  #
  #   }
  # }
  #
  # samples_subset <- beta_ord_output_reparam
  # # samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  # samples_subset <- samples_subset[,idxcov + 0:(d - 1)*ncov_ord,drop=F]

  beta_ord_output <- returnOrdinationCovariatesOutput(fitmodel)

  if(is.null(covName)){
    stop("No name provided")
  }

  ordCovNames <- colnames(fitmodel$X_ord)
  idxcov <- which(ordCovNames == covName)

  if(length(idxcov) == 0){
    stop("Covariate name not found. If you are using a categorical covariates,
         the name might have changed to code the level. Use
         colnames(fitmodel$X_ord) to find the new names")
  }

  d <- fitmodel$infos$d

  if(is.null(idx_factors)){
    idx_factors <- 1:d
  }

  samples_subset <- matrix(beta_ord_output[,idxcov, idx_factors],
                           dim(beta_ord_output)[1], length(idx_factors))

  # beta_ord_output[,idxcov,]
  # samples_subset <- samples_subset[,idx_factors,drop=F]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Factor = 1:d) %>%
    filter(Factor %in% idx_factors)

  plot_covs <- data_plot %>%
    ggplot(aes(x = Factor,
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Factors") +
    # ylim(c(0,1)) +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + geom_hline(aes(yintercept = 0), color = "red")

  plot_covs

}

#' returnDetectionCovariates
#'
#' Detection covariate coefficients.
#'
#' @details
#' Returns the detection covariates coefficients posterior sample
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnDetectionCovariates <- function(fitmodel){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_theta <- fitmodel$infos$ncov_theta
  speciesNames <- fitmodel$infos$speciesNames
  detCovNames <- colnames(fitmodel$X_theta)
  # idxcov <- which(occCovNames == covName)

  # samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  # samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_psi]

  beta_theta_output0 <-
    matrix_of_draws[,grepl("beta_theta\\[", colnames(matrix_of_draws))]

  niter <- nrow(beta_theta_output0)

  beta_theta_output <- array(NA, dim = c(niter, ncov_theta, S))
  for(iter in 1:niter){
    beta_theta_output[iter,,] <- matrix(beta_theta_output0[iter,], ncov_theta, S, byrow = F)
  }

  dimnames(beta_theta_output)[[2]] <- detCovNames
  dimnames(beta_theta_output)[[3]] <- speciesNames

  beta_theta_output

  # data_plot <- apply(samples_subset, 2, function(x) {
  #   quantile(x, probs = c(0.025, 0.975))
  # }) %>%
  #   t %>%
  #   as.data.frame %>%
  #   mutate(Species = speciesNames) %>%
  #   mutate(speciesOrder = order(`2.5%`)) %>%
  #   filter(Species %in% speciesNames[idx_species])
  #
  # orderSpecies <- order(data_plot$`2.5%`)
  #
  # plot_occcovs <- data_plot %>%
  #   ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
  #              ymin = `2.5%`,
  #              ymax = `97.5%`)) + geom_errorbar() +
  #   xlab("Species") +
  #   # ylim(c(0,1)) +
  #   ggtitle(covName) +
  #   theme_bw() +
  #   theme(
  #     axis.text = element_text(angle = 90,
  #                              size = 8),
  #     plot.title = element_text(hjust = .5,
  #                               size = 15)
  #   ) + geom_hline(aes(yintercept = 0), color = "red")
  #
  # plot_occcovs

}

#' plotDetectionCovariates
#'
#' Detection covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the detection covariates coefficients
#'
#' @param fitmodel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotDetectionCovariates <- function(fitmodel,
                                    covName = NULL,
                                    idx_species = NULL
){

  beta_theta_output <- returnDetectionCovariates(fitmodel)

  if(is.null(covName)){
    stop("No name provided")
  }

  # matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_theta <- fitmodel$infos$ncov_theta
  speciesNames <- fitmodel$infos$speciesNames
  detCovNames <- colnames(fitmodel$X_theta)
  idxcov <- which(detCovNames == covName)

  if(length(idxcov) == 0){
    stop("Covariate name not found. If you are using a categorical covariates,
         the name might have changed to code the level. Use
         colnames(fitmodel$X_det) to find the new names")
  }

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  samples_subset <- matrix(beta_theta_output[,idxcov, idx_species],
                           dim(beta_theta_output)[1], length(idx_species))

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames[idx_species]) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  plot_occcovs <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + geom_hline(aes(yintercept = 0), color = "red")

  plot_occcovs

}


plotSpeciesRates <- function(data_plot,
                             orderSpecies,
                             subset){

  data_plot %>%
    filter(Species %in% speciesNames[orderSpecies[subset]]) %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Collection rates") +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    )

}


#' returnOccupancyProbs
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Returns the 95% credible interval of the baseline occupancy rates
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' returnOccupancyRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnOccupancyRates <- function(fitmodel){


  matrix_of_draws <- fitmodel$matrix_of_draws

  # conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

  S <- fitmodel$infos$S
  speciesNames <- fitmodel$infos$speciesNames
  n <- length(fitmodel$infos$siteNames)

  logit_psi_samples <- matrix_of_draws[,grepl("logit_psi\\[", colnames(matrix_of_draws))]

  niter <- nrow(logit_psi_samples)

  logit_psi_samples_array <- array(logistic(logit_psi_samples), dim = c(niter, n, S))

  # psi_quantiles <- apply(logit_psi_samples_array, c(2,3), function(x){
  #   quantile(logistic(x), probs = conflevels)
  # })

  # dimnames(psi_quantiles)[[2]] <- fitmodel$infos$siteNames
  # dimnames(psi_quantiles)[[3]] <- fitmodel$infos$speciesNames

  logit_psi_samples_array

  # psi_quantiles

}

#' plotOccupancyRates
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline occupancy rates
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotSpeciesRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyRates <- function(fitmodel,
                               idx_species = NULL){

  psi_output <- computeOccupancyProbs(fitmodel)
  # matrix_of_draws <- fitmodel$matrix_of_draws
  #
  # beta0_psi_output <-
  #   matrix_of_draws[,grepl("beta0_psi\\[", colnames(matrix_of_draws))]
  # U_output <-
  #   matrix_of_draws[,grepl("U\\[", colnames(matrix_of_draws))]
  # L_output <-
  #   matrix_of_draws[,grepl("LL\\[", colnames(matrix_of_draws))]
  #
  X_psi <- fitmodel$X_psi
  ncov_psi <- ncol(X_psi)

  niter <- dim(psi_output)[1]
  S <- fitmodel$infos$S
  d <- fitmodel$infos$d
  n <- length(fitmodel$infos$siteNames)
  speciesNames <- fitmodel$infos$speciesNames
  #
  # beta0psiUL_output <- array(NA, dim = c(niter, n, S))
  #
  # for (iter in 1:niter) {
  #   beta0psiUL_output[iter,,] <-
  #     matrix(beta0_psi_output[iter,], n, S, byrow = T) +
  #     X_psi %*% matrix(beta_psi_output[iter,], ncov_psi, S) +
  #     matrix(U_output[iter,], n, d, byrow = F) %*% matrix(L_output[iter,], d, S)
  # }


  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  # UL_mean <- apply(beta0psiUL_output, 3, mean)


  psi_mean_output <- apply(psi_output, c(1,3), function(x){
    mean(x)
  })

  data_plot <- apply(psi_mean_output, 2, function(x) {
    # quantile(logistic(x), probs = c(0.025, 0.975))
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)


  plot_occupancyrates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Baseline Occupancy rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("") +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  plot_occupancyrates



}

#' plotCollectionRates
#'
#' Baseline detection rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline detection rates
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotSpeciesRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCollectionRates <- function(fitmodel,
                                idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  ncov_theta <- fitmodel$infos$ncov_theta
  speciesNames <- fitmodel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "beta_theta"

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  samples_subset <- samples_subset[,1 + 0:(S-1)*ncov_theta]


  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(logistic(x), probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  plot_collectionrates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Collection rates") +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  plot_collectionrates

}

#' plotFPTPStage2Rates
#'
#' True and false positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the true and false positives at the lab stage
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotFPTPStage2Rates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotFPTPStage2Rates <- function(fitmodel,
                               idx_species = NULL,
                               primerName = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  # ncov_theta <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  primerNames <- fitmodel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  if(is.null(primerName)){
    primerName <- primerNames[1]
  }

  p_output <- matrix_of_draws[,grepl("p\\[", colnames(matrix_of_draws))]
  q_output <- matrix_of_draws[,grepl("q\\[", colnames(matrix_of_draws))]

  data_plot_p <- apply(p_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    rename(p1 = `2.5%`,
           p2 = `97.5%`)

  data_plot_q <- apply(q_output, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    rename(q1 = `2.5%`,
           q2 = `97.5%`)

  texts <- rownames(data_plot_p)
  idx_speciesprimer <- stringr::str_match(texts, "\\[(\\d+),(\\d+)\\]")

  data_plot <- cbind(data_plot_p, data_plot_q) %>%
    mutate(Species = as.numeric(idx_speciesprimer[,3]),
           Primer = as.numeric(idx_speciesprimer[,2])) %>%
    mutate(Species = speciesNames[Species],
           Primer = primerNames[Primer]) %>%
    mutate(speciesOrder = order(p1)) %>%
    filter(Species %in% speciesNames[idx_species]) %>%
    filter(Primer == primerName)

  # orderSpecies <- order(data_plot$`2.5%`[data_plot$Primer == data_plot$Primer[1]])

  detectionRates <- data_plot %>%
    ggplot()  +
    geom_errorbar(aes(x = factor(Species, level = speciesNames[speciesOrder]),
                      # factor(Species, level = speciesNames[orderSpecies]),
                      ymin = p1,
                      ymax = p2,
                      color = "TP rate"), position = position_dodge(width = .15), # Use the SAME width as geom_col
                  width = .5) +
    geom_errorbar(aes(x = factor(Species, level = speciesNames[speciesOrder]),
                      # factor(Species, level = speciesNames[orderSpecies]),
                      ymin = q1,
                      ymax = q2,
                      color = "FP rate"), position = position_dodge(width = .15), # Use the SAME width as geom_col
                  width = .5) +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Detection rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("Detection probability") +
    scale_color_manual(
      name = "Colour",
      values = c("TP rate" = "blue", "FP rate" = "red")
    ) +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  detectionRates

}

#' plotDetectionRates
#'
#' True positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the true positives at the lab stage
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotDetectionRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotDetectionRates <- function(fitmodel,
                               idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  # ncov_theta <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  primerNames <- fitmodel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "p\\["

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  texts <- rownames(data_plot)
  idx_speciesprimer <- stringr::str_match(texts, "\\[(\\d+),(\\d+)\\]")

  data_plot <- data_plot %>%
    mutate(Species = as.numeric(idx_speciesprimer[,3]),
           Primer = as.numeric(idx_speciesprimer[,2])) %>%
    mutate(Species = speciesNames[Species],
           Primer = primerNames[Primer]) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  # orderSpecies <- order(data_plot$`2.5%`[data_plot$Primer == data_plot$Primer[1]])

  detectionRates <- data_plot %>%
    ggplot(aes(x =
                 factor(Species, level = speciesNames),
                 # factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`,
               color = factor(Primer),
               group = factor(Primer))) +
    geom_errorbar(position = position_dodge(width = .15), # Use the SAME width as geom_col
                  width = .5) +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Detection rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("p") +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  detectionRates

}

#' plotStage1FPRates
#'
#' False positives rates at the field stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives rate at the field stage
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotStage1FPRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotStage1FPRates <- function(fitmodel,
                              idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  # ncov_theta <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  primerNames <- fitmodel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "theta0\\["

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  data_plot <- data_plot %>%
    mutate(Species = speciesNames) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`)

  detectionRates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNames[orderSpecies]),
               ymin = `2.5%`,
               ymax = `97.5%`)) +
    geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Stage 1 FP rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("") +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  detectionRates

}

#' plotStage2FPRates
#'
#' False positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives at the lab stage
#'
#' @param fitmodel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotStage2FPRates(fitmodel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotStage2FPRates <- function(fitmodel,
                              idx_species = NULL){

  matrix_of_draws <- fitmodel$matrix_of_draws

  S <- fitmodel$infos$S
  # ncov_theta <- fitmodel$infos$ncov_psi
  speciesNames <- fitmodel$infos$speciesNames
  primerNames <- fitmodel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  param <- "q\\["

  samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame

  texts <- rownames(data_plot)
  idx_speciesprimer <- stringr::str_match(texts, "\\[(\\d+),(\\d+)\\]")

  data_plot <- data_plot %>%
    mutate(Species = as.numeric(idx_speciesprimer[,3]),
           Primer = as.numeric(idx_speciesprimer[,2])) %>%
    mutate(Species = speciesNames[Species],
           Primer = primerNames[Primer]) %>%
    mutate(speciesOrder = order(`2.5%`)) %>%
    filter(Species %in% speciesNames[idx_species])

  orderSpecies <- order(data_plot$`2.5%`[data_plot$Primer == data_plot$Primer[1]])

  detectionRates <- data_plot %>%
    ggplot(aes(x =
                 # factor(Species, level = speciesNames[orderSpecies]),
                 factor(Species, level = speciesNames),
               ymin = `2.5%`,
               ymax = `97.5%`,
               color = factor(Primer))) +
    geom_errorbar(position = position_dodge(width = .15), # Use the SAME width as geom_col
                  width = .5) +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Stage 2 FP rates") +
    theme_bw() +
    # ylim(c(0,1)) +
    ylab("q") +
    theme(
      axis.text = element_text(angle = 0,
                               size = 8),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + coord_flip()

  detectionRates

}

#' plotReadIntensity
#'
#' Plot the reads distribution under the true positives and false positives
#'
#' @details
#' Plots the 95% credible interval of density of true positives and false positives
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotReadIntensity(fitmodel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotReadIntensity <- function(fitmodel){

  matrix_of_draws <- fitmodel$matrix_of_draws

  mu1_output <-
    matrix_of_draws[,grepl("mu1", colnames(matrix_of_draws))]
  sigma0_output <-
    matrix_of_draws[,grepl("sigma0", colnames(matrix_of_draws))]
  sigma1_output <-
    matrix_of_draws[,grepl("sigma1", colnames(matrix_of_draws))]

  niter <- length(mu1_output)

  # x_grid <- seq(1, fitmodel$infos$maxexplogy1, by = 5)
  x_grid <- exp(seq(log(1), log(fitmodel$infos$maxexplogy1), length.out = 250))

  # seq(1, fitmodel$infos$maxexplogy1, by = 5)

  x <- log(x_grid + 1)

  densities_plot_pos <- matrix(NA, length(x_grid), niter)
  densities_plot_neg <- matrix(NA, length(x_grid), niter)

  for (iter in 1:niter) {
    mu1 <- mu1_output[iter]
    sigma1 <- sigma1_output[iter]
    sigma0 <- sigma0_output[iter]

    densities_plot_pos[,iter] <- dnorm(x, mean = mu1, sd = sigma1)
    densities_plot_neg[,iter] <- dnorm(x, mean = 0, sd = sigma0)
  }

  densities_plot_pos_quantiles <-
    apply(densities_plot_pos, 1,
          function(x) {quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))}) %>% t %>%
    as.data.frame %>%
    mutate(x = x_grid,
           Type = "True Positive")

  densities_plot_neg_quantiles <-
    apply(densities_plot_neg, 1,
          function(x) {quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))}) %>% t %>%
    as.data.frame %>%
    mutate(x = x_grid,
           Type = "False Positives")

  densities_plot_quantiles <-
    rbind(densities_plot_pos_quantiles,
          densities_plot_neg_quantiles)

  # ggplot() +
  #   geom_ribbon(data = densities_plot_pos_quantiles,
  #               aes(x = x_grid,
  #                   ymax = `97.5%`,
  #                   ymin = `2.5%`)) +
  #   geom_line(data = densities_plot_pos_quantiles,
  #             aes(x = x_grid,
  #                 y = `50%`))
  #
  # df0 <- as_tibble(densities_plot_neg) %>%
  #   mutate(x = x_grid) %>%
  #   pivot_longer(cols = -x, names_to = "iter", values_to = "density") %>%
  #   mutate(Type = "False positives")
  #
  # df1 <- as_tibble(densities_plot_pos) %>%
  #   mutate(x = x_grid) %>%
  #   pivot_longer(cols = -x, names_to = "iter", values_to = "density") %>%
  #   mutate(Type = "True positives")
  #
  # # Combine into one data frame
  # df_combined <- bind_rows(df0, df1) %>%
  #   mutate(iter = as.numeric(gsub("V", "", iter)))  # Clean iteration labels

  x_grid_breaks <- round( c(0, 10, 20,
                            x_grid[seq(1, length(x_grid), by = 10)] - 1))

  # ggplot(df_combined, aes(x = x, y = density, group = interaction(iter, Type), color = Type)) +
  #   geom_line(alpha = 0.1, aes(color = Type)) +
  #   scale_color_manual(values = c("False positives" = "blue", "True positives" = "red")) +
  #   labs(title = "Reads distributions",
  #        x = "x", y = "Density") +
  #   theme_minimal() +
  #   scale_x_continuous(
  #     name = "Number of reads",
  #     breaks = x_grid_breaks,  # Custom breaks (e.g., e⁻², e⁻¹, e⁰, e¹, ...)
  #     # labels = function(x) sprintf("%.2f", exp(x) - 1)  # Format labels
  #     trans = "log"
  #   ) +
  #   theme(axis.text.x = element_text(angle = 90),
  #         plot.title = element_text(hjust = 0.5,
  #                                   size = 16,
  #                                   face = "bold"))


  ggplot() +
    geom_ribbon(data = densities_plot_quantiles,
                aes(x = x,
                    ymax = `97.5%`,
                    ymin = `2.5%`,
                    fill = Type),
                alpha = .5) +
    geom_ribbon(data = densities_plot_quantiles,
                aes(x = x,
                    ymax = `90%`,
                    ymin = `10%`,
                    fill = Type),
                alpha = .75) +
    theme_minimal() +
    scale_x_continuous(
      name = "Number of reads",
      breaks = x_grid_breaks,  # Custom breaks (e.g., e⁻², e⁻¹, e⁰, e¹, ...)
      # labels = function(x) sprintf("%.2f", exp(x) - 1)  # Format labels
      trans = "log"
    ) +
    theme(axis.text.x = element_text(angle = 90,
                                     size = 12),
          axis.text.y = element_text(angle = 90,
                                     size = 12),
          plot.title = element_text(hjust = 0.5,
                                    size = 16,
                                    face = "bold"))

}

generateCorrelationMatrixOutput <- function(fitmodel,
                                            idx_species = NULL){

  matrix_of_draws <- fitmodel$results_output

  L_output <- matrix_of_draws$LL_output
  beta0_psi_output <- matrix_of_draws$beta_psi_output
  S <- fitmodel$infos$S
  d <- fitmodel$infos$d

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  niter <- nrow(L_output)

  Lambda_output <- array(NA, dim = c(niter, S, S))

  for (iter in 1:niter) {
    L_output_current <- matrix(L_output[iter,], S, d, byrow = T)
    beta0psi_output_current <- matrix(beta0_psi_output[iter,], S, 1, byrow = T)
    L_all_output_current <- cbind(L_output_current, beta0psi_output_current)

    Lambda_output[iter,,] <- cov2cor(L_all_output_current %*% t(L_all_output_current))
  }

  Lambda_output[,idx_species, idx_species]

}




#' plotCorrelationMatrix
#'
#' Plot the correlation matrix
#'
#' @details
#' Plots the 95% credible interval of density of true positives and false positives
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotCorrelationMatrix(fitmodel)
#' }
#'
#' @export
#' @import dplyr
#' @importFrom ggcorrplot ggcorrplot
#'
plotCorrelationMatrix <- function(fitmodel,
                                  idx_species = NULL){

  Lambda_output <- generateCorrelationMatrixOutput(fitmodel, idx_species)

  Lambda_quantiles <- apply(Lambda_output, c(2,3),
                            function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})

  ggcorrplot::ggcorrplot(Lambda_quantiles[2,,], method = "square", type = "lower",
                         lab = F, lab_size = 3,
                         colors = c("blue", "white", "red"),
                         title = "Covariance Matrix (as Correlation)") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16,
                                    face = "bold"))


}

#' plotSigElementsCorMatrix
#'
#' Plot the significant element of the Correlation matrix
#'
#' @details
#' Plots the 95% credible interval of density of true positives and false positives
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotSigElementsCorMatrix(fitmodel)
#' }
#'
#' @export
#' @import dplyr
#' @importFrom ggcorrplot ggcorrplot
#'
plotSigElementsCorMatrix <- function(fitmodel,
                                     idx_species = NULL){

  Lambda_output <- generateCorrelationMatrixOutput(fitmodel, idx_species)

  S <- fitmodel$infos$S

  Lambda_quantiles <- apply(Lambda_output, c(2,3),
                            function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  Lambda_indexes_possign <- which(Lambda_quantiles[1,,] > 0, arr.ind = T)

  Lambda_indexes_negsign <- which(Lambda_quantiles[3,,] < 0, arr.ind = T)

  Lambda_sign <- matrix(0, length(idx_species), length(idx_species))
  Lambda_sign[Lambda_indexes_possign] <- 1
  Lambda_sign[Lambda_indexes_negsign] <- -1

  rownames(Lambda_sign) <- fitmodel$infos$speciesNames[idx_species]
  colnames(Lambda_sign) <- rownames(Lambda_sign)


  #
  #   ggcorrplot2(Lambda_sign, method = "square", type = "lower",
  #              lab = F, lab_size = 3,
  #              colors = c("blue", "white", "red"),
  #              title = "Significant correlations") +
  #     theme(plot.title = element_text(hjust = 0.5,
  #                                     size = 16,
  #                                     face = "bold"))


  ggcorrplot(
    # Lambda_quantiles[2,,],
    Lambda_sign,
    method = "square", type = "lower",
    lab = F, lab_size = 3, insig = "blank",
    colors = c("blue", "white", "red"),
    title = "Significant correlations") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16,
                                    face = "bold")) +
    theme(legend.position = "none") #+ theme(
  # axis.text.x = element_blank(),
  # axis.text.y = element_blank()
  # )


}

#' computeOccupancyProbs
#'
#' Computes the quantiles of the occupancy probability
#'
#' @details
#' Compute the credible interval of the occupancy probability
#'
#' @param fitmodel Output from the function runOccPlus
#'
#' @return An array with the quantiles
#'
#' @examples
#' \dontrun{
#' computeOccupancyProbs(fitmodel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
computeOccupancyProbs <- function(fitmodel#,
                                  # confidence = .95
                                  ){


  matrix_of_draws <- fitmodel$matrix_of_draws

  X_psi <- fitmodel$X_psi
  ncov_psi <- ncol(X_psi)

  beta0_psi_output <-
    matrix_of_draws[,grepl("beta0_psi\\[", colnames(matrix_of_draws))]
  beta_psi_output <-
    matrix_of_draws[,grepl("beta_psi\\[", colnames(matrix_of_draws))]
  U_output <-
    matrix_of_draws[,grepl("U\\[", colnames(matrix_of_draws))]
  L_output <-
    matrix_of_draws[,grepl("LL\\[", colnames(matrix_of_draws))]

  niter <- nrow(beta0_psi_output)
  S <- fitmodel$infos$S
  d <- fitmodel$infos$d
  n <- length(fitmodel$infos$siteNames)
  speciesNames <- fitmodel$infos$speciesNames

  psi_output <- array(NA, dim = c(niter, n, S))

  dimnames(psi_output)[[2]] <- fitmodel$infos$siteNames
  dimnames(psi_output)[[3]] <- speciesNames

  for (iter in 1:niter) {
    psi_output[iter,,] <-
      logistic(
        matrix(beta0_psi_output[iter,], n, S, byrow = T) +
          X_psi %*% matrix(beta_psi_output[iter,], ncov_psi, S) +
          matrix(U_output[iter,], n, d, byrow = F) %*% matrix(L_output[iter,], d, S)
      )
  }

  psi_output


}

# to write
plotFactorScores <- function(fitmodel,
                             idx_factor = c(1,2)){

  d <- fitmodel$d

  if(d == 0){

  }

}
