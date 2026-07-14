#' thinOutput
#'
#' Thins the MCMC output.
#'
#' @details
#' Return the same fitModel output but thinned
#'
#' @param fitModel Output from the function runOccPlus
#' @param fitModel Number of iterations to thin
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
thinOutput <- function(fitModel, thin = 5){

  niter <- dim(fitModel$results_output$beta_ord_output)[3]
  idx_thinned <- seq(1, niter, by = 5)

  results_output <- fitModel$results_output

  results_output_thinned <- lapply(1:length(results_output), function(i){

    nameElem <- names(results_output)[i]

    x <- results_output[[i]]

    if(nameElem == "z_output"){

      if(length(dim(x))== 4){

        x[,,idx_thinned,,drop=F]

      } else if(length(dim(x))== 2){

        x

      } else {

        print("Dimension not recognised")

      }

    } else {

      if(length(dim(x)) == 4){

        x[,,idx_thinned,,drop=F]

      } else if (length(dim(x)) == 3) {

        x[,idx_thinned,,drop=F]

      } else if (length(dim(x)) == 2) {

        x[idx_thinned,,drop=F]

      } else {

        print("Dimension not recognised")

      }

    }

  })

  names(results_output_thinned) <- names(results_output)

  fitModel$results_output <- results_output_thinned

  fitModel

}

logit <- function(x){
  log(x / (1 - x))
}

logistic <- function(x) 1 / (1 + exp(-x))

# TRAITS ---------

#' returnTraitsCoeff
#'
#' Traits covariate coefficients.
#'
#' @details
#' Returns the traits covariates coefficients posterior sample
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnTraitsCoeff <- function(fitModel){

  # S <- fitModel$infos$S
  # g <-
  # ncov_psi <- fitModel$infos$ncov_psi
  # # speciesNames <- fitModel$infos$speciesNames
  occCovNames <- colnames(fitModel$X_psi)
  traitNames <- colnames(fitModel$Tr)

  traitsCoeffOutput <- fitModel$results_output$jsdm_output$G_output
  traitsCoeffOutput <- apply(traitsCoeffOutput, c(1,2), c)

  niter <- dim(traitsCoeffOutput)[1]
  dimnames(traitsCoeffOutput)[[2]] <- traitNames
  dimnames(traitsCoeffOutput)[[3]] <- occCovNames

  traitsCoeffOutput <- aperm(traitsCoeffOutput, c(1,3,2))

  traitsCoeffOutput

}

#' plotTraitsCoefficients
#'
#' Traits covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitModel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_traits Indexes of the traits to be plotted (leave out to plot all the traits).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotTraitsCoefficients <- function(fitModel,
                                    covName = NULL,
                                    idx_traits = NULL
){

  traits_output <- returnTraitsCoeff(fitModel)

  if(is.null(covName)){
    stop("No name provided")
  }

  plotCoefficient(traits_output, covName, idx_traits) + xlab("Traits")

}

# OCCUPANCY COVARIATES --------

#' returnOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Returns the occupancy covariates coefficients posterior sample
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnOccupancyCovariates <- function(fitModel){

  occCovNames <- colnames(fitModel$X_psi)
  speciesNames <- fitModel$infos$speciesNames

  occCoeffOutput <- fitModel$results_output$jsdm_output$B_output
  occCoeffOutput <- apply(occCoeffOutput, c(1,2), c)

  niter <- dim(occCoeffOutput)[1]
  dimnames(occCoeffOutput)[[2]] <- occCovNames
  dimnames(occCoeffOutput)[[3]] <- speciesNames

  occCoeffOutput

}

#' plotOccupancyCovariates
#'
#' Occupancy covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the occupancy covariates coefficients
#'
#' @param fitModel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotOccupancyCovariates <- function(fitModel,
                                    covName = NULL,
                                    idx_species = NULL
){

  occCov_output <- returnOccupancyCovariates(fitModel)

  if(is.null(covName)){
    stop("No name provided")
  }

  plotCoefficient(occCov_output, covName, idx_species) + xlab("Species")

}


#' returnBaselineOccupancyRates
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Returns the 95% credible interval of the baseline occupancy rates
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' returnBaselineOccupancyRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnBaselineOccupancyRates <- function(fitModel){

  speciesNames <- fitModel$infos$speciesNames

  logisticB0_output <- logistic(fitModel$results_output$jsdm_output$B0_output)
  dimnames(logisticB0_output)[[1]] <- speciesNames

  logisticB0_output
}

#' plotBaselineOccupancyRates
#'
#' Baseline occupancy rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline occupancy rates
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#' @param confidence Confidence level of the estimate,  default to .95
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotBaselineOccupancyRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotBaselineOccupancyRates <- function(fitModel,
                               idx_species = NULL,
                               confidence = .95){

  confInt <- c((1 - confidence) / 2, (1 + confidence) / 2)

  psi0_output <- returnBaselineOccupancyRates(fitModel)

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  psi0_output_vec <- apply(psi0_output, 2, c)

  data_plot <- apply(psi0_output_vec, 1, function(x) {
    quantile(x, probs = confInt)
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = dimnames(psi0_output)[[1]])

  colnames(data_plot)[1:2] <- c("Min","Max")

  data_plot <- data_plot %>%
    mutate(speciesOrder = order(Min)) %>%
    filter(Species %in% speciesNames[idx_species])

  speciesNameOrdered <- speciesNames[order(data_plot$Min)]

  plot_occupancyrates <- data_plot %>%
    ggplot(aes(x =  factor(Species, level = speciesNameOrdered),
               ymin = Min,
               ymax = Max)) + geom_errorbar() +
    xlab("Species") +
    # ylim(c(0,1)) +
    ggtitle("Baseline Occupancy rates") +
    theme_bw() +
    ylim(c(0,1)) +
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


# COLLECTION COVARIATES --------

#' returnCollectionCovariates
#'
#' Collection covariate coefficients.
#'
#' @details
#' Returns the collection covariates coefficients posterior sample
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnCollectionCovariates <- function(fitModel){

  matrix_of_draws <- fitModel$results_output

  S <- fitModel$infos$S
  ncov_theta <- fitModel$infos$ncov_theta
  speciesNames <- fitModel$infos$speciesNames
  collcovNames <- colnames(fitModel$X_theta)
  # idxcov <- which(occCovNames == covName)

  # samples_subset <- matrix_of_draws[,grepl(param, colnames(matrix_of_draws))]
  # samples_subset <- samples_subset[,idxcov + 0:(S - 1)*ncov_psi]

  beta_theta_output <- matrix_of_draws$beta_theta_output

  beta_theta_output <- apply(beta_theta_output, c(1,2), c)
  # beta_theta_output <- aperm(beta_theta_output, c(2,3,1))

  dimnames(beta_theta_output)[[2]] <- collcovNames
  dimnames(beta_theta_output)[[3]] <- speciesNames

  beta_theta_output

}

#' plotCollectionCovariates
#'
#' Collection covariate coefficients.
#'
#' @details
#' Plots the 95% credible interval of the collection covariates coefficients
#'
#' @param fitModel Output from the function runOccPlus
#' @param covName Name of the covariate to be plotted (same name as in data$info)
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCollectionCovariates <- function(fitModel,
                                    covName = NULL,
                                    idx_species = NULL
){

  collCov_output <- returnCollectionCovariates(fitModel)

  if(is.null(covName)){
    stop("No name provided")
  }

  plotCoefficient(collCov_output, covName, idx_species) + xlab("Species")

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



#' plotCollectionRates
#'
#' Baseline collection rate for each species.
#'
#' @details
#' Plots the 95% credible interval of the baseline collection rates
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return The credible interval plot
#'
#' @examples
#' \dontrun{
#' plotSpeciesRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCollectionRates <- function(fitModel,
                                idx_species = NULL){

  S <- fitModel$infos$S
  ncov_theta <- fitModel$infos$ncov_theta
  speciesNames <- fitModel$infos$speciesNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  samples_subset <- fitModel$results_output$beta_theta_output[1,,,]
  samples_subset <- apply(samples_subset, 1, c)

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(logistic(x), probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Species = speciesNames) #%>%
    # mutate(speciesOrder = order(`2.5%`)) %>%
    # filter(Species %in% speciesNames[idx_species])

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

# SECOND STAGE RATES ----

#' plotFPTPStage2Rates
#'
#' True and false positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the true and false positives at the lab stage
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotFPTPStage2Rates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotFPTPStage2Rates <- function(fitModel,
                                idx_species = NULL,
                                primerName = NULL){

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  if(is.null(primerName)){
    primerName <- primerNames[1]
  }

  p_output <- fitModel$results_output$p_output
  q_output <- fitModel$results_output$q_output

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
    mutate(Species = speciesNames) %>%
    filter(Species %in% speciesNames[idx_species]) %>%
    mutate(speciesOrder = order(p1))

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
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotDetectionRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotDetectionRates <- function(fitModel,
                               idx_species = NULL){

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  p_output <- fitModel$results_output$p_output

  data_plot_array <- apply(p_output, c(1,2), function(x) {
    quantile(x, probs = c(0.025, 0.975))
  })

  dimnames(data_plot_array)[[2]] <- primerNames

  data_plot <- data.frame(
    lower = as.vector(data_plot_array[1,,]),
    upper = as.vector(data_plot_array[2,,]),
    Primer = rep(primerNames, times = dim(data_plot_array)[3]),
    Species = rep(speciesNames, each = dim(data_plot_array)[2])
  ) %>%
    filter(Species %in% speciesNames[idx_species]) %>%
    group_by(Species) %>%
    mutate(mean_lower = mean(lower)) %>%
    ungroup() %>%
    mutate(Species = reorder(factor(Species), mean_lower))

  plotDetectionRates <- ggplot() +
    geom_errorbar(data = data_plot, aes(x = Species,
                                 ymin = lower, ymax = upper, color = Primer)) +
    labs(
      x = "Species",
      y = "p",
      color = "Primer"
    ) +
    theme_bw() + coord_flip()

  plotDetectionRates

}

#' plotStage1FPRates
#'
#' False positives rates at the field stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives rate at the field stage
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotStage1FPRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotStage1FPRates <- function(fitModel,
                              idx_species = NULL){

  S <- fitModel$infos$S
  # ncov_theta <- fitModel$infos$ncov_psi
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  samples_subset <- fitModel$results_output$theta0_output
  samples_subset <- apply(samples_subset, 1, c)

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

  collectionRates <- data_plot %>%
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

  collectionRates

}

#' plotStage2FPRates
#'
#' False positives rates at the lab stage for each species.
#'
#' @details
#' Plots the 95% credible interval of the false positives at the lab stage
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotStage2FPRates(fitModel, idx_species = 1:5)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2
#'
plotStage2FPRates <- function(fitModel,
                              idx_species = NULL){

  S <- fitModel$infos$S
  # ncov_theta <- fitModel$infos$ncov_psi
  speciesNames <- fitModel$infos$speciesNames
  primerNames <- fitModel$infos$primerNames

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  samples_subset <- fitModel$results_output$q_output
  samples_subset <- apply(samples_subset, c(1,2), c)

  data_plot <- apply(samples_subset, 3, function(x) {
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
    ggplot(aes(x = factor(Species, level =  speciesNames[orderSpecies]),
                 # factor(Species, level = speciesNames[orderSpecies]),
                 # factor(Species, level = speciesNames),
               ymin = `2.5%`,
               ymax = `97.5%`#,
               # color = factor(Primer))
    )
    ) +
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

# CORRELATION MATRIX -----

#' plotResidualCorrelationMatrix
#'
#' Plot the residual correlation matrix
#'
#' @details
#' Plots the posterior median of the correlation matrix,
#' with nonsignificant correlation marked with an X
#'
#' @param fitModel Output from the function runOccPlus
#' @param idx_species Indexes of the species to be plotted (leave out to plot all the species).
#' @param showSignificance Should an X be shown for non significant elements?
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotResidualCorrelationMatrix(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @importFrom ggcorrplot ggcorrplot
#'
plotResidualCorrelationMatrix <- function(fitModel,
                                            idx_species = NULL,
                                            showSignificance = T,
                                            confidence = .95){

  L_output <- fitModel$results_output$jsdm_output$L_output
  speciesNames <- fitModel$infos$speciesNames

  plotCorrelationMatrix(L_output,
                        idx_species,
                        speciesNames,
                        showSignificance,
                        confidence)

}

# PREDICTIONS --------

#' computePredictiveOccupancyProbs
#'
#' Computes the quantiles of the predictive occupancy probability
#'
#' @details
#' Compute the credible interval of the occupancy probability
#'
#' @param fitModel Output from the function runOccPlus
#' @param X_psi Occupancy covariates matrix for the new locations
#' @param X_ord Ordination covariates matrix for the new locations
#' @param summarised Should the output be return in the form of quantiles? Set to TRUE if the number of sites is very large
#' @param confidence If quantiles are returned, the confidence level of the quantiles.
#'
#' @return An array of size (,sites,species) with either the quantiles or the iterations in the first dimension
#'
#' @examples
#' \dontrun{
#' computePredictiveOccupancyProbs(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
computePredictiveOccupancyProbs <- function(fitModel,
                                            X_psi,
                                            X_s,
                                            summarised = F,
                                            confidence = .95
){


  X_psi <- as.matrix(X_psi)

  S <- fitModel$infos$S
  speciesNames <- fitModel$infos$speciesNames
  n <- nrow(X_psi)

  if(is.null(X_psi)) {
    X_psi <- fitModel$X_psi
  }

  if(is.null(X_s)) {
    X_s <- fitModel$X_s
  }

  #

  beta_psi_output <- fitModel$results_output$beta_psi_output
  beta_ord_output <- fitModel$results_output$beta_ord_output
  LL_output <- fitModel$results_output$LL_output

  nchain <- dim(beta_psi_output)[4]
  niter <- dim(beta_psi_output)[3]

  if(!summarised){

    psi_output <- array(NA, dim = c(nchain * niter, n, S))
    for (chain in 1:nchain) {
      for (iter in 1:niter) {
        psi_output[iter + (chain - 1)*niter,,] <-
          logistic(
            computePsiE(X_psi, beta_psi_output[,,iter,chain], X_ord,
                        beta_ord_output[,,iter,chain],
                        LL_output[,,iter,chain])
          )
      }
    }

  } else {

    conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

    beta_ord_output <- aperm(apply(beta_ord_output, c(1,2), c), c(2,3,1))
    beta_psi_output <- aperm(apply(beta_psi_output, c(1,2), c), c(2,3,1))
    LL_output <- aperm(apply(LL_output, c(1,2), c), c(2,3,1))

    # niter <- dim(beta_ord_output)[3]

    psi_output <- computePsiOutput(
      X_psi,
      beta_psi_output,
      X_ord,
      beta_ord_output,
      LL_output,
      conflevels)

    # psi_output <- array(NA, dim = c(3, n, S))
    # for (i in 1:n) {
    #   for (j in 1:S) {
    #     mcmc_output <- rep(NA, niter)
    #     for (iter in 1:niter) {
    #       mcmc_output[iter] <- logistic(
    #         computePsiE(X_psi[i,,drop=F], beta_psi_output[,j,iter],
    #                     X_ord[i,,drop=F],
    #                     beta_ord_output[,,iter],
    #                     LL_output[,j,iter])
    #       )
    #
    #     }
    #
    #     psi_output[,i,j] <- quantile(mcmc_output, conflevels)
    #
    #   }
    # }

  }

  psi_output

  # for (iter in 1:niter) {
  #   psi_output[iter,,] <-
  #     logistic(
  #       matrix(beta0_psi_output[iter,], n, S, byrow = T) +
  #         X_psi %*% matrix(beta_psi_output[iter,], ncov_psi, S) +
  #         matrix(U_output[iter,], n, n_factors, byrow = F) %*% matrix(L_output[iter,], n_factors, S)
  #     )
  # }
  #
  # psi_output


}


#' computeConditionalOccupancyProbs
#'
#' Computes the posterior mean of the conditional occupancy probability
#'
#' @details
#' Computes the posterior mean of the conditional occupancy probability
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A matrix of size (site X species) with the posterior menan occupancy at
#' each site for each species.
#'
#' @examples
#' \dontrun{
#' computeConditionalOccupancyProbs(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
computeConditionalOccupancyProbs <- function(fitModel){

  z_output <- fitModel$results_output$z_output

  z_mean <- apply(z_output, c(1,2), mean)

  rownames(z_mean) <- fitModel$infos$siteNames
  colnames(z_mean) <- fitModel$infos$speciesNames

  z_mean

}

# OTHER ----------

#' returnLatentPresences
#'
#' Compute the latent presences
#'
#' @details
#' Compute the latent presences
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A matrix
#'
#' @examples
#' \dontrun{
#' returnLatentPresences(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnLatentPresences <- function(fitModel, idx_species = 1){

  z_output <- fitModel$results_output$z_output
  w_output <- fitModel$results_output$w_output

  z_ouput_vec <- apply(z_output, c(1,2), c)
  z_output_species <- z_ouput_vec[,,idx_species]
  z_mean <- apply(z_output_species, 2, mean)

  w_ouput_vec <- apply(w_output, c(1,2), c)
  w_output_species <- w_ouput_vec[,,idx_species]
  w_mean <- apply(w_output_species, 2, mean)

  speciesNames <- fitModel$infos$speciesNames

  returnVariancePartitioningMatrix(varPart_output, speciesNames)

}


#' returnVariancePartitioning
#'
#' Compute the variance partitioning for each species
#'
#' @details
#' Compute the variance partitioning for each species
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A matrix of size (species x 4)
#'
#' @examples
#' \dontrun{
#' returnVariancePartitioning(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
returnVariancePartitioning <- function(fitModel){

  varPart_output <- fitModel$results_output$jsdm_output$varPart_output

  speciesNames <- fitModel$infos$speciesNames

  returnVariancePartitioningMatrix(varPart_output, speciesNames)

}

#' plotVariancePartitioning
#'
#' Plot the variance partitioning for each species
#'
#' @details
#' Plot the variance partitioning for each species
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotVariancePartitioning(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @import ggtern
#'
plotVariancePartitioning <- function(fitModel){

  varPart_output <- fitModel$results_output$jsdm_output$varPart_output

  speciesNames <- fitModel$infos$speciesNames

  plotVarPart(varPart_output, speciesNames)

}

#' plotReadIntensity
#'
#' Plot the reads distribution under the true positives and false positives
#'
#' @details
#' Plots the 95% credible interval of density of true positives and false positives
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' plotReadIntensity(fitModel)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotReadIntensity <- function(fitModel){

  mu1_output <- fitModel$results_output$mu1_output
  mu0_output <- fitModel$results_output$mu0_output
  sigma1_output <- fitModel$results_output$sigma1_output
  sigma0_output <- fitModel$results_output$sigma0_output

  niter <- length(mu1_output)

  # x_grid <- seq(1, fitModel$infos$maxexplogy1, by = 5)
  x_grid <- exp(seq(log(1), log(fitModel$infos$maxexplogy1), length.out = 250))

  # seq(1, fitModel$infos$maxexplogy1, by = 5)

  x <- log(x_grid + 1)

  densities_plot_pos <- matrix(NA, length(x_grid), niter)
  densities_plot_neg <- matrix(NA, length(x_grid), niter)

  for (iter in 1:niter) {

    mu1 <- mu1_output[iter]
    mu0 <- mu0_output[iter]
    sigma1 <- sigma1_output[iter]
    sigma0 <- sigma0_output[iter]

    densities_plot_pos[,iter] <- dnorm(x, mean = mu1, sd = sigma1)
    densities_plot_neg[,iter] <- dnorm(x, mean = mu0, sd = sigma0)
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

#' plotCumulativeSpeciesDetections
#'
#' Cumulative species detection plot
#'
#' @details
#' Plot the credible interval of the overall species detections for several
#' number of PCR, conditional on species presence in the sample
#'
#' @param fitModel Output from the function runOccPlus
#'
#' @return A plot with the credible interval of cumulative detections
#'
#' @examples
#' \dontrun{
#' plotCumulativeSpeciesDetections(fitModel, K = 5)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#'
plotCumulativeSpeciesDetections <- function(fitModel, K, primer = 0, alpha = .95){

  p_output <- fitModel$results_output$p_output

  ab_p <- apply(p_output, c(1,2), function(x){
    x <- as.vector(x)
    mean_x <- mean(x)
    var_x <- var(x)

    alpha <- mean_x * ((mean_x * (1 - mean_x) / var_x) - 1)
    beta <- (1 - mean_x) * ((mean_x * (1 - mean_x) / var_x) - 1)
    c(alpha, beta)
  })

  speciesDetected <- computeSpeciesDetected(ab_p, K, primer, alpha)

  ggplot(data = NULL, aes(x = 1:K,
                          ymin = speciesDetected[1,],
                          ymax = speciesDetected[2,])) +
    geom_errorbar() + theme_bw() +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 15, face = "bold")) +
    scale_x_continuous(breaks = 1:K, name = "Number of technical replicates") +
    scale_y_continuous(name = "Species Detected")

}

computeSpeciesDetected <- function(ab_p, K, primer, alpha){

  S <- dim(ab_p)[3]

  B <- 100

  P_all <- dim(ab_p)[2]

  if(primer == 0){
    idx_primer <- 1:P_all
  } else {
    idx_primer <- primer
  }

  P <- length(idx_primer)

  matrixDetectedPCR <- sapply(1:B, function(b){

    detectionMatrix <- sapply(1:S, function(s){

      ps <- rbeta(P, ab_p[1,idx_primer,s], ab_p[2,idx_primer,s])
      sapply(1:K, function(k){

        detectionsAcrossPrimers <-
          sapply(1:P, function(p){
            rbinom(1, 1, ps[p])
          })

        as.numeric(any(detectionsAcrossPrimers) > 0)

      })
    })

    detectionsPerPcr <- sapply(1:K, function(k){
      sum(apply(detectionMatrix[1:k,,drop=F], 2, function(x){ any(x > 0)}))
    })

    detectionsPerPcr
  })

  apply(matrixDetectedPCR, 1, function(x){
    quantile(x, probs = c((1 - alpha)/2, (1 + alpha)/2))
  })

}

## TODO
plotOrdination <- function(fitModel,
                           idx_factor = c(1,2)){

  n_factors <- fitModel$infos$n_factors

  if(n_factors== 0){
    stop("No factor used")
  }

  if(n_factors > 2){
    print("More than 2 factors present, the ordination plot will use the first
          two factors only")
  }

  U_output <- fitModel$results_output$U_output

  plot_data <- as.data.frame.table(U_output) %>%
    rename(Chain = Var1, Iter = Var2, Obs = Var3, Dim = Var4) %>%
    mutate(Obs = as.numeric(Obs)) %>%
    group_by(Obs, Dim) %>%
    summarise(
      mean_val = mean(Freq),
      lower    = quantile(Freq, 0.025),
      upper    = quantile(Freq, 0.975),
      .groups  = "drop"
    ) %>%
    # Pivot wider so Dim 1 and Dim 2 are in separate columns for 2D plotting
    pivot_wider(
      names_from = Dim,
      values_from = c(mean_val, lower, upper),
      names_sep = "_d"
    )

  # --- 3. Plot with ggplot2 ---
  ggplot(plot_data, aes(x = mean_val_d1, y = mean_val_d2)) +
    # Horizontal error bars (Uncertainty in Dimension 1)
    geom_errorbarh(aes(xmin = lower_d1, xmax = upper_d1), color = "gray60", alpha = 0.7) +
    # Vertical error bars (Uncertainty in Dimension 2)
    geom_errorbar(aes(ymin = lower_d2, ymax = upper_d2), color = "gray60", alpha = 0.7) +
    # Central estimate points
    geom_point(color = "firebrick", size = 2) +
    labs(
      title = "Observation Estimates with 95% Credible Intervals",
      x = "Dimension 1",
      y = "Dimension 2"
    ) +
    theme_minimal()

}

plotOccupancyStates <- function(fitModel){

  speciesNames <- fitModel$infos$speciesNames
  siteNames <- fitModel$infos$siteNames

  z_output <- fitModel$results_output$z_output
  z_mean <- apply(z_output, c(1,2), mean)
  dimnames(z_mean) <- list(siteNames, speciesNames)

  observedOccupancies <- cbind(data_info, OTU) %>%
    group_by(Site, Sample) %>%
    summarise(
      across(starts_with("OTU_"), function(x) sum(x > 0))
    ) %>%
    group_by(Site) %>%
    summarise(
      across(starts_with("OTU_"), function(x) round(mean(x > 0), 2))
    )  %>% as.data.frame %>% dplyr::select(-Site) %>% as.matrix
  dimnames(observedOccupancies) <- list(siteNames, speciesNames)

  df_colors  <- as.data.frame.table(z_mean) %>%
    rename(Site = Var1, Species = Var2, Occupancy = Freq)

  df_numbers <- as.data.frame.table(observedOccupancies) %>%
    rename(Site = Var1, Species = Var2, Frequency = Freq)

  plot_data <- left_join(df_colors, df_numbers, by = c("Site", "Species")) #%>%

  ggplot(plot_data, aes(x = Species, y = Site)) +
    geom_tile(aes(fill = Occupancy), color = "white", linewidth = 0.5) +
    geom_text(aes(label = Frequency), color = "black", fontface = "bold", size = 4) +
    scale_fill_viridis_c(option = "plasma", name = "Occupancy") +
    # scale_y_reverse() +
    # scale_x_continuous(breaks = 1:5, position = "top") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 12, face = "bold")
    )

}

# DEPRECATED


# plotCorrelationMatrix <- function(fitModel,
#                                   idx_species = NULL){
#
#   Lambda_output <- generateCorrelationMatrixOutput(fitModel, idx_species)
#
#   Lambda_quantiles <- apply(Lambda_output, c(2,3),
#                             function(x){quantile(x, probs = c(0.025, 0.5, 0.975))})
#
#   ggcorrplot::ggcorrplot(Lambda_quantiles[2,,],
#                          method = "square",
#                          lab = F, lab_size = 3,
#                          colors = c("blue", "white", "red"),
#                          title = "Covariance Matrix (as Correlation)") +
#     theme(plot.title = element_text(hjust = 0.5,
#                                     size = 16,
#                                     face = "bold"))
#
#
# }


