#' runOccPlus
#'
#' Run the OccPlus model.
#'
#' @details
#' Fit the model described in Diana et al..
#'
#' @param data A list with two elements:
#' \describe{
#'   \item{info}{A data.frame of N rows (where N is the total number of samples
#'   analyzed). Further information in the vignette.}
#'   \item{OTU}{A matrix of dimension (N x S), where S is the number of species,
#'   with the number of reads of each species in each sample.}
#' }
#' @param n_factors Number of factors.
#' @param threshold (Optional) threshold used to truncated the reads to binary data.
#' Data greater than or equal to the threshold will be considered as detection.
#' If this is set to 0 (which is the default), occPlus estimates two modes as described
#' in the paper.
#' @param occCovariates vector of the name of the covariates for the occupancy probabilities.
#' Names should match the column name in data$info.
#' @param ordCovariates vector of the name of the covariates for the occupancy probabilities.
#' Names should match the column name in data$info.
#' @param collCovariates vector of the name of the covariates for the collection probabilities.
#' Names should match the column name in data$info.
#'
#' @return Description of the return value (e.g., a list, data frame, or numeric output).
#'
#' @examples
#' \dontrun{
#' # Example usage
#' fitmodel  <- runOccPlus(data,
#' n_factors = 2,
#' occCovariates = c("X_psi.1","X_psi.2"),
#' ordCovariates = c("X_ord.1","X_ord.2"),
#' collCovariates = c("X_theta"))
#' }
#'
#' @export
#' @import dplyr
#'
runOccPlus <- function(data,
                       n_factors,
                       threshold = 0,
                       occCovariates = c(),
                       ordCovariates = c(),
                       collCovariates = c(),
                       MCMCparams = list(nchain = 2,
                                         nburn = 5000,
                                         niter = 5000,
                                         nthin = 1)){

  {
    # n_factors = 2
    # threshold = 1
    # occCovariates = c()
    # ordCovariates = c("season")
    # collCovariates = c("genetic_predator")
    # MCMCparams = list(nchain = 1,
    #                   nburn = 5000,
    #                   niter = 5000,
    #                   thin = 1)
  }

  data_info <- as.data.frame(data$info)
  OTU <- data$OTU

  # data checks
  {
    if(!all(c(occCovariates, ordCovariates, collCovariates) %in% colnames(data$info))){
      stop("Covariate names provided not in data$info")
    }

    if(any(is.na(data_info$Site)) | any(is.na(data_info$Sample)) | any(is.na(data_info$Primer))){
      stop("NA in Site, Sample or Primer columns")
    }

    if(n_factors > ncol(OTU)){
      print("More species than factors. The number of factors will be capped to the
            number of species")
      n_factors = ncol(OTU)
    }
  }

  # sort the data
  {
    data_info <- data_info %>%
      dplyr::arrange(Site, Sample, Primer)
  }

  # samples per site
  {
    M_df <- data_info %>%
      dplyr::group_by(Site, Sample) %>%
      slice(1) %>%
      dplyr::group_by(Site) %>%
      dplyr::summarise(M = n())

    M <- M_df$M
    names(M) <- M_df$Site

    n <- length(M)

    sumM <- c(0, cumsum(M)[-n])

    siteNames <- unique(data_info$Site)

    data_info_sample <- as.numeric(factor(data_info$Sample, levels = unique(data_info$Sample)))

  }

  # marker per samples
  {
    # number of markers

    L_df <- data_info %>%
      dplyr::group_by(Site, Sample, Primer) %>%
      dplyr::slice(1) %>%
      dplyr::group_by(Site, Sample) %>%
      dplyr::summarise(L = n())

    L <- L_df$L
    names(L) <- L_df$Sample

    # number of observations
    {
      L_all <- data_info %>%
        dplyr::group_by(Site, Sample, Primer) %>%
        dplyr::slice(1) %>%
        dplyr::group_by(Site, Sample) %>%
        dplyr::summarise(L_m = n())
      L_all <- L_all$L_m

      sumL <- c(0, cumsum(L_all)[-length(L_all)])

    }

    primerNames <- unique(data_info$Primer)
    # sumL <- c(0, cumsum(L)[-length(L)])

  }

  # pcr per marker
  {
    data_K <- data_info %>%
      dplyr::group_by(Site, Sample, Primer) %>%
      dplyr::summarise(K = n(),
                       dplyr::across(contains("Species"),function(x){sum(x > 0)})
      ) %>%
      dplyr::ungroup()#%>%
    # ungroup() %>%
    # mutate(across(contains("Species"), function(x){sum(x > 0)}))

    K <- data_K$K

    numL <- length(K)

    sumK <- c(0, cumsum(K)[-numL])

    maxL = max(L)

    idx_z <- rep(1:n, M)
    idx_k <- as.numeric(as.factor(data_info_sample))

    primerIdx <- as.numeric(as.factor(data_info$Primer))
  }

  y <- OTU

  logy1 <- log(OTU + 1)

  # data infos
  {
    n <- length(M)
    S <- ncol(logy1)
    N <- sum(M)
    N2 <- numL
    N3 <- nrow(logy1)

    speciesNames <- colnames(data$OTU)
    if(is.null(speciesNames)){
      speciesNames <- 1:S
    }

  }

  # create delta
  {
    {
      M_marker_df <- data_info %>%
        dplyr::group_by(Site, Sample) %>%
        dplyr::summarise(M = n())

      M_marker <- M_marker_df$M
      # names(M) <- M_df$Site

      sumM_marker <- c(0, cumsum(M_marker)[-N])

    }

    # delta <- matrix(NA, N, S)
    #
    # for (s in 1:S) {
    #
    #   for (i in 1:n) {
    #
    #     for (m in 1:M[i]) {
    #
    #       delta[sumM[i] + m,s] <-
    #         as.numeric(all(OTU[sumM_marker[sumM[i] + m] + 1:M_marker[sumM[i] + m], s] == 0))
    #
    #     }
    #
    #   }
    #
    # }
    #
    # delta[is.na(delta)] <- 0

  }

  # create covariates matrix
  {

    if(is.null(dim(OTU))){
      S <- 1
    } else {
      S <- ncol(OTU)
    }

    if(length(ordCovariates) > 0){

      X_ord <- data_info %>%
        dplyr::group_by(Site) %>%
        dplyr::summarise(across(all_of(ordCovariates),
                                function(x) {x[1]}))

      sitesNames <- X_ord$Site

      X_ord <- X_ord %>%
        dplyr::select(-Site) %>%
        dplyr::mutate_if(is.numeric, scale) %>%
        dplyr::mutate(dplyr::across(tidyselect::where(~ !is.numeric(.x)), as.factor)) %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
        stats::model.matrix(~., .)

      X_ord <- X_ord[,-1,drop=F]

      rownames(X_ord) <- sitesNames

    } else {
      X_ord <- matrix(0, n, 0)
    }

    if(length(occCovariates) > 0){

      X_psi <- data_info %>%
        dplyr::group_by(Site) %>%
        dplyr::summarise(dplyr::across(all_of(occCovariates),
                                       function(x) {x[1]})) %>%
        dplyr::select(-Site) %>%
        dplyr::mutate_if(is.numeric, scale) %>%
        dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor)) %>%
        model.matrix(~., .)

      # X_psi <- X_psi[,-1,drop=F]

    } else {

      X_psi <- matrix(1, n, 1)

    }

    if(length(collCovariates) > 0){

      X_theta <- data_info %>%
        dplyr::group_by(Sample) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(collCovariates),
                                       function(x) {x[1]})) %>%
        dplyr::select(-Sample) %>%
        dplyr::mutate_if(is.numeric, scale) %>%
        dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor)) %>%
       model.matrix(~., .)

    } else {

      X_theta <- matrix(1, N, 1)

    }

    ncov_psi <- ncol(X_psi)
    ncov_theta <- ncol(X_theta)
    ncov_ord <- ncol(X_ord)

  }

  # check for nas
  {
    logy_na <- is.na(logy1)
    mode(logy_na) <- "integer"
    logy1[is.na(logy1)] <- -1
  }

  # priors
  {
    prior_beta_psi <- 0
    prior_beta_psi_sd <- 1
    prior_beta_theta <- 0
    prior_beta_theta_sd <- 1
    a_theta0 <- 1
    b_theta0 <- 20
    a_p <- 5
    b_p <- 1
    a_q <- 1
    b_q <- 20
    a_sigma0 <- 1
    b_sigma0 <- 5
    a_sigma1 <- 1
    b_sigma1 <- 1

    mu_mu1 <- 3
    sd_mu1 <- 3
    mu_mu0 <- 1
    sd_mu0 <- .5

    b_betatheta <- rep(1, ncov_theta)
    B_betatheta <- diag(1, nrow = ncov_theta)

    b_betatheta[1] <- prior_beta_theta
    B_betatheta[1,1] <- prior_beta_theta_sd

  }

  y <- OTU

  if(threshold != 0){

    y[OTU >= threshold] <- 1
    y[OTU < threshold] <- 0

  }

  # run MCMC

  print("Running MCMC")

  nchain <- MCMCparams$nchain
  nburn <- MCMCparams$nburn
  niter <- MCMCparams$niter

  # chain output
  {
    beta_ord_output <- array(NA, dim = c(ncov_ord, n_factors, niter, nchain))
    beta_psi_output <- array(NA, dim = c(ncov_psi, S, niter, nchain))
    beta_theta_output <- array(NA, dim = c(ncov_theta, S, niter, nchain))
    LL_output <- array(NA, dim = c(n_factors, S, niter, nchain))
    E_output <- array(NA, dim = c(n, n_factors, niter, nchain))
    U_output <- array(NA, dim = c(n, n_factors, niter, nchain))
    UL_output <- array(NA, dim = c(n, S, niter, nchain))
    z_output <- array(NA, dim = c(n, S, niter, nchain))
    p_output <- array(NA, dim = c(maxL, S, niter, nchain))
    q_output <- array(NA, dim = c(maxL, S, niter, nchain))
    theta0_output <- array(NA, dim = c(S, niter, nchain))
    mu1_output <- array(NA, dim = c(niter, nchain))
    sigma1_output <- array(NA, dim = c(niter, nchain))
    mu0_output <- array(NA, dim = c(niter, nchain))
    sigma0_output <- array(NA, dim = c(niter, nchain))
  }

  for (chain in 1:nchain) {

    # chain output
    {
      beta_ord_output_chain <- array(NA, dim = c(ncov_ord, n_factors, niter))
      beta_psi_output_chain <- array(NA, dim = c(ncov_psi, S, niter))
      beta_theta_output_chain <- array(NA, dim = c(ncov_theta, S, niter))
      LL_output_chain <- array(NA, dim = c(n_factors, S, niter))
      E_output_chain <- array(NA, dim = c(n, n_factors, niter))
      U_output_chain <- array(NA, dim = c(n, n_factors, niter))
      UL_output_chain <- array(NA, dim = c(n, S, niter))
      z_output_chain <- array(NA, dim = c(n, S, niter))
      p_output_chain <- array(NA, dim = c(maxL, S, niter))
      q_output_chain <- array(NA, dim = c(maxL, S, niter))
      theta0_output_chain <- array(NA, dim = c(S, niter))
      mu1_output_chain <- rep(NA, niter)
      sigma1_output_chain <- rep(NA, niter)
      mu0_output_chain <- rep(NA, niter)
      sigma0_output_chain <- rep(NA, niter)
    }

    # starting values
    {
      if(threshold == 0){

        trueStartingPoint <- F

        # true starting points
        if (trueStartingPoint){

          z <- z_true
          w <- delta
          theta <- theta_true
          theta0 <- theta0_true
          beta_psi <- beta_psi_true
          beta_ord <- beta_ord_true
          beta_theta <- beta_theta_true
          U <- U_true
          E <- E_true
          LL <- L_true
          p <- p_true
          q <- q_true
          c_imk <- cimk_true
          mu1 <- mu1_true
          sigma1 <- sd1_true
          mu0 <- mu0_true
          sigma0 <- sd0_true

        } else {

          c_imk <- y > 0

          w <- matrix(NA, N, S)

          for (s in 1:S) {
            for (i in 1:n) {
              for (m in 1:M[i]) {
                idx_im1 <- sumK[sumL[sumM[i] + m] + 1] + 1
                idx_im2 <- sumK[sumL[sumM[i] + m] + L[sumM[i] + m]] +
                  K[sumL[sumM[i] + m] + L[sumM[i] + m]]
                w[sumM[i] + m,s] <- as.numeric(any(y[idx_im1:idx_im2,s] > 0))
              }
            }
          }

          z <- matrix(NA, n, S)

          for (s in 1:S) {
            for (i in 1:n) {
              idx_i1 <- sumM[i] + 1
              idx_i2 <- sumM[i] + M[i]

              z[i,s] <- as.numeric(any(w[idx_i1:idx_i2, s] > 0))
            }
          }



          mu1 <- 7
          sigma1 <- 3
          mu0 <- 0
          sigma0 <- 1
        }



      } else {

        c_imk <- y > 0

        w <- matrix(NA, N, S)

        for (s in 1:S) {
          for (i in 1:n) {
            for (m in 1:M[i]) {
              idx_im1 <- sumK[sumL[sumM[i] + m] + 1] + 1
              idx_im2 <- sumK[sumL[sumM[i] + m] + L[sumM[i] + m]] +
                K[sumL[sumM[i] + m] + L[sumM[i] + m]]
              w[sumM[i] + m,s] <- as.numeric(any(y[idx_im1:idx_im2,s] > 0))
            }
          }
        }

        z <- matrix(NA, n, S)

        for (s in 1:S) {
          for (i in 1:n) {
            idx_i1 <- sumM[i] + 1
            idx_i2 <- sumM[i] + M[i]

            z[i,s] <- as.numeric(any(w[idx_i1:idx_i2, s] > 0))
          }
        }


      }

      beta_psi <- matrix(0, ncov_psi, S)
      beta_theta <- matrix(0, ncov_theta, S)
      beta_ord <- matrix(0, ncov_ord, n_factors)
      E <- matrix(0, n, n_factors)
      LL <- matrix(1, n_factors, S)


      theta <- computeTheta(X_theta, beta_theta)

      p <- matrix(.9, maxL, S)
      q <- matrix(.05, maxL, S)
      theta0 <- rep(.05, S)


      U <- computeU(X_ord, beta_ord, E)
      psi <- computePsi(X_psi, beta_psi, U, LL)

    }

    for (iter in 1:(nburn + niter * nthin)) {

      if(iter %% 10 == 0){

        if(iter > nburn){
          print(paste0("Chain ", chain, " - Iteration ",iter - nburn))
        } else {
          print(paste0("Chain ", chain, " - Burn in Iteration ",iter))
        }

      }

      # sample z

      z <- sample_z_cpp(w, psi, theta, theta0, M, sumM)

      # sample psi

      {
        list_betapsiLL <- sample_psivars(z, X_psi, beta_psi, X_ord, beta_ord,
                                         E, LL, prior_beta_psi, prior_beta_psi_sd)
        beta_psi <- list_betapsiLL$beta_psi
        beta_ord <- list_betapsiLL$beta_ord
        LL <- list_betapsiLL$LL
        E <- list_betapsiLL$E
        U <- computeU(X_ord, beta_ord, E)
        psi <- computePsi(X_psi, beta_psi, U, LL)
      }

      # sample w
      if(threshold > 0){

        w <- sample_w_cim_cipp(y, theta, theta0, p, q,
                               M, K, sumL, sumM, sumK, maxL, z)

      } else {

        w <- sample_w_cpp(logy1, mu0, sigma0, mu1, sigma1, theta, theta0, p, q,
                          M, K, sumL, sumM, sumK, maxL, z)


      }

      # sample cimk

      if (threshold > 0) {

        w_all <- w[idx_k,,drop=F]

        # faster way to assign c_imk to 1 if logy1 > 0 for w_all = 1 and to 2
        # if log1 > 0 when w_all = 0
        c_imk <- (logy1 > 0) * (2 - (w_all == 1))

      }

      # sample theta
      beta_theta <- sample_betatheta_cpp(w, z, beta_theta, idx_z, X_theta,
                                     b_betatheta, B_betatheta)
      theta <- computeTheta(X_theta, beta_theta)

      # sample pq
      # list_pq <- sample_pq(y, w, primerIdx, idx_k, maxL, a_p, b_p, a_q, b_q)
      list_pq <- sample_pq_cpp(c_imk, w, primerIdx, idx_k, maxL, a_p, b_p, a_q, b_q)
      (p <- list_pq$p)
      (q <- list_pq$q)

      # sample theta0
      theta0 <- sample_theta0(z, w, idx_z, a_theta0, b_theta0)

      if(threshold == 0){

        # sample mu1sigma1

        list_mu1sigma1 <- sample_musigma(sigma1,
                                         logy1, c_imk,
                                         mu_mu1, sd_mu1,
                                         a_sigma1, b_sigma1, tpfp = T)
        mu1 <- list_mu1sigma1$mu1
        sigma1 <- list_mu1sigma1$sigma1

        # sample sigma0
        sigma0 <- sample_sigma0(logy1, c_imk, a_sigma0, b_sigma0)

      }

      {

        if(iter > nburn & (iter - nburn) %% nthin == 0){
          currentIter <- (iter - nburn) / nthin
          beta_psi_output_chain[,,currentIter] <- beta_psi
          beta_ord_output_chain[,,currentIter] <- beta_ord
          beta_theta_output_chain[,,currentIter] <- beta_theta
          LL_output_chain[,,currentIter] <- LL
          U_output_chain[,,currentIter] <- U
          E_output_chain[,,currentIter] <- E
          p_output_chain[,,currentIter] <- p
          q_output_chain[,,currentIter] <- q
          theta0_output_chain[,currentIter] <- theta0
          z_output_chain[,,currentIter] <- z

          if(threshold == 0){
            mu1_output_chain[currentIter] <- mu1
            sigma1_output_chain[currentIter] <- sigma1
            mu0_output_chain[currentIter] <- mu0
            sigma0_output_chain[currentIter] <- sigma0
          }

        }

      }

      beta_psi_output[,,,chain] <- beta_psi_output_chain
      beta_ord_output[,,,chain] <- beta_ord_output_chain
      beta_theta_output[,,,chain] <- beta_theta_output_chain
      LL_output[,,,chain] <- LL_output_chain
      U_output[,,,chain] <- U_output_chain
      E_output[,,,chain] <- E_output_chain
      p_output[,,,chain] <- p_output_chain
      q_output[,,,chain] <- q_output_chain
      theta0_output[,,chain] <- theta0_output_chain
      z_output[,,,chain] <- z_output_chain
      mu1_output[,chain] <- mu1_output_chain
      sigma1_output[,chain] <- sigma1_output_chain
      mu0_output[,chain] <- mu0_output_chain
      sigma0_output[,chain] <- sigma0_output_chain

    }

  }

  results_output <- list(
    "beta_ord_output" = beta_ord_output,
    "beta_psi_output" = beta_psi_output,
    "beta_theta_output" = beta_theta_output,
    "LL_output" = LL_output,
    "E_output" = E_output,
    "UL_output" = UL_output,
    "U_output" = U_output,
    "z_output" = z_output, # TODO: summarise this object
    "p_output" = p_output,
    "q_output" = q_output,
    "theta0_output" = theta0_output,
    "mu1_output" = mu1_output,
    "sigma1_output" = sigma1_output,
    "mu0_output" = mu0_output,
    "sigma0_output" = sigma0_output
    )

  infos <- list(
    "S" = S,
    "M" = M,
    "n" = n,
    "K" = K,
    "data_info" = data_info,
    "speciesNames" = speciesNames,
    "primerNames" = primerNames,
    "siteNames" = siteNames,
    "n_factors" = n_factors,
    "ncov_theta" = ncov_theta,
    "ncov_psi" = ncov_psi,
    "ncov_ord" = ncov_ord,
    "maxexplogy1" = max(exp(logy1), na.rm = T)
  )

  list(
    "results_output" = results_output,
    "infos" = infos,
    "X_ord" = X_ord,
    "X_theta" = X_theta,
    "X_psi" = X_psi)

}
