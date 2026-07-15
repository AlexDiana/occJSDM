get_param <- function(params, key, default = 0) {
  if (!is.null(params[[key]])) {
    return(params[[key]])
  } else {
    return(default)
  }
}

process_covariates <- function(data_info, covariates, group_by_col, n_obs,
                               remove_intercept = FALSE) {

  # If no covariates provided, return default matrix
  if (length(covariates) == 0) {
    if (!remove_intercept) {
      return(matrix(1, n_obs, 1))
    } else {
      return(matrix(0, n_obs, 0))
    }
  }

  # Process the covariates
  result <- data_info %>%
    dplyr::group_by(!!rlang::sym(group_by_col)) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(covariates), ~ dplyr::first(.x))) %>%
    dplyr::select(-dplyr::all_of(group_by_col)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), scale)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
    stats::model.matrix(~., .)

  # Remove intercept if requested (first column is always intercept)
  if (remove_intercept) {
    result <- result[, -1, drop = FALSE]
  }

  return(result)
}

createDataIdx <- function(n, M, P, K){

  N <- sum(M)
  N2 <- P * N
  N3 <- sum(K)

  sumM <- c(0, cumsum(M)[-n])
  sumP <- c(0, cumsum(rep(P, N))[-N])
  sumK <- c(0, cumsum(K)[-N2])

  idx_z_w <- rep(NA, N)
  idx_z_p <- rep(NA, N2)
  idx_w_p <- rep(NA, N2)
  idx_z_k <- rep(NA, N3)
  idx_p_k <- rep(NA, N3)
  idx_w_k <- rep(NA, N3)

  idx_z <- 1
  idx_w <- 1
  idx_p <- 1
  idx_k <- 1

  for (i in 1:n) {
    for (m in 1:M[i]) {
      idx_z_w[idx_w] <- i
      for (p in 1:P) {
        idx_z_p[idx_p] <- i

        for (k in 1:K[sumP[sumM[i] + m] + p]) {

          idx_z_k[idx_k] <- idx_z
          idx_w_k[idx_k] <- idx_w
          idx_p_k[idx_k] <- p

          idx_k <- idx_k + 1
        }
        idx_z_p[idx_p] <- idx_z
        idx_w_p[idx_p] <- idx_w

        idx_p <- idx_p + 1
      }
      idx_w <- idx_w + 1
    }
    idx_z <- idx_z + 1
  }

  list("idx_z_w" = idx_z_w,
       "idx_z_p" = idx_z_p,
       "idx_w_p" = idx_w_p,
       "idx_p_k" = idx_p_k,
       "idx_w_k" = idx_w_k,
       "idx_z_k" = idx_z_k)

}

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
                       listParams = list(),
                       threshold = 1,
                       occCovariates = c(),
                       collCovariates = c(),
                       spatCovariates = c(),
                       traitsMatrix = NULL,
                       MCMCparams = list(nchain = 2,
                                         nburn = 5000,
                                         niter = 5000,
                                         nthin = 1),
                       summarisedLatentPresences = T,
                       listPriors = list()){

  {
    # data <- data
    # listParams = list(n_factors = 2)
    # threshold = 1
    # occCovariates = c("X_psi.EnvCov.1","X_psi.EnvCov.2")
    # collCovariates = c("X_theta.1","X_theta.2")
    # spatCovariates = NULL#c("latitude","longitude")
    # traitsMatrix = NULL
    # MCMCparams = list(nchain = 2,
    #                   nburn = 5000,
    #                   niter = 5000)
    # summarisedLatentPresences <- T
    # listPriors = list()
  }

  data_info <- as.data.frame(data$info)
  OTU <- data$OTU

  # data checks
  {
    if(is.null(occCovariates)) occCovariates <- c()
    if(is.null(collCovariates)) collCovariates <- c()
    if(is.null(spatCovariates)) spatCovariates <- c()

    if(!all(c(occCovariates, collCovariates, spatCovariates) %in% colnames(data$info))){
      stop("Covariate names provided not in data$info")
    }

    if(any(is.na(data_info$Site)) |
       any(is.na(data_info$Sample)) |
       any(is.na(data_info$Primer))){
      stop("NA in Site, Sample or Primer columns")
    }

  }

  # clean the data
  {
    data_info <- data_info %>%
      dplyr::arrange(Site, Sample, Primer)

    # samples per site
    {
      M_df <- data_info %>%
        dplyr::group_by(Site, Sample) %>%
        slice(1) %>%
        dplyr::group_by(Site) %>%
        dplyr::summarise(M = n(),
                         .groups = "keep")

      M <- M_df$M
      names(M) <- M_df$Site
      siteNames <- M_df$Site

      n <- length(M)
      N <- sum(M)

      sumM <- c(0, cumsum(M)[-n])

      # siteNames <- unique(data_info$Site)

      # data_info_sample <- as.numeric(factor(data_info$Sample, levels = unique(data_info$Sample)))

      }

    # marker per samples
    {
      # number of markers

      P_df <- data_info %>%
        dplyr::group_by(Site, Sample, Primer) %>%
        dplyr::slice(1) %>%
        dplyr::group_by(Site, Sample) %>%
        dplyr::summarise(P = n(),
                         .groups = "keep")

      P <- P_df$P # this in theory allows for different primer per sample, but in practice later we don't
      names(P) <- P_df$Sample
      sumP <- c(0, cumsum(rep(P, N))[-N])
      maxP <- P[1]

      primerNames <- unique(data_info$Primer)

      N2 <- maxP * N

    }

    # pcr per marker
    {
      data_K <- data_info %>%
        dplyr::group_by(Site, Sample, Primer) %>%
        dplyr::summarise(K = n(),
                         dplyr::across(contains("Species"),function(x){sum(x > 0)}),
                         .groups = "keep"
        ) %>%
        dplyr::ungroup()

      K <- data_K$K
      sumK <- c(0, cumsum(K)[-N2])

      N3 <- sum(K)

    }

    list_idx <- createDataIdx(n, M, maxP, K)
    idx_z_w <- list_idx$idx_z_w
    idx_z_k <- list_idx$idx_z_k
    idx_w_p <- list_idx$idx_w_p
    idx_z_p <- list_idx$idx_z_p
    idx_w_k <- list_idx$idx_w_k
    idx_p_k <- list_idx$idx_p_k

    # read OTU data
    {
      y <- OTU
      S <- ncol(OTU)

      # truncate data
      if(threshold != 0){

        y[OTU >= threshold] <- 1
        y[OTU < threshold] <- 0

      }

      # check for nas
      {
        y_NA <- is.na(y)
        mode(y_NA) <- "integer"
        y_NA[is.na(y_NA)] <- -1
      }

      speciesNames <- colnames(data$OTU)
      if(is.null(speciesNames)){
        speciesNames <- 1:S
      }

    }

    if(nrow(y) != N3) stop("Number of rows in data$OTU different from what obtained from data$info")

  }

  # create delta
  {
    {
      # M_marker_df <- data_info %>%
      #   dplyr::group_by(Site, Sample) %>%
      #   dplyr::summarise(M = n(),
      #                    .groups = "keep")
      #
      # M_marker <- M_marker_df$M
      # # names(M) <- M_df$Site
      #
      # sumM_marker <- c(0, cumsum(M_marker)[-N])

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

    # For occupancy covariates (group by Site, includes intercept)
    {
      X_psi <- process_covariates(data_info, occCovariates, "Site", n,
                                  remove_intercept = TRUE)

      list_Xpsi_standardised <- standardiseCovMatrix(X_psi)
      X_psi <- list_Xpsi_standardised$X
    }

    # For collection covariates (group by Sample, includes intercept)
    {
      X_theta <- process_covariates(data_info, collCovariates, "Sample", N,
                                    remove_intercept = FALSE)

    }

    # For the spatial field
    {
      Xs <- process_covariates(data_info, spatCovariates, "Site", n,
                               remove_intercept = TRUE)
      list_Xs_standardised <- standardiseCovMatrix(Xs)
      Xs <- list_Xs_standardised$X
    }

    # For the traits
    {
      if(!is.null(traitsMatrix)){
        speciesNamesInTraitsMatrix <- rownames(traitsMatrix)
        if(!all(speciesNames %in% speciesNamesInTraitsMatrix)){
          stop("Species names in OTU not present in traits matrix")
        }

        idx_speciesNames <- match(speciesNames, speciesNamesInTraitsMatrix)
        Tr <- traitsMatrix
        Tr <- Tr[idx_speciesNames,]
        Tr <- as.matrix(Tr)
        traitsNames <- colnames(Tr)
      } else {
        Tr <- matrix(NA, S, 0)
      }

    }

    if(F){
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
    }

    ncov_psi <- ncol(X_psi)
    ncov_theta <- ncol(X_theta)
    g <- ncol(Tr)

  }

  # set unknonwn parameters
  {
    d <- get_param(listParams, "n_factors")

    if(d > ncol(OTU)){
      print("More species than factors. The number of factors will be capped to the
            number of species")
      d <- ncol(OTU)
    }

    gt <- get_param(listParams, "n_lattrait", 2)
    if(ncol(Xs) > 0){
      ns <- nrow(unique(Xs))
      ps <- get_param(listParams, "n_supportpoints", getDefaultSupportPoints(ns))

    } else {
      ps <- 0
    }


  }

  # priors
  {
    prior_beta_psi <- ifelse(is.null(listPriors$prior_beta_psi), 0,  listPriors$prior_beta_psi)
    prior_beta_psi_sd <- ifelse(is.null(listPriors$prior_beta_psi_sd), 1, listPriors$prior_beta_psi_sd)
    prior_beta_theta <- 0
    prior_beta_theta_sd <- 1
    a_theta0 <- ifelse(is.null(listPriors$a_theta0), 1, listPriors$a_theta0)
    b_theta0 <- ifelse(is.null(listPriors$b_theta0), 30, listPriors$b_theta0)
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

  # run MCMC

  print("Running MCMC")

  # chain parameters
  {
    nchain <- MCMCparams$nchain
    nburn <- MCMCparams$nburn
    niter <- MCMCparams$niter
    nthin <- ifelse(is.null(MCMCparams$nthin),1,MCMCparams$nthin)
  }

  # chain output
  {
    beta_theta_output <- array(NA, dim = c(ncov_theta, S, niter, nchain))
    p_output <- array(NA, dim = c(maxP, S, niter, nchain))
    q_output <- array(NA, dim = c(maxP, S, niter, nchain))
    theta0_output <- array(NA, dim = c(S, niter, nchain))
    mu1_output <- array(NA, dim = c(niter, nchain))
    sigma1_output <- array(NA, dim = c(niter, nchain))
    mu0_output <- array(NA, dim = c(niter, nchain))
    sigma0_output <- array(NA, dim = c(niter, nchain))

    # jsdm params

    {
      B0_output <- array(NA, dim = c(S, niter, nchain))
      B_output <- array(NA, dim = c(ncov_psi, S, niter, nchain))
      G_output <- array(NA, dim = c(g, ncov_psi, niter, nchain))
      A_output <- array(NA, dim = c(S, gt, niter, nchain))
      C_output <- array(NA, dim = c(gt, ncov_psi, niter, nchain))
      Bs_output <- array(NA, dim = c(ps, S, niter, nchain))
      Gs_output <- array(NA, dim = c(g, ps, niter, nchain))
      As_output <- array(NA, dim = c(S, gt, niter, nchain))
      Cs_output <- array(NA, dim = c(gt, ps, niter, nchain))
      U_output  <- array(NA, dim = c(n, d, niter, nchain))
      L_output  <- array(NA, dim = c(d, S, niter, nchain))
      tau_output <- array(NA, dim = c(S, niter, nchain))
      sigmab_output <- array(NA, dim = c(niter, nchain))
      sigmabs_output <- array(NA, dim = c(niter, nchain))
      varPart_output <- array(NA, dim = c(S, 4, niter, nchain))
    }

    if(summarisedLatentPresences){
      w_output_mean <- matrix(0, N, S)
      z_output_mean <- matrix(0, n, S)
    } else {
      w_output <- array(NA, dim = c(N, S, niter, nchain))
      z_output <- array(NA, dim = c(n, S, niter, nchain))
    }

  }

  # precompute spatial quantities
  {
    ps <- ifelse(ncol(Xs) > 0, getDefaultSupportPoints(n), 0)

    # Spatial covariates matrix
    list_Xs <- computeSpatialSummaries(Xs, ps, maxPoints = 5)
    Xs_centers <- list_Xs$Xs_centers
    Xs_index <- list_Xs$Xs_index
    X_s_centers <- list_Xs$X_s_centers
    X_tilde <- list_Xs$X_tilde
    X_s <- list_Xs$X_s

    length_grid_ls <- 10
    l_s_grid <- seq(0.01, 0.3, length.out = length_grid_ls)
    list_SoRSummaries <- precomputeSORmatrices(l_s_grid, list_Xs)
  }

  # precompute jsdm params
  {
    list_data <- list(
      "X" = X_psi,
      "Tr" = Tr)

    a_sigmab <- .01; b_sigmab <- .01
    a_sigmabs <- .01; b_sigmabs <- .01
    a_tau <- 5; b_tau <- 5
    a_l_s <- 1; b_l_s <- 1

    list_priors <- list(
      "a_sigmab" = a_sigmab,
      "b_sigmab" = b_sigmab,
      "a_sigmabs" = a_sigmabs,
      "b_sigmabs" = b_sigmabs,
      "a_tau" = a_tau,
      "b_tau" = b_tau,
      "a_l_s" = a_l_s,
      "b_l_s" = b_l_s
    )
  }

  for (chain in 1:nchain) {

    # chain output
    {
      beta_theta_output_chain <- array(NA, dim = c(ncov_theta, S, niter))
      p_output_chain <- array(NA, dim = c(maxP, S, niter))
      q_output_chain <- array(NA, dim = c(maxP, S, niter))
      theta0_output_chain <- array(NA, dim = c(S, niter))
      mu1_output_chain <- rep(NA, niter)
      sigma1_output_chain <- rep(NA, niter)
      mu0_output_chain <- rep(NA, niter)
      sigma0_output_chain <- rep(NA, niter)

      # jsdm params
      {
        B0_output_chain <- array(NA, dim = c(S, niter))
        B_output_chain <- array(NA, dim = c(ncov_psi, S, niter))
        G_output_chain <- array(NA, dim = c(g, ncov_psi, niter))
        A_output_chain <- array(NA, dim = c(S, gt, niter))
        C_output_chain <- array(NA, dim = c(gt, ncov_psi, niter))
        Bs_output_chain <- array(NA, dim = c(ps, S, niter))
        Gs_output_chain <- array(NA, dim = c(g, ps, niter))
        As_output_chain <- array(NA, dim = c(S, gt, niter))
        Cs_output_chain <- array(NA, dim = c(gt, ps, niter))
        U_output_chain <- array(NA, dim = c(n, d, niter))
        L_output_chain <- array(NA, dim = c(d, S, niter))
        sigmab_output_chain <- rep(NA, niter)
        sigmabs_output_chain <- rep(NA, niter)
        tau_output_chain <- matrix(NA, S, niter)
        varPart_output_chain <-  array(NA, dim = c(S, 4, niter))
      }

      if(!summarisedLatentPresences){
        z_output_chain <- array(NA, dim = c(n, S, niter))
        w_output_chain <- array(NA, dim = c(N, S, niter))
      }

    }

    # starting values
    {
      c_imk <- y > 0

      w <- matrix(NA, N, S)

      for (s in 1:S) {
        for (i in 1:n) {
          for (m in 1:M[i]) {
            idx_im1 <- sumK[sumP[sumM[i] + m] + 1] + 1
            idx_im2 <- sumK[sumP[sumM[i] + m] + P[sumM[i] + m]] +
              K[sumP[sumM[i] + m] + P[sumM[i] + m]]
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

      beta_theta <- matrix(0, ncov_theta, S)
      theta <- computeTheta(X_theta, beta_theta)

      p <- matrix(.9, maxP, S)
      q <- matrix(.05, maxP, S)
      theta0 <- rep(.05, S)

      # jsdmParams
      {
        B0 <- rep(0, S)
        B <- matrix(0, ncov_psi, S)
        Bs <- matrix(0, ps, S)
        L <- matrix(1, d, S)
        diag(L) <- 1
        L[lower.tri(L)] <- 0
        G <- matrix(0, g, ncov_psi)
        C <- matrix(1, gt, ncov_psi)
        diag(C) <- 1
        C[lower.tri(C)] <- 0
        A <- matrix(0, S, gt)
        Gs <- matrix(0, g, ps)
        Cs <- matrix(1, gt, ps)
        diag(Cs) <- 1
        Cs[lower.tri(Cs)] <- 0
        As <- matrix(0, S, gt)
        U <- matrix(0, n, d)
        sigma_b <- 1
        sigma_bs <- .001
        idx_ls <- 3 # dim(list_SoRSummaries$Ks_all)[3][5]
        tau <- rep(1, S)

        Bt <- t(B) - computeBtcoef(G, Tr, A, C, matrix(0, S, ncov_psi))
        Bst <- t(Bs) - computeBtcoef(Gs, Tr, As, Cs, matrix(0, S, ps))

        Ks <- list_SoRSummaries$Ks_all[,,idx_ls]

        list_jSDMparams <- list(
          "B0" = B0,
          "B" = B,
          "G" = G,
          "A" = A,
          "C" = C,
          "Bt" = Bt,
          "Bs" = Bs,
          "Gs" = Gs,
          "As" = As,
          "Cs" = Cs,
          "Bst" = Bst,
          "U" = U,
          "L" = L,
          "sigma_b" = sigma_b,
          "sigma_bs" = sigma_bs,
          "idx_ls" = idx_ls,
          "tau" = tau
        )

        list_psiCoef <- computePsiCoef(
          X_psi, Ks, list_Xs$Xs_centers, Tr,
          B0, G, A, C, Bt,
          Gs, As, Cs, Bst,
          U, L)

      }

      psi <- logistic(list_psiCoef$eta)

    }

    for (iter in 1:(nburn + niter * nthin)) {

      if(iter > nburn){
        currentIter <- (iter - nburn) / nthin
        if(currentIter %% 100 == 0){
          print(paste0("Chain ", chain, " - Iteration ",currentIter))
        }
      } else {
        if(iter %% 100 == 0){
          print(paste0("Chain ", chain, " - Burn in Iteration ",iter))
        }
      }

      # sample z

      z <- sample_z_cpp(w, psi, theta, theta0, M, sumM)

      # sample psi
      {
        list_data$z <- z

        list_jSDMparams <- update_jSDMcoef(
          list_data,
          list_jSDMparams,
          list_priors,
          list_Xs,
          list_SoRSummaries,
          model = "binary"
        )
        psi <- list_jSDMparams$psi
      }

      # sample psi - old

      if (F) {
        list_betapsiLL <- sample_psivars(z, X_psi, beta_psi, X_ord, beta_ord,
                                         E, Es, LL, prior_beta_psi, prior_beta_psi_sd)
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
                               M, K, sumP, sumM, sumK, maxP, z)

      } else {

        w <- sample_w_cpp(logy1, mu0, sigma0, mu1, sigma1, theta, theta0, p, q,
                          M, K, sumP, sumM, sumK, maxP, z)

      }

      # sample cimk

      if (threshold > 0) {

        w_all <- w[idx_w_k,,drop=F]

        # faster way to assign c_imk to 1 if logy1 > 0 for w_all = 1 and to 2
        # if log1 > 0 when w_all = 0
        c_imk <- (y > 0) * (2 - (w_all == 1))

      }

      # sample theta
      beta_theta <- sample_betatheta_cpp(w, z, beta_theta, idx_z_w, X_theta,
                                         b_betatheta, B_betatheta)
      theta <- computeTheta(X_theta, beta_theta)

      # sample pq
      list_pq <- sample_pq_cpp(c_imk, w, idx_p_k, idx_w_k, maxP, a_p, b_p, a_q, b_q)
      p <- list_pq$p
      q <- list_pq$q

      # sample theta0
      theta0 <- sample_theta0(z, w, idx_z_w, a_theta0, b_theta0)

      {

        if(iter > nburn & (iter - nburn) %% nthin == 0){
          currentIter <- (iter - nburn) / nthin

          beta_theta_output_chain[,,currentIter] <- beta_theta
          p_output_chain[,,currentIter] <- p
          q_output_chain[,,currentIter] <- q
          theta0_output_chain[,currentIter] <- theta0

          # save jsdm params
          {
            B0_output_chain[,currentIter] <- list_jSDMparams$B0
            B_output_chain[,,currentIter] <- list_jSDMparams$B
            G_output_chain[,,currentIter] <- list_jSDMparams$G
            A_output_chain[,,currentIter] <- list_jSDMparams$A
            C_output_chain[,,currentIter] <- list_jSDMparams$C
            Bs_output_chain[,,currentIter] <- list_jSDMparams$Bs
            Gs_output_chain[,,currentIter] <- list_jSDMparams$Gs
            As_output_chain[,,currentIter] <- list_jSDMparams$As
            Cs_output_chain[,,currentIter] <- list_jSDMparams$Cs
            U_output_chain[,,currentIter] <- list_jSDMparams$U
            L_output_chain[,,currentIter] <- list_jSDMparams$L
            sigmab_output_chain[currentIter] <- list_jSDMparams$sigma_b
            sigmabs_output_chain[currentIter] <- list_jSDMparams$sigma_bs

            varPart_output_chain[,,currentIter] <-
              as.matrix(list_jSDMparams$variancePartitioning)
          }

          if(summarisedLatentPresences){
            z_output_mean <- z_output_mean +
              (1 / (niter) * nchain) * z
            w_output_mean <- w_output_mean +
              (1 / (niter) * nchain) * w
          } else {
            z_output_chain[,,currentIter] <- z
            w_output_chain[,,currentIter] <- w
          }

          if(threshold == 0){
            mu1_output_chain[currentIter] <- mu1
            sigma1_output_chain[currentIter] <- sigma1
            mu0_output_chain[currentIter] <- mu0
            sigma0_output_chain[currentIter] <- sigma0
          }

        }

      }

    }

    beta_theta_output[,,,chain] <- beta_theta_output_chain
    p_output[,,,chain] <- p_output_chain
    q_output[,,,chain] <- q_output_chain
    theta0_output[,,chain] <- theta0_output_chain
    mu1_output[,chain] <- mu1_output_chain
    sigma1_output[,chain] <- sigma1_output_chain
    mu0_output[,chain] <- mu0_output_chain
    sigma0_output[,chain] <- sigma0_output_chain

    if(!summarisedLatentPresences){
      z_output[,,,chain] <- z_output_chain
      w_output[,,,chain] <- w_output_chain
    }

    # save jsdm params
    {
      B0_output[,,chain] <- B0_output_chain
      B_output[,,,chain] <- B_output_chain
      G_output[,,,chain] <- G_output_chain
      A_output[,,,chain] <- A_output_chain
      C_output[,,,chain] <- C_output_chain
      Bs_output[,,,chain] <- Bs_output_chain
      Gs_output[,,,chain] <- Gs_output_chain
      As_output[,,,chain] <- As_output_chain
      Cs_output[,,,chain] <- Cs_output_chain
      U_output[,,,chain] <- U_output_chain
      L_output[,,,chain] <- L_output_chain
      sigmab_output[,chain] <- sigmab_output_chain
      sigmabs_output[,chain] <- sigmabs_output_chain
      varPart_output[,,,chain] <- varPart_output_chain
      tau_output[,,chain] <- tau_output_chain
    }

  }

  # reparametrise factor coefficients
  {
    if(d > 0){

      list_UL_output_reparm <- reparamFactorModel(U_output, L_output)
      U_output <- list_UL_output_reparm$U_output
      L_output <- list_UL_output_reparm$L_output

    }

    if(gt > 0){

      list_AC_output_reparm <- reparamFactorModel(A_output, C_output)
      A_output <- list_AC_output_reparm$U_output
      C_output <- list_AC_output_reparm$L_output

    }

  }

  jsdm_results_output <- list(
    "B0_output" = B0_output,
    "B_output" = B_output,
    "G_output" = G_output,
    "A_output" =  A_output,
    "C_output" = C_output,
    "Bs_output" = Bs_output,
    "Gs_output" = Gs_output,
    "As_output" = As_output,
    "Cs_output" = Cs_output,
    "U_output" = U_output,
    "L_output" = L_output,
    "sigmab_output" = sigmab_output,
    "sigmabs_output" = sigmabs_output,
    "varPart_output" = varPart_output,
    "tau_output" = tau_output
  )

  results_output <- list(
    "jsdm_output" = jsdm_results_output,
    "beta_theta_output" = beta_theta_output,
    "p_output" = p_output,
    "q_output" = q_output,
    "theta0_output" = theta0_output,
    "mu1_output" = mu1_output,
    "sigma1_output" = sigma1_output,
    "mu0_output" = mu0_output,
    "sigma0_output" = sigma0_output
  )

  if(summarisedLatentPresences){
    results_output$z_output <- z_output_mean
    results_output$w_output <- w_output_mean
  } else {
    results_output$z_output <- z_output
    results_output$w_output <- w_output
  }

  minESS <- computeMinESS(results_output)

  if(minESS < 50) {
    print(paste0("Minimum ffective sample size equal to ", minESS,", please rerun with more iterations"))
  }

  infos <- list(
    "S" = S,
    "g" = g,
    "M" = M,
    "n" = n,
    "K" = K,
    "n_factors" = d,
    "list_idx" = list_idx,
    "data_info" = data_info,
    "OTU" = OTU,
    "speciesNames" = speciesNames,
    "primerNames" = primerNames,
    "siteNames" = siteNames,
    "ncov_theta" = ncov_theta,
    "ncov_psi" = ncov_psi,
    "Xpsi_standardised" = list_Xpsi_standardised,
    "Xs_standardised" = list_Xs_standardised

  )

  list(
    "results_output" = results_output,
    "infos" = infos,
    "Tr" = Tr,
    "X_theta" = X_theta,
    "X_s" = X_s,
    "X_psi" = X_psi)

}
