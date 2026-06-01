
runMCMCOccPlus <- function(data,
                           d,
                           threshold = F,
                           occCovariates = c(),
                           ordCovariates = c(),
                           detCovariates = c(),
                           MCMCparams = list(nchain = 2,
                                             nburn = 1000,
                                             niter = 1000)){

  data_info <- as.data.frame(data$info)
  OTU <- data$OTU

  # data checks
  {
    if(!all(c(occCovariates, ordCovariates, detCovariates) %in% colnames(data$info))){
      stop("Covariate names provided not in data$info")
    }

    if(any(is.na(data_info$Site)) | any(is.na(data_info$Sample)) | any(is.na(data_info$Primer))){
      stop("NA in Site, Sample or Primer columns")
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

    data_info_sample <- as.numeric(as.factor(data_info$Sample))

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
  }

  maxL = max(L)

  idx_z <- rep(1:n, M)
  idx_k <- as.numeric(as.factor(data_info_sample))

  primerIdx <- as.numeric(as.factor(data_info$Primer))

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

      X_psi <- X_psi[,-1,drop=F]

    } else {

      X_psi <- matrix(0, n, 0)

    }

    if(length(detCovariates) > 0){

      X_theta <- data_info %>%
        dplyr::group_by(Sample) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(detCovariates),
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
    a_sigma0 <- 10
    b_sigma0 <- 10
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

  if(threshold != F){

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
    beta_ord_output <- array(NA, dim = c(ncov_ord, d, niter, nchain))
    beta_psi_output <- array(NA, dim = c(ncov_psi, S, niter, nchain))
    beta_theta_output <- array(NA, dim = c(ncov_theta, S, niter, nchain))
    LL_output <- array(NA, dim = c(d, S, niter, nchain))
    E_output <- array(NA, dim = c(n, d, niter, nchain))
    U_output <- array(NA, dim = c(n, d, niter, nchain))
    UL_output <- array(NA, dim = c(n, S, niter, nchain))
    z_output <- array(NA, dim = c(n, S, niter, nchain))
    p_output <- array(NA, dim = c(maxL, S, niter, nchain))
    q_output <- array(NA, dim = c(maxL, S, niter, nchain))
    theta0_output <- array(NA, dim = c(S, niter, nchain))
  }

  for (chain in 1:nchain) {

    # chain output
    {
      beta_ord_output_chain <- array(NA, dim = c(ncov_ord, d, niter))
      beta_psi_output_chain <- array(NA, dim = c(ncov_psi, S, niter))
      beta_theta_output_chain <- array(NA, dim = c(ncov_theta, S, niter))
      LL_output_chain <- array(NA, dim = c(d, S, niter))
      E_output_chain <- array(NA, dim = c(n, d, niter))
      U_output_chain <- array(NA, dim = c(n, d, niter))
      UL_output_chain <- array(NA, dim = c(n, S, niter))
      z_output_chain <- array(NA, dim = c(n, S, niter))
      p_output_chain <- array(NA, dim = c(maxL, S, niter))
      q_output_chain <- array(NA, dim = c(maxL, S, niter))
      theta0_output_chain <- array(NA, dim = c(S, niter))
    }

    # starting values
    {
      if(!threshold){


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

          beta_psi <- matrix(0, ncov_psi, S)
          beta_theta <- matrix(0, ncov_theta, S)
          beta_ord <- matrix(0, ncov_ord, d)
          E <- matrix(0, n, d)
          LL <- matrix(1, d, S)

          U <- computeU(X_ord, beta_ord, E)
          psi <- computePsi(X_psi, beta_psi, U, LL)
          theta <- computeTheta(X_theta, beta_theta)

          p <- matrix(.9, maxL, S)
          q <- matrix(.05, maxL, S)
          theta0 <- rep(.05, S)

          mu1 <- 7
          sigma1 <- 3
          mu0 <- 1
          sigma0 <- 1
        }

        U <- computeU(X_ord, beta_ord, E)
        psi <- computePsi(X_psi, beta_psi, U, LL)
      } else {


      }
    }

    for (iter in 1:(niter + nburn)) {

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
      w <- sample_w_cpp(logy1, mu0, sigma0, mu1, sigma1, theta, theta0, p, q,
                        M, K, sumL, sumM, sumK, maxL, z)

      # sample cimk

      {
        w_all <- w[idx_k,,drop=F]

        c_imk <- ifelse(w_all == 1 & logy1 > 0, 1,
                        ifelse(w_all == 0 & logy1 > 0, 2, 0))

        # c_imk <- sample_cimk(logy1, mu1, sigma1, pi0, sigma0,
        # p, q, idx_k, primerIdx)

        # c_imk <- sample_cimk_cpp(y, logy1, w, mu1, sigma1, pi0, sigma0,
        #                          p, q, idx_k, primerIdx, maxL)
      }

      # sample theta
      beta_theta <- sample_betatheta(w, beta_theta, idx_z, X_theta,
                                     b_betatheta, B_betatheta)
      theta <- computeTheta(X_theta, beta_theta)

      # sample pq
      # list_pq <- sample_pq(y, w, primerIdx, idx_k, maxL, a_p, b_p, a_q, b_q)
      list_pq <- sample_pq_cpp(c_imk, w, primerIdx, idx_k, maxL, a_p, b_p, a_q, b_q)
      (p <- list_pq$p)
      (q <- list_pq$q)

      # sample theta0
      theta0 <- sample_theta0(z, w, idx_z, a_theta0, b_theta0)

      # sample mu1sigma1
      list_mu1sigma1 <- sample_musigma(sigma1,
                                       logy1, c_imk,
                                       mu_mu1, sd_mu1,
                                       a_sigma1, b_sigma1, tpfp = T)
      mu1 <- list_mu1sigma1$mu1
      sigma1 <- list_mu1sigma1$sigma1

      # sample mu0sigma0
      list_mu1sigma1 <- sample_musigma(sigma0,
                                       logy1, c_imk,
                                       mu_mu0, sd_mu0,
                                       a_sigma0, b_sigma0, tpfp = F)
      mu0 <- list_mu1sigma1$mu1
      sigma0 <- list_mu1sigma1$sigma1

      {

        if(iter > nburn){
          currentIter <- iter - nburn
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

    }

  }

  results_output <- list(
    "beta_ord_output" = beta_ord_output,
    "beta_psi_output" = beta_psi_output,
    "beta_theta_output" = beta_theta_output,
    "LL_output" = LL_output,
    "E_output" = E_output,
    "UL_output" = UL_output,
    "z_output" = z_output,
    "p_output" = p_output,
    "q_output" = q_output,
    "theta0_output" = theta0_output)

  infos <- list(
    "S" = S,
    "speciesNames" = speciesNames,
    "primerNames" = primerNames,
    "siteNames" = siteNames,
    "d" = d,
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
