
logistic <- function(x){
  1 / (1 + exp(-x))
}

reparamFactorModel <- function(U_output, L_output){

  L_output_reparam <- L_output
  U_output_reparam <- U_output

  d <- dim(L_output)[1]

  for (chain in 1:nchain) {

    for (iter in 1:niter) {

      if(d == 1){

        L1 <- L_output[1,1,iter,chain]
        L_output_reparam[1,,iter,chain] <- L_output[1,,iter,chain] / L1
        U_output_reparam[,1,iter,chain] <- U_output[,1,iter,chain] * L1

      } else {

        L_current <- L_output[,,iter,chain]
        U_current <- U_output[,,iter,chain]

        qr_decomp <- qr(L_current)
        Q_current <- qr.Q(qr_decomp)
        R_current <- qr.R(qr_decomp)

        Q2 <- Q_current %*% diag(diag(R_current), nrow = d)
        invQ2 <- diag(1 / diag(R_current), nrow = d) %*% t(Q_current)

        L_new <- invQ2 %*% L_current
        U_new <- U_current %*% Q2

        L_output_reparam[,,iter,chain] <- L_new
        U_output_reparam[,,iter,chain] <- U_new
      }

    }
  }

  list("L_output" = L_output_reparam,
       "U_output" = U_output_reparam)
}

buildGrid <- function(XY_sp, gridStep){

  x_grid <- seq(min(XY_sp[,1]) - (1.5) * gridStep,
                max(XY_sp[,1]) + (1.5) * gridStep, by = gridStep)
  y_grid <- seq(min(XY_sp[,2]) - (1.5) * gridStep,
                max(XY_sp[,2]) + (1.5) * gridStep, by = gridStep)

  pointInGrid <- matrix(T, nrow = length(x_grid), ncol = length(y_grid))

  for (i in 2:(length(x_grid) - 1)) {

    for (j in 2:(length(y_grid) - 1)) {

      isAnyPointInBandRight <- isPointInBandRight(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandLeft <- isPointInBandLeft(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandUp <- isPointInBandUp(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandDown <- isPointInBandDown(XY_sp, x_grid, y_grid, i - 1, j - 1)

      if(!isAnyPointInBandRight | !isAnyPointInBandLeft | !isAnyPointInBandUp | !isAnyPointInBandDown){
        pointInGrid[i,j] <- F
      }

    }

  }

  pointInGrid <- pointInGrid[-c(1,nrow(pointInGrid)),]
  pointInGrid <- pointInGrid[,-c(1,ncol(pointInGrid))]
  x_grid <- x_grid[-c(1,length(x_grid))]
  y_grid <- y_grid[-c(1,length(y_grid))]

  allPoints <- cbind(expand.grid(x_grid, y_grid), as.vector((pointInGrid)))
  allPoints <- allPoints[allPoints[,3],-3]

  allPoints
}

computeSpatialSummaries <- function(Xs, ps, maxPoints){

  if(ps > 0){

    # isolate unique locations and assign indexes to sites
    uniqueXs <- which(!duplicated(Xs))
    X_s <- Xs[uniqueXs,]

    # indexes assigning original locations (Xs) to new locations (X_s)
    Xs_index <- match(
      do.call(paste, as.data.frame(Xs)),
      do.call(paste, as.data.frame(X_s))
    )

    # location of support points
    # X_tilde <- as.matrix(buildGrid(X_s, gridStep = .4))
    list_kmeans <- kmeans(X_s, centers = ps)
    X_tilde <- list_kmeans$centers

    {
      # ggplot() +
      #   geom_point(data = NULL, aes(x = X_tilde[,1],
      #                               y = X_tilde[,2], color = "red")) +
      #   geom_point(data = NULL, aes(x = X_s[,1],
      #                               y = X_s[,2]))

    }

    # distance from support points
    X_s_Xtilde_dist <- t(apply(X_s, 1, function(x){
      apply(X_tilde, 1, function(y){
        (x[1] - y[1])^2 + (x[2] - y[2])^2
      })
    }))

    # indexes of closest support points to the unique points
    X_s_centers <- t(apply(X_s_Xtilde_dist, 1, function(x){
      order(x)[1:maxPoints]
    }))

    # closest support points to the original locations
    Xs_centers <- X_s_centers[Xs_index,]

  } else {

    Xs_index <- rep(NA, n)
    X_tilde <- matrix(NA, ps, 0)
    X_s_centers <- matrix(NA, n, 0)
    Xs_centers <- matrix(NA, n, 0)
    X_s <- matrix(NA, nrow(Xs), 2)

  }

  list(
    "Xs_index" = Xs_index, # indexes matching original locations to new locations
    "X_tilde" = X_tilde, # location of support points
    "X_s_centers" = X_s_centers, # indexes to match new locations to the coefficients
    "Xs_centers" = Xs_centers, # indexes to match original locations to the coefficients
    "X_s" = X_s) # new unique locations

}

computeSORmatrix <- function(l_s, X_tilde, X_s, Xs_index, X_s_centers){

  ps <- nrow(X_tilde)

  if(ps > 0){

    K_uu <- K2(X_tilde, X_tilde, 1, l_s) + diag(exp(-10), nrow = nrow(X_tilde))
    L_Kmm <- FastGP::rcppeigen_get_chol(K_uu)
    invL_Kmm <- FastGP::rcppeigen_invert_matrix(L_Kmm)
    K_staru <- K2(X_s, X_tilde, 1, l_s)
    KnmLmt <- K_staru %*% t(invL_Kmm)
    Ks <- t(sapply(1:nrow(KnmLmt), function(i){ KnmLmt[i,X_s_centers[i,]]}))
    Ks <- Ks[Xs_index,]

    logDetKuu <- sum(log(FastGP::rcppeigen_get_diag(K_uu))) * 2
    sq_term <- FastGP::rcppeigen_get_chol(FastGP::rcppeigen_invert_matrix(K_uu))

  } else {

    Ks <- matrix(NA, nrow = length(Xs_index), 0)
    logDetKuu <- NULL
    sq_term <- NULL

  }

  list("Ks" = Ks,
       "logDetKuu" = logDetKuu,
       "sq_term" = sq_term)

}

precomputeSORmatrices <- function(l_s_grid, list_Xs){


  length_grid_ls <- length(l_s_grid)

  X_tilde <- list_Xs$X_tilde
  X_s <- list_Xs$X_s
  Xs_index <- list_Xs$Xs_index
  X_s_centers <- list_Xs$X_s_centers

  X_centers <- nrow(X_tilde)
  maxPoints <- ncol(X_s_centers)
  n <- length(Xs_index)

  Ks_all <- array(NA, dim = c(n, maxPoints, length_grid_ls))
  logDetKuu_grid <- rep(NA, length_grid_ls)
  Lm1_grid <- array(NA, c(X_centers, X_centers, length_grid_ls))

  if(X_centers > 0){
    for (j in 1:length_grid_ls) {

      if(T){
        print(paste0("Precomputing covariance matrix ",j," out of ",length_grid_ls))
      }

      l_s_current <- l_s_grid[j]

      list_SoRelem <- computeSORmatrix(l_s_current, X_tilde, X_s, Xs_index, X_s_centers)

      Ks_all[,,j] <- list_SoRelem$Ks
      logDetKuu_grid[j] <- list_SoRelem$logDetKuu
      Lm1_grid[,,j] <- list_SoRelem$sq_term

    }
  }

  list("Ks_all" = Ks_all,
       "logDetKuu_grid" = logDetKuu_grid,
       "Lm1_grid" = Lm1_grid,
       "l_s_grid" = l_s_grid)

}

createSplinesObjects <- function(X, df){

  # list_ns <- list()

  list_ns <- lapply(seq_len(ncol(X)), function(j) {
    ns_j <- ns(X[, j], df = 5)
    ns_j
    # list_ns[[j]] <-
    # Zj <- ns_j
    # colnames(Zj) <- paste0(colnames(X)[j], " - s", seq_len(ncol(Zj)))
    # Zj
  })

  # Z <- do.call(cbind, Z_list)

  # Z

  list_ns

}

createSplinesMatrixSingleCov <- function(n_s, X_cov){

  Zj <- predict(n_s, X_cov)
  colnames(Zj) <- paste0(colnames(X_cov), " - s", seq_len(ncol(Zj)))

  Zj
}

createSplinesMatrix <- function(list_ns, X_new){

  Z_list <- lapply(seq_len(ncol(X_new)), function(j) {
    # Zj <- predict(list_ns[[j]], X_new[,j])
    # colnames(Zj) <- paste0(colnames(X)[j], " - s", seq_len(ncol(Zj)))
    Zj <- createSplinesMatrixSingleCov(list_ns[[j]], X_new[,j,drop=F])
    Zj
  })

  Z <- do.call(cbind, Z_list)

  Z

}

# DATA SIMULATION -------

sampleEffects <- function(n){

  sample((-1):1, n, replace = T)

}

simulateData <- function(
    n, S, p, g, gt, d, tau, ds,
    sigma_b, sigma_bs, sigma_ts, sigma_h, l_s,
    useSpatField, usingSplines, model){

  # simulate data
  {
    speciesNames <- as.character(1:S)

    # Fixed effects matrix
    X <- matrix(rnorm(n * p), n, p)
    if(p > 0) colnames(X) <- paste("EnvCov", seq_len(p))
    if(usingSplines){
      X0 <- X
      list_ns <- createSplinesObjects(X, df= 3)
      X <- createSplinesMatrix(list_ns, X0)
      p0 <- p
      p <- ncol(X)
    }

    # Traits matrix
    Tr <- matrix(rnorm(S * g), S, g)

    # Spatial locations
    ns <- n
    Xs <- matrix(runif(ns * 2), ns, 2)
    if(ns < n){
      Xs <- Xs[sample(1:ns, n, replace = T),]
    }

    ps <- 0
  }

  # params
  {

    # Traits covariate coefficients
    G <- matrix(sampleEffects(p * g), g, p)

    # Unobserved traits
    A <- matrix(rnorm(gt * S), S, gt)

    # Unobserved response to traits
    C <- matrix(sampleEffects(p * gt), gt, p)
    diag(C) <- 1
    C[lower.tri(C)] <- 0

    # Residual environmental covariates variation
    Bt <- matrix(rnorm(p * S, sd = sigma_b), S, p)

    # Environmental covariates coefficient
    B <- t(computeBtcoef(G, Tr, A, C, Bt))

    # Spatial response to traits
    Gs <- matrix(0, g, ps)

    # Unobserved spatial traits
    As <- matrix(0, S, gt)

    # Unobserved response to spatial traits
    Cs <- matrix(0, gt, ps)
    diag(Cs) <- 1
    Cs[lower.tri(Cs)] <- 0

    # Residual spatial variation
    Bst <- matrix(0, S, ps)

    # Spatial covariates coefficient
    Bs <- t(computeBtcoef(Gs, Tr, As, Cs, Bst))

    #  factor scores
    U <- matrix(rnorm(n * d, sd = sigma_h), n, d)

    # factor loadings
    L <- matrix(sampleEffects(d * S), d, S)
    diag(L) <- 1
    L[lower.tri(L)] <- 0

  }

  # generate spatial field
  {
    if(useSpatField){

      K_mat <- K2(Xs, Xs, 1, l_s) + diag(0.001, nrow = n)

      Lambda <- matrix(rnorm(ds * S), ds, S)
      Sigma <- t(Lambda) %*% Lambda + diag(0.001, nrow = S)

      LU <- chol(K_mat)
      LV <- chol(Sigma)

      Z <- matrix(rnorm(n * S), n, S)

      spatField <- t(LU) %*% Z %*% LV

    }

    Xs_centers <- matrix(NA, 0, 0)
    Ks <- matrix(NA, 0, 0)

  }

  # simulate observations
  list_eta <- computePsiCoef(
    X, Ks, Xs_centers, Tr,
    G, A, C, Bt,
    Gs, As, Cs, Bst,
    U, L)
  eta <- list_eta$eta
  XB <- list_eta$XB
  SE <- list_eta$SE
  UL <- list_eta$UL

  if(ps > 0) eta <- eta + spatField

  varPart <- computeVariancePartitioning(XB, SE + spatField, UL)

  # outcomes
  if(model == "continuous"){

    z <- t(sapply(1:n, function(i){
      sapply(1:S, function(j){
        rnorm(1, eta[i,j], tau[j])
      })
    }))

  } else if(model == "binary"){

    psi <- logistic(eta)

    z <- t(sapply(1:n, function(i){
      sapply(1:S, function(j){
        rbinom(1, 1, psi[i,j])
      })
    }))

  } else if (model == "count"){

    mu <- exp(eta)

    z <- t(sapply(1:n, function(i){
      sapply(1:S, function(j){
        rpois(1, lambda = mu[i,j])
      })
    }))

  }

  # save true values
  {
    data <- list(
    "z" = z,
    "X" = X,
    "Tr" = Tr,
    "Xs" = Xs
    )

    trueParams <- list(
      "B" = B,
      "Bt" = Bt,
      "G" = G,
      "A" = A,
      "C" = C,
      "SE" = SE,
      "U" = U,
      "L" = L,
      "sigma_b" = sigma_b,
      "sigma_bs" = sigma_bs,
      "idx_ls" = idx_ls,
      "eta" = eta,
      "tau" = tau,
      "varPart" = varPart
    )
  }

  list("data" = data,
       "trueParams" = trueParams)
}

# MCMC FUN ---------------------

computeBtcoef <- function(G, Tr, A, C, Btilde){

  Tr %*% G + A %*% C + Btilde

}

computePsiCoef <- function(
    X, Ks, Xs_centers, Tr,
    G, A, C, Bt,
    Gs, As, Cs, Bst,
    H, L){

  ps <- ncol(Bst)

  B <- t(computeBtcoef(G, Tr, A, C, Bt))

  XB_prod <- X %*% B

  if(ps > 0){
    # Bs <- computeBscoef(Tr, Gs, As, Cs)
    Bs <- t(computeBtcoef(Gs, Tr, As, Cs, Bst))
    XsBs_prod <- KsBproduct(Ks, Bs, Xs_centers)
  } else {
    XsBs_prod <- matrix(0, nrow(XB_prod), ncol(XB_prod))
  }

  HL <- H %*% L

  eta <- XB_prod + XsBs_prod + HL

  list("eta" = eta,
       "XB" = XB_prod,
       "SE" = XsBs_prod,
       "UL" = HL)

}

# sample residual variance of env. covariates coefficients
sample_sigmab <- function(B, Tr, G, A, C, a_sigmab, b_sigmab){

  p <- nrow(B)
  S <- ncol(B)

  Bt <- t(B) - computeBtcoef(G, Tr, A, C, matrix(0, S, p))

  sumsq <- sum(Bt^2)
  n_samples <- p * S

  sqrt(
    rinvgamma_cpp(a_sigmab + (n_samples / 2), b_sigmab + (sumsq / 2))
  )

}

# sample variances of responses
sample_tau <- function(z, psiCoef, a_tau, b_tau){

  n <- nrow(z)
  S <- ncol(z)

  sumsqs <- colSums((z - psiCoef)^2)

  tau <- sapply(1:S, function(s){

    sqrt(
      rinvgamma_cpp(a_sigmab + (n / 2), b_sigmab + (sumsqs[s] / 2))
    )
  })

  tau
}

# sample the fixed effects and the factor loadings
sample_BCsL <- function(k, X, H, G, Tr,
                      A, C, sigma_b,
                      Ks, Xs_centers, Gs, As,
                      Omega) {

  p <- ncol(X)
  d <- ncol(U)
  S <- ncol(Omega)
  ps <- nrow(As)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))

  B <- matrix(0, p, S)
  Cs <- matrix(0, gt, S)
  L <- matrix(0, d, S)

  # no spatial factor loadings if there is no spatial effect
  ps <- ifelse(ps > 0, ps, 0)

  if(p + gts + d > 0){

    for (s in 1:S) {

      k_current <- k[,s]

      XU <- cbind(X, U)

      b_current <- c(M_B[,s], rep(0, ps), rep(0, d))
      B_current <- diag(1, nrow = p + gts + d)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      # if(ps == 0){
      #   BL <- sampleB(XU, B_current, b_current, Omega[,s], k_current)
      # } else {
        invB <- diag(1 / diag(B_current), nrow = p + ps + d)

        BCsL <- sampleB_SoR(XU, B_current, b_current, k_current,
                            Omega[,s], Xs_centers, Ks, ps)
        if(gts > 0){
          Cs[seq_len(gt),s] <- BCsL[p + d + seq_len(gts)]
        }

      # }
      B[seq_len(p),s] <- BCsL[seq_len(p)]
      L[seq_len(d),s] <- BCsL[p + seq_len(d)]


    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  Bs <- computeBscoef(Tr, Gs, As, Cs)

  list("B" = B,
       "L" = L,
       "Bt" = Bt,
       "Bs" = Bs,
       "Cs" = Cs)
}

# sample the fixed effects, spatial fixed effects and the factor loadings
sample_BBsL <- function(k, X, Tr, U,
                        G, A, C, sigma_b,
                        Gs, As, Cs, sigma_bs,
                        Ks, Xs_centers,
                        Omega) {

  p <- ncol(X)
  ps <- ncol(Cs)
  d <- ncol(U)
  S <- ncol(Omega)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))
  M_Bs <- t(computeBtcoef(Gs, Tr, As, Cs, matrix(0, S, ps)))

  B <- matrix(0, p, S)
  Bs <- matrix(0, ps, S)
  L <- matrix(0, d, S)

  if(p + ps + d > 0){

    for (s in 1:S) {

      if(model == "continuous"){
        k_current <- k[,s] * Omega[,s]
      } else if(model == "binary"){
        k_current <- k[,s]
      }


      XU <- cbind(X, U)

      b_current <- c(M_B[,s], rep(0, d), rep(0, ps))
      B_current <- diag(1, nrow = p + d + ps)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      invB_current <- diag(1 / diag(B_current), nrow = p + d + ps)

      BBsL <- sampleB_SoR(XU, invB_current, b_current, k_current,
                          Omega[,s], Xs_centers, Ks, ps)

      B[seq_len(p),s] <- BBsL[seq_len(p)]
      L[seq_len(d),s] <- BBsL[p + seq_len(d)]
      Bs[seq_len(ps),s] <- BBsL[p + d + seq_len(ps)]


    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  Bts <- t(Bs) - Tr %*% Gs - As %*% Cs

  list("B" = B,
       "Bt" = Bt,
       "L" = L,
       "Bs" = Bs,
       "Bts" = Bts)
}

# sample traits (observed and unobserved) response to covariates
sample_GC <- function(B, Tr, A, sigma_b){

  p <- nrow(B)
  S <- ncol(B)
  g <- ncol(Tr)
  gt <- ncol(A)

  G <- matrix(0, g, p)
  C <- matrix(0, gt, p)

  if(p > 0){

    tB <- t(B)

    B_prior <- diag(2, nrow = g + gt)
    b_prior <- rep(0, g + gt)

    for (k in 1:p) {

      B_current <- tB[,k]

      if(g + gt > 0){

        TA <- cbind(Tr, A)

        b_prior <- rep(0, g + gt)
        B_prior <- diag(1, nrow = g + gt)

        GC <- sampleBuniv(TA, B_prior, b_prior, B_current, sigma_b)
        G[seq_len(g),k] <- GC[seq_len(g)]
        C[seq_len(gt),k] <- GC[g + seq_len(gt)]

      }

    }

  }



  Bt <- t(B) - Tr %*% G - A %*% C

  list("G" = G,
       "C" = C,
       "Bt" = Bt)
}

# sample unobserved traits
sample_A <- function(B, C, Tr, G, sigma_b){

  p <- nrow(B)
  gt <- ncol(A)
  S <- ncol(B)

  Btilde <- B - t(Tr %*% G)

  A <- matrix(0, S, gt)

  if(p > 0 & gt > 0){

    B_current <- diag(1, nrow = gt)
    b_current <- rep(0, gt)

    for (s in 1:S) {

      A[s,] <- sampleBuniv(t(C), B_current, b_current, Btilde[,s], sigma_b)

    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("A" = A,
       "Bt" = Bt)

}

# sample factor scores
sample_U <- function(k, L, X, B, Omega, model){

  d <- nrow(L)
  n <- nrow(k)

  U <- matrix(NA, n, d)

  if(d > 0){

    if(model == "continuous"){
      k_new <- k - (X %*% B)
    } else if(model == "binary"){
      k_new <- k - Omega * (X %*% B)
    }

    B_current <- diag(1, nrow = d)
    b_current <- rep(0, d)

    for (i in 1:n) {

      U[i,] <- sampleB(t(L), B_current, b_current, Omega[i,], k_new[i,])

    }
  }

  U

}


loglik_spatialEffect <- function(KsBs_s, Lm1, logdet, sigma_s){

  xsq <- t(KsBs_s) %*% Lm1

  n_p <- nrow(Lm1)

  loglikelihood <- - n_p / 2 * (log(2 * pi) + log(sigma_s^2)) - .5 * logdet -
    (1/2) * (1 / sigma_s^2) * (xsq %*% t(xsq))

  loglikelihood
}

# sample scale parameter of spatial field
sample_ls <- function(idx_ls, KsBs, list_SoRSummaries,
                      a_l_s, b_l_s, sigma_s){

  if(!is.null(list_SoRSummaries)){

    l_s_grid <- list_SoRSummaries$l_s_grid
    ldet_grid <- list_SoRSummaries$logDetKuu_grid
    Lm1_grid <- list_SoRSummaries$Lm1_grid

    if(idx_ls == 1){
      idx_ls_star <- 2
    } else if(idx_ls == length(l_s_grid)){
      idx_ls_star <- length(l_s_grid) - 1
    } else {
      idx_ls_star <- ifelse(runif(1) < .5, idx_ls - 1, idx_ls + 1)
    }

    # current point
    l_s_current <- l_s_grid[idx_ls]

    loglikelihood_current <- sum(
      sapply(1:S, function(s){
        loglik_spatialEffect(KsBs[,s], Lm1_grid[,,idx_ls], ldet_grid[idx_ls], sigma_s, n_data)
      })
    )

    logPrior <- dgamma(l_s_current, a_l_s, b_l_s, log = T)

    logposterior_current <- logPrior + loglikelihood

    # proposed point
    l_s_star <- l_s_grid[idx_ls_star]

    loglikelihood_star <- sum(
      sapply(1:S, function(s){
        loglik_spatialEffect(KsBs[,s], sq_grid[,,idx_ls_star], ldet_grid[idx_ls_star], sigma_s, n_data)
      })
    )

    logPrior <- dgamma(l_s_star, a_l_s, b_l_s, log = T)

    logposterior_star <- logPrior + loglikelihood

    if(runif(1) < exp(logposterior_star - logposterior_current)){
      idx_ls <- idx_ls_star
    }
  }

  idx_ls
}

update_jSDMcoef <- function(list_data,
                            list_params,
                            list_priors,
                            list_Xs,
                            list_SoRSummaries,
                            model){

  # read data
  {
    z <- list_data$z
    X <- list_data$X
    Tr <- list_data$Tr
  }

  # read params
  {
    B <- list_params$B
    G <- list_params$G
    A <- list_params$A
    C <- list_params$C
    Bt <- list_params$Bt
    Bs <- list_params$B
    Gs <- list_params$G
    As <- list_params$As
    Cs <- list_params$Cs
    Bst <- list_params$Bst
    U <- list_params$U
    L <- list_params$L
    sigma_b <- list_params$sigma_b
    sigma_bs <- list_params$sigma_bs
    idx_ls <- list_params$idx_ls
    tau <- list_params$tau

    l_s <- list_SoRSummaries$l_s_grid[idx_ls]
    Ks <- list_SoRSummaries$Ks_all[,,idx_ls]
  }

  # read priors
  {
    a_tau <- list_priors$a_tau
    b_tau <- list_priors$b_tau
    a_sigmab <- list_priors$a_sigmab
    b_sigmab <- list_priors$b_sigmab
    a_sigmabs <- list_priors$a_sigmabs
    b_sigmabs <- list_priors$b_sigmabs
    a_l_s <- list_priors$a_l_s
    b_l_s <- list_priors$b_l_s
  }

  # read state variables
  {
    ps <- ncol(Bs)
  }

  # define transformed target variables
  if(model == "continuous"){
    k <- z
  } else if(model == "binary"){
    k <- z - .5
  } else if(model == "count"){
    # k <- k / w
  }

  # compute linear predictor
  list_psiCoef <- computePsiCoef(
    X, Ks, list_Xs$Xs_centers, Tr,
    G, A, C, Bt,
    Gs, As, Cs, Bst,
    U, L)
  psiCoef <- list_psiCoef$eta
  XB <- list_psiCoef$XB
  SE <- list_psiCoef$SE
  UL <- list_psiCoef$UL
  variancePartitioning <- computeVariancePartitioning(XB, SE, UL)

  # sample variance of continuous output
  if(model == "continuous"){
    tau <- sample_tau(z, psiCoef, a_tau, b_tau)
  }

  # sample Omega
  if(model == "continuous"){
    Omega <- matrix(1 / tau^2, nrow(z), ncol(z), byrow = T)
  } else if(model == "binary"){
    Omega <- samplePGvariables(psiCoef)
  }

  # sample fixed effects, spatial trait loadings and factor loadings
  list_BBsL <- sample_BBsL(k, X, Tr, U,
                           G, A, C, sigma_b,
                           Gs, As, Cs, sigma_bs,
                           Ks, list_XsXs_centers,
                           Omega)
  B <- list_BBsL$B
  Bt <- list_BBsL$Bt
  Bs <- list_BBsL$Bs
  Bst <- list_BBsL$Bts
  L <- list_BBsL$L

  # update variance of residuals of environmental covariates
  sigma_b <- sample_sigmab(B, Tr, G, A, C, a_sigmab, b_sigmab)
  if(ps > 0){
    sigma_bs <- sample_sigmab(Bs, Tr, Gs, As, Cs, a_sigmabs, b_sigmabs)
  }

  # sample response to traits (observed and unobsered)
  list_GC <- sample_GC(B, Tr, A, sigma_b)
  G <- list_GC$G
  C <- list_GC$C
  Bt <- list_GC$Bt

  # sample unobserved traits
  list_A <- sample_A(B, C, Tr, G, sigma_b)
  A <- list_A$A
  Bt <- list_A$Bt

  # sample spatial response to traits (observed and unobsered)
  if(ps > 0){
    list_GCs <- sample_GC(Bs, Tr, As, sigma_bs)
    Gs <- list_GCs$G
    Cs <- list_GCs$C
    Bst <- list_GCs$Bt
  }

  # sample unobserved spatial traits
  if(ps > 0){
    list_As <- sample_A(Bs, Cs, Tr, Gs, sigma_bs)
    As <- list_As$A
    Bst <- list_As$Bt
  }

  # sample factor scores
  U <- sample_U_cpp(k, L, X, B, Omega, model)

  # sample spatial field scale
  if(ps > 0){
    idx_ls <- sample_ls(idx_ls, l_s_grid,
                        SE, list_SoRSummaries,
                        a_l_s, b_l_s, sigma_s = 1)
    l_s <- l_s_grid[idx_ls]
    Ks <- list_SoRSummaries$Ks_all[,,idx_ls]
  }

  # output params

  list_params <- list(
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
   "tau" = tau,
   "variancePartitioning" = variancePartitioning
  )

  list_params

}

# OUTPUT -------------

computePredictiveProbs <- function(X_new){



}

computeVariancePartitioning <- function(XB, SE, UL){

  V <- function(eta)
    apply(logistic(eta), 2, var)

  S <- ncol(XB)

  # Variances for every subset
  V0   <- rep(0, S)
  VE   <- V(XB)
  VS   <- V(SE)
  VF   <- V(UL)
  VES  <- V(XB + SE)
  VEF  <- V(XB + UL)
  VSF  <- V(SE + UL)
  VESF <- V(XB + SE + UL)

  # Shapley contributions

  CE <-
    (VE - V0 +
       (VES - VS) +
       (VEF - VF) +
       (VESF - VSF)) / 4

  CS <-
    (VS - V0 +
       (VES - VE) +
       (VSF - VF) +
       (VESF - VEF)) / 4

  CF <-
    (VF - V0 +
       (VEF - VE) +
       (VSF - VS) +
       (VESF - VES)) / 4

  Total <- CE + CS + CF

  out <- data.frame(
    Environmental = CE / Total,
    Spatial = CS / Total,
    Biotic = CF / Total,
    Total = Total
  )

  out

}

plotVariancePartitioning <- function(varPart_output, speciesNames){

  varPart_output_vec <- apply(varPart_output, c(1,2), c)
  varPart_output_mean <- apply(varPart_output_vec, c(2,3), mean)
  varPart_output_sd <- apply(varPart_output_vec, 2, function(x){
    sqrt(sum(apply(x[,1:3], 2, function(y){
      var(y)
    })))
  })

  vp <- data.frame(
    Species = speciesNames,
    Env = varPart_output_mean[,1],
    Spatial = varPart_output_mean[,2],
    Biotic = varPart_output_mean[,3],
    StDev = varPart_output_sd
  )

  ggtern(vp,
         aes(x = Env,
             y = Spatial,
             z = Biotic,
             size = StDev)) +
    geom_point(alpha = 0.7) +
    labs(
      T = "Spatial",
      L = "Environment",
      R = "Biotic"
    ) +
    theme_bw() +
    theme_showarrows() +
    theme(
      legend.position = "none",
      tern.axis.title.T = element_text(size = 13),
      tern.axis.title.L = element_text(size = 13),
      tern.axis.title.R = element_text(size = 13)
    )

}

plotCovariateTrend <- function(idxCov, idxSpecies, X0, B, list_ns){

  covNames <- colnames(X0)

  gridVals <- seq(
    min(X0[,idxCov]),
    max(X0[,idxCov]),
    length.out = 100)

  X_new <- createSplinesMatrixSingleCov(list_ns[[idxCov]],
                                        matrix(gridVals, length(gridVals), 1))

  df <- ncol(X_new)

  idxCoeffs <- (idxCov - 1) * df + 1:df
  covEffect <- as.data.frame(as.matrix(X_new) %*% B[idxCoeffs, idxSpecies] )
  covEffect$x <- gridVals

  df_long <- pivot_longer(
    covEffect,
    cols = -x,
    names_to = "Species",
    values_to = "Value"
  )

  ggplot(df_long, aes(x = x, y = Value, colour = Species)) +
    geom_line() +
    theme_bw() + xlab(covNames[idxCov]) + ylab("Covariate Effect")
}

returnSpatialEffectMean <- function(Bs_output, Ks){

  Bs_output_vec <- apply(Bs_output, c(1,2), c)
  Bs_output_vec <- aperm(Bs_output_vec, c(2,3,1))

  spatEffect_mean <- spatEffectMeanCpp(Bs_output_vec, Ks, Xs_centers)

  spatEffect_mean

}

plotSpatialEffect <- function(spatEffect_output, Xs, idx_species = 1){

  spatEffect_output <- returnSpatialEffect(Bs_output, Ks)

  spatEffect_median <- apply(spatEffect_output, c(2,3),
                             function(x) {quantile(x, probs = .5)})

  ggplot(data = NULL, aes(x = Xs[,1],
                          y = Xs[,2],
                          color = spatEffect_median[,idx_species])) + geom_point()

}

# DEPRECATED ---------

# sample traits (observed and unobserved) response to covariates
sample_GC_fixed <- function(B, Tr, A, sigma_b){

  p <- nrow(B)
  S <- ncol(B)
  g <- ncol(Tr)
  gt <- ncol(A)

  tB <- t(B)

  G <- matrix(0, g, p)
  C <- matrix(0, gt, p)

  B_prior <- diag(2, nrow = g + gt)
  b_prior <- rep(0, g + gt)

  for (k in 1:p) {

    elemZero <- max(gt - k, 0)
    elem1 <- ifelse(k <= gt, 1, 0)
    elemNonZero <- gt - elemZero - elem1

    if(elem1 > 0){

      B_current <- tB[,k] - A[,k]

    } else {

      B_current <- tB[,k]

    }

    if(g + elemNonZero > 0){

      TA <- cbind(Tr,
                  matrix(A[,seq_len(elemNonZero)], S, elemNonZero))

      b_prior <- rep(0, g + elemNonZero)
      B_prior <- diag(1, nrow = g + elemNonZero)

      GC <- sampleBuniv(TA, B_prior, b_prior, B_current, sigma_b)
      G[seq_len(g),k] <- GC[seq_len(g)]
      C[seq_len(elemNonZero),k] <- GC[g + seq_len(elemNonZero)]

    }

    if(elem1 > 0) C[k,k] <- 1

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("G" = G,
       "C" = C,
       "Bt" = Bt)
}

# multivariate sample of a matrix normal regression
sampleB_m <- function(k, X, eta, Omega, B, b){

  S <- ncol(Omega)
  B_output <- matrix(NA, length(b), S)

  for (s in 1:S) {

    k_current <- k[,s] - Omega[,s] * eta[,s]

    B_output[,s] <- sampleB(X, B, b, Omega[,s], k_current)

  }

  B_output
}

# sample the ordination covariate coefficients
sample_Br <- function(k, L, Xr, Utilde, eta, Omega){

  ncov_ord <- ncol(Xr)
  d <- nrow(L)
  S <- ncol(L)

  Omega_vec <- as.vector(Omega)

  ktilde <- k - (eta + Utilde %*% L) * Omega

  X <- kronecker(t(L), Xr)

  if(ncov_ord > 0 & d > 0){
    Br <- sample_beta_cpp(X,
                          diag(100, nrow = ncol(X)),
                          rep(0, ncol(X)),
                          Omega_vec,
                          ktilde)

  }

  matrix(Br, ncov_ord, d, byrow = F)
}


# sample the fixed effects and the factor loadings
sample_BL_fixed <- function(k, X, H, G, Tr,
                            A, C, sigma_b,
                            Omega) {

  p <- ncol(X)
  d <- ncol(U)
  S <- ncol(Omega)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))

  B <- matrix(0, p, S)
  L <- matrix(0, d, S)

  for (s in 1:S) {

    # elemZero <- max(d - s, 0)
    # elem1 <- ifelse(s <= d, 1, 0)
    # elemNonZero <- d - elemZero - elem1
    #
    # if(elem1 > 0){
    #
    #   k_current <- k[,s] - Omega[,s] * U[,s]
    #
    # } else {
    #
    k_current <- k[,s]
    #
    # }

    elemNonZero <- d

    if(p + elemNonZero > 0){

      XU <- cbind(X, matrix(U[,seq_len(elemNonZero)], n, elemNonZero))

      b_current <- c(M_B[,s], rep(0, elemNonZero))
      B_current <- diag(1, nrow = p + elemNonZero)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      BL <- sampleB(XU, B_current, b_current, Omega[,s], k_current)
      B[seq_len(p),s] <- BL[seq_len(p)]
      L[seq_len(elemNonZero),s] <- BL[p + seq_len(elemNonZero)]

    }

    # if(elem1 > 0) L[s,s] <- 1

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("B" = B,
       "L" = L,
       "Bt" = Bt)
}

