# Historical helpers retired on 7 October 2026 from R/jsdmfun.R.
# Original source: 0c5b1ff19c4de43be329fbd06730aee24548ca60.
# Preserved verbatim for reference; not loaded or compiled by the package.


createSplinesObjects <- function(X, df){

  list_ns <- lapply(seq_len(ncol(X)), function(j) {
    ns_j <- bs(X[, j], df = 5, intercept = F)
    ns_j
  })

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


# SPATIAL FUNCTIONS -----------

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


# sample size parameter of responses
sample_rnb <- function(z, eta, tune_sd = 5){

  n <- nrow(z)
  S <- ncol(z)

  mu <- exp(eta)

  rnb <- sapply(1:S, function(s){

    r_current <- rnb[s]

    r_star <- exp(rnorm(1, mean = log(r_current), sd = tune_sd))

    ll_star <- sum(dnbinom(z[,s], size = r_star, mu = mu[,s], log = TRUE))
    ll_current <- sum(dnbinom(z[,s], size = r_current, mu = mu[,s], log = TRUE))

    lp_star <- 0#dgamma(r_star, shape = prior_shape, rate = prior_rate, log = TRUE)
    lp_current <- 0#dgamma(r_current, shape = prior_shape, rate = prior_rate, log = TRUE)

    jacobian <- log(r_star) - log(r_current)

    log_alpha <- (ll_star - ll_current) + (lp_star - lp_current) + jacobian

    # 6. Accept or reject
    if (log(runif(1)) < log_alpha) {
      return(r_star)
    } else {
      return(r_current)
    }

  })

  rnb
}


computeVariancePartitioning_R2 <- function(XB, SE, UL, y, model){

  S <- ncol(XB)

  if(model == "binary"){
    # Variances for every subset
    R20   <- rep(0, S)
    R2E   <- pseudo_R2(XB, y)
    R2S   <- pseudo_R2(SE, y)
    R2F   <- pseudo_R2(UL, y)
    R2ES  <- pseudo_R2(XB + SE, y)
    R2EF  <- pseudo_R2(XB + UL, y)
    R2SF  <- pseudo_R2(SE + UL, y)
    R2ESF <- pseudo_R2(XB + SE + UL, y)

    # Shapley contributions

    CE <-
      (R2E - R20 +
         (R2ES - R2S) +
         (R2EF - R2F) +
         (R2ESF - R2SF)) / 4

    CS <-
      (R2S - R20 +
         (R2ES - R2E) +
         (R2SF - R2F) +
         (R2ESF - R2EF)) / 4

    CF <-
      (R2F - R20 +
         (R2EF - R2E) +
         (R2SF - R2S) +
         (R2ESF - R2ES)) / 4

    Total <- CE + CS + CF
  } else {
    CE <- NA
    CS <- NA
    CF <- NA
    Total <- NA
  }

  out <- data.frame(
    Environmental = CE / Total,
    Spatial = CS / Total,
    Biotic = CF / Total,
    Total = Total
  )

  out

}


computeVariancePartitioning_linvars <- function(XB, SE, UL){

  S <- ncol(XB)

  out <- data.frame(
    Environmental = rep(NA, S),
    Spatial = rep(NA, S),
    Biotic = rep(NA, S),
    Total = rep(NA, S)
  )

  for (s in 1:S) {

    alpha_12 <- var(XB[,s]) / (var(XB[,s]) + var(SE[,s]))
    alpha_23 <- var(SE[,s]) / (var(UL[,s]) + var(SE[,s]))
    alpha_13 <- var(XB[,s]) / (var(XB[,s]) + var(UL[,s]))

    total_var <- var(XB[,s]+SE[,s]+UL[,s])
    contrib <- c(var(XB[,s])+ alpha_12 * cov(XB[,s],SE[,s])+ alpha_13 * cov(XB[,s],UL[,s]),
                 var(SE[,s])+ (1 - alpha_12) * cov(XB[,s],SE[,s])+ alpha_23 * cov(SE[,s],UL[,s]),
                 var(UL[,s])+ (1 - alpha_13) *cov(XB[,s],UL[,s])+ (1 - alpha_23) * cov(SE[,s],UL[,s]))
    percent <- contrib / total_var

    out$Environmental[s] <- percent[1]
    out$Spatial[s] <- percent[2]
    out$Biotic[s] <- percent[3]
    out$Total[s] <- total_var

  }

  out

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
