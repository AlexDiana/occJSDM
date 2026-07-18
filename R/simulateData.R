

logit <- function(x) log(x / (1-x))

logistic <- function(x) 1 / (1 + exp(-x))

#' simulateOccPlusData
#'
#' Simulate data
#'
#' @param list_datasettings List of data dimension settings (n, S, g, M, P, K,
#'   ncov_psi, ncov_theta).
#' @param list_params List of detection/amplification parameters (p, q, theta0,
#'   theta_baseline). Optionally also `mu1`, `sigma1`, `mu0`, `sigma0`
#'   controlling the simulated read-count intensities (see Details); if
#'   omitted these default to `mu1 = 5`, `sigma1 = 1`, `mu0 = 1.5`,
#'   `sigma0 = 1`.
#' @param list_jsdmParams List of JSDM parameters (gt, d, ds, sigma_b, sigma_bs,
#'   sigma_ts, sigma_h, l_s).
#' @param useSpatField Logical. If `TRUE`, simulated occupancy incorporates a
#'   spatially autocorrelated random field over the site coordinates `Xs`
#'   (using `ds` latent spatial factors and length scale `l_s`). If `FALSE`
#'   (the default), site coordinates are still simulated and returned but do
#'   not influence occupancy.
#'
#' @details
#' Simulate data. `data_list$y` contains simulated read counts (matching the
#' `threshold = 0` mode of `runOccPlus()`), rather than binary detections:
#' for a true detection, `log(y + 1) ~ Normal(mu1, sigma1)`; for a
#' false-positive/contamination detection, `log(y + 1) ~ Normal(mu0, sigma0)`;
#' otherwise `y = 0`. The realised counts are `round(exp(log(y + 1)) - 1)`.
#'
#' @return Description of the return value (e.g., a list, data frame, or numeric output).
#'
#'
#' @export
#' @import tidyverse
#'
simulateOccPlusData <- function(list_datasettings,
                         list_params,
                         list_jsdmParams,
                         useSpatField = FALSE){
  # read data settings
  {
    n <- list_datasettings$n
    S <- list_datasettings$S
    g <- list_datasettings$g
    M <- list_datasettings$M
    P <- list_datasettings$P
    K <- list_datasettings$K
    ncov_psi <- list_datasettings$ncov_psi
    ncov_theta <- list_datasettings$ncov_theta
  }

  # read jsdm params
  {
    gt <- list_jsdmParams$gt
    d <- list_jsdmParams$d
    ds <- list_jsdmParams$ds
    sigma_b <- list_jsdmParams$sigma_b
    sigma_bs <- list_jsdmParams$sigma_bs
    sigma_ts <- list_jsdmParams$sigma_ts
    sigma_h <- list_jsdmParams$sigma_h
    sigma_s <- list_jsdmParams$sigma_s
    l_s <- list_jsdmParams$l_s
  }

  # read param settings
  {
    p <- list_params$p
    q <- list_params$q
    theta0 <- list_params$theta0
    theta_baseline <- list_params$theta_baseline
  }


  list_simjSDMData <- simulateData(
    n, S, ncov_psi,
    g, gt, d, tau = NULL, ds, n,
    sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
    useSpatField = useSpatField, usingSplines = F, model = "binary")

  # data
  {
    z <- list_simjSDMData$data$z
    X_psi <- list_simjSDMData$data$X
    Tr <- list_simjSDMData$data$Tr
    Xs <- list_simjSDMData$data$Xs
  }

  # create data structure
  {
    N <- sum(M)
    N2 <- P * N
    N3 <- sum(K)

    sumM <- c(0, cumsum(M)[-n])
    sumP <- c(0, cumsum(rep(P, N))[-N])
    sumK <- c(0, cumsum(K)[-N2])

    list_idx <- createDataIdx(n, M, P, K)
    idx_z_w <- list_idx$idx_z_w
    idx_z_k <- list_idx$idx_z_k
    idx_w_p <- list_idx$idx_w_p
    idx_z_p <- list_idx$idx_z_p
    idx_w_k <- list_idx$idx_w_k
    idx_p_k <- list_idx$idx_p_k
  }

  X_theta <- cbind(1, matrix(rnorm(N * (ncov_theta - 1)), N, ncov_theta - 1))

  beta_theta_true <- matrix(sample(c(-1,1,0), ncov_theta * S, replace = T), ncov_theta, S)
  beta_theta_true[1,] <- logit(theta_baseline)
  theta_true <- logistic(X_theta %*% beta_theta_true)

  theta0_true <- theta0
  p_true <- p
  q_true <- q

  z_rep <- z[idx_z_w,]
  y <- matrix(NA, N3, S)

  w <- sapply(1:S, function(s){
    sapply(1:N, function(i){
      if(z_rep[i,s] == 1){
        rbinom(1, 1, theta_true[i,s])
      } else {
        rbinom(1, 1, theta0_true[s])
      }
    })
  })

  # L_rep <- rep(1:N, L)
  # LK_rep <- rep(L_rep, K)
  w_rep <- w[idx_w_k,]

  # idx_z <- rep(1:n, M)
  # idx_k <- as.numeric(as.factor(data_info_sample))

  # primer_rep <- rep(rep(1:L[1], each = K[1]), times = N) # only right in some cases

  cimk_true <- sapply(1:S, function(s){
    sapply(1:N3, function(i){
      if(w_rep[i,s] == 1){
        rbinom(1, 1, p_true[idx_p_k[i],s])
      } else {
        2 * rbinom(1, 1, q_true[idx_p_k[i],s])
      }
    })
  })

  # Read-count model (matches the threshold = 0 mode of runOccPlus()): given
  # a true detection (cimk_true == 1), log(y + 1) ~ Normal(mu1, sigma1); given
  # a false-positive/contamination detection (cimk_true == 2),
  # log(y + 1) ~ Normal(mu0, sigma0). No reads (cimk_true == 0) gives y = 0.
  mu1_true <- get_param(list_params, "mu1", 5)
  sigma1_true <- get_param(list_params, "sigma1", 1)
  mu0_true <- get_param(list_params, "mu0", 1.5)
  sigma0_true <- get_param(list_params, "sigma0", 1)

  logy1 <- matrix(0, N3, S)
  logy1[cimk_true == 1] <- rnorm(sum(cimk_true == 1), mu1_true, sigma1_true)
  logy1[cimk_true == 2] <- rnorm(sum(cimk_true == 2), mu0_true, sigma0_true)

  y <- round(exp(logy1) - 1)
  y[y < 0] <- 0

  speciesNames <- paste0("OTU_", 1:S)

  colnames(y) <- speciesNames
  rownames(Tr) <- speciesNames
  colnames(Tr) <- paste0("Trait_", seq_len(ncol(Tr)))

  true_params <- list(
    "beta_theta_true" = beta_theta_true,
    "jsdmParams_true" = list_simjSDMData$trueParams,
    "z_true" = z,
    "w_true" = w,
    "theta_true" = theta_true,
    "p_true" = p_true,
    "q_true" = q_true,
    "mu1_true" = mu1_true,
    "sigma1_true" = sigma1_true,
    "mu0_true" = mu0_true,
    "sigma0_true" = sigma0_true
  )

  data_list <- list(
    "X_psi" = X_psi,
    "X_theta" = X_theta,
    "Xs" = Xs,
    "Tr" = Tr,
    "y" = y
  )

  list(true_params = true_params,
       data_list = data_list)
}

#' toRunOccPlusFormat
#'
#' Convert the output of \code{simulateOccPlusData()} into the
#' \code{info}/\code{OTU} list format expected by \code{runOccPlus()}
#' (the same shape as \code{sampledata}).
#'
#' @param sim Output of \code{simulateOccPlusData()}, i.e. a list with
#'   elements \code{true_params} and \code{data_list}.
#' @param n Number of sites (must match \code{list_datasettings$n} used to
#'   generate \code{sim}).
#' @param M Integer vector of length \code{n}: number of samples per site
#'   (must match \code{list_datasettings$M}).
#' @param P Integer: number of primers per sample (must match
#'   \code{list_datasettings$P}).
#' @param K Integer vector of length \code{P * sum(M)}: number of PCR
#'   replicates per sample/primer combination (must match
#'   \code{list_datasettings$K}).
#' @param drop_theta_intercept Logical. \code{simulateOccPlusData()} builds
#'   \code{X_theta} with an intercept column of 1s in the first position.
#'   \code{runOccPlus()} adds its own intercept internally via
#'   \code{process_covariates()}, so by default (\code{TRUE}) that column is
#'   dropped from the returned \code{info} data.frame.
#'
#' @details
#' \code{simulateOccPlusData()} returns covariates at the site level
#' (\code{X_psi}) and sample level (\code{X_theta}), while \code{runOccPlus()}
#' expects a single \code{info} data.frame with one row per PCR replicate and
#' \code{Site}/\code{Sample}/\code{Primer} id columns. This function expands
#' \code{X_psi}/\code{X_theta} to the PCR-replicate level (using the same
#' indexing as \code{createDataIdx()}) and combines them with the id columns
#' and \code{sim$data_list$y} (renamed to \code{OTU}).
#'
#' Note \code{n}, \code{M}, \code{P}, \code{K} are not stored in \code{sim}
#' itself, so they must be supplied here matching the \code{list_datasettings}
#' originally passed to \code{simulateOccPlusData()}.
#'
#' @return A list with elements \code{info} (data.frame with \code{Site},
#'   \code{Sample}, \code{Primer} id columns plus expanded \code{X_psi.*},
#'   \code{X_theta.*}, and \code{Xs.*} covariates, one row per PCR
#'   replicate), \code{OTU} (the simulated detection matrix, same number of
#'   rows as \code{info}), \code{spatCovariates} (character vector of the
#'   \code{Xs.*} column names in \code{info}, ready to pass as the
#'   \code{spatCovariates} argument of \code{runOccPlus()}), and
#'   \code{traitsMatrix} (the \code{S} x \code{g} trait matrix
#'   \code{sim$data_list$Tr}, with row names matching the species names in
#'   \code{OTU}'s column names, ready to pass as the \code{traitsMatrix}
#'   argument of \code{runOccPlus()}).
#'
#' @export
toRunOccPlusFormat <- function(sim, n, M, P, K, drop_theta_intercept = TRUE) {

  X_psi <- sim$data_list$X_psi
  X_theta <- sim$data_list$X_theta
  Xs <- sim$data_list$Xs
  Tr <- sim$data_list$Tr
  y <- sim$data_list$y

  idx <- createDataIdx(n, M, P, K)

  N3 <- length(idx$idx_z_k)
  if (nrow(y) != N3) {
    stop("nrow(sim$data_list$y) does not match n/M/P/K -- ",
         "these must be the same list_datasettings used to create `sim`.")
  }

  if (drop_theta_intercept) {
    X_theta <- X_theta[, -1, drop = FALSE]
  }

  X_psi_rep <- X_psi[idx$idx_z_k, , drop = FALSE]
  X_theta_rep <- X_theta[idx$idx_w_k, , drop = FALSE]
  Xs_rep <- Xs[idx$idx_z_k, , drop = FALSE]

  colnames(X_psi_rep) <- if (ncol(X_psi_rep) == 1) {
    "X_psi"
  } else {
    paste0("X_psi.", seq_len(ncol(X_psi_rep)))
  }

  colnames(X_theta_rep) <- if (ncol(X_theta_rep) == 1) {
    "X_theta"
  } else {
    paste0("X_theta.", seq_len(ncol(X_theta_rep)))
  }

  colnames(Xs_rep) <- if (ncol(Xs_rep) == 1) {
    "Xs"
  } else {
    paste0("Xs.", seq_len(ncol(Xs_rep)))
  }

  info <- data.frame(
    Site = idx$idx_z_k,
    Sample = idx$idx_w_k,
    Primer = idx$idx_p_k,
    X_psi_rep,
    X_theta_rep,
    Xs_rep,
    check.names = FALSE
  )

  list(
    info = info,
    OTU = y,
    spatCovariates = colnames(Xs_rep),
    traitsMatrix = Tr
  )
}



logit <- function(x) log(x / (1-x))

logistic <- function(x) 1 / (1 + exp(-x))

#' simulateOccJSDMData
#'
#' Simulate data
#'
#' @param list_datasettings List of data dimension settings (n, S, g, M, P, K,
#'   ncov_psi, ncov_theta).
#' @param list_params List of detection/amplification parameters (p, q, theta0,
#'   theta_baseline). Optionally also `mu1`, `sigma1`, `mu0`, `sigma0`
#'   controlling the simulated read-count intensities (see Details); if
#'   omitted these default to `mu1 = 5`, `sigma1 = 1`, `mu0 = 1.5`,
#'   `sigma0 = 1`.
#' @param list_jsdmParams List of JSDM parameters (gt, d, ds, sigma_b, sigma_bs,
#'   sigma_ts, sigma_h, l_s).
#' @param useSpatField Logical. If `TRUE`, simulated occupancy incorporates a
#'   spatially autocorrelated random field over the site coordinates `Xs`
#'   (using `ds` latent spatial factors and length scale `l_s`). If `FALSE`
#'   (the default), site coordinates are still simulated and returned but do
#'   not influence occupancy.
#' @param model Supported data formats are "binary","continous","occupancy"
#'  and "two-stage"
#'
#' @details
#' Simulate data. `data_list$y` contains simulated read counts (matching the
#' `threshold = 0` mode of `runOccPlus()`), rather than binary detections:
#' for a true detection, `log(y + 1) ~ Normal(mu1, sigma1)`; for a
#' false-positive/contamination detection, `log(y + 1) ~ Normal(mu0, sigma0)`;
#' otherwise `y = 0`. The realised counts are `round(exp(log(y + 1)) - 1)`.
#'
#' @return Description of the return value (e.g., a list, data frame, or numeric output).
#'
#'
#' @export
#' @import tidyverse
#'
simulateOccJSDMData <- function(list_datasettings,
                                list_params,
                                list_jsdmParams,
                                model){

  if(!(model %in% c("binary","occupancy","continuous","two_stage"))){
    stop("The model used is not supported")
  }

  # read data settings
  {
    n <- list_datasettings$n
    S <- list_datasettings$S
    g <- list_datasettings$g
    M <- list_datasettings$M
    P <- list_datasettings$P
    K <- list_datasettings$K
    ncov_psi <- list_datasettings$ncov_psi
    ncov_theta <- list_datasettings$ncov_theta
  }

  # read jsdm params
  {
    gt <- list_jsdmParams$gt
    d <- list_jsdmParams$d
    ds <- list_jsdmParams$ds
    sigma_b <- list_jsdmParams$sigma_b
    sigma_bs <- list_jsdmParams$sigma_bs
    sigma_ts <- list_jsdmParams$sigma_ts
    sigma_h <- list_jsdmParams$sigma_h
    sigma_s <- list_jsdmParams$sigma_s
    l_s <- list_jsdmParams$l_s
    tau <- list_jsdmParams$tau
    useSpatField <- list_jsdmParams$useSpatField
  }

  # read param settings
  {
    p <- list_params$p
    q <- list_params$q
    theta0 <- list_params$theta0
    theta_baseline <- list_params$theta_baseline
  }

  if(model %in% c("binary","occupancy","two_stage")){
    jSDMsimModel <- "binary"
  } else if(model == "continuous"){
    jSDMsimModel <- "continuous"
  }

  list_simjSDMData <- simulateData(
    n, S, ncov_psi,
    g, gt, d, tau, ds, n,
    sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
    useSpatField = useSpatField, usingSplines = F, model = jSDMsimModel)

  # data
  {
    z <- list_simjSDMData$data$z
    X_psi <- list_simjSDMData$data$X
    Tr <- list_simjSDMData$data$Tr
    Xs <- list_simjSDMData$data$Xs
  }

  # create data structure (for occupancy or two-stage model)
  {
    N <- sum(M)

    if(model %in% c("occupancy","two_stage")){
      N2 <- P * N
      N3 <- sum(K)

      sumM <- c(0, cumsum(M)[-n])
      sumP <- c(0, cumsum(rep(P, N))[-N])
      sumK <- c(0, cumsum(K)[-N2])

      list_idx <- createDataIdx(n, M, P, K, model == "two_stage")
      idx_z_w <- list_idx$idx_z_w
      idx_z_k <- list_idx$idx_z_k
      idx_w_p <- list_idx$idx_w_p
      idx_z_p <- list_idx$idx_z_p
      idx_w_k <- list_idx$idx_w_k
      idx_p_k <- list_idx$idx_p_k
    }

  }

  if(model %in% c("occupancy","two_stage")){
    X_theta <- cbind(1, matrix(rnorm(N * (ncov_theta - 1)), N, ncov_theta - 1))

    beta_theta_true <- matrix(sample(c(-1,1,0), ncov_theta * S, replace = T), ncov_theta, S)
    beta_theta_true[1,] <- logit(theta_baseline)
    theta_true <- logistic(X_theta %*% beta_theta_true)

    theta0_true <- theta0

    z_rep <- z[idx_z_w,]
    w <- sapply(1:S, function(s){
      sapply(1:N, function(i){
        if(z_rep[i,s] == 1){
          rbinom(1, 1, theta_true[i,s])
        } else {
          rbinom(1, 1, theta0_true[s])
        }
      })
    })

    if(model == "two_stage"){
      p_true <- p
      q_true <- q


      y <- matrix(NA, N3, S)
      w_rep <- w[idx_w_k,]

      cimk_true <- sapply(1:S, function(s){
        sapply(1:N3, function(i){
          if(w_rep[i,s] == 1){
            rbinom(1, 1, p_true[idx_p_k[i],s])
          } else {
            2 * rbinom(1, 1, q_true[idx_p_k[i],s])
          }
        })
      })

      # Read-count model (matches the threshold = 0 mode of runOccPlus()): given
      # a true detection (cimk_true == 1), log(y + 1) ~ Normal(mu1, sigma1); given
      # a false-positive/contamination detection (cimk_true == 2),
      # log(y + 1) ~ Normal(mu0, sigma0). No reads (cimk_true == 0) gives y = 0.
      mu1_true <- get_param(list_params, "mu1", 5)
      sigma1_true <- get_param(list_params, "sigma1", 1)
      mu0_true <- get_param(list_params, "mu0", 1.5)
      sigma0_true <- get_param(list_params, "sigma0", 1)

      logy1 <- matrix(0, N3, S)
      logy1[cimk_true == 1] <- rnorm(sum(cimk_true == 1), mu1_true, sigma1_true)
      logy1[cimk_true == 2] <- rnorm(sum(cimk_true == 2), mu0_true, sigma0_true)

      y <- round(exp(logy1) - 1)
      y[y < 0] <- 0
    }

  }

  if(model %in% c("continuous","binary")){
    y <- z
  } else if(model == "occupancy"){
    y <- w
  } else if(model == "two_stage"){
    y <- y
  }

  speciesNames <- paste0("OTU_", 1:S)

  colnames(y) <- speciesNames
  rownames(Tr) <- speciesNames
  colnames(Tr) <- paste0("Trait_", seq_len(ncol(Tr)))

  # data_list <- list(
  #   "X_psi" = X_psi,
  #   "X_theta" = X_theta,
  #   "Xs" = Xs,
  #   "Tr" = Tr,
  #   "y" = y
  # )
  #
  # data_list <- list_simulatedData$data_list
  # X_psi <- data_list$X_psi
  # X_theta <- data_list$X_theta
  # Tr <- data_list$Tr
  # y <- data_list$y

  if(model %in% c("occupancy","two_stage")){
    list_idx <- createDataIdx(n, M, P, K, model == "two_stage")
    idx_z_k <- list_idx$idx_z_k
    idx_w_k <- list_idx$idx_w_k
    idx_p_k <- list_idx$idx_p_k
  }

  if(model %in% c("binary","continuous")){
    data_info <- data.frame(
      X_psi = X_psi,
      Xs = Xs
    )
  } else if(model %in% c("occupancy")){
    data_info <- data.frame(
      Site = idx_z_k,
      X_psi = X_psi[idx_z_k,],
      Xs = Xs,
      X_theta = X_theta[,-1]
    )
  } else if (model == "two_stage"){
    data_info <- data.frame(
      Site = idx_z_k,
      Sample = idx_w_k,
      Primer = idx_p_k,
      X_psi = X_psi[idx_z_k,],
      Xs = Xs,
      X_theta = X_theta[idx_w_k,-1]
    )
  }

  data <- list(info = data_info,
               OTU = y,
               traits = Tr)

  true_params <- list(
    "jsdmParams_true" = list_simjSDMData$trueParams
  )

  if(model %in% c("occupancy","two_stage")){
    true_params$beta_theta_true <- beta_theta_true
    true_params$z_true <- z
  }

  if(model == "two_stage"){
    true_params$w_true <- w
    true_params$p_true <- p_true
    true_params$q_true <- q_true
  }

  list(true_params = true_params,
       data_list = data)
}
