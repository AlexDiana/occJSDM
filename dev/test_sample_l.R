n <- 500
S <- 10
p <- 10 # environmental covariates
g <- 6 # observed traits
gt <- 2 # unobserved traits
d <- 0 # number of factors
ds <- 2 # latent spatial factors

model <- "continuous" # "binary" #  "count"
usingSplines <- T

# variances
{
  if(model == "continuous"){
    tau <- rep(0.01, S)# rgamma(S, 5, 5)
    rnb <- NULL
  } else if(model == "count"){
    rnb <- rpois(S, lambda = 50)
    tau <- NULL
  } else {
    tau <- NULL
    rnb <- NULL
  }

  # variation of residual environmental covariates
  sigma_b <- .01

  # variation of spatial traits
  sigma_ts <- .0001

  # variation of residual spatial field
  sigma_bs <- .0001

  # variation of factor scores
  sigma_h <- 1

  # variation of spatial field
  sigma_s <- .5

  # spatial field scale
  length_grid_ls <- 10
  l_s_grid <- seq(0.01, 0.4, length.out = length_grid_ls)
  idx_ls <- 10
  l_s <- l_s_grid[idx_ls]

}

length_grid_ls <- 20
l_s_grid <- seq(0.01, 0.4, length.out = length_grid_ls)
idx_ls <- 3
l_s <- l_s_grid[idx_ls]
l_s_true <- l_s

a_l_s <- 1
b_l_s <- 1

list_simData <- simulateData(
  n, S, p,
  g, gt, d, tau, rnb, ds, n,
  sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
  useSpatField = T, usingSplines, model)

# data
{
  z <- list_simData$data$z
  X <- list_simData$data$X
  Tr <- list_simData$data$Tr
  Xs <- list_simData$data$Xs
}

ps <- getDefaultSupportPoints(n)

# Spatial covariates matrix
list_Xs <- computeSpatialSummaries(Xs, ps, maxPoints = 5)
Xs_centers <- list_Xs$Xs_centers
Xs_index <- list_Xs$Xs_index
X_s_centers <- list_Xs$X_s_centers
X_tilde <- list_Xs$X_tilde
X_s <- list_Xs$X_s

list_SoRSummaries <- precomputeSORmatrices(l_s_grid, list_Xs)

# simulate spatial field
K_mat <- K2(Xs, Xs, sigma_s, l_s) + diag(10^(-5), nrow = ns)
LU <- chol(K_mat)
Z <- rnorm(n)
SE <- matrix(t(LU) %*% Z, n, S, byrow = F)

niter <- 1000

for (iter in 1:niter) {
  print(idx_ls)

  #
  if(F){
    S <- ncol(SE)

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
        loglik_spatialEffect(SE[,s], Lm1_grid[,,idx_ls], ldet_grid[idx_ls], sigma_s)
      })
    )

    logPrior_current <- dgamma(l_s_current, a_l_s, b_l_s, log = T)

    logposterior_current <- logPrior_current + loglikelihood_current

    # proposed point
    l_s_star <- l_s_grid[idx_ls_star]

    loglikelihood_star <- sum(
      sapply(1:S, function(s){
        loglik_spatialEffect(SE[,s], Lm1_grid[,,idx_ls_star], ldet_grid[idx_ls_star], sigma_s)
      })
    )

    logPrior_star <- dgamma(l_s_star, a_l_s, b_l_s, log = T)

    logposterior_star <- logPrior_star + loglikelihood_star

    if(runif(1) < exp(logposterior_star - logposterior_current)){
      idx_ls <- idx_ls_star
    }
  }

  idx_ls <- sample_ls(idx_ls, SE,
                      list_SoRSummaries,
                      a_l_s, b_l_s, sigma_s = sigma_s)
  l_s <- list_SoRSummaries$l_s_grid[idx_ls]
  Ks <- list_SoRSummaries$Ks_all[,,idx_ls]
}

