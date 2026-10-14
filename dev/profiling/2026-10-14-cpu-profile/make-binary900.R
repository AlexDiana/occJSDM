# Usage: Rscript make-binary900.R <out.rds>  (OCCJSDM_PKG names the checkout to load)
# One-stage binary spatial dataset: 900 sites, 12 species, spatial field with range 0.08.
suppressMessages(pkgload::load_all(Sys.getenv("OCCJSDM_PKG", unset = "."), quiet = TRUE))
set.seed(5162)
n <- 900; S <- 12; g <- 3; gt <- 2; d <- 2; ds <- 2; ncov_psi <- 2; ncov_theta <- 2
M <- rep(1, n); P <- 1; N <- sum(M); K <- rep(1, P * N)
list_datasettings <- list(n = n, S = S, g = g, M = M, P = P, K = K, ncov_psi = ncov_psi, ncov_theta = ncov_theta)
list_jsdmParams <- list(gt = gt, d = d, ds = ds, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                        sigma_h = 1, sigma_s = 1, l_s = 0.08, useSpatField = TRUE)
list_params <- list(p = matrix(0.5, P, S), q = matrix(0.03, P, S), theta0 = rep(0.05, S), theta_baseline = rep(0.3, S))
sim <- simulateOccJSDMData(list_datasettings, list_params, list_jsdmParams, model = "binary")
data <- sim$data_list
str(data, max.level = 1); print(names(data$info)); print(colMeans(data$OTU))
saveRDS(data, commandArgs(trailingOnly = TRUE)[1])
