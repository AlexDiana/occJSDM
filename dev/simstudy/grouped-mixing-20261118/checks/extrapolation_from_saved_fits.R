# Do chains that agree on the fitted effect at the training sites also agree, and do as well, at new sites?
# Uses saved per-chain posterior means of B and the true B; regenerates the training embedding as repro.R does.
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
root <- "dev/simstudy/grouped-mixing-20261118/results"
sim_X <- function(seed) { n <- 256L; S <- 48L; K <- 64L; set.seed(seed)
  ds <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n), ncov_psi = 5L, ncov_theta = 1L)
  jp <- list(gt = 4L, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5, sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, S), useSpatField = FALSE, embedding = list(n_emb = K, n_blends = 3L, n_active = 10L))
  sim <- simulateOccJSDMData(ds, list(), jp, model = "binary"); as.matrix(sim$data_list$info[, sim$data_list$embeddingVars]) }
rms <- function(x) sqrt(mean(x^2))
fits <- list(c("attempt2/proto", "lowrank_n256_S48_K64_gt4_b5000_s1", 1), c("attempt1/proto", "lowrank_n256_S48_K64_gt4_b5000_s1", 1),
             c("attempt1/base", "lowrank_n256_S48_K64_gt4_b5000_s1", 1), c("attempt1/proto", "lowrank_n256_S48_K64_gt4_b1500_s2", 2), c("attempt2/proto", "lowrank_n256_S48_K64_gt4_b1500_s2", 2))
for (f in fits) {
  x <- readRDS(file.path(root, f[1], paste0(f[2], "_summary.rds"))); X <- sim_X(as.integer(f[3])); emb <- x$emb_rows
  Bt <- x$true_B[emb, ]; B1 <- x$B_chain_mean[emb, , 1]; B2 <- x$B_chain_mean[emb, , 2]
  ev <- eigen(cov(X), symmetric = TRUE); set.seed(99); m <- 2000L
  Z <- matrix(rnorm(m * 64), m, 64)
  sites <- list(training = X,
                "new, like training" = Z %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors),
                "new, weak directions x3" = Z %*% diag(sqrt(pmax(ev$values, 0)) * ifelse(seq_along(ev$values) > 8, 3, 1)) %*% t(ev$vectors),
                "new, all directions equal" = Z)
  cat(sprintf("\n%s %s | log sigma_emb by chain %s | embedding variance in its first 8 directions %.2f\n", f[1], f[2], paste(round(x$log_sigma_emb_chain_mean, 2), collapse = ", "), sum(ev$values[1:8]) / sum(ev$values)))
  for (nm in names(sites)) { Xn <- sites[[nm]]; T <- Xn %*% Bt; E1 <- Xn %*% B1; E2 <- Xn %*% B2
    cat(sprintf("  %-26s size of true effect %5.2f | chains differ by %4.0f%% of it | error of chain 1 %4.0f%%, chain 2 %4.0f%% of it\n", nm, rms(T), 100 * rms(E1 - E2) / rms(T), 100 * rms(E1 - T) / rms(T), 100 * rms(E2 - T) / rms(T))) }
  cat(sprintf("  coefficients themselves: chains differ by %.0f%% of the true coefficients' size; error %0.f%% and %.0f%%\n", 100 * rms(B1 - B2) / rms(Bt), 100 * rms(B1 - Bt) / rms(Bt), 100 * rms(B2 - Bt) / rms(Bt)))
}
