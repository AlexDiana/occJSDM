source("slab_fns.R")
N_DRAWS <- 1500L
K <- 64L
set.seed(77L)
X <- embedding_columns(N_SITES, K)
cat("K = 64, approach A with the axis scale rescaled, S = 100\n")
for (zeta in c(2, 3)) for (mult in c(0.7, 2, 4, 8)) {
  m <- run_A(X, mult / K, zeta, 100L, sweeps = N_DRAWS)
  summ(sprintf("A: zeta = %g, a = %g / K", zeta, mult), m)
}
cat("B with the leftover capped, c = 2, checking it does not depend on a rerun\n")
summ("B+: c = 2", run_B(X, 0.7 / K, 2, TRUE))
summ("B+: c = 2.5", run_B(X, 0.7 / K, 2.5, TRUE))
