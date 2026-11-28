# Scratch prior-only comparison of ways to regularise the axis scales of the
# grouped prior. Not package code. One species' coefficients on K embedding
# columns: beta_k = mu_k + a' c_k + e_k, no measured traits.
#
#   current : lambda_k ~ C+(a), sigma ~ C+(1); sds lambda, lambda/sqrt(q), sigma*lambda
#   A (product form, Nishimura and Suchard): the current prior times an
#       independent Normal(0, zeta^2) factor on mu_k and e_ks and
#       Normal(0, zeta^2/q) on c_jk. The scales' implied prior then depends
#       on the number of species S, so (sigma, lambda) are sampled by
#       alternating exact one-dimensional grid draws.
#   B (Piironen and Vehtari): lambda_k ~ C+(a) kept, effective scale
#       lambda c / sqrt(c^2 + lambda^2).
#   C : half-t(3) on the axis scales at 0.6 / sqrt(K), no slab.

args <- commandArgs(trailingOnly = TRUE)
N_SITES <- 250L
N_DRAWS <- 2000L
q <- 8L

half_cauchy <- function(n, scale) abs(rcauchy(n, 0, scale))

embedding_columns <- function(n, K) {
  Fx <- min(8L, K)
  W <- matrix(rnorm(K * Fx), K, Fx) %*% diag(1 / seq_len(Fx), Fx)
  scale(matrix(rnorm(n * Fx), n, Fx) %*% t(W) + 0.3 * matrix(rnorm(n * K), n, K))
}

contribution <- function(X, sd_mu, sd_c, sd_e) {
  K <- ncol(X)
  mu <- rnorm(K, 0, sd_mu)
  C <- matrix(rnorm(q * K, 0, rep(sd_c, each = q)), q, K)
  a <- rnorm(q)
  e <- rnorm(K, 0, sd_e)
  beta <- mu + as.vector(a %*% C) + e
  c(sd = sd(X %*% beta), maxb = max(abs(beta)))
}

summ <- function(label, m) {
  s <- quantile(m["sd", ], c(0.05, 0.5, 0.95, 0.99), names = FALSE)
  b <- quantile(m["maxb", ], c(0.5, 0.95), names = FALSE)
  cat(sprintf("  %-34s sd: 5%% %5.2f  med %5.2f  95%% %6.2f  99%% %7.1f | max|b|: med %5.2f  95%% %6.2f\n",
              label, s[1], s[2], s[3], s[4], b[1], b[2]))
}

run_current <- function(X, a) vapply(seq_len(N_DRAWS), function(i) {
  l <- half_cauchy(ncol(X), a); s <- half_cauchy(1, 1)
  contribution(X, l, l / sqrt(q), s * l)
}, numeric(2))

run_B <- function(X, a, c, cap_sigma = FALSE) vapply(seq_len(N_DRAWS), function(i) {
  l <- half_cauchy(ncol(X), a); s <- half_cauchy(1, 1)
  lt <- l * c / sqrt(c^2 + l^2)
  sd_e <- if (cap_sigma) 1 / sqrt(1 / (s * l)^2 + 1 / c^2) else s * lt
  contribution(X, lt, lt / sqrt(q), sd_e)
}, numeric(2))

run_C <- function(X, a) vapply(seq_len(N_DRAWS), function(i) {
  l <- a * abs(rt(ncol(X), 3)); s <- half_cauchy(1, 1)
  contribution(X, l, l / sqrt(q), s * l)
}, numeric(2))

# Approach A. Implied prior of the scales with S species:
#   p(sigma, lambda) prop. to C+(sigma; 1) prod_k C+(lambda_k; a)
#     (lambda_k^2 + zeta^2)^(-(1 + q) / 2) (sigma^2 lambda_k^2 + zeta^2)^(-S / 2)
run_A <- function(X, a, zeta, S, sweeps = N_DRAWS, thin = 5L, burn = 200L) {
  K <- ncol(X)
  grid <- seq(log(1e-7), log(1e5), length.out = 3000L)
  lg <- exp(grid)
  draw_grid <- function(logdens) {          # rows: variables, cols: grid
    logdens <- logdens - apply(logdens, 1, max)
    p <- exp(logdens); cs <- t(apply(p, 1, cumsum)); cs <- cs / cs[, ncol(cs)]
    u <- runif(nrow(cs))
    idx <- pmin(rowSums(cs < u) + 1L, ncol(cs))
    exp(grid[idx] + runif(length(idx), -0.5, 0.5) * (grid[2] - grid[1]))
  }
  lam <- rep(a, K); sig <- 1
  out <- matrix(NA_real_, 2, sweeps)
  scales <- matrix(NA_real_, 3, sweeps)
  it <- 0L
  for (t in seq_len(burn + sweeps * thin)) {
    # lambda | sigma, all K at once on a common grid (density in log lambda)
    ld <- grid - log1p((lg / a)^2) - (1 + q) / 2 * log(lg^2 + zeta^2) -
      S / 2 * log(sig^2 * lg^2 + zeta^2)
    lam <- draw_grid(matrix(ld, K, length(grid), byrow = TRUE))
    # sigma | lambda
    ls <- grid - log1p(lg^2) -
      S / 2 * colSums(log(outer(lam^2, lg^2) + zeta^2))
    sig <- draw_grid(matrix(ls, 1))
    if (t > burn && (t - burn) %% thin == 0L) {
      it <- it + 1L
      eff <- function(s) 1 / sqrt(1 / s^2 + 1 / zeta^2)
      out[, it] <- contribution(X, eff(lam), eff(lam) / sqrt(q), eff(sig * lam))
      scales[, it] <- c(sig, median(lam), max(lam))
    }
  }
  rownames(out) <- c("sd", "maxb")
  attr(out, "scales") <- scales
  out
}

for (K in c(64L, 128L)) {
  set.seed(20261125L + K)
  X <- embedding_columns(N_SITES, K)
  a0 <- 0.7 / K
  cat(sprintf("\nK = %d, q = %d, axis scale a = 0.7 / K unless stated\n", K, q))
  summ("current", run_current(X, a0))
  summ("C: half-t(3), 0.6 / sqrt(K)", run_C(X, 0.6 / sqrt(K)))
  for (c in c(1, 2, 3)) {
    summ(sprintf("B: cap c = %g", c), run_B(X, a0, c))
    summ(sprintf("B: cap c = %g, leftover capped too", c), run_B(X, a0, c, TRUE))
  }
  for (S in c(30L, 100L, 1000L)) for (zeta in c(1, 2, 3)) {
    m <- run_A(X, a0, zeta, S)
    summ(sprintf("A: zeta = %g, S = %d", zeta, S), m)
    sc <- attr(m, "scales")
    cat(sprintf("      sigma_emb med %.2f 95%% %.2f | max lambda med %.3f 95%% %.3f\n",
                median(sc[1, ]), quantile(sc[1, ], 0.95), median(sc[3, ]),
                quantile(sc[3, ], 0.95)))
  }
}
