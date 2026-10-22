# Kernel, factorisation and quadrature utilities for the recoverability pilot.
sq_exp <- function(a, b, range) {
  exp(-(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2) / (2 * range^2))
}

half_correlation <- function(range) range * sqrt(2 * log(2))

# Square-root factor of a covariance matrix. Cholesky with jitter when it succeeds and the
# condition number is at most `cond_ceiling`; otherwise a truncated eigendecomposition,
# zero-padded so that `lower` is always n by n and can feed the elliptical-slice samplers.
stable_factor <- function(K, jitter = 1e-8, rel_tol = 1e-10, cond_ceiling = 1e12) {
  K <- (K + t(K)) / 2; n <- nrow(K)
  values <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  condition <- if (min(values) > 0) max(values) / min(values) else Inf
  ch <- tryCatch(chol(K + diag(jitter, n)), error = function(e) NULL)
  if (!is.null(ch) && condition <= cond_ceiling)
    return(list(lower = t(ch), route = "cholesky", condition = condition, rank = n, error = 0))
  ev <- eigen(K, symmetric = TRUE)
  keep <- ev$values > rel_tol * max(ev$values)
  L <- ev$vectors[, keep, drop = FALSE] %*% diag(sqrt(ev$values[keep]), sum(keep))
  list(lower = cbind(L, matrix(0, n, n - sum(keep))), route = "eigen", condition = condition,
       rank = sum(keep), error = max(abs(K - L %*% t(L))))
}

# Conditional distribution of the lattice field given the site field: mean map A and a
# truncated-eigen factor of the conditional covariance. Never forms an explicit inverse.
lattice_conditional <- function(xy_sites, xy_lattice, range, jitter = 1e-8, rel_tol = 1e-10) {
  Kss <- sq_exp(xy_sites, xy_sites, range) + diag(jitter, nrow(xy_sites))
  Kls <- sq_exp(xy_lattice, xy_sites, range)
  ch <- tryCatch(chol(Kss), error = function(e) NULL)
  if (!is.null(ch)) {
    A <- t(backsolve(ch, forwardsolve(t(ch), t(Kls)))); route <- "cholesky"
  } else {
    ev <- eigen(Kss, symmetric = TRUE); keep <- ev$values > rel_tol * max(ev$values)
    V <- ev$vectors[, keep, drop = FALSE]
    A <- Kls %*% V %*% diag(1 / ev$values[keep], sum(keep)) %*% t(V); route <- "eigen"
  }
  C <- sq_exp(xy_lattice, xy_lattice, range) - A %*% t(Kls); C <- (C + t(C)) / 2
  ev <- eigen(C, symmetric = TRUE); keep <- ev$values > rel_tol * max(ev$values)
  L <- ev$vectors[, keep, drop = FALSE] %*% diag(sqrt(ev$values[keep]), sum(keep))
  list(A = A, lower = L, rank = sum(keep), min_eigen = min(ev$values),
       error = max(abs(C - L %*% t(L))), route = route)
}

# Gauss-Hermite nodes and weights (physicists' convention) by Golub-Welsch.
gauss_hermite <- function(n) {
  J <- matrix(0, n, n); off <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1, ]^2)
}
