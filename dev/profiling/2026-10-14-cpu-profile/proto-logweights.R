# Prototype of spatial_range_logweights() with (1) the symmetric Gram product
# crossprod(sqrt(w) * Z) in place of crossprod(Z, w * Z), (2) the prior
# precision matrix, sqrt(Omega) and the fixed design block hoisted out of the
# grid and species loops, (3) backsolve(C, h, transpose = TRUE) in place of
# forwardsolve(t(C), h), and (4) no scatter into H when the basis is full.
# Same algebra, same grouped dispatch; the first three calls are checked
# against the original and the largest discrepancy is printed.
ns <- asNamespace("occJSDM")
orig <- get("spatial_range_logweights", envir = ns)
proto <- function(X, U, M_B, M_Bs, sigma_b, sigma_bs, kappa, Omega, Xs_centers,
                  list_SoRSummaries, a_l_s, b_l_s, location = NULL, sigma_b0 = 1) {
  n <- nrow(X); p <- ncol(X); d <- ncol(U); ps <- nrow(M_Bs); S <- ncol(Omega)
  prior_precision <- c(1/sigma_b0^2, rep(1/sigma_b^2, p), rep(1, d), rep(1/sigma_bs^2, ps))
  prior_means <- rbind(rep(0, S), M_B, matrix(0, d, S), M_Bs)
  prior_linear <- prior_precision * prior_means
  constant_omega <- vapply(seq_len(S), function(s) all(Omega[, s] == Omega[1, s]), logical(1))
  if (!is.null(location)) {
    if (length(location) != n || anyNA(location))
      stop("location must contain one non-missing group label per site row")
    if (anyDuplicated(location) && !all(constant_omega)) {
      return(spatial_range_logweights_grouped(X, U, M_B, M_Bs, sigma_b, sigma_bs,
        kappa, Omega, Xs_centers, list_SoRSummaries, a_l_s, b_l_s, location, sigma_b0))
    }
  }
  grid <- list_SoRSummaries$l_s_grid
  P <- diag(prior_precision, length(prior_precision))
  sqrt_omega <- sqrt(Omega)
  A <- cbind(1, X, U)
  full_basis <- ncol(Xs_centers) == ps && all(Xs_centers == rep(seq_len(ps), each = n))
  vapply(seq_along(grid), function(j) {
    Ks <- matrix(list_SoRSummaries$Ks_all[, , j], nrow = n)
    if (full_basis) {
      H <- Ks
    } else {
      H <- matrix(0, n, ps)
      H[cbind(rep(seq_len(n), ncol(Xs_centers)), as.vector(Xs_centers))] <- as.vector(Ks)
    }
    Z <- cbind(A, H)
    gram <- if (any(constant_omega)) crossprod(Z) else NULL
    linear <- crossprod(Z, kappa) + prior_linear
    score <- vapply(seq_len(S), function(s) {
      Q <- if (constant_omega[s]) Omega[1, s] * gram else crossprod(sqrt_omega[, s] * Z)
      C <- chol(Q + P)
      v <- backsolve(C, linear[, s], transpose = TRUE)
      .5 * sum(v^2) - sum(log(diag(C)))
    }, numeric(1))
    sum(score) + dgamma(grid[j], a_l_s, b_l_s, log = TRUE)
  }, numeric(1))
}
environment(proto) <- ns
calls <- 0L
checked <- function(...) {
  calls <<- calls + 1L
  new <- proto(...)
  if (calls <= 3L) {
    old <- orig(...)
    cat(sprintf("proto check %d: max |diff| = %.3e, max |rel diff| = %.3e, same argmax = %s\n",
                calls, max(abs(new - old)), max(abs(new - old) / pmax(abs(old), 1)),
                identical(which.max(new), which.max(old))))
  }
  new
}
environment(checked) <- ns
if (bindingIsLocked("spatial_range_logweights", ns)) unlockBinding("spatial_range_logweights", ns)
assign("spatial_range_logweights", checked, envir = ns)
cat("prototype range step installed\n")
