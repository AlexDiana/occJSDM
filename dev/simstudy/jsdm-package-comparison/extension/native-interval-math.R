# Study-only native Wald intervals. Coefficient blocks are on the fitted scale.
native_coefficient_intervals <- function(beta, covariance, centre, spread) {
  k <- nrow(beta); species <- ncol(beta)
  stopifnot(length(covariance) == species, k == length(spread) + 1L)
  transform <- diag(k)
  transform[-1, -1] <- diag(1/spread, length(spread))
  transform[1, -1] <- -centre/spread
  estimate <- as.vector(raw_coefficients(beta, centre, spread))
  se <- rep(NA_real_, length(estimate))
  valid <- logical(species)
  for (s in seq_len(species)) {
    v <- covariance[[s]]
    if (!is.matrix(v) || !identical(dim(v), c(k, k)) || any(!is.finite(v))) next
    if (max(abs(v-t(v))) > 1e-8 * max(1, max(abs(v)))) next
    if (min(eigen(v, symmetric = TRUE, only.values = TRUE)$values) <= 0) next
    raw <- transform %*% v %*% t(transform)
    if (any(diag(raw) <= 0)) next
    se[(s-1L)*k + seq_len(k)] <- sqrt(diag(raw))
    valid[s] <- TRUE
  }
  data.frame(estimate = estimate, se = se,
    lower = estimate - qnorm(.975)*se, upper = estimate + qnorm(.975)*se,
    covariance_valid = rep(valid, each = k))
}
