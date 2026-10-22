# Accurate one-species scale profile. Rotate the two standard-normal factors
# so the focal species depends on the first factor only. Integrate its steep
# transition adaptively, retaining Gaussian quadrature for the other factor.
profile_likelihood <- function(parameters, x, y, species, scale, inner_nodes = 61L) {
  loading <- parameters$loading
  direction <- loading[, species] / sqrt(sum(loading[, species]^2))
  rotation <- rbind(direction, c(-direction[2], direction[1]))
  loading <- rotation %*% loading
  stopifnot(abs(loading[2, species]) < 1e-12)
  loading[2, species] <- 0
  beta <- parameters$beta
  beta[, species] <- beta[, species] * scale
  loading[, species] <- loading[, species] * scale
  mu <- cbind(1, as.matrix(x)) %*% beta
  rule <- normal_rule(inner_nodes)
  hidden_2 <- outer(rule$nodes, loading[2, ])
  values <- errors <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    threshold <- -mu[i, species] / loading[1, species]
    transition <- threshold + c(-20, -5, 0, 5, 20) / loading[1, species]
    split <- sort(unique(c(-9, transition[transition > -9 & transition < 9], 9)))
    integrand <- function(z) vapply(z, function(z1) {
      eta <- sweep(hidden_2, 2, mu[i, ] + z1*loading[1, ], "+")
      p <- .999999*plogis(eta) + .0000005
      likelihood <- exp(as.vector(log(p) %*% y[i, ] + log1p(-p) %*% (1-y[i, ])))
      sum(rule$weights * likelihood) * dnorm(z1)
    }, numeric(1))
    terms <- lapply(seq_len(length(split)-1L), function(k) {
      integrate(integrand, split[k], split[k+1L], rel.tol = 1e-8,
                abs.tol = 1e-12, subdivisions = 500L)
    })
    values[i] <- sum(vapply(terms, `[[`, numeric(1), "value"))
    errors[i] <- sum(vapply(terms, `[[`, numeric(1), "abs.error"))
  }
  list(loglik = sum(log(values)), site_likelihood = values,
       estimated_loglik_error = sum(errors / values),
       parameters = list(beta = beta, loading = loading), inner_nodes = inner_nodes)
}
