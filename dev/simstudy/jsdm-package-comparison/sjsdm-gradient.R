# Deterministic diagnostic of sjSDM's unpenalised, two-factor logit likelihood.
# Parameter order: R-column-major beta (predictors by species), then loading.
# Uses sjSDM's exact probability stabilisation, not a modified statistical model.
make_joint_gradient <- function(x, y, nodes = 61L) {
  design <- cbind(1, as.matrix(x))
  species <- ncol(y)
  covariates <- ncol(design)
  rule <- normal_rule(nodes)
  grid <- as.matrix(expand.grid(rule$nodes, rule$nodes))
  weight <- as.vector(outer(rule$weights, rule$weights))
  keep <- weight > 0
  grid <- grid[keep, , drop = FALSE]
  log_weight <- log(weight[keep])
  function(theta) {
    stopifnot(length(theta) == (covariates + 2L) * species)
    beta <- matrix(theta[seq_len(covariates * species)], covariates, species)
    loading <- matrix(theta[-seq_len(covariates * species)], 2, species)
    mu <- design %*% beta
    hidden <- grid %*% loading
    gradient_mu <- mu * 0
    gradient_loading <- loading * 0
    loglik <- 0
    for (i in seq_len(nrow(y))) {
      eta <- sweep(hidden, 2, mu[i, ], "+")
      raw <- plogis(eta)
      probability <- 0.999999 * raw + 0.0000005
      site_log_weight <- as.vector(log(probability) %*% y[i, ] +
                                     log1p(-probability) %*% (1-y[i, ])) + log_weight
      maximum <- max(site_log_weight)
      posterior <- exp(site_log_weight - maximum)
      loglik <- loglik + maximum + log(sum(posterior))
      posterior <- posterior / sum(posterior)
      derivative <- (sweep(1/probability, 2, y[i, ], "*") -
                       sweep(1/(1-probability), 2, 1-y[i, ], "*")) *
        (0.999999 * raw * (1-raw))
      weighted_derivative <- derivative * posterior
      gradient_mu[i, ] <- colSums(weighted_derivative)
      gradient_loading <- gradient_loading + crossprod(grid, weighted_derivative)
    }
    list(value = -loglik,
         gradient = -c(crossprod(design, gradient_mu), gradient_loading))
  }
}
