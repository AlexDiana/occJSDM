# Development helpers for the new-site teaching example, not a package API.
# Standard-normal Gaussian quadrature, obtained from the Jacobi eigenproblem.
normal_rule <- function(n) {
  stopifnot(n >= 3L)
  jacobi <- matrix(0, n, n)
  for (k in seq_len(n - 1L)) jacobi[k, k + 1L] <- jacobi[k + 1L, k] <- sqrt(k)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1, ]^2)
}

normal_logistic_mean <- function(location, spread, rule = normal_rule(31L)) {
  stopifnot(length(spread) %in% c(1L, length(location)), all(spread >= 0))
  result <- location * 0
  for (k in seq_along(rule$nodes)) {
    result <- result + rule$weights[k] * plogis(location + spread * rule$nodes[k])
  }
  result
}

normal_logistic_matrix <- function(location, spread, rule) {
  stopifnot(is.matrix(location), length(spread) == ncol(location), all(spread >= 0))
  result <- location * 0
  for (k in seq_along(rule$nodes)) {
    result <- result + rule$weights[k] *
      plogis(sweep(location, 2, spread * rule$nodes[k], "+"))
  }
  result
}

new_site_draw_parameters <- function(fit, species) {
  j <- fit$results_output$jsdm_output
  d <- fit$infos$n_factors
  loadings <- matrix(j$L_output[, species, , ], nrow = d)
  list(
    intercept = as.vector(j$B0_output[species, , ]),
    slopes = matrix(j$B_output[, species, , ], nrow = ncol(fit$X_psi)),
    spread = as.vector(j$sigmah_output) * sqrt(colSums(loadings^2)),
    iterations = dim(j$B0_output)[2], chains = dim(j$B0_output)[3]
  )
}

# Select accuracy from a declared grid covering each fit's residual SD tail.
# This changes only numerical integration accuracy, never data or fitted draws.
checked_normal_rule <- function(max_spread) {
  grid <- expand.grid(location = seq(-30, 30, by = .5),
                      spread = seq(0, max_spread, length.out = 9))
  for (nodes in c(31L, 61L, 121L)) {
    rule <- normal_rule(nodes)
    reference <- normal_rule(nodes * 2L - 1L)
    error <- max(abs(normal_logistic_mean(grid$location, grid$spread, rule) -
                     normal_logistic_mean(grid$location, grid$spread, reference)))
    if (error < 1e-5) return(list(rule = rule, nodes = nodes,
      grid_max_difference = error, max_spread = max_spread))
  }
  stop("Quadrature did not meet its predeclared accuracy check")
}

predict_marginal_teaching <- function(fit, raw_covariates, label) {
  scaling <- fit$infos$list_X_psi_mat
  stopifnot(identical(colnames(raw_covariates), scaling$names_df), fit$infos$ps == 0L)
  x <- sweep(sweep(raw_covariates, 2, scaling$mean_df, "-"), 2, scaling$sd_df, "/")
  parameters <- lapply(seq_along(fit$infos$speciesNames), function(s) new_site_draw_parameters(fit, s))
  accuracy <- checked_normal_rule(max(vapply(parameters, function(p) max(p$spread), numeric(1))))
  tables <- list()
  diagnostic_rows <- list()
  for (s in seq_along(parameters)) {
    p <- parameters[[s]]
    community_mean <- numeric(length(p$intercept))
    rows <- list()
    blocks <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / 50))
    for (b in seq_along(blocks)) {
      ii <- blocks[[b]]
      eta <- sweep(x[ii, , drop = FALSE] %*% p$slopes, 2, p$intercept, "+")
      probabilities <- normal_logistic_matrix(eta, p$spread, accuracy$rule)
      intervals <- t(apply(probabilities, 1, quantile, probs = c(.025, .975), names = FALSE))
      rows[[b]] <- data.frame(
        arm = label, Site = rownames(x)[ii], species = fit$infos$speciesNames[s],
        estimate = rowMeans(probabilities), lower = intervals[, 1], upper = intervals[, 2])
      community_mean <- community_mean + colSums(probabilities) / nrow(x)
    }
    tables[[s]] <- do.call(rbind, rows)
    draws <- matrix(community_mean, p$iterations, p$chains)
    diagnostic_rows[[s]] <- data.frame(arm = label, species = fit$infos$speciesNames[s],
      rhat = posterior::rhat(draws), ess = posterior::ess_mean(draws),
      mcse = posterior::mcse_mean(draws))
    cat(label, fit$infos$speciesNames[s], "all posterior draws integrated\n")
  }
  list(cells = do.call(rbind, tables), diagnostics = do.call(rbind, diagnostic_rows),
       quadrature = accuracy[c("nodes", "grid_max_difference", "max_spread")])
}
