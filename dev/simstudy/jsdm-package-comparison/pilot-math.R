# Study-only probability integration and extraction, not package API changes.
normal_rule <- function(n) {
  j <- matrix(0, n, n)
  for (k in seq_len(n - 1L)) j[k, k + 1L] <- j[k + 1L, k] <- sqrt(k)
  e <- eigen(j, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1, ]^2)
}

logistic_normal <- function(mu, sd, nodes = 31L) {
  rule <- normal_rule(nodes)
  out <- mu * 0
  for (k in seq_along(rule$nodes)) {
    out <- out + rule$weights[k] * plogis(mu + sd * rule$nodes[k])
  }
  out
}

# Use normal-density integration at ordinary SDs, where the probability-domain
# identity can miss rare-event mass at its endpoints. At very large SDs that
# identity is smoother than the near-step in the normal-density integrand.
logistic_normal_adaptive <- function(mu, sd, tolerance = 1e-9) {
  stopifnot(length(sd) %in% c(1L, length(mu)))
  sd <- rep_len(sd, length(mu))
  vapply(seq_along(mu), function(i) {
    if (sd[i] < 0.05) return(logistic_normal(mu[i], sd[i], 31L))
    if (sd[i] >= 10) {
      return(integrate(function(u) pnorm((mu[i] - qlogis(u)) / sd[i]),
                       lower = 0, upper = 1, rel.tol = tolerance,
                       abs.tol = tolerance, subdivisions = 2000L)$value)
    }
    integrate(function(z) plogis(mu[i] + sd[i] * z) * dnorm(z),
              lower = -Inf, upper = Inf, rel.tol = tolerance,
              abs.tol = tolerance, subdivisions = 2000L)$value
  }, numeric(1))
}

fixed_grid <- function() {
  data.frame(environment_1 = c(0, -1, 1, 0, 0),
             environment_2 = c(0, 0, 0, -1, 1))
}

point_parameters <- function(fit, package) {
  if (package == "gllvm") {
    list(beta = rbind(fit$params$beta0, t(fit$params$Xcoef)),
         loading = t(sweep(fit$params$theta, 2, fit$params$sigma.lv, "*")))
  } else {
    list(beta = t(fit$weights[[1]]), loading = t(fit$sigma))
  }
}

point_marginal <- function(parameters, x) {
  mu <- cbind(1, as.matrix(x)) %*% parameters$beta
  sd <- sqrt(colSums(parameters$loading^2))
  answer <- mu * 0
  for (s in seq_len(ncol(mu))) {
    answer[, s] <- logistic_normal_adaptive(mu[, s], sd[s])
  }
  answer
}

# Integrate the whole site's Bernoulli likelihood jointly over two factors.
# This conditions on training data only, never on known generating truth.
joint_integration <- function(parameters, x, y, nodes = 61L, clamp = FALSE) {
  rule <- normal_rule(nodes)
  grid <- as.matrix(expand.grid(rule$nodes, rule$nodes))
  weights <- as.vector(outer(rule$weights, rule$weights))
  keep <- weights > 0
  grid <- grid[keep, , drop = FALSE]
  weights <- weights[keep]
  hidden <- grid %*% parameters$loading
  mu <- cbind(1, as.matrix(x)) %*% parameters$beta
  conditional <- mu * 0
  site_loglik <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    eta <- sweep(hidden, 2, mu[i, ], "+")
    probabilities <- plogis(eta)
    if (clamp) probabilities <- probabilities * 0.999999 + 0.0000005
    if (clamp) {
      log_probability <- log(probabilities)
      log_absence <- log1p(-probabilities)
    } else {
      log_probability <- plogis(eta, log.p = TRUE)
      log_absence <- plogis(-eta, log.p = TRUE)
    }
    log_weight <- as.vector(log_probability %*% y[i, ] +
                            log_absence %*% (1 - y[i, ])) + log(weights)
    maximum <- max(log_weight)
    posterior_weight <- exp(log_weight - maximum)
    site_loglik[i] <- maximum + log(sum(posterior_weight))
    posterior_weight <- posterior_weight / sum(posterior_weight)
    conditional[i, ] <- colSums(probabilities * posterior_weight)
  }
  list(loglik = sum(site_loglik), site_loglik = site_loglik,
       conditional = conditional, nodes = nodes)
}

bayesian_parameters <- function(fit, package) {
  if (package == "occJSDM") {
    j <- fit$results_output$jsdm_output
    n <- dim(j$B0_output)[2]
    chains <- dim(j$B0_output)[3]
    beta <- array(NA_real_, c(3, 10, n, chains))
    beta[1, , , ] <- j$B0_output
    beta[2:3, , , ] <- j$B_output
    loading <- j$L_output
    for (c in seq_len(chains)) for (i in seq_len(n)) {
      loading[, , i, c] <- loading[, , i, c] * j$sigmah_output[i, c]
    }
    score <- j$U_output
    for (c in seq_len(chains)) for (i in seq_len(n)) {
      score[, , i, c] <- score[, , i, c] / j$sigmah_output[i, c]
    }
    stopifnot(dim(j$A_output)[2] == 0L,
              max(abs(fit$infos$list_X_psi_mat$mean_df)) < 1e-12,
              max(abs(fit$infos$list_X_psi_mat$sd_df - 1)) < 1e-12)
    link <- "logit"
  } else {
    chains <- length(fit$postList)
    n <- length(fit$postList[[1]])
    beta <- array(NA_real_, c(3, 10, n, chains))
    loading <- array(NA_real_, c(2, 10, n, chains))
    score <- array(NA_real_, c(100, 2, n, chains))
    for (c in seq_len(chains)) for (i in seq_len(n)) {
      d <- fit$postList[[c]][[i]]
      beta[, , i, c] <- d$Beta
      loading[, , i, c] <- d$Lambda[[1]]
      score[, , i, c] <- d$Eta[[1]][fit$Pi[, 1], ]
      stopifnot(all(d$sigma == 1), nrow(d$Lambda[[1]]) == 2)
    }
    stopifnot(max(abs(fit$XScalePar[1, ])) == 0,
              all(fit$XScalePar[2, ] == 1))
    link <- "probit"
  }
  list(beta = beta, loading = loading, score = score,
       n = n, chains = chains, link = link)
}

bayesian_marginal <- function(p, x, nodes = 31L) {
  result <- array(NA_real_, c(nrow(x), 10, p$n, p$chains))
  for (s in 1:10) {
    beta <- matrix(p$beta[, s, , ], nrow = 3)
    mu <- cbind(1, as.matrix(x)) %*% beta
    variance <- colSums(matrix(p$loading[, s, , ], nrow = 2)^2)
    result[, s, , ] <- if (p$link == "probit") {
      pnorm(sweep(mu, 2, sqrt(1 + variance), "/"))
    } else {
      logistic_normal(mu, rep(sqrt(variance), each = nrow(mu)), nodes)
    }
  }
  result
}

bayesian_training <- function(p, x) {
  result <- array(NA_real_, c(nrow(x), 10, p$n, p$chains))
  inverse_link <- if (p$link == "probit") pnorm else plogis
  for (c in seq_len(p$chains)) for (i in seq_len(p$n)) {
    result[, , i, c] <- inverse_link(cbind(1, as.matrix(x)) %*% p$beta[, , i, c] +
                                     p$score[, , i, c] %*% p$loading[, , i, c])
  }
  result
}

mcmc_diagnostics <- function(draws, labels, block, n, chains) {
  draws <- array(draws, c(length(labels), n, chains))
  rows <- lapply(seq_along(labels), function(k) {
    m <- matrix(draws[k, , ], n, chains)
    data.frame(block = block, variable = labels[k], mean = mean(m),
               rhat = posterior::rhat(m), bulk_ess = posterior::ess_bulk(m),
               tail_ess = posterior::ess_tail(m), mcse = posterior::mcse_mean(m))
  })
  do.call(rbind, rows)
}
