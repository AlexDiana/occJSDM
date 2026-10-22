# Source the already verified pilot-math.R first. This extension avoids its
# pilot-specific dimensions and adds native trait coefficient extraction.
extract_bayesian <- function(fit, package) {
  if (package == "occJSDM") {
    j <- fit$results_output$jsdm_output
    species <- dim(j$B0_output)[1]; draws <- dim(j$B0_output)[2]; chains <- dim(j$B0_output)[3]
    slopes <- dim(j$B_output)[1]
    beta <- array(NA_real_, c(slopes + 1L, species, draws, chains))
    beta[1, , , ] <- j$B0_output
    beta[-1, , , ] <- j$B_output
    loading <- j$L_output
    for (chain in seq_len(chains)) for (draw in seq_len(draws)) {
      loading[, , draw, chain] <- loading[, , draw, chain] * j$sigmah_output[draw, chain]
    }
    stopifnot(dim(j$A_output)[2] == 0L,
      max(abs(fit$infos$list_X_psi_mat$mean_df)) < 1e-10,
      max(abs(fit$infos$list_X_psi_mat$sd_df - 1)) < 1e-10)
    trait <- j$G_output
    link <- "logit"
  } else {
    chains <- length(fit$postList); draws <- length(fit$postList[[1]])
    shape <- dim(fit$postList[[1]][[1]]$Beta)
    beta <- array(NA_real_, c(shape, draws, chains))
    loading <- array(NA_real_, c(2, shape[2], draws, chains))
    gamma_shape <- dim(fit$postList[[1]][[1]]$Gamma)
    trait <- array(NA_real_, c(gamma_shape, draws, chains))
    for (chain in seq_len(chains)) for (draw in seq_len(draws)) {
      d <- fit$postList[[chain]][[draw]]
      beta[, , draw, chain] <- d$Beta
      loading[, , draw, chain] <- d$Lambda[[1]]
      trait[, , draw, chain] <- d$Gamma
      stopifnot(all(d$sigma == 1), nrow(d$Lambda[[1]]) == 2)
    }
    link <- "probit"
  }
  list(beta = beta, loading = loading, trait = trait, link = link,
       n = draws, chains = chains)
}

posterior_marginal <- function(p, x, nodes = 31L, retain_draws = FALSE, draw_indices = NULL) {
  shape <- dim(p$beta)
  if (is.null(draw_indices)) draw_indices <- seq_len(shape[3])
  design <- cbind(1, as.matrix(x))
  if (retain_draws) output <- array(NA_real_, c(nrow(x), shape[2], length(draw_indices), shape[4]))
  else output <- matrix(NA_real_, nrow(x), shape[2])
  for (species in seq_len(shape[2])) {
    beta <- matrix(p$beta[, species, draw_indices, , drop = FALSE], nrow = shape[1])
    mu <- design %*% beta
    variance <- colSums(matrix(p$loading[, species, draw_indices, , drop = FALSE], nrow = 2)^2)
    probability <- if (p$link == "probit") pnorm(sweep(mu, 2, sqrt(1 + variance), "/")) else
      logistic_normal(mu, rep(sqrt(variance), each = nrow(mu)), nodes)
    if (retain_draws) output[, species, , ] <- probability else output[, species] <- rowMeans(probability)
  }
  output
}

extract_gllvm <- function(fit, input, trait_model) {
  if (!trait_model) return(point_parameters(fit, "gllvm"))
  env <- colnames(input$x)
  beta <- matrix(rep(fit$params$B[env], ncol(input$y)), nrow = length(env))
  if (!is.null(input$traits)) beta <- beta + t(as.matrix(fit$TR) %*% t(fit$fourth.corner))
  beta <- beta + fit$params$Br
  p <- list(beta = rbind(fit$params$beta0, beta),
            loading = t(sweep(fit$params$theta, 2, fit$params$sigma.lv, "*")))
  native <- predict(fit, newX = input$x, type = "link", level = 0)
  stopifnot(max(abs(cbind(1, as.matrix(input$x)) %*% p$beta - native)) < 1e-7)
  p
}

parameter_diagnostics <- function(p, input) {
  shape <- dim(p$beta); species <- shape[2]
  d <- mcmc_diagnostics(p$beta, paste0("beta_", seq_len(prod(shape[1:2]))), "coefficient", p$n, p$chains)
  pairs <- which(lower.tri(matrix(0, species, species), diag = TRUE), arr.ind = TRUE)
  covariance <- array(NA_real_, c(nrow(pairs), p$n, p$chains))
  for (i in seq_len(nrow(pairs))) {
    covariance[i, , ] <- apply(p$loading[, pairs[i, 1], , , drop = FALSE] *
      p$loading[, pairs[i, 2], , , drop = FALSE], 3:4, sum)
  }
  d <- rbind(d, mcmc_diagnostics(covariance, paste0("covariance_", seq_len(nrow(pairs))), "covariance", p$n, p$chains))
  grid <- input$x[c(1, 21, 41, 61, 81), , drop = FALSE]
  prediction <- posterior_marginal(p, grid, nodes = 61L, retain_draws = TRUE)
  coarser <- posterior_marginal(p, grid, nodes = 31L)
  integration_error <- max(abs(apply(prediction, 1:2, mean) - coarser))
  d <- rbind(d, mcmc_diagnostics(prediction, paste0("probability_", seq_len(5 * species)), "new_site_probability", p$n, p$chains))
  if (!is.null(input$traits)) {
    d <- rbind(d, mcmc_diagnostics(p$trait, paste0("trait_", seq_len(prod(dim(p$trait)[1:2]))), "trait", p$n, p$chains))
  }
  d$pass <- with(d, is.finite(rhat) & rhat <= 1.01 & bulk_ess >= 400 & tail_ess >= 400)
  list(table = d, passed = all(d$pass) && integration_error < .0001,
       integration_error = integration_error,
       max_rhat = max(d$rhat), min_bulk_ess = min(d$bulk_ess), min_tail_ess = min(d$tail_ess))
}

checked_posterior_mean <- function(p, x, grid, reference_grid, tolerance = .001) {
  # Choose posterior subsetting using a fixed training grid, before truth is read.
  # All global draws remain archived. Check the retained mean against the full
  # posterior mean already computed by the fit diagnostics.
  n <- dim(p$beta)[3]
  sizes <- unique(pmin(n, c(1000L, 2000L, 4000L, 8000L, n)))
  for (size in sizes) {
    idx <- unique(round(seq(1, n, length.out = size)))
    grid_mean <- posterior_marginal(p, grid, nodes = 61L, draw_indices = idx)
    change <- max(abs(grid_mean - reference_grid))
    if (change < tolerance || size == n) break
  }
  prediction <- posterior_marginal(p, x, nodes = 31L, draw_indices = idx)
  finer <- posterior_marginal(p, x, nodes = 61L, draw_indices = idx)
  integration_error <- max(abs(prediction - finer)); prediction <- finer; nodes <- 61L
  if (integration_error > .0001) {
    finer <- posterior_marginal(p, x, nodes = 121L, draw_indices = idx)
    integration_error <- max(abs(prediction - finer)); prediction <- finer; nodes <- 121L
  }
  stopifnot(integration_error <= .0001, change <= tolerance)
  list(mean = prediction, draws_per_chain = length(idx), nodes = nodes,
       subset_grid_change = change, integration_error = integration_error,
       uses_all_draws = length(idx) == n)
}
