# Independent ecological simulator. Fitters receive make_input() only.
# Fixed seeds and coefficients are declared before any fit or truth score.
simulate_community <- function(replicate, scenario) {
  stopifnot(replicate %in% 1:10, scenario %in% c("baseline", "rare", "correlated", "curved", "traits"))
  set.seed(26093000L + replicate)
  n <- 600L
  species <- if (scenario == "traits") 30L else 10L
  ids <- c(sprintf("train_%03d", 1:300), sprintf("test_%03d", 1:300))
  raw_x <- matrix(rnorm(n * 2), n, 2, dimnames = list(ids, c("environment_1", "environment_2")))
  hidden <- matrix(rnorm(n * 2), n, 2)
  if (scenario == "correlated") {
    raw_x[, 2] <- .85 * raw_x[, 1] + sqrt(1 - .85^2) * raw_x[, 2]
  }
  # Generate all 30 species first so the trait experiment has nested species sets.
  all_traits <- matrix(rnorm(60), 30, 2,
    dimnames = list(sprintf("species_%02d", 1:30), c("drought_tolerance", "irrelevant_trait")))
  all_intercepts <- rnorm(30, 0, .9)
  all_deviations <- matrix(rnorm(60, 0, .35), 30, 2)
  all_loadings <- matrix(rnorm(60, 0, .6), 2, 30)
  trait_effect <- matrix(c(1, 0, -.6, 0), 2, 2,
    dimnames = list(colnames(all_traits), colnames(raw_x)))
  if (scenario == "traits") {
    beta <- rbind(all_intercepts, t(all_traits %*% trait_effect + all_deviations))
    loading <- all_loadings
  } else {
    beta <- rbind(seq(-1.8, 1.8, length.out = 10),
      c(1.2, .8, -1, -.7, 1, -1.1, .6, -.8, 1, -1.2),
      c(-.6, .9, .7, -1, .4, .8, -.9, .6, -.7, .5))
    # Communities vary independently, while scenarios within a replicate are paired.
    beta <- beta + matrix(rnorm(30, 0, .15), 3, 10)
    loading <- rbind(c(.9, .7, -.8, -.6, .1, .2, .8, -.9, .5, -.6),
                     c(.2, -.3, .2, -.4, .9, -.8, .6, -.5, -.8, .7))
    loading <- loading + matrix(rnorm(20, 0, .1), 2, 10)
    if (scenario == "rare") beta[1, 1:3] <- -4
  }
  species_ids <- sprintf("species_%02d", seq_len(species))
  dimnames(beta) <- list(c("intercept", colnames(raw_x)), species_ids)
  colnames(loading) <- species_ids
  curvature <- rep(0, species)
  if (scenario == "curved") curvature[seq(1, 9, 2)] <- -.9
  eta <- cbind(1, raw_x) %*% beta + outer(raw_x[, 1]^2, curvature) + hidden %*% loading
  probability <- plogis(eta)
  uniforms <- matrix(runif(n * species), n, species, dimnames = list(ids, species_ids))
  y <- (uniforms < probability) * 1
  list(replicate = replicate, scenario = scenario, raw_x = raw_x, hidden = hidden,
       beta = beta, loading = loading, curvature = curvature,
       traits = as.data.frame(all_traits[seq_len(species), , drop = FALSE]),
       trait_effect = trait_effect, probability = probability, uniforms = uniforms, y = y)
}

basis <- function(raw_x, response) {
  x <- as.data.frame(raw_x)
  if (response == "quadratic") x$environment_1_squared <- x$environment_1^2
  x
}

make_input <- function(community, n_sites, n_species, response, use_traits) {
  stopifnot(n_sites %in% c(100, 300), n_species <= ncol(community$y),
            response %in% c("linear", "quadratic"))
  raw <- basis(community$raw_x, response)
  centre <- colMeans(raw[seq_len(n_sites), , drop = FALSE])
  spread <- apply(raw[seq_len(n_sites), , drop = FALSE], 2, sd)
  scaled <- as.data.frame(sweep(sweep(as.matrix(raw), 2, centre, "-"), 2, spread, "/"))
  list(x = scaled[seq_len(n_sites), , drop = FALSE],
       y = community$y[seq_len(n_sites), seq_len(n_species), drop = FALSE],
       test_x = scaled[301:600, , drop = FALSE],
       traits = if (use_traits) community$traits[seq_len(n_species), , drop = FALSE] else NULL,
       centre = centre, spread = spread, response = response)
}

scaled_truth <- function(community, input) {
  beta <- community$beta[, colnames(input$y), drop = FALSE]
  if (input$response == "quadratic") beta <- rbind(beta, community$curvature[seq_len(ncol(beta))])
  beta[1, ] <- beta[1, ] + input$centre %*% beta[-1, , drop = FALSE]
  beta[-1, ] <- beta[-1, , drop = FALSE] * input$spread
  beta
}

marginal_truth <- function(community, raw_x) {
  mu <- cbind(1, as.matrix(raw_x)) %*% community$beta + outer(raw_x[, 1]^2, community$curvature)
  sd <- sqrt(colSums(community$loading^2))
  out <- mu * 0
  for (s in seq_len(ncol(out))) {
    out[, s] <- vapply(mu[, s], function(m) {
      if (sd[s] == 0) return(plogis(m))
      integrate(function(z) plogis(m + sd[s] * z) * dnorm(z), -Inf, Inf,
                rel.tol = 1e-9)$value
    }, numeric(1))
  }
  out
}

error_summary <- function(truth, estimate) {
  stopifnot(length(truth) == length(estimate), all(is.finite(estimate)))
  data.frame(n = length(truth), bias_pp = 100 * mean(estimate - truth),
             mae_pp = 100 * mean(abs(estimate - truth)))
}
