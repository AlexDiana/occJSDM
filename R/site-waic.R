# Observed-data likelihood for a complete survey at an independent site.
# Helpers are intentionally separate from the MCMC sampler.
site_waic_logadd <- function(a, b) {
  high <- pmax(a, b)
  value <- high + log1p(exp(-abs(a - b)))
  value[is.infinite(high) & high < 0] <- -Inf
  value
}

site_waic_log_sigmoid <- function(x) {
  -pmax(-x, 0) - log1p(exp(-abs(x)))
}

site_waic_rule <- function(order, factors, max_nodes) {
  if (order^factors > max_nodes) {
    stop("Joint quadrature exceeds max_nodes; reduce the factor count or use held-out-site validation.", call. = FALSE)
  }
  jacobi <- matrix(0, order, order)
  for (k in seq_len(order - 1L)) {
    jacobi[k, k + 1L] <- jacobi[k + 1L, k] <- sqrt(k)
  }
  decomposition <- eigen(jacobi, symmetric = TRUE)
  weights <- decomposition$vectors[1, ]^2
  # Some extremely small tail weights round to zero in the eigensolver.
  keep <- which(weights > 0)
  grid <- as.matrix(expand.grid(rep(list(keep), factors)))
  list(nodes = matrix(decomposition$values[grid], nrow(grid), factors),
       log_weights = rowSums(matrix(log(weights[grid]), nrow(grid), factors)))
}

site_waic_array <- function(x, leading, iterations, chains, name) {
  expected <- c(leading, iterations, chains)
  if (!identical(as.integer(dim(x)), as.integer(expected)) || any(!is.finite(x))) {
    stop("Missing, nonfinite or incompatible posterior array: ", name, ".", call. = FALSE)
  }
  array(x, c(leading, iterations * chains))
}

site_waic_index <- function(x, size, upper, name) {
  if (length(x) != size || anyNA(x) || any(!is.finite(x)) ||
      any(x != floor(x) | x < 1 | x > upper)) {
    stop("Invalid fitted observation mapping: ", name, ".", call. = FALSE)
  }
  as.integer(x)
}

#' Observed-data WAIC for independent sites
#'
#' Compare non-spatial models of the same community observations, with a whole
#' site as the prediction unit. This is post-processing: it does not refit the
#' model or change occupancy estimates.
#'
#' @param fitModel Output from [runOccJSDM()]. Only non-spatial `binary`,
#'   `occupancy` and `two_stage` models are currently supported.
#' @param threshold Detection threshold used when fitting. New fits save it;
#'   older occupancy/two-stage fits require it explicitly. A different threshold
#'   cannot be used to rescore a new fit.
#' @param draws Optional unique indices of retained posterior draws, ordered
#'   by iteration within chain. By default use all draws. At least two are needed.
#'   Subsetting is for computational checks; assess MCMC precision before using
#'   a subset to select a model.
#' @param quadrature Increasing vector of Gaussian quadrature orders per factor.
#'   Each draw is checked against the next order until the largest absolute
#'   change in any site's log likelihood is below `tolerance`. At least two
#'   orders are required. Increase this sequence to check numerical stability.
#'   The default node budget permits checking up to three factors; difficult
#'   three-factor fits may exhaust it when a finer grid is needed.
#' @param tolerance Positive tolerance for successive site log likelihoods.
#'   Agreement of grids is a numerical diagnostic, not a rigorous error bound.
#' @param max_nodes Maximum tensor-grid nodes allowed per integration rule.
#'   Cost increases rapidly with the number of factors. If integration fails,
#'   the function stops instead of returning an unchecked score.
#'
#' @details
#' Occupancy and collection states are summed out exactly. Shared Gaussian
#' site factors are integrated jointly, after multiplying conditional species
#' likelihoods, using the fitted `sigmah_output` as their standard deviation.
#' Fitted site scores and sampled latent states are never treated as data.
#' Missing PCRs contribute no likelihood term. Sites with no observed responses
#' are omitted and named in the result.
#'
#' For each site, WAIC is `-2 * (log posterior mean likelihood - variance
#' of log likelihood)`. Lower WAIC indicates better estimated predictive accuracy
#' for the same observations and sampling design. Pointwise penalties above 0.4
#' trigger a warning: WAIC can then be unreliable, and actual held-out-site
#' validation is preferable. A corrected likelihood does not ensure a reliable
#' factor-count choice or recovery of the generating count.
#'
#' The reported standard error describes variation across independent sites,
#' not MCMC or integration error. Check chain convergence and repeat with more
#' posterior draws and finer quadrature before interpreting small differences.
#' Spatial models require a separate prediction target and are not scored here.
#'
#' @return A list of class `occJSDM_site_waic` with `WAIC`, `SE`, `elpd_waic`,
#'   `p_waic`, a `pointwise` table, a draws-by-site `log_lik` matrix, `draw_ids`,
#'   `integration` diagnostics, `omitted_sites`, and `observations` for checking
#'   that comparisons use identical data and grouping.
#' @seealso [extractWAIC()], [compareSiteWAIC()]
#' @references
#' Merkle, Furr and Rabe-Hesketh (2019), \doi{10.1007/s11336-019-09679-0}.
#' Vehtari, Gelman and Gabry (2017), \doi{10.1007/s11222-016-9696-4}.
#' @examples
#' \dontrun{
#' score <- computeSiteWAIC(fitModel)
#' score$pointwise
#' # For an old saved fit originally fitted at threshold = 1:
#' score <- computeSiteWAIC(oldFit, threshold = 1)
#' }
#' @md
#' @export
computeSiteWAIC <- function(fitModel, threshold = NULL, draws = NULL,
                            quadrature = c(15L, 25L, 41L, 65L),
                            tolerance = 1e-4, max_nodes = 200000L) {
  info <- fitModel$infos
  if (length(info$ps) != 1L || is.na(info$ps) || info$ps != 0) {
    stop("Site WAIC currently supports non-spatial fits only.", call. = FALSE)
  }
  if (!info$model %in% c("binary", "occupancy", "two_stage")) {
    stop("Site WAIC supports binary, occupancy and two_stage models only.", call. = FALSE)
  }
  if (!is.numeric(quadrature) || length(quadrature) < 2L ||
      any(!is.finite(quadrature)) || any(quadrature < 2 | quadrature != floor(quadrature)) ||
      any(diff(quadrature) <= 0)) {
    stop("quadrature must contain at least two increasing integer orders >= 2.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      length(max_nodes) != 1L || !is.finite(max_nodes) || max_nodes < 1) {
    stop("tolerance and max_nodes must be positive finite numbers.", call. = FALSE)
  }
  y <- as.matrix(info$OTU)
  n <- info$n
  S <- info$S
  factors <- info$n_factors
  if (!is.numeric(y) || ncol(y) != S || any(!is.finite(y) & !is.na(y))) {
    stop("Invalid observed response matrix.", call. = FALSE)
  }
  if (length(factors) != 1L || !is.finite(factors) || factors < 0 || factors != floor(factors)) {
    stop("Invalid fitted number of factors.", call. = FALSE)
  }
  if (info$model != "binary") {
    saved_threshold <- info$threshold
    if (is.null(threshold)) threshold <- saved_threshold
    if (is.null(threshold)) {
      stop("This old fit did not save its threshold; supply the threshold used to fit it.", call. = FALSE)
    }
    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
      stop("threshold must be a positive finite number.", call. = FALSE)
    }
    if (!is.null(saved_threshold) && threshold != saved_threshold) {
      stop("threshold must match the threshold used to fit the model.", call. = FALSE)
    }
    y <- 1 * (y >= threshold)
  } else if (any(!is.na(y) & !y %in% c(0, 1))) {
    stop("Binary observations must be zero or one.", call. = FALSE)
  }
  x <- as.matrix(fitModel$X_psi)
  if (nrow(x) != n || any(!is.finite(x))) stop("Invalid fitted X_psi.", call. = FALSE)
  output <- fitModel$results_output
  jsdm <- output$jsdm_output
  dims <- dim(jsdm$B0_output)
  if (length(dims) != 3L) stop("Missing posterior intercept draws.", call. = FALSE)
  iterations <- dims[2]
  chains <- dims[3]
  total <- iterations * chains
  if (is.null(draws)) draws <- seq_len(total)
  if (length(draws) < 2L) stop("Use at least two posterior draws.", call. = FALSE)
  if (!is.numeric(draws) || any(!is.finite(draws)) || any(draws != floor(draws)) ||
      any(draws < 1 | draws > total) || anyDuplicated(draws)) {
    stop("draws must be unique retained-draw indices.", call. = FALSE)
  }
  get_array <- function(a, leading, name) site_waic_array(a, leading, iterations, chains, name)
  intercept <- get_array(jsdm$B0_output, S, "B0_output")
  slopes <- get_array(jsdm$B_output, c(ncol(x), S), "B_output")
  loadings <- get_array(jsdm$L_output, c(factors, S), "L_output")
  sigma <- get_array(jsdm$sigmah_output, integer(), "sigmah_output")
  if (any(sigma < 0)) stop("Negative factor standard deviation.", call. = FALSE)

  idx <- info$list_idx
  if (info$model == "binary") {
    if (nrow(y) != n) stop("Binary observations must have one row per site.", call. = FALSE)
    site_index <- seq_len(n)
  } else {
    xt <- as.matrix(fitModel$X_theta)
    samples <- nrow(xt)
    if (any(!is.finite(xt))) stop("Invalid fitted X_theta.", call. = FALSE)
    sample_site <- site_waic_index(idx$idx_z_w, samples, n, "idx_z_w")
    if (!identical(sort(unique(sample_site)), seq_len(n))) stop("Missing site mapping.", call. = FALSE)
    beta_theta <- get_array(output$beta_theta_output, c(ncol(xt), S), "beta_theta_output")
    theta0 <- get_array(output$theta0_output, S, "theta0_output")
    if (any(theta0 < 0 | theta0 > 1)) stop("Invalid theta0 probabilities.", call. = FALSE)
    if (info$model == "two_stage") {
      primers <- length(info$primerNames)
      sample_index <- site_waic_index(idx$idx_w_k, nrow(y), samples, "idx_w_k")
      primer_index <- site_waic_index(idx$idx_p_k, nrow(y), primers, "idx_p_k")
      if (!identical(sort(unique(sample_index)), seq_len(samples))) stop("Missing sample mapping.", call. = FALSE)
      p <- get_array(output$p_output, c(primers, S), "p_output")
      q <- get_array(output$q_output, c(primers, S), "q_output")
      if (any(p < 0 | p > 1) || any(q < 0 | q > 1)) stop("Invalid PCR probabilities.", call. = FALSE)
      site_index <- sample_site[sample_index]
    } else {
      if (nrow(y) != samples) stop("Occupancy observations must match samples.", call. = FALSE)
      sample_index <- seq_len(samples)
      primer_index <- NULL
      site_index <- sample_site
    }
  }
  observed <- rowsum(rowSums(!is.na(y)), site_index, reorder = TRUE)[, 1] > 0
  if (!any(observed)) stop("No observed responses to score.", call. = FALSE)
  # Binary fits historically label siteNames as 1:n even when data_info
  # records real site identities. Preserve those identities for paired scores.
  site_names <- if (info$model == "binary" && !is.null(info$data_info$Site)) {
    as.character(info$data_info$Site)
  } else {
    as.character(info$siteNames)
  }
  if (length(site_names) != n || anyNA(site_names) || anyDuplicated(site_names)) {
    stop("Site names must uniquely identify the fitted sites.", call. = FALSE)
  }
  log_lik <- matrix(NA_real_, length(draws), sum(observed), dimnames = list(NULL, site_names[observed]))
  orders <- integer(length(draws))
  differences <- numeric(length(draws))
  rules <- vector("list", length(quadrature))
  log_bernoulli <- function(prob) {
    value <- ifelse(y == 1, log(prob), log1p(-prob))
    value[is.na(y)] <- 0
    value
  }
  for (r in seq_along(draws)) {
    draw <- draws[r]
    eta <- sweep(x %*% matrix(slopes[, , draw], ncol(x), S), 2, intercept[, draw], "+")
    if (info$model == "binary") {
      g0 <- ifelse(is.na(y) | y == 0, 0, -Inf)
      g1 <- ifelse(is.na(y) | y == 1, 0, -Inf)
    } else {
      if (info$model == "two_stage") {
        a <- rowsum(log_bernoulli(matrix(p[, , draw], primers, S)[primer_index, , drop = FALSE]), sample_index, reorder = TRUE)
        b <- rowsum(log_bernoulli(matrix(q[, , draw], primers, S)[primer_index, , drop = FALSE]), sample_index, reorder = TRUE)
      } else {
        a <- ifelse(is.na(y) | y == 1, 0, -Inf)
        b <- ifelse(is.na(y) | y == 0, 0, -Inf)
      }
      collection_eta <- xt %*% matrix(beta_theta[, , draw], ncol(xt), S)
      fp <- matrix(theta0[, draw], samples, S, byrow = TRUE)
      g0 <- rowsum(site_waic_logadd(log(fp) + a, log1p(-fp) + b), sample_site, reorder = TRUE)
      g1 <- rowsum(site_waic_logadd(site_waic_log_sigmoid(collection_eta) + a,
                                  site_waic_log_sigmoid(-collection_eta) + b), sample_site, reorder = TRUE)
    }
    scaled_loadings <- matrix(loadings[, , draw], factors, S) * sigma[draw]
    evaluate <- function(rule) site_loglik_cpp(eta[observed, , drop = FALSE],
      g0[observed, , drop = FALSE], g1[observed, , drop = FALSE],
      scaled_loadings, rule$nodes, rule$log_weights)
    if (factors == 0L) {
      value <- evaluate(list(nodes = matrix(numeric(), 1, 0), log_weights = 0))
    } else {
      previous <- NULL
      converged <- FALSE
      for (k in seq_along(quadrature)) {
        if (is.null(rules[[k]])) rules[[k]] <- site_waic_rule(quadrature[k], factors, max_nodes)
        value <- evaluate(rules[[k]])
        if (!is.null(previous) && all(is.finite(value)) && all(is.finite(previous))) {
          difference <- max(abs(value - previous))
          if (difference < tolerance) {
            orders[r] <- quadrature[k]
            differences[r] <- difference
            converged <- TRUE
            break
          }
        }
        previous <- value
      }
      if (!converged) stop("Joint quadrature did not converge at posterior draw ", draw,
        "; increase quadrature orders or use held-out-site validation.", call. = FALSE)
    }
    if (any(!is.finite(value))) stop("Nonfinite site log likelihood at posterior draw ", draw, ".", call. = FALSE)
    log_lik[r, ] <- value
  }
  maxima <- apply(log_lik, 2, max)
  lppd <- maxima + log(colMeans(exp(sweep(log_lik, 2, maxima, "-"))))
  penalty <- apply(log_lik, 2, stats::var)
  pointwise <- data.frame(site = colnames(log_lik), lppd = unname(lppd),
    p_waic = unname(penalty), elpd_waic = unname(lppd - penalty),
    waic = unname(-2 * (lppd - penalty)))
  if (any(penalty > 0.4)) warning(sum(penalty > 0.4), " of ", length(penalty),
    " sites have p_waic > 0.4; WAIC may be unreliable. Check held-out-site predictions.", call. = FALSE)
  observations <- list(y = y, site_index = site_index, sites = site_names,
    species = as.character(info$speciesNames), threshold = threshold)
  if (info$model != "binary") {
    observations$sample_index <- sample_index
    observations$primer <- if (info$model == "two_stage") as.character(info$primerNames[primer_index]) else NULL
  }
  structure(list(WAIC = sum(pointwise$waic),
    SE = sqrt(nrow(pointwise) * stats::var(pointwise$waic)),
    elpd_waic = sum(pointwise$elpd_waic), p_waic = sum(penalty), pointwise = pointwise,
    log_lik = log_lik, draw_ids = data.frame(iteration = (draws - 1L) %% iterations + 1L,
      chain = (draws - 1L) %/% iterations + 1L),
    integration = list(orders = orders, max_log_lik_change = max(differences),
      tolerance = tolerance, quadrature = quadrature),
    omitted_sites = site_names[!observed], observations = observations), class = "occJSDM_site_waic")
}

#' Compare two site-level WAIC results
#'
#' @param first,second Results from [computeSiteWAIC()] for identical observed
#'   data, site grouping and response order. Model covariates and factor counts
#'   may differ. Reordering data or changing thresholds is rejected.
#' @return A list with `difference` (first minus second; positive favours the
#'   second), its paired `SE`, and a `pointwise` table. This does not include
#'   numerical integration or MCMC error. A small difference relative to its
#'   SE is weak evidence for choosing one model.
#' @md
#' @export
compareSiteWAIC <- function(first, second) {
  if (!inherits(first, "occJSDM_site_waic") || !inherits(second, "occJSDM_site_waic")) {
    stop("Supply two computeSiteWAIC results.", call. = FALSE)
  }
  if (!identical(first$observations, second$observations) ||
      !identical(first$pointwise$site, second$pointwise$site)) {
    stop("Comparisons require the same scored observations, order, threshold and site grouping.", call. = FALSE)
  }
  difference <- first$pointwise$waic - second$pointwise$waic
  list(difference = sum(difference), SE = sqrt(length(difference) * stats::var(difference)),
       pointwise = data.frame(site = first$pointwise$site, difference = difference))
}
