# Study-only summaries. Independent simulated communities receive equal weight.
calibration_summary <- function(rows) {
  keep <- is.finite(rows$estimate) & is.finite(rows$truth)
  d <- rows[keep, , drop = FALSE]
  valid <- is.finite(d$lower) & is.finite(d$upper) & d$lower <= d$upper
  d$error <- d$estimate - d$truth
  d$covered <- ifelse(valid, d$lower <= d$truth & d$truth <= d$upper, NA_real_)
  d$width <- ifelse(valid, d$upper - d$lower, NA_real_)
  average <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  mcse <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x))
  }
  groups <- split(d, d$replicate)
  means <- function(f) vapply(groups, f, numeric(1))
  bias <- means(function(z) mean(z$error))
  mse <- means(function(z) mean(z$error^2))
  coverage <- means(function(z) average(z$covered))
  widths <- means(function(z) average(z$width))
  data.frame(communities = length(groups), elements = nrow(d),
    bias = average(bias), bias_mcse = mcse(bias), rmse = sqrt(average(mse)),
    intervals = sum(valid), interval_communities = sum(is.finite(coverage)),
    coverage = average(coverage), coverage_mcse = mcse(coverage), width = average(widths),
    missing_intervals = sum(!is.finite(d$lower) | !is.finite(d$upper)),
    invalid_intervals = sum(is.finite(d$lower) & is.finite(d$upper) & d$lower > d$upper))
}

raw_coefficients <- function(beta, centre, spread) {
  shape <- dim(beta)
  stopifnot(length(shape) %in% c(2L, 4L), shape[1] == length(spread) + 1L,
            length(centre) == length(spread), all(is.finite(spread)), all(spread > 0))
  b <- matrix(beta, nrow = shape[1])
  b[-1, ] <- b[-1, , drop = FALSE] / spread
  b[1, ] <- b[1, ] - as.vector(centre %*% b[-1, , drop = FALSE])
  array(b, dim = shape, dimnames = dimnames(beta))
}

probability_intervals <- function(parameters, x, tolerance = 1e-4) {
  # Use every saved global draw. Mean-stability checks from the old scorer
  # do not establish stability of tail quantiles, so check all three anew.
  summarise <- function(nodes) {
    draws <- posterior_marginal(parameters, x, nodes = nodes, retain_draws = TRUE)
    answer <- data.frame(estimate = as.vector(apply(draws, 1:2, mean)),
      lower = as.vector(apply(draws, 1:2, quantile, probs = .025)),
      upper = as.vector(apply(draws, 1:2, quantile, probs = .975)))
    stopifnot(all(is.finite(as.matrix(answer))), all(as.matrix(answer) >= 0),
              all(as.matrix(answer) <= 1))
    answer
  }
  if (parameters$link == "probit") {
    result <- summarise(31L)
    attr(result, "integration") <- list(method = "analytic probit-normal", nodes = NA_integer_, error = 0)
    return(result)
  }
  coarse <- summarise(31L)
  for (nodes in c(61L, 121L, 241L)) {
    fine <- summarise(nodes)
    error <- max(abs(as.matrix(fine) - as.matrix(coarse)))
    if (error <= tolerance) {
      attr(fine, "integration") <- list(method = "normal quadrature; mean and both endpoints checked",
                                        nodes = nodes, error = error)
      return(fine)
    }
    coarse <- fine
  }
  stop("Saved-draw probability integration did not meet the endpoint tolerance")
}

calibration_group_summary <- function(rows) {
  keys <- c("package", "scenario", "n_sites", "n_species", "response", "use_traits", "target", "term")
  collect <- function(d, population) {
    if (!nrow(d)) return(NULL)
    id <- do.call(paste, c(d[keys], sep = "|"))
    do.call(rbind, lapply(split(d, id), function(z) {
      cbind(z[1, keys, drop = FALSE], population = population,
            flagged_communities = length(unique(z$replicate[!z$diagnostic_pass])),
            calibration_summary(z))
    }))
  }
  rbind(collect(rows, "All scored fits"),
        collect(rows[rows$diagnostic_pass, , drop = FALSE], "Passed diagnostics only"))
}
