# The collapsed shift of the community mean mu_b0 and every species' occupancy
# intercept B0, for the experimental b0_prior = "community" in occupancy and
# two-stage fits (dev/simstudy/mu-b0-drift-20261018/REPORT.md).
#
# Each species' B0 is sharp given the latent occupancy states z but broad
# otherwise, so the direction in which all the intercepts move together is
# identified only through z, and the Gibbs sampler walks it one flipped site
# at a time. This step proposes a common shift delta ~ Normal(0, sd^2) of
# every B0 and of mu_b0 and accepts it by Metropolis under the first-stage
# likelihood with z integrated out, site by site and species by species:
# psi' L1 + (1 - psi') L0 over psi L1 + (1 - psi) L0, where psi' =
# logistic(eta + delta), L1 is the probability of the site's field-sample
# outcomes w given z = 1 (Bernoulli in theta) and L0 given z = 0 (Bernoulli
# in theta0), the quantities sample_z_cpp_parallel() forms for the z draw.
# The deviations B0 - mu_b0 are unchanged, so the Normal(mu_b0, sigma_b0^2)
# prior terms of the B0 cancel and only the hyperprior of mu_b0 enters. It
# runs at the top of the iteration, immediately before the z draw, which then
# comes from its full conditional at the shifted predictor; that fresh draw of
# the collapsed variable is what makes the step valid, and the Polya-Gamma
# variables need nothing because update_jSDMcoef() draws them afresh from eta.
# Cost O(n S) per proposal.

# Settings, read from MCMCparams: NULL (nothing recorded, nothing run) unless
# the community prior is on in an occupancy or two-stage fit, where z is
# latent. b0_shift = FALSE switches the step off; b0_shift_sd fixes the
# proposal sd (NULL = adapt it during burn-in); b0_shift_proposals is the
# number of proposals per iteration.
read_b0_shift <- function(MCMCparams, intercept_prior, model) {
  keys <- c("b0_shift", "b0_shift_sd", "b0_shift_proposals")
  given <- keys[keys %in% names(MCMCparams)]
  active <- identical(intercept_prior$type, "community") &&
    model %in% c("occupancy", "two_stage")
  if (!active) {
    if (length(given))
      stop(given[1], " requires b0_prior = 'community' in an occupancy or two-stage fit",
           call. = FALSE)
    return(NULL)
  }
  enabled <- get_param(MCMCparams, "b0_shift", TRUE)
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled))
    stop("b0_shift must be TRUE or FALSE", call. = FALSE)
  if (!enabled) {
    if (any(c("b0_shift_sd", "b0_shift_proposals") %in% given))
      stop("b0_shift_sd and b0_shift_proposals require b0_shift = TRUE", call. = FALSE)
    return(list(enabled = FALSE))
  }
  sd <- MCMCparams[["b0_shift_sd"]]
  if (!is.null(sd) && (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0))
    stop("b0_shift_sd must be a finite positive number", call. = FALSE)
  proposals <- get_param(MCMCparams, "b0_shift_proposals", 5)
  if (!is.numeric(proposals) || length(proposals) != 1L || !is.finite(proposals) ||
      proposals < 1 || proposals != floor(proposals))
    stop("b0_shift_proposals must be a single whole number of at least 1", call. = FALSE)
  list(enabled = TRUE, sd = if (is.null(sd)) NULL else unname(sd),
       proposals = as.integer(proposals))
}

# The starting proposal sd when it is adapted: about the scale of mu_b0's
# posterior sd, which shrinks with the number of species (0.15 to 0.2 at 100
# species, 0.45 to 0.55 at 10, in the simulation studies).
b0_shift_initial_sd <- function(S) 2 / sqrt(S)

# The acceptance rate a one-dimensional random walk is tuned to.
b0_shift_target_acceptance <- 0.44

# Per site and species, the log-probability of the site's field-sample
# outcomes given z = 1 (l1, Bernoulli in theta) and z = 0 (l0, Bernoulli in
# theta0), summed over the site's field samples exactly as the z draw forms
# them. w and theta are N by S over field samples, theta0 is length S, and
# idx_z_w maps each field sample to its site.
b0_shift_pair_loglik <- function(w, theta, theta0, idx_z_w, n) {
  one <- w == 1
  lt <- log1p(-theta)
  lt[one] <- log(theta[one])
  lt0 <- matrix(log1p(-theta0), nrow(w), ncol(w), byrow = TRUE)
  lt0[one] <- log(theta0)[col(w)[one]]
  l1 <- matrix(0, n, ncol(w))
  l0 <- matrix(0, n, ncol(w))
  s1 <- rowsum(lt, idx_z_w, reorder = TRUE)
  s0 <- rowsum(lt0, idx_z_w, reorder = TRUE)
  rows <- as.integer(rownames(s1))
  l1[rows, ] <- s1
  l0[rows, ] <- s0
  list(l1 = l1, l0 = l0)
}

# The first-stage log-likelihood with z integrated out, summed over sites and
# species: log(psi L1 + (1 - psi) L0) at psi = logistic(eta), formed on the
# log scale so that extreme predictors neither overflow nor lose the small term.
b0_shift_loglik <- function(eta, l1, l0) {
  a <- stats::plogis(eta, log.p = TRUE) + l1
  b <- stats::plogis(-eta, log.p = TRUE) + l0
  m <- pmax(a, b)
  sum(m + log1p(exp(-abs(a - b))))
}

# The log acceptance ratio of shifting every B0 and mu_b0 by delta, given the
# current z-integrated log-likelihood (recomputed when not supplied), and the
# proposed log-likelihood for the caller to keep if it accepts.
b0_shift_logratio <- function(delta, eta, l1, l0, mu_b0, hyper,
                              current = b0_shift_loglik(eta, l1, l0)) {
  proposed <- b0_shift_loglik(eta + delta, l1, l0)
  list(logratio = proposed - current +
         stats::dnorm(mu_b0 + delta, hyper$mu_b0_mean, hyper$mu_b0_sd, log = TRUE) -
         stats::dnorm(mu_b0, hyper$mu_b0_mean, hyper$mu_b0_sd, log = TRUE),
       proposed = proposed)
}

# One iteration's proposals. Each draws delta from R's RNG, then a uniform,
# and accepts when log(u) < the log ratio. With adapt = TRUE (burn-in, sd not
# set by the user) the log of sd moves by (acceptance probability - target)
# over the square root of the running proposal count t, a Robbins-Monro
# recursion that settles at the target rate; it is kept within loose bounds.
# Returns the shifted B0, mu_b0 and eta, the number accepted, and the new sd
# and t.
b0_shift_step <- function(B0, mu_b0, eta, l1, l0, hyper, sd, proposals, adapt, t) {
  current <- b0_shift_loglik(eta, l1, l0)
  accepted <- 0L
  for (k in seq_len(proposals)) {
    delta <- stats::rnorm(1, 0, sd)
    u <- stats::runif(1)
    r <- b0_shift_logratio(delta, eta, l1, l0, mu_b0, hyper, current = current)
    lr <- r$logratio
    if (is.finite(lr) && log(u) < lr) {
      B0 <- B0 + delta
      mu_b0 <- mu_b0 + delta
      eta <- eta + delta
      current <- r$proposed
      accepted <- accepted + 1L
    }
    if (adapt) {
      t <- t + 1L
      alpha <- if (is.finite(lr)) min(1, exp(lr)) else 0
      sd <- exp(log(sd) + (alpha - b0_shift_target_acceptance) / sqrt(t))
      sd <- min(max(sd, 1e-3), 1e2)
    }
  }
  list(B0 = B0, mu_b0 = mu_b0, eta = eta, accepted = accepted, sd = sd, t = t)
}

# What infos$b0_shift records: NULL when the settings are NULL, so default
# fits record nothing new; the switched-off list as it is; otherwise the
# settings with each chain's final proposal sd and its acceptance rates over
# the kept iterations and over burn-in (NaN when there was no burn-in).
b0_shift_record <- function(settings, parts) {
  if (is.null(settings)) return(NULL)
  if (!settings$enabled) return(settings)
  rate <- function(phase) vapply(parts, function(p) {
    p[[paste0("accepted_", phase)]] / p[[paste0("proposals_", phase)]]
  }, numeric(1), USE.NAMES = FALSE)
  list(enabled = TRUE, proposals = settings$proposals, adapted = is.null(settings$sd),
       sd = vapply(parts, function(p) p$sd, numeric(1), USE.NAMES = FALSE),
       acceptance = rate("kept"), acceptance_burnin = rate("burnin"))
}
