# Pieces of the count model that live in R: the effort offset and the update
# of the negative-binomial size, one size per species
# (dev/superpowers/specs/2026-11-12-count-model-design.md, "How it is fitted").
#
# A count y has mean E * mu and variance mean + mean^2 / r, where E is the
# sample's known effort and r is the species' size. Given the log mean, the
# sizes of different species are independent, and none has a standard full
# conditional, so each is updated by a random-walk Metropolis step on log(r):
# propose log(r) + step * Normal(0, 1) and accept under the species' summed
# negative-binomial log density plus the prior. The prior is a half-normal on
# 1 / sqrt(r), which pulls towards the Poisson (1 / sqrt(r) = 0) when the data
# say little, as they do for a rare species. The step size of each species is
# adapted during burn-in, as in b0_shift_step().
#
# This is stage 1 of the count model: nothing here is called by the sampler
# loop yet. Stage 2 calls count_log_mean() to form the log mean of every
# sample and count_size_step() once per iteration.

# The acceptance rate a one-dimensional random walk is tuned to.
count_size_target_acceptance <- 0.44

# The starting proposal standard deviation on the log size.
count_size_initial_step <- 0.5

# The size is kept inside these bounds; a proposal outside is rejected.
count_size_bounds <- c(1e-6, 1e6)

# The log mean of every sample and species: the site's predictor plus the
# sample's log effort. eta is sites by species, idx_site maps each sample to
# its site, and log_effort has one value per sample.
count_log_mean <- function(eta, idx_site, log_effort) {
  if (!is.numeric(idx_site) || anyNA(idx_site) || any(idx_site != floor(idx_site)) ||
      any(idx_site < 1) || any(idx_site > nrow(eta)))
    stop("idx_site must index the rows of eta", call. = FALSE)
  if (!is.numeric(log_effort) || length(log_effort) != length(idx_site) ||
      !all(is.finite(log_effort)))
    stop("log_effort must be finite, with one value per sample", call. = FALSE)
  eta[idx_site, , drop = FALSE] + log_effort
}

# The log prior density of log(r), up to a constant, when 1 / sqrt(r) is
# half-normal with standard deviation prior_sd. The second term is the
# Jacobian of the change from 1 / sqrt(r) to log(r).
count_size_logprior <- function(log_r, prior_sd) {
  -1 / (2 * prior_sd^2 * exp(log_r)) - log_r / 2
}

# One Metropolis proposal for every species' size. y and log_mean are samples
# by species; include is a logical matrix of the same shape, TRUE where the
# cell counts towards the likelihood (every cell in the plain count JSDM,
# occupied cells in the shared-predictor model). A missing count is ignored.
# Each call draws S normals and then S uniforms from R's RNG and nothing
# else, so a run is reproducible from set.seed(). With adapt = TRUE (burn-in)
# the log of each step moves by (acceptance probability - target) over the
# square root of the running count of adapting calls t, the Robbins-Monro
# recursion of b0_shift_step(), within loose bounds. Returns the sizes,
# which proposals were accepted, and the new steps and t.
count_size_step <- function(r, y, log_mean, include, prior_sd, step, adapt, t) {
  S <- ncol(y)
  if (length(r) != S || length(step) != S)
    stop("one size per species", call. = FALSE)
  if (!identical(dim(y), dim(log_mean)) || !identical(dim(y), dim(include)))
    stop("y, log_mean and include must have the same dimensions", call. = FALSE)

  z <- stats::rnorm(S)
  u <- stats::runif(S)

  log_r <- log(r)
  log_proposed <- log_r + step * z
  proposed <- exp(log_proposed)
  inside <- is.finite(proposed) &
    proposed >= count_size_bounds[1] & proposed <= count_size_bounds[2]

  # A proposal outside the bounds is rejected without being evaluated, so the
  # current size stands in for it in the density call.
  used <- include & !is.na(y)
  species_loglik <- function(size) {
    ll <- nbLogPmf_cpp(y, log_mean, size)
    ll[!used] <- 0
    colSums(ll)
  }
  logratio <- species_loglik(ifelse(inside, proposed, r)) - species_loglik(r) +
    count_size_logprior(log_proposed, prior_sd) - count_size_logprior(log_r, prior_sd)
  logratio[!inside] <- -Inf

  accepted <- is.finite(logratio) & log(u) < logratio
  r[accepted] <- proposed[accepted]

  if (adapt) {
    t <- t + 1L
    alpha <- ifelse(is.finite(logratio), pmin(1, exp(logratio)), 0)
    step <- exp(log(step) + (alpha - count_size_target_acceptance) / sqrt(t))
    step <- pmin(pmax(step, 1e-3), 1e2)
  }
  list(r = r, accepted = accepted, step = step, t = t)
}
