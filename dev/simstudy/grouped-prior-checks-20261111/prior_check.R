# Prior-only check of the grouped covariate prior's default scales.
#
# Question: under the default scales, how large is the embedding's prior
# contribution to the linear predictor? The design asks for a prior median
# between 1 and 2 for the standard deviation across sites of that contribution
# (dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md, section
# "The model, part 1: covariate side", the paragraph on what the defaults are
# meant to do; the outcome is in the same section under "What the prior check
# found").
#
# Run from the package root:
#   Rscript dev/simstudy/grouped-prior-checks-20261111/prior_check.R
#
# First run on 12 November 2026 under the default axis_scale of the time,
# 1 / sqrt(K), which failed the check. The default was changed to 0.7 / K the
# same day and the script was run again. Both outputs are in the README.
#
# Nothing is fitted. Each draw is one draw from the prior of one species'
# coefficients on K embedding columns, with no measured traits, so the trait
# matrix of the covariate side is the single column of ones (the community
# row of G). Where the definitions come from, by function and not by line
# number, which goes stale:
#
#   - Default scales and c_base: read_coef_prior() in R/coef-prior.R.
#     sigma_emb_scale defaults to 1, axis_scale to 0.7 / number of embedding
#     columns, and c_base is 1 / sqrt(gt) on embedding columns. The function
#     is called here, not copied, so the check follows the package's defaults.
#   - How the scales combine: coef_prior_sds() in R/coef-prior.R, also called
#     here. For embedding column k the leftover has standard deviation
#     sigma_emb * lambda[k], the entries of G have lambda[k], and the entries
#     of C have c_base[k] * lambda[k].
#   - The scales' own priors: half-Cauchy with scale axis_scale for each
#     lambda[k] and with scale sigma_emb_scale for sigma_emb, as updated by
#     sample_coef_scales() in R/coef-prior.R.
#   - G and C are Normal with mean 0 and those standard deviations: the prior
#     that sample_GC() in R/jsdmfun.R builds for each covariate.
#   - A is standard Normal: the prior that sample_A() in R/jsdmfun.R uses.
#   - A species' coefficients are its trait row times G, plus its row of A
#     times C, plus its leftover: computeBtcoef() in R/jsdmfun.R.
#   - The site by column matrix is the simulator's factor model for an
#     embedding block, in simulateData() in R/jsdmfun.R (the branch that runs
#     when no embedding$X is supplied): min(8, K) shared factors
#     with weights that decay as 1 / factor number, plus independent noise of
#     standard deviation 0.3, then each column standardised.

suppressMessages(devtools::load_all(".", quiet = TRUE))

N_SITES <- 250L
N_DRAWS <- 2000L
TARGET <- c(1, 2)
TARGET_MEDIAN <- 1.5

half_cauchy <- function(n, scale) abs(stats::rcauchy(n, 0, scale))

# The simulator's embedding columns (simulateData() in R/jsdmfun.R).
embedding_columns <- function(n, K) {
  Fx <- min(8L, K)
  W <- matrix(rnorm(K * Fx), K, Fx) %*% diag(1 / seq_len(Fx), Fx)
  scale(matrix(rnorm(n * Fx), n, Fx) %*% t(W) +
          0.3 * matrix(rnorm(n * K), n, K))[, , drop = FALSE]
}

# The package's prior for K embedding columns and gt latent traits. `priors`
# is what a user would pass as listPriors; empty gives the defaults.
grouped_prior <- function(K, gt, priors = list()) {
  emb <- paste0("X_emb.", seq_len(K))
  columns <- data.frame(column = emb, variable = emb, kind = "numeric")
  read_coef_prior(priors, list(embeddingVars = emb), S = 2L, gt = gt,
                  columns = columns)
}

# Standard deviation across sites of X %*% beta for n_draws prior draws of
# beta. With fixed_scales = TRUE the scales are held at the scales of their
# priors (lambda = axis_scale, sigma_emb = sigma_emb_scale), which is the case
# the design's "about 1.7" refers to; otherwise they are drawn.
prior_contribution_sd <- function(X, cp, gt, n_draws = N_DRAWS,
                                  fixed_scales = FALSE) {
  K <- ncol(X)
  vapply(seq_len(n_draws), function(i) {
    if (fixed_scales) {
      lambda <- rep(cp$axis_scale, K)
      sigma_emb <- cp$sigma_emb_scale
    } else {
      lambda <- half_cauchy(K, cp$axis_scale)
      sigma_emb <- half_cauchy(1L, cp$sigma_emb_scale)
    }
    # sigma_b is not used: every column here is an embedding column.
    sds <- coef_prior_sds(cp, sigma_b = 1, sigma_emb = sigma_emb,
                          lambda = lambda, g_b = 1L, p = K)
    G <- matrix(rnorm(K, 0, sds$g[1, ]), 1L, K)            # community row
    C <- matrix(rnorm(gt * K, 0, rep(sds$c, each = gt)), gt, K)
    A <- matrix(rnorm(gt), 1L, gt)                         # one species
    leftover <- matrix(rnorm(K, 0, sds$b), 1L, K)
    beta <- computeBtcoef(G, matrix(1, 1L, 1L), A, C, leftover)
    stats::sd(X %*% t(beta))
  }, numeric(1))
}

summarise <- function(x) {
  q <- stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  c(q05 = q[1], median = q[2], q95 = q[3])
}

show <- function(label, x) {
  s <- summarise(x)
  cat(sprintf("  %-46s median %8.3f   5%% %8.3f   95%% %10.3f\n",
              label, s["median"], s["q05"], s["q95"]))
  invisible(s)
}

settings <- list(list(K = 64L, gt = 8L, seed = 20261111L),
                 list(K = 128L, gt = 8L, seed = 20261112L))

cat("Prior check of the grouped covariate prior\n")
cat("R ", as.character(getRversion()), ", ", N_DRAWS, " prior draws, ",
    N_SITES, " sites, target median ", TARGET[1], " to ", TARGET[2], "\n",
    sep = "")

for (s in settings) {
  set.seed(s$seed)
  X <- embedding_columns(N_SITES, s$K)
  cp <- grouped_prior(s$K, s$gt)
  cat(sprintf("\nK = %d embedding columns, gt = %d\n", s$K, s$gt))
  cat(sprintf("  defaults: axis_scale = %.5f (%.3g / K), sigma_emb_scale = %g, c_base = %.5f\n",
              cp$axis_scale, cp$axis_scale * s$K, cp$sigma_emb_scale, cp$c_base[1]))

  drawn <- show("default scales, scales drawn (the check)",
                prior_contribution_sd(X, cp, s$gt))
  inside <- drawn["median"] >= TARGET[1] && drawn["median"] <= TARGET[2]
  cat("  median inside ", TARGET[1], " to ", TARGET[2], ": ",
      if (inside) "YES" else "NO", "\n", sep = "")

  # Diagnostics, to say where the size comes from.
  show("scales held at their prior scales",
       prior_contribution_sd(X, cp, s$gt, fixed_scales = TRUE))
  X_indep <- scale(matrix(rnorm(N_SITES * s$K), N_SITES, s$K))
  show("independent columns, scales held",
       prior_contribution_sd(X_indep, cp, s$gt, fixed_scales = TRUE))
  show("independent columns, scales drawn",
       prior_contribution_sd(X_indep, cp, s$gt))

  # Every coefficient is lambda[k] times something that does not depend on
  # axis_scale, so the contribution is proportional to axis_scale and the
  # value that puts the median at 1.5 follows by proportion. Confirmed with
  # fresh draws.
  axis_new <- cp$axis_scale * TARGET_MEDIAN / drawn["median"]
  cat(sprintf("  axis_scale for a median of %.1f: %.5f = %.3f / K = %.3f / sqrt(K)\n",
              TARGET_MEDIAN, axis_new, axis_new * s$K, axis_new * sqrt(s$K)))
  show(sprintf("axis_scale = %.5f, fresh draws", axis_new),
       prior_contribution_sd(X, grouped_prior(s$K, s$gt, list(axis_scale = axis_new)),
                             s$gt))

  # The other scale a user can set. It reaches only the leftover.
  show("sigma_emb_scale = 0.01, axis_scale default",
       prior_contribution_sd(X, grouped_prior(s$K, s$gt, list(sigma_emb_scale = 0.01)),
                             s$gt))
}
