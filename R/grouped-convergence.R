# A convergence check for fits made with the grouped prior and an embedding.
#
# R-hat on single embedding coefficients is the wrong test for these fits.
# Single blends and single embedding coefficients are not pinned down, and
# their draws are heavy-tailed, so it fails fits that are fine. And when every
# chain ends with the leftovers carrying the species' responses, the chains
# agree with each other, so it passes a fit that found no gradients. The
# checks here ask instead which state each chain is in, whether the measured
# coefficients and the scales have converged, whether the chains found the
# same set of gradients, and how far apart their fitted effects are. Design:
# dev/superpowers/specs/2026-11-19-grouped-convergence-check-design.md.
# Evidence: dev/simstudy/grouped-mixing-20261118/README.md.

# The limits of the checks. `state`: a chain whose median sigma_emb is at
# least this has leftovers as large as the shared part. `rhat`: the usual
# limit. `gradients`: the least share of one chain's effect that must lie
# within another's gradients. `disagreement`: the most a chain's effect may
# differ from the pooled one, as a share of its size. `small_effect`: below
# this spread across sites, on the scale of the linear predictor, the
# embedding's effect is too small for the last two checks to mean anything.
# `hold`: the gradients are as many leading directions as hold this share of
# the pooled effect. The limits for the gradients and the disagreement were
# set against 18 saved fits whose state was already known, which is the only
# evidence there was and is circular to that extent.
GROUPED_LIMITS <- list(state = 1, rhat = 1.1, gradients = 0.90,
                       disagreement = 0.15, small_effect = 0.5, hold = 0.90)

# The two checks of the embedding's combined effect. `E` is a sites by species
# by chains array: each chain's posterior mean of the embedding columns times
# their coefficients. Both checks compare chain means, so they say whether the
# chains reached the same answer, not whether each chain explored its own
# uncertainty.
#
# The gradients of a chain are the leading directions, across sites, of its
# effect; their number is the smallest that holds `limits$hold` of the effect
# pooled over chains. `agreement_pairs[a, b]` is the share of chain a's effect
# that lies within chain b's gradients, which does not depend on how either
# chain splits its gradients into single blends. `disagreement_by_chain` is
# how far each chain's effect is from the pooled one, as a share of the pooled
# effect's size (root mean squares). Not assessed, with the reason in
# `why_not`, with fewer than two chains or when the effect is small.
grouped_effect_checks <- function(E, limits = GROUPED_LIMITS) {
  nchain <- dim(E)[3]
  pooled <- apply(E, c(1, 2), mean)
  size <- stats::median(apply(pooled, 2, stats::sd))
  out <- list(assessed = FALSE, why_not = NA_character_, size = size,
              n_gradients = NA_integer_, agreement = NA_real_,
              agreement_pairs = matrix(NA_real_, nchain, nchain),
              disagreement = NA_real_,
              disagreement_by_chain = rep(NA_real_, nchain))
  if (nchain < 2L) {
    out$why_not <- "needs two chains or more"
    return(out)
  }
  if (!is.finite(size) || size < limits$small_effect) {
    out$why_not <- "the embedding's effect is small"
    return(out)
  }
  d <- svd(pooled, nu = 0, nv = 0)$d
  r <- which(cumsum(d^2) / sum(d^2) >= limits$hold)[1]
  gradients <- lapply(seq_len(nchain), function(k) svd(E[, , k], nu = r, nv = 0)$u)
  for (a in seq_len(nchain)) for (b in seq_len(nchain)) if (a != b) {
    out$agreement_pairs[a, b] <-
      sum(crossprod(gradients[[b]], E[, , a])^2) / sum(E[, , a]^2)
  }
  rms <- function(x) sqrt(mean(x^2))
  out$disagreement_by_chain <- vapply(seq_len(nchain), function(k)
    rms(E[, , k] - pooled), numeric(1)) / rms(pooled)
  out$assessed <- TRUE
  out$n_gradients <- as.integer(r)
  out$agreement <- min(out$agreement_pairs, na.rm = TRUE)
  out$disagreement <- mean(out$disagreement_by_chain)
  out
}
