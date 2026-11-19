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

# All four checks, from the pieces of a fit made with embeddingVars.
# `B_output` is design columns by species by iterations by chains;
# `sigmab_output` and `sigmaemb_output` are iterations by chains;
# `lambda_output` is embedding columns by iterations by chains; `X` is the
# design matrix and `is_emb` marks its embedding columns. Returns the verdict
# ("ok", or "review" when any check is flagged), `checks` (one row per check),
# `chains` (one row per chain) and `effect` (grouped_effect_checks()).
#
# 1. state: the number of chains whose median sigma_emb is at least
#    limits$state, the chains in which the leftovers are at least as large as
#    the shared part. Flagged when some chains are and some are not, and when
#    all are: the second can be right, and its meaning says so.
# 2. rhat_measured and rhat_scales: the package's R-hat (blockRhatESS()) on
#    the measured rows of B, and on the logs of sigma_emb, of the axis scales
#    and, when there are measured columns, of sigma_b. The largest of each.
# 3. gradients and 4. disagreement: grouped_effect_checks() on each chain's
#    posterior mean of the embedding's combined effect at the fitted sites.
# 5. disagreement_beyond: how far each chain's posterior mean of the embedding
#    coefficients is from the pooled one, as a share of their size. Chains
#    whose effects agree at the fitted sites can still differ in directions
#    of the embedding those sites barely vary in, and so differ at sites
#    unlike them: this is the disagreement at sites that vary equally in every
#    embedding direction. Reported, never flagged: there is no evidence yet
#    for a limit.
grouped_convergence <- function(B_output, sigmab_output, sigmaemb_output,
                                lambda_output, X, is_emb,
                                limits = GROUPED_LIMITS) {
  d <- dim(B_output)
  S <- d[2]; niter <- d[3]; nchain <- d[4]
  n_meas <- sum(!is_emb)
  sigmaemb <- matrix(sigmaemb_output, niter, nchain)

  # 1. state
  med <- apply(sigmaemb, 2, stats::median)
  dominate <- med >= limits$state
  chains <- data.frame(chain = seq_len(nchain), sigma_emb_median = med,
                       state = ifelse(dominate, "leftovers dominate", "shared structure"),
                       stringsAsFactors = FALSE)
  n_dom <- sum(dominate)
  state_meaning <- if (n_dom == 0L) {
    "in every chain the shared part is larger than the leftovers"
  } else if (n_dom < nchain) {
    paste0("the chains are in different states: in ", n_dom, " of ", nchain,
           " the leftovers are at least as large as the shared part, so the ",
           "chains describe the species' responses to the embedding differently")
  } else {
    paste0("in every chain the leftovers are at least as large as the shared ",
           "part; right if species respond independently to the embedding, ",
           "wrong if the chains missed shared gradients, which is the likely ",
           "reading with fewer than 5,000 burn-in iterations or a model type ",
           "other than binary")
  }

  # 2. R-hat on the measured coefficients and on the logs of the scales
  max_rhat <- function(arr, n_params) {
    block <- blockRhatESS(arr, n_params, niter, nchain)
    r <- block$rhat[block$variable]
    if (length(r) == 0L || all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)
  }
  rhat_measured <- NA_real_
  rhat_scales <- NA_real_
  if (nchain > 1L) {
    if (n_meas > 0L)
      rhat_measured <- max_rhat(B_output[!is_emb, , , , drop = FALSE], n_meas * S)
    n_emb <- sum(is_emb)
    scales <- array(log(lambda_output), c(n_emb, niter, nchain))
    scales <- array(c(aperm(scales, c(2, 3, 1)), log(sigmaemb),
                      if (n_meas > 0L) log(matrix(sigmab_output, niter, nchain))),
                    c(niter, nchain, n_emb + 1L + (n_meas > 0L)))
    scales <- aperm(scales, c(3, 1, 2))
    rhat_scales <- max_rhat(scales, dim(scales)[1])
  }
  rhat_meaning <- function(value, what, none) {
    if (nchain < 2L) return("not assessed: needs two chains or more")
    if (is.na(value)) return(none)
    if (value > limits$rhat)
      return(paste0("the chains disagree about ", what, ": run longer chains"))
    paste0("the chains agree about ", what)
  }

  # 3 and 4. the embedding's combined effect, from each chain's mean of B
  Xe <- X[, is_emb, drop = FALSE]
  B_mean <- vapply(seq_len(nchain), function(ch)
    rowMeans(array(B_output[is_emb, , , ch], c(sum(is_emb), S, niter)), dims = 2),
    matrix(0, sum(is_emb), S))
  B_mean <- array(B_mean, c(sum(is_emb), S, nchain))
  E <- vapply(seq_len(nchain), function(ch) Xe %*% B_mean[, , ch],
              matrix(0, nrow(X), S))
  effect <- grouped_effect_checks(array(E, c(nrow(X), S, nchain)), limits)
  not_assessed <- paste0("not assessed: ", effect$why_not)

  # 5. the coefficients themselves
  beyond <- NA_real_
  if (effect$assessed) {
    B_pooled <- apply(B_mean, c(1, 2), mean)
    rms <- function(x) sqrt(mean(x^2))
    beyond <- mean(vapply(seq_len(nchain), function(ch)
      rms(B_mean[, , ch] - B_pooled), numeric(1))) / rms(B_pooled)
  }

  flag <- function(x) !is.na(x) && x
  checks <- data.frame(
    check = c("state", "rhat_measured", "rhat_scales", "gradients", "disagreement",
              "disagreement_beyond"),
    value = c(n_dom, rhat_measured, rhat_scales, effect$agreement, effect$disagreement,
              beyond),
    limit = c(limits$state, limits$rhat, limits$rhat, limits$gradients,
              limits$disagreement, NA_real_),
    flagged = c(n_dom > 0L,
                flag(rhat_measured > limits$rhat), flag(rhat_scales > limits$rhat),
                flag(effect$agreement < limits$gradients),
                flag(effect$disagreement > limits$disagreement), FALSE),
    meaning = c(
      state_meaning,
      rhat_meaning(rhat_measured, "the measured coefficients",
                   "not assessed: the fit has no measured covariate"),
      rhat_meaning(rhat_scales, "the scales of the prior", "not assessed"),
      if (!effect$assessed) not_assessed else if (effect$agreement < limits$gradients)
        paste0("the chains found different gradients at the fitted sites: at ",
               "least one chain's fitted effect of the embedding lies partly ",
               "outside another's ", effect$n_gradients, " gradients")
      else paste0("the chains found the same ", effect$n_gradients,
                  " gradients at the fitted sites"),
      if (!effect$assessed) not_assessed else if (effect$disagreement > limits$disagreement)
        paste0("at the fitted sites the chains' effects of the embedding ",
               "differ by more than 15 percent of the effect's size")
      else "at the fitted sites the chains' effects of the embedding are close",
      if (!effect$assessed) not_assessed else
        paste0("the chains' embedding coefficients differ by this share of ",
               "their size, which is how far the chains differ at sites ",
               "unlike the fitted ones (sites that vary equally in every ",
               "direction of the embedding); reported, not part of the verdict")),
    stringsAsFactors = FALSE)

  list(verdict = if (any(checks$flagged)) "review" else "ok",
       checks = checks, chains = chains, effect = effect)
}

# Prints grouped_convergence()'s result with message(): one line per check,
# the verdict, and the meaning of every flagged check on a line starting
# "  ! ". runOccJSDM() calls it at the end of a fit with embedding columns.
print_grouped_convergence <- function(x) {
  v <- function(name) x$checks$value[x$checks$check == name]
  pct <- function(z) sprintf("%.0f%%", 100 * z)
  num <- function(z) if (is.na(z)) "not assessed" else sprintf("%.2f", z)
  nchain <- nrow(x$chains)
  message("\nGrouped prior: convergence checks")
  message("--------------------------------------------------")
  message(sprintf("  State        | median sigma_emb by chain: %s | leftovers dominate in %d of %d chains",
                  paste(sprintf("%.2f", x$chains$sigma_emb_median), collapse = ", "),
                  v("state"), nchain))
  message(sprintf("  R-hat        | measured coefficients, max: %s | scales (log), max: %s",
                  num(v("rhat_measured")), num(v("rhat_scales"))))
  if (x$effect$assessed) {
    message(sprintf("  Gradients    | %d found | agreement between chains at the fitted sites: %.2f (limit %.2f)",
                    x$effect$n_gradients, v("gradients"), GROUPED_LIMITS$gradients))
    message(sprintf("  Disagreement | at the fitted sites: %s of the effect's size (limit %s) | beyond them: %s",
                    pct(v("disagreement")), pct(GROUPED_LIMITS$disagreement), pct(v("disagreement_beyond"))))
  } else {
    message(sprintf("  Gradients    | not assessed: %s", x$effect$why_not))
    message(sprintf("  Disagreement | not assessed: %s", x$effect$why_not))
  }
  message(sprintf("  Verdict      | %s", x$verdict))
  for (i in which(x$checks$flagged)) message(sprintf("  ! %s", x$checks$meaning[i]))
  message("\n==================================================")
  invisible(NULL)
}

#' Convergence checks for a fit made with the grouped prior and an embedding
#'
#' The usual check, R-hat on every coefficient, misleads for a fit made with
#' `embeddingVars` (see \link[=grouped-prior]{grouped-prior}). Single blends
#' and single embedding coefficients are not pinned down, so chains that found
#' the same gradients can still differ on them; and chains that all ended with
#' the leftovers carrying the species' responses agree with each other while
#' reporting no gradients at all. This function asks four questions that do
#' not have those faults. [runOccJSDM()] prints its result at the end of such
#' a fit.
#'
#' @details
#' The four checks.
#' * `state`: which state each chain is in. A chain whose median embedding
#'   leftover scale is 1 or more has leftovers at least as large as the shared
#'   part (the community average, the traits and the blends). `value` is the
#'   number of such chains. It is flagged when the chains are in different
#'   states, and when every chain is in that state. The second can be right,
#'   if species respond independently to the embedding; it is wrong if the
#'   chains missed shared gradients.
#' * `rhat_measured` and `rhat_scales`: the largest R-hat among the measured
#'   coefficients, and among the logs of the scales of the prior. Flagged
#'   over 1.1. These are the only checks here that are an R-hat.
#' * `gradients`: whether the chains found the same set of gradients. For
#'   each chain, the posterior mean of the embedding's combined effect on the
#'   linear predictor (the embedding columns times their coefficients, summed,
#'   for every site and species) is taken, and its leading directions across
#'   sites are that chain's gradients. Their number is the smallest that holds
#'   90 percent of the effect pooled over chains. `value` is the smallest
#'   share, over pairs of chains, of one chain's effect that lies within the
#'   other's gradients. It does not depend on how a chain splits its gradients
#'   into single blends. Flagged under 0.90.
#' * `disagreement`: how far each chain's combined effect is from the pooled
#'   one, as a share of the effect's size, averaged over chains. Flagged over
#'   0.15.
#' * `disagreement_beyond`: how far each chain's embedding coefficients are
#'   from the pooled ones, as a share of their size. Reported and never
#'   flagged.
#'
#' The last three compare the chains' posterior means. They say whether the
#' chains reached the same answer and how far apart they are. They do not say
#' whether each chain explored its own uncertainty, which is what R-hat is
#' for. They are not assessed with a single chain, or when the embedding's
#' effect is small (a spread across sites under 0.5 on the scale of the linear
#' predictor, for the typical species), as when no species responds to the
#' embedding.
#'
#' `gradients` and `disagreement` are about the fitted sites, and hold for
#' new sites that resemble them. They understate how far the chains differ at
#' sites unlike the fitted ones, such as parts of a map whose habitat was not
#' sampled. The data say little about directions of the embedding in which
#' the fitted sites barely vary, so each chain fills those in its own way,
#' and that only shows where a site's embedding moves in those directions. In
#' one simulated fit whose chains differed by 6 percent of the effect's size
#' at the fitted sites, they differed by 16 percent at sites varying equally
#' in every direction, and the error against the truth rose from 27 to 42
#' percent. `disagreement_beyond` is that second kind of difference. It has no
#' limit yet. A fit that is "ok" here is not thereby safe to extrapolate.
#'
#' A flagged fit is marked "review" and not "failed" because the first check
#' can flag a fit that is right. The limits of the last two checks were set
#' against a small number of simulated fits (256 sites, 48 species, 64
#' embedding columns) and may change.
#'
#' @param fitModel Output from the function runOccJSDM, made with
#'   `embeddingVars`. Any model type.
#'
#' @return A list with four elements.
#' * `verdict`: `"ok"`, or `"review"` when any check is flagged.
#' * `checks`: a data frame with one row per check and columns `check`,
#'   `value`, `limit`, `flagged` and `meaning`, a sentence on what the value
#'   says about the fit.
#' * `chains`: a data frame with one row per chain and columns `chain`,
#'   `sigma_emb_median` and `state` (`"shared structure"` or
#'   `"leftovers dominate"`).
#' * `effect`: what lies behind the last two checks, including `n_gradients`,
#'   `size` (the spread of the pooled effect), `agreement_pairs` (entry a, b
#'   is the share of chain a's effect within chain b's gradients) and
#'   `disagreement_by_chain`.
#'
#' @seealso \link[=grouped-prior]{grouped-prior} for the prior that these fits
#'   use, and [computeDiagnostics()] for R-hat and effective sample size by
#'   block of parameters.
#'
#' @md
#' @export
returnGroupedConvergence <- function(fitModel) {
  jo <- fitModel$results_output$jsdm_output
  if (!is_grouped_fit(fitModel) || is.null(jo$lambda_output)) {
    stop("returnGroupedConvergence() needs a fit with embeddingVars")
  }
  grouped_convergence(jo$B_output, jo$sigmab_output, jo$sigmaemb_output,
                      jo$lambda_output, fitModel$X_psi,
                      fitModel$infos$coef_prior$is_emb)
}
