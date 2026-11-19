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
# limit. `gradients`: the least agreement on the gradients between two chains
# (see grouped_effect_checks()). `disagreement`: the most two chains' effects
# may differ, as a share of the effect's size. `small_effect`: below this
# spread across sites, on the scale of the linear predictor, the embedding's
# effect is too small for the last two checks to mean anything. `hold`: the
# gradients are as many leading directions as hold this share of the pooled
# effect. The limits for the gradients and the disagreement were looked at
# against 18 saved fits whose state was already known, which is the only
# evidence there was and is circular to that extent. Among them the chains of
# the fits with nothing known against them differ by 0.06 to 0.30, and those
# of the two fits with chains known to be apart by 0.34 and 0.68, so the limit
# of 0.30 has almost no margin; and none of the 18 has chains with clearly
# different gradients (the lowest agreement is 0.93), so the limit of 0.90
# rests on made-up effects only
# (dev/simstudy/grouped-mixing-20261118/README.md).
GROUPED_LIMITS <- list(state = 1, rhat = 1.1, gradients = 0.90,
                       disagreement = 0.30, small_effect = 0.5, hold = 0.90)

# The two checks of the embedding's combined effect. `E` is a sites by species
# by chains array: each chain's posterior mean of the embedding columns times
# their coefficients. Both checks compare chain means, so they say whether the
# chains reached the same answer, not whether each chain explored its own
# uncertainty.
#
# The gradients of a chain are the leading directions, across sites, of its
# effect; their number is the smallest that holds `limits$hold` of the effect
# pooled over chains. `agreement_pairs[a, b]` is how much of chain a's effect
# lies within chain b's gradients, as a share of how much lies within chain
# a's own. It is 1 when the two chains have the same gradients, whatever share
# of the effect those hold, and cannot exceed 1; it does not depend on how
# either chain splits its gradients into single blends. `disagreement_pairs[a,
# b]` is how far apart the two chains' effects are, as a share of the pooled
# effect's size (root mean squares). `agreement` is the smallest of the first
# and `disagreement` the largest of the second, so one chain apart from many
# that agree counts in full. Not assessed, with the reason in `why_not`, with
# fewer than two chains, when the effect has missing or non-finite values
# (`missing`), or when the effect is small.
grouped_effect_checks <- function(E, limits = GROUPED_LIMITS) {
  nchain <- dim(E)[3]
  blank <- matrix(NA_real_, nchain, nchain)
  out <- list(assessed = FALSE, why_not = NA_character_, missing = FALSE,
              size = NA_real_, n_gradients = NA_integer_, agreement = NA_real_,
              agreement_pairs = blank, disagreement = NA_real_,
              disagreement_pairs = blank)
  if (!all(is.finite(E))) {
    out$missing <- TRUE
    out$why_not <- "the embedding's effect has missing or non-finite values"
    return(out)
  }
  pooled <- rowMeans(E, dims = 2)
  out$size <- stats::median(apply(pooled, 2, stats::sd))
  if (nchain < 2L) {
    out$why_not <- "needs two chains or more"
    return(out)
  }
  if (!is.finite(out$size) || out$size < limits$small_effect) {
    out$why_not <- "the embedding's effect is small"
    return(out)
  }
  d <- svd(pooled, nu = 0, nv = 0)$d
  r <- which(cumsum(d^2) / sum(d^2) >= limits$hold)[1]
  gradients <- vector("list", nchain)
  own <- numeric(nchain)                 # what each chain's own gradients hold
  for (k in seq_len(nchain)) {
    s <- svd(E[, , k], nu = r, nv = 0)
    gradients[[k]] <- s$u
    own[k] <- sum(s$d[seq_len(r)]^2)
  }
  rms <- function(x) sqrt(mean(x^2))
  size_pooled <- rms(pooled)
  for (a in seq_len(nchain)) for (b in seq_len(nchain)) if (a != b) {
    out$agreement_pairs[a, b] <- sum(crossprod(gradients[[b]], E[, , a])^2) / own[a]
    out$disagreement_pairs[a, b] <- rms(E[, , a] - E[, , b]) / size_pooled
  }
  out$assessed <- TRUE
  out$n_gradients <- as.integer(r)
  out$agreement <- min(out$agreement_pairs, na.rm = TRUE)
  out$disagreement <- max(out$disagreement_pairs, na.rm = TRUE)
  out
}

# All the checks, from the pieces of a fit made with embeddingVars.
# `B_output` is design columns by species by iterations by chains;
# `sigmab_output` and `sigmaemb_output` are iterations by chains;
# `lambda_output` is embedding columns by iterations by chains; `X` is the
# design matrix and `is_emb` marks its embedding columns. Returns the verdict
# ("ok", or "review" when any check is flagged), `checks` (one row per check),
# `chains` (one row per chain) and `effect` (grouped_effect_checks()).
#
# 1. state: the number of chains whose median sigma_emb is at least
#    limits$state, the chains in which the leftovers are at least as large as
#    the shared part. Flagged when it is over 0: when some chains are and some
#    are not, and when all are. The second can be right, and its meaning says
#    so.
# 2. rhat_measured and rhat_scales: the package's R-hat (blockRhatESS()) on
#    the measured rows of B, and on the logs of sigma_emb and, when there are
#    measured columns, of sigma_b. The largest of each.
#    rhat_axis: the largest R-hat among the logs of the axis scales, with how
#    many are over the limit. Reported, never part of the verdict: a single
#    axis scale, like a single embedding coefficient, is not pinned down, and
#    in full-length fits whose chains agree on everything else a few of the
#    64 fail R-hat (README, 19 November 2026). Flagged only when its draws
#    have missing values.
# 3. gradients and 4. disagreement: grouped_effect_checks() on each chain's
#    posterior mean of the embedding's combined effect at the fitted sites.
# 5. disagreement_beyond: the largest difference between two chains' posterior
#    means of the embedding coefficients, as a share of their size. Chains
#    whose effects agree at the fitted sites can still differ in directions
#    of the embedding those sites barely vary in, and so differ at sites
#    unlike them: this is the disagreement at sites that vary equally in every
#    embedding direction. Reported, never flagged: there is no evidence yet
#    for a limit.
#
# B_output is the largest object a fit holds, and none of it is copied: its
# chain means and R-hat come from one pass over the whole array, and the
# measured and embedding rows are told apart afterwards. Draws that are
# missing or not finite flag the check they belong to; they are never an
# error, because this runs at the end of a fit that must still be returned.
grouped_convergence <- function(B_output, sigmab_output, sigmaemb_output,
                                lambda_output, X, is_emb,
                                limits = GROUPED_LIMITS) {
  d <- dim(B_output)
  p <- d[1]; S <- d[2]; niter <- d[3]; nchain <- d[4]
  is_emb <- as.logical(is_emb)
  if (length(is_emb) != p) {
    stop("is_emb must have one entry per row of B_output (", p, "), not ",
         length(is_emb))
  }
  n_meas <- sum(!is_emb); n_emb <- sum(is_emb)
  sigmaemb <- matrix(sigmaemb_output, niter, nchain)

  # 1. state
  med <- apply(sigmaemb, 2, stats::median)
  known <- is.finite(med)
  dominate <- known & med >= limits$state
  chains <- data.frame(chain = seq_len(nchain), sigma_emb_median = med,
                       state = ifelse(!known, "not known",
                                      ifelse(dominate, "leftovers dominate",
                                             "shared structure")),
                       stringsAsFactors = FALSE)
  n_dom <- if (all(known)) sum(dominate) else NA_real_
  state_meaning <- if (is.na(n_dom)) {
    paste0("the state of ", sum(!known), " of ", nchain, " chains is not ",
           "known: the draws of sigma_emb have missing or non-finite values")
  } else if (n_dom == 0L) {
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

  # One pass over B: every coefficient's mean in each chain, and its R-hat.
  # Coefficient i of species j is entry i + (j - 1) * p.
  block <- blockRhatESS(B_output, p * S, niter, nchain)
  emb_coef <- rep(is_emb, times = S)

  # 2. R-hat on the measured coefficients and on the logs of the scales
  max_rhat <- function(block, keep = TRUE) {
    r <- block$rhat[block$variable & keep]
    if (length(r) == 0L || all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)
  }
  rhat_measured <- NA_real_
  rhat_scales <- NA_real_
  rhat_axis <- NA_real_
  axis_over <- NA_integer_
  missing_measured <- n_meas > 0L && anyNA(block$var[!emb_coef])
  # sigma_emb and, with measured columns, sigma_b: one row each, on the log scale
  shared <- array(c(log(sigmaemb), if (n_meas > 0L) log(matrix(sigmab_output, niter, nchain))),
                  c(niter, nchain, 1L + (n_meas > 0L)))
  shared <- aperm(shared, c(3, 1, 2))
  axis <- array(log(lambda_output), c(n_emb, niter, nchain))
  missing_scales <- !all(is.finite(shared))
  missing_axis <- !all(is.finite(axis))
  if (nchain > 1L) {
    if (n_meas > 0L && !missing_measured) rhat_measured <- max_rhat(block, !emb_coef)
    if (!missing_scales)
      rhat_scales <- max_rhat(blockRhatESS(shared, dim(shared)[1], niter, nchain))
    if (!missing_axis) {
      axis_block <- blockRhatESS(axis, n_emb, niter, nchain)
      rhat_axis <- max_rhat(axis_block)
      axis_over <- sum(axis_block$rhat[axis_block$variable] > limits$rhat, na.rm = TRUE)
    }
  }
  rhat_meaning <- function(value, what, none, missing) {
    if (missing)
      return(paste0("the draws of ", what, " have missing or non-finite values"))
    if (nchain < 2L) return("not assessed: needs two chains or more")
    if (is.na(value)) return(none)
    if (value > limits$rhat)
      return(paste0("the chains disagree about ", what, ": run longer chains"))
    paste0("the chains agree about ", what)
  }

  # 3 and 4. the embedding's combined effect, from each chain's mean of B
  Xe <- X[, is_emb, drop = FALSE]
  B_mean <- array(block$chain_mean[emb_coef, , drop = FALSE], c(n_emb, S, nchain))
  E <- vapply(seq_len(nchain), function(ch) Xe %*% B_mean[, , ch],
              matrix(0, nrow(X), S))
  effect <- grouped_effect_checks(array(E, c(nrow(X), S, nchain)), limits)
  not_assessed <- paste0("not assessed: ", effect$why_not)

  # 5. the coefficients themselves
  beyond <- NA_real_
  if (effect$assessed) {
    rms <- function(x) sqrt(mean(x^2))
    apart <- 0
    for (a in seq_len(nchain - 1L)) for (b in (a + 1L):nchain)
      apart <- max(apart, rms(B_mean[, , a] - B_mean[, , b]))
    beyond <- apart / rms(rowMeans(B_mean, dims = 2))
  }

  flag <- function(x) !is.na(x) && x
  percent <- function(x) paste0(format(100 * x), " percent")
  checks <- data.frame(
    check = c("state", "rhat_measured", "rhat_scales", "rhat_axis", "gradients",
              "disagreement", "disagreement_beyond"),
    value = c(n_dom, rhat_measured, rhat_scales, rhat_axis, effect$agreement,
              effect$disagreement, beyond),
    limit = c(0, limits$rhat, limits$rhat, NA_real_, limits$gradients,
              limits$disagreement, NA_real_),
    flagged = c(is.na(n_dom) || n_dom > 0L,
                missing_measured || flag(rhat_measured > limits$rhat),
                missing_scales || flag(rhat_scales > limits$rhat),
                missing_axis,
                effect$missing || flag(effect$agreement < limits$gradients),
                flag(effect$disagreement > limits$disagreement), FALSE),
    meaning = c(
      state_meaning,
      rhat_meaning(rhat_measured, "the measured coefficients",
                   "not assessed: the fit has no measured covariate", missing_measured),
      rhat_meaning(rhat_scales,
                   if (n_meas > 0L) "sigma_emb and sigma_b" else "sigma_emb",
                   "not assessed", missing_scales),
      if (missing_axis) "the draws of the axis scales have missing or non-finite values"
      else if (nchain < 2L) "not assessed: needs two chains or more"
      else paste0(axis_over, " of ", n_emb, " axis scales have an R-hat over ",
                  format(limits$rhat), "; a single axis scale is not pinned ",
                  "down, so this is reported, not part of the verdict"),
      if (!effect$assessed) not_assessed else if (effect$agreement < limits$gradients)
        paste0("the chains found different gradients at the fitted sites: at ",
               "least one chain's fitted effect of the embedding lies partly ",
               "outside another's ", effect$n_gradients, " gradients")
      else paste0("the chains found the same ", effect$n_gradients,
                  " gradients at the fitted sites"),
      if (!effect$assessed) not_assessed else if (effect$disagreement > limits$disagreement)
        paste0("at the fitted sites two chains' effects of the embedding ",
               "differ by more than ", percent(limits$disagreement),
               " of the effect's size")
      else "at the fitted sites the chains' effects of the embedding are close",
      if (!effect$assessed) not_assessed else
        paste0("the largest difference between two chains' embedding ",
               "coefficients, as a share of their size, which is how far the ",
               "chains differ at sites unlike the fitted ones (sites that vary ",
               "equally in every direction of the embedding); reported, not ",
               "part of the verdict")),
    stringsAsFactors = FALSE)

  list(verdict = if (any(checks$flagged)) "review" else "ok",
       checks = checks, chains = chains, effect = effect)
}

# The step at the end of a fit made with embeddingVars: the checks, printed,
# and one warning when the verdict is "review". `jsdm_output` is the fit's
# jsdm_output (B_output, sigmab_output, sigmaemb_output, lambda_output). The
# sampler has finished and the fit must be returned whatever happens here, so
# an error in the checks or in printing them becomes a warning. Returns the
# checks, or NULL when they could not be computed, invisibly.
report_grouped_convergence <- function(jsdm_output, X, is_emb) {
  check <- tryCatch({
    x <- grouped_convergence(jsdm_output$B_output, jsdm_output$sigmab_output,
                             jsdm_output$sigmaemb_output, jsdm_output$lambda_output,
                             X, is_emb)
    print_grouped_convergence(x)
    x
  }, error = function(e) {
    warning("The grouped prior's convergence checks could not be computed (",
            conditionMessage(e), "). The fit is returned: try ",
            "returnGroupedConvergence() on it.", call. = FALSE)
    NULL
  })
  if (!is.null(check) && check$verdict == "review") {
    warning("The grouped prior's convergence checks flag this fit for review: ",
            "see returnGroupedConvergence().", call. = FALSE)
  }
  invisible(check)
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
  message(sprintf("  State        | median sigma_emb by chain: %s | leftovers dominate in %s of %d chains",
                  paste(sprintf("%.2f", x$chains$sigma_emb_median), collapse = ", "),
                  if (is.na(v("state"))) "an unknown number" else format(v("state")), nchain))
  message(sprintf("  R-hat        | measured coefficients, max: %s | sigma_emb and sigma_b (log), max: %s",
                  num(v("rhat_measured")), num(v("rhat_scales"))))
  axis_said <- x$checks$meaning[x$checks$check == "rhat_axis"]
  message(sprintf("               | axis scales (log): %s",
                  if (is.na(v("rhat_axis"))) axis_said else
                    sprintf("%s over %.2f, max: %.2f, not part of the verdict",
                            sub(" axis scales.*$", "", axis_said), GROUPED_LIMITS$rhat, v("rhat_axis"))))
  if (x$effect$assessed) {
    message(sprintf("  Gradients    | %d found | agreement between chains at the fitted sites: %.2f (limit %.2f)",
                    x$effect$n_gradients, v("gradients"), GROUPED_LIMITS$gradients))
    message(sprintf("  Disagreement | largest between two chains, at the fitted sites: %s of the effect's size (limit %s) | beyond them: %s",
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
#' reporting no gradients at all. This function asks questions that do not
#' have those faults. [runOccJSDM()] prints its result at the end of such a
#' fit.
#'
#' @details
#' The checks, one row each in `checks`.
#' * `state`: which state each chain is in. A chain whose median embedding
#'   leftover scale is 1 or more has leftovers at least as large as the shared
#'   part (the community average, the traits and the blends). `value` is the
#'   number of such chains, and it is flagged when it is over 0: when the
#'   chains are in different states, and when every chain is in that state.
#'   The second can be right, if species respond independently to the
#'   embedding; it is wrong if the chains missed shared gradients.
#' * `rhat_measured` and `rhat_scales`: the largest R-hat among the measured
#'   coefficients, and among the logs of the two scales shared by all columns
#'   (the embedding leftover scale and, with measured covariates, the
#'   measured-group scale). Flagged over 1.1.
#' * `rhat_axis`: the largest R-hat among the logs of the axis scales, with
#'   the number over 1.1 in `meaning`. Reported and not part of the verdict.
#'   A single axis scale, like a single embedding coefficient, is not pinned
#'   down, and a few of them fail R-hat in fits whose chains agree on
#'   everything else.
#' * `gradients`: whether the chains found the same set of gradients. For
#'   each chain, the posterior mean of the embedding's combined effect on the
#'   linear predictor (the embedding columns times their coefficients, summed,
#'   for every site and species) is taken, and its leading directions across
#'   sites are that chain's gradients. Their number is the smallest that holds
#'   90 percent of the effect pooled over chains. For a pair of chains, the
#'   agreement is how much of one chain's effect lies within the other's
#'   gradients, as a share of how much lies within its own: 1 when the two
#'   have the same gradients, and less the further apart they are. It does not
#'   depend on how a chain splits its gradients into single blends. `value` is
#'   the smallest over pairs. Flagged under 0.90.
#' * `disagreement`: the largest difference between two chains' combined
#'   effects, as a share of the effect's size (root mean squares). One chain
#'   apart from several that agree counts in full. Flagged over 0.30.
#' * `disagreement_beyond`: the largest difference between two chains'
#'   embedding coefficients, as a share of their size. Reported and never
#'   flagged.
#'
#' The three R-hat rows are the only ones that are an R-hat. The last three
#' rows compare the chains' posterior means. They say whether the
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
#' one simulated fit whose two chains differed by 6 percent of the size of
#' the true effect at the fitted sites, they differed by 16 percent at sites
#' varying equally in every direction, and the error against the truth rose
#' from 27 to 42 percent. `disagreement_beyond` is that second kind of
#' difference. It has no limit yet. A fit that is "ok" here is not thereby
#' safe to extrapolate.
#'
#' Draws that are missing or not finite flag the check they belong to, with
#' that as the stated meaning.
#'
#' A flagged fit is marked "review" and not "failed" because the first check
#' can flag a fit that is right.
#'
#' @section How far the limits have been tested:
#' The limits of `gradients` and `disagreement` were chosen on simulated fits
#' whose answer was already known, all of one size (256 sites, 48 species, 64
#' embedding columns, 4 chains of 5,000 burn-in and 5,000 kept iterations).
#' They have not been tested on fits they were not chosen on, and may change.
#' Among 18 such fits, the 16 with nothing known against them had a
#' disagreement of 0.06 to 0.30 and the two with chains known to be apart
#' 0.34 and 0.68, so the limit of 0.30 has almost no margin. None of the 18
#' had chains with clearly different gradients (the agreement was 0.93 to
#' 1.00), so the limit of 0.90 rests on made-up effects alone.
#'
#' The whole function has been run on four full-length fits of that size, two
#' of them repeats of fits among the 18. Data with no response to the
#' embedding: "ok". Data with two shared gradients, chains in agreement: "ok",
#' with 2 of the 64 axis scales over 1.1, which is why `rhat_axis` does not
#' decide the verdict. Data with three shared gradients: "review", on the
#' R-hat of the measured coefficients (1.16). Data in which species respond
#' independently to the embedding: "review" on `state`, every chain's leftover
#' scale being near 8, although that fit recovered the truth as well as the
#' others. That last flag is the one that can be right to set aside.
#'
#' @param fitModel Output from the function runOccJSDM, made with
#'   `embeddingVars`. Any model type.
#'
#' @return A list with four elements.
#' * `verdict`: `"ok"`, or `"review"` when any check is flagged.
#' * `checks`: a data frame with one row per check (seven) and columns `check`,
#'   `value`, `limit`, `flagged` and `meaning`, a sentence on what the value
#'   says about the fit.
#' * `chains`: a data frame with one row per chain and columns `chain`,
#'   `sigma_emb_median` and `state` (`"shared structure"` or
#'   `"leftovers dominate"`).
#' * `effect`: what lies behind `gradients` and `disagreement`, including
#'   `n_gradients`, `size` (the spread of the pooled effect),
#'   `agreement_pairs` (entry a, b is how much of chain a's effect lies within
#'   chain b's gradients, as a share of how much lies within chain a's own)
#'   and `disagreement_pairs` (entry a, b is the difference between the two
#'   chains' effects, as a share of the effect's size).
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
