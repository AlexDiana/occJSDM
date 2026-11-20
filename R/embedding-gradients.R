# Gradients of the species' whole response to the embedding, and their
# comparison with covariates known for the same sites. Design, with the
# mathematics: dev/superpowers/specs/2026-11-20-embedding-gradients-design.md.
#
# The effect is E = X_emb %*% B_emb, sites by species, from the posterior mean
# of the embedding coefficients: each chain's, and their mean over chains
# (pooled). These are the arrays the convergence check builds
# (embedding_effect() in grouped-convergence.R), so the number of gradients and
# the agreement between chains are the check's own.

# The number of singular values `d` (largest first) of an n by S matrix that
# are above rounding error.
effect_rank <- function(d, n, S) {
  if (length(d) == 0L || !is.finite(d[1]) || d[1] <= 0) return(0L)
  sum(d > max(n, S) * d[1] * .Machine$double.eps)
}

# The first `r` gradients of one effect matrix `E1` (sites by species), where
# E1 = X_emb %*% B1 and `B1` is embedding columns by species. Site scores have
# mean square 1, so a species' loading is the change in its linear predictor
# per unit of the gradient, and scores %*% t(loadings) is the best rank-r
# approximation of E1. `weights` gives the scores as X_emb %*% weights. Each
# gradient is signed so that its weight of largest absolute value is positive.
# `share` has one row per direction up to the rank. NULL when the rank of E1
# is below `r`.
gradient_parts <- function(E1, B1, r) {
  n <- nrow(E1); S <- ncol(E1)
  s <- svd(E1, nu = r, nv = r)
  rank <- effect_rank(s$d, n, S)
  if (rank < r) return(NULL)
  keep <- seq_len(r)
  d <- s$d[keep]
  V <- s$v[, keep, drop = FALSE]
  scores <- sqrt(n) * s$u[, keep, drop = FALSE]
  loadings <- V %*% diag(d, nrow = r) / sqrt(n)
  weights <- B1 %*% V %*% diag(1 / d, nrow = r) * sqrt(n)
  for (k in keep) {
    if (weights[which.max(abs(weights[, k])), k] < 0) {
      scores[, k] <- -scores[, k]
      loadings[, k] <- -loadings[, k]
      weights[, k] <- -weights[, k]
    }
  }
  held <- s$d[seq_len(rank)]^2 / sum(s$d^2)
  share <- data.frame(gradient = paste0("gradient", seq_len(rank)), share = held,
                      cumulative = cumsum(held), stringsAsFactors = FALSE)
  list(site_scores = scores, species_loadings = loadings, weights = weights,
       share = share)
}

#' returnEmbeddingGradients
#'
#' The leading gradients of the species' whole response to the embedding.
#'
#' @details
#' The embedding columns have no meaning one by one. What can be read is the
#' small set of gradients that the species share in their responses to them.
#' This function takes each species' whole response to the embedding at the
#' fitted sites (the embedding columns times the posterior mean of their
#' coefficients, shared part and leftover together) and returns its leading
#' directions across sites: a score for every site, a loading for every
#' species, and the share of the response that each gradient holds.
#'
#' The gradients are defined only up to rotation: any rotation of them
#' describes the same fit. Read them together, as a plane (two gradients) or a
#' block (three), never one by one. [compareEmbeddingGradients()] compares them
#' with covariates in a way that does not depend on the rotation.
#'
#' [returnEmbeddingBlends()] works from the shared part of the responses only.
#' How a response is divided between the shared part and the leftover is poorly
#' determined, and at weak signal a fit can put most of it in the leftovers
#' although the species share gradients; the leading directions of the whole
#' response are far more stable. This is the function to use for
#' interpretation.
#'
#' With `n_gradients = NULL` the number of gradients is the smallest that holds
#' 90 percent of the response, the rule [returnGroupedConvergence()] uses. When
#' that function assesses the embedding's effect, its number of gradients and
#' its agreement between chains are the ones returned here. When it does not,
#' because the effect is small, this function still returns the gradients, with
#' a warning; with one chain the agreement is `NA`.
#'
#' @param fitModel Output from the function runOccJSDM, made with
#' `embeddingVars`
#' @param n_gradients Number of gradients to return. Leave as `NULL` for the
#' 90 percent rule, or give one whole number, at most the rank of the
#' embedding's effect.
#'
#' @return A list.
#' * `n_gradients`: the number of gradients, and `rule`: `"hold 0.90"` when the
#'   rule chose it and `"set"` when you did.
#' * `site_scores`: sites by gradients, with the site names as row names. Each
#'   column has mean 0 and mean square 1 over the fitted sites.
#' * `species_loadings`: species by gradients. A loading is the change in that
#'   species' linear predictor per unit (one standard deviation across the
#'   fitted sites) of the gradient. `site_scores %*% t(species_loadings)` is
#'   the part of the embedding's effect that the gradients hold.
#' * `weights`: embedding columns by gradients. The embedding columns of the
#'   standardised design matrix times `weights` give `site_scores`, so the
#'   weights can score sites that were not fitted.
#' * `share`: a data frame with columns `gradient`, `share` and `cumulative`:
#'   the share of the response each direction holds, for every direction the
#'   response has and not only those returned, so you can see what one more
#'   gradient would add.
#' * `size`: the size of the embedding's effect, the median over species of its
#'   standard deviation across sites.
#' * `agreement`: whether the chains found the same gradients, as in
#'   [returnGroupedConvergence()] but for the number of gradients returned
#'   here: 1 when they did, less the further apart they are. `agreement_pairs`
#'   holds the value for each pair of chains.
#' * `by_chain`: one entry per chain, with that chain's own `site_scores`,
#'   `species_loadings` and `share`. Single gradients cannot be matched from
#'   one chain to another; only the sets can be compared. An entry is `NULL`
#'   when that chain's response has fewer directions than `n_gradients`.
#'
#' @seealso \link[=grouped-prior]{grouped-prior} for the prior that these fits
#' use and for how to interpret the gradients.
#'
#' @md
#' @export
returnEmbeddingGradients <- function(fitModel, n_gradients = NULL) {

  is_emb <- as.logical(fitModel$infos$coef_prior$is_emb)
  if (!is_grouped_fit(fitModel) || !any(is_emb)) {
    stop("returnEmbeddingGradients() needs a fit with embeddingVars")
  }
  if (!is.null(n_gradients) &&
      (!is.numeric(n_gradients) || length(n_gradients) != 1L || is.na(n_gradients) ||
       n_gradients < 1 || n_gradients != round(n_gradients))) {
    stop("n_gradients must be a single positive whole number")
  }

  B_output <- fitModel$results_output$jsdm_output$B_output
  dB <- dim(B_output)
  p <- dB[1]; S <- dB[2]; niter <- dB[3]; nchain <- dB[4]
  X <- fitModel$X_psi
  n <- nrow(X); K <- sum(is_emb)

  # Chain means in one pass, as the convergence check takes them: no part of
  # B_output is copied.
  chain_mean <- blockRhatESS(B_output, p * S, niter, nchain)$chain_mean
  parts <- embedding_effect(chain_mean, X, is_emb, S)
  E <- parts$E
  if (!all(is.finite(E))) {
    stop("The embedding's effect has missing or non-finite values, so its ",
         "gradients cannot be found.")
  }
  pooled <- rowMeans(E, dims = 2)
  B_pooled <- rowMeans(parts$B_mean, dims = 2)

  d <- svd(pooled, nu = 0, nv = 0)$d
  rank <- effect_rank(d, n, S)
  if (rank < 1L) {
    stop("The embedding's effect is zero at every fitted site, so it has no gradients.")
  }
  if (is.null(n_gradients)) {
    r <- which(cumsum(d^2) / sum(d^2) >= GROUPED_LIMITS$hold)[1]
    rule <- paste0("hold ", format(GROUPED_LIMITS$hold, nsmall = 2))
  } else {
    r <- n_gradients
    rule <- "set"
    if (r > rank) {
      stop("n_gradients cannot exceed the rank of the embedding's effect, which is ", rank)
    }
  }
  r <- as.integer(r)

  site_names <- as.character(fitModel$infos$siteNames)
  if (length(site_names) != n) site_names <- as.character(seq_len(n))
  species_names <- as.character(fitModel$infos$speciesNames)
  if (length(species_names) != S) species_names <- as.character(seq_len(S))
  gradient_names <- paste0("gradient", seq_len(r))
  named <- function(g) {
    dimnames(g$site_scores) <- list(site_names, gradient_names)
    dimnames(g$species_loadings) <- list(species_names, gradient_names)
    dimnames(g$weights) <- list(colnames(X)[is_emb], gradient_names)
    g
  }

  out <- named(gradient_parts(pooled, B_pooled, r))

  size <- stats::median(apply(pooled, 2, stats::sd))
  agreement <- NA_real_
  agreement_pairs <- matrix(NA_real_, nchain, nchain)
  if (nchain > 1L) {
    agreement_pairs <- gradient_agreement(E, r)
    agreement <- min(agreement_pairs, na.rm = TRUE)
  }

  by_chain <- vector("list", nchain)
  for (ch in seq_len(nchain)) {
    g <- gradient_parts(matrix(E[, , ch], n, S), matrix(parts$B_mean[, , ch], K, S), r)
    if (!is.null(g)) {
      g <- named(g)
      by_chain[[ch]] <- g[c("site_scores", "species_loadings", "share")]
    }
  }
  short <- which(vapply(by_chain, is.null, NA))

  if (is.finite(size) && size < GROUPED_LIMITS$small_effect) {
    warning("The embedding's effect is small (its size is ", format(size, digits = 2),
            ", under ", GROUPED_LIMITS$small_effect, "), so its gradients are ",
            "mostly noise.", call. = FALSE)
  }
  if (!is.na(agreement) && agreement < GROUPED_LIMITS$gradients) {
    warning("The chains found different gradients: with ", r, " gradient",
            if (r > 1L) "s", " the agreement between chains is ",
            format(agreement, digits = 2), ", under the limit of ",
            GROUPED_LIMITS$gradients, ". See returnGroupedConvergence().",
            call. = FALSE)
  }
  if (length(short) > 0L) {
    warning("The response in chain", if (length(short) > 1L) "s", " ",
            paste(short, collapse = ", "), " has a rank below the ", r,
            " gradients asked for, so by_chain holds nothing for ",
            if (length(short) > 1L) "them" else "it", ".", call. = FALSE)
  }

  list(n_gradients = r, rule = rule,
       site_scores = out$site_scores, species_loadings = out$species_loadings,
       weights = out$weights, share = out$share, size = size,
       agreement = agreement, agreement_pairs = agreement_pairs,
       by_chain = by_chain)
}
