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

# ---- comparison with covariates ---------------------------------------------

# Residual and total sums of squares of each column of `Y` regressed on the
# columns of `X` and an intercept, by QR, so collinear columns of `X` are
# allowed. `terms` is the rank of `X` beyond the intercept.
regress_out <- function(Y, X = NULL) {
  Y <- as.matrix(Y)
  Q <- qr(cbind(rep(1, nrow(Y)), X))
  resid <- qr.resid(Q, Y)
  list(resid = resid, rss = colSums(resid^2),
       tss = colSums(sweep(Y, 2, colMeans(Y))^2), terms = Q$rank - 1L)
}

adjusted_share <- function(share, m, terms) 1 - (1 - share) * (m - 1) / (m - terms - 1)

# Matches the rows of `covariates` to the sites of the gradients, by name, and
# cleans the covariates. Returns `rows` (the positions, among the fitted sites
# and in their order, of the sites used), `cov` (the covariates of those
# sites: numeric columns and factors, none constant) and `sites` (how many
# were used and which were left out, by reason). Says what it leaves out.
match_gradient_sites <- function(fit_sites, covariates, site_col, r) {
  if (!is.data.frame(covariates)) {
    stop("covariates must be a data frame with one row per site")
  }
  if (!is.character(site_col) || length(site_col) != 1L || !site_col %in% names(covariates)) {
    stop("site_col must name the column of covariates that holds the site names")
  }
  sites <- as.character(covariates[[site_col]])
  if (anyNA(sites)) stop("The site column of covariates has missing values")
  if (anyDuplicated(sites) > 0L) {
    stop("These sites appear more than once in covariates: ",
         paste(unique(sites[duplicated(sites)]), collapse = ", "))
  }
  cov <- as.data.frame(covariates)[setdiff(names(covariates), site_col)]
  if (ncol(cov) == 0L) stop("covariates has no column other than the site column")

  count <- function(k, one, many) paste(k, if (k == 1L) one else many)
  not_in_fit <- sites[!sites %in% fit_sites]
  if (length(not_in_fit) > 0L) {
    message(count(length(not_in_fit), "row of covariates is for a site that is not in the fit and is ignored",
                  "rows of covariates are for sites that are not in the fit and are ignored"))
  }
  idx <- match(fit_sites, sites)
  present <- !is.na(idx)
  if (!any(present)) {
    stop("No site in covariates matches a site of the fit. The fit's sites are ",
         "named in the row names of the site scores (infos$siteNames of the fit).")
  }
  not_in_covariates <- fit_sites[!present]
  if (length(not_in_covariates) > 0L) {
    message(count(length(not_in_covariates), "fitted site is not in covariates and is left out",
                  "fitted sites are not in covariates and are left out"))
  }
  rows <- which(present)
  cov <- cov[idx[present], , drop = FALSE]

  for (v in names(cov)) {
    x <- cov[[v]]
    if (is.logical(x)) cov[[v]] <- as.numeric(x)
    else if (is.character(x)) cov[[v]] <- factor(x)
    else if (!is.numeric(x) && !is.factor(x)) {
      stop("Covariate \"", v, "\" must be numeric, logical, a factor or character")
    }
  }
  complete <- stats::complete.cases(cov)
  missing_values <- fit_sites[rows[!complete]]
  if (length(missing_values) > 0L) {
    message(count(length(missing_values), "site has missing covariate values and is left out",
                  "sites have missing covariate values and are left out"))
  }
  rows <- rows[complete]
  cov <- cov[complete, , drop = FALSE]
  m <- length(rows)
  few <- function() {
    stop("There are too few sites to compare: ", m, " with covariate values, for ",
         r, " gradient", if (r > 1L) "s", " and the covariates given")
  }
  if (m < r + 2L) few()
  for (v in names(cov)) {
    if (is.factor(cov[[v]])) cov[[v]] <- droplevels(cov[[v]])
    if (length(unique(cov[[v]])) < 2L) {
      stop("Covariate \"", v, "\" is constant over the sites used, so it cannot be compared")
    }
  }
  rownames(cov) <- NULL
  list(rows = rows, cov = cov, few = few,
       sites = list(n_used = m, not_in_fit = not_in_fit,
                    not_in_covariates = not_in_covariates,
                    missing_values = missing_values))
}

# The covariates as columns to be predicted in check 1: a numeric covariate
# as it is, a two-level one as the indicator of its second level, and one with
# more levels as an indicator per level, named "<covariate>: <level>".
covariate_targets <- function(cov) {
  out <- list()
  for (v in names(cov)) {
    x <- cov[[v]]
    if (!is.factor(x)) {
      out[[v]] <- x
    } else if (nlevels(x) == 2L) {
      out[[v]] <- as.numeric(x == levels(x)[2])
    } else {
      for (lev in levels(x)) out[[paste0(v, ": ", lev)]] <- as.numeric(x == lev)
    }
  }
  out
}

#' compareEmbeddingGradients
#'
#' Compares the gradients of the embedding response with covariates known for
#' the same sites.
#'
#' @details
#' The gradients from [returnEmbeddingGradients()] are defined only up to
#' rotation, so correlating one gradient with one covariate gives an answer
#' that depends on an arbitrary choice. This function compares the gradients
#' as a set, in three ways.
#'
#' * Covariate from gradients: how well do the site scores of the first
#'   gradients, taken together, predict each covariate?
#' * Gradients from covariates: how much of the gradients do all the covariates
#'   together account for, and how much does each account for alone and
#'   uniquely? "Alone" is with that covariate only; "unique" is what is lost
#'   when it is taken away and the others are left, so it is what that
#'   covariate accounts for that the others cannot.
#' * Unmeasured gradients: the leading directions of what the covariates leave
#'   unexplained. This is habitat variation that the embedding records and the
#'   species respond to, but that the covariates do not describe.
#'
#' The main use is with covariates that were not in the model, to suggest what
#' the gradients may stand for. A covariate that was in the model can still
#' match the gradients, because the embedding columns are correlated with it:
#' that shows overlap between the covariate and the embedding, not a new
#' driver. To ask what a new covariate adds beyond those in the model, give
#' both and read the new one's `unique` value.
#'
#' Sites are matched by name, never by row order. Rows of `covariates` for
#' sites that are not in the fit, fitted sites that are not in `covariates`,
#' and sites with a missing value in any covariate are left out, with a
#' message saying how many.
#'
#' The function gives no test and no interval. Sites near one another are not
#' independent, so a test would overstate the evidence, and a match is a
#' hypothesis about what drives the species, not evidence of it: anything
#' correlated with the covariate across your sites would match as well.
#'
#' @param gradients The list returned by [returnEmbeddingGradients()]
#' @param covariates A data frame with one row per site: a column of site
#' names, and one column for each covariate. Numeric and logical columns are
#' used as numbers; factor and character columns are categorical.
#' @param site_col Name of the column of `covariates` that holds the site
#' names, which must be those of the fit (`infos$siteNames`)
#' @param squared If `TRUE`, each numeric covariate also enters the second and
#' third comparisons as its square, for species that respond to a covariate
#' with a hump. The square counts as part of that covariate.
#'
#' @return A list.
#' * `covariate_from_gradients`: a data frame with columns `covariate`, `chain`
#'   (`"pooled"`, or the number of a chain), `n_gradients`, `r_squared` and
#'   `adj_r_squared`: how well the first `n_gradients` site scores together
#'   predict the covariate. A categorical covariate with two levels is its
#'   second level against its first; one with more levels has a row set for
#'   each level, named `"covariate: level"`.
#' * `gradients_from_covariates`: a list of three data frames. `overall` has
#'   columns `chain`, `explained` and `adj_explained`: the share of the
#'   gradients that all the covariates together account for. `by_covariate`
#'   has columns `covariate`, `alone` and `unique`. `by_gradient` has columns
#'   `gradient`, `share`, `r_squared` and `adj_r_squared`, where `share` is the
#'   gradient's share of the response at the sites used; this last table
#'   depends on the rotation and is for description only.
#' * `unmeasured`: a list with `site_scores` (the sites used by unmeasured
#'   gradients), `species_loadings` and `share` (a data frame with columns
#'   `gradient`, `share` and `cumulative`, as shares of the gradients at the
#'   sites used). The shares sum to what the covariates leave unexplained.
#'   There can be fewer unmeasured gradients than gradients, and none when the
#'   covariates account for everything.
#' * `sites`: `n_used`, and the names of the sites left out: `not_in_fit`,
#'   `not_in_covariates` and `missing_values`.
#'
#' @seealso \link[=grouped-prior]{grouped-prior} for how to interpret the
#' gradients.
#'
#' @md
#' @export
compareEmbeddingGradients <- function(gradients, covariates, site_col = "Site",
                                      squared = FALSE) {

  if (!is.list(gradients) ||
      !all(c("n_gradients", "site_scores", "species_loadings") %in% names(gradients)) ||
      is.null(rownames(gradients$site_scores))) {
    stop("gradients must be the list returned by returnEmbeddingGradients()")
  }
  Z <- as.matrix(gradients$site_scores)
  r <- ncol(Z)
  matched <- match_gradient_sites(rownames(Z), covariates, site_col, r)
  rows <- matched$rows
  cov <- matched$cov
  m <- length(rows)

  # The covariates as regressors: `F_cov` without the intercept (regress_out()
  # adds it), and `owner`, the covariate each column belongs to.
  design <- stats::model.matrix(~ ., data = cov)
  F_cov <- design[, -1, drop = FALSE]
  owner <- attr(design, "assign")[-1]
  if (isTRUE(squared)) {
    for (j in seq_along(cov)) {
      x <- cov[[j]]
      if (is.numeric(x) && length(unique(x)) > 2L) {
        F_cov <- cbind(F_cov, (x - mean(x))^2)
        owner <- c(owner, j)
      }
    }
  }
  terms <- qr(cbind(1, F_cov))$rank - 1L
  if (m < terms + 3L) matched$few()

  # The gradient sets compared: pooled, then every chain that has one.
  sets <- list(pooled = list(Z = Z[rows, , drop = FALSE],
                             L = as.matrix(gradients$species_loadings)))
  for (ch in seq_along(gradients$by_chain)) {
    g <- gradients$by_chain[[ch]]
    if (!is.null(g)) {
      sets[[as.character(ch)]] <- list(Z = as.matrix(g$site_scores)[rows, , drop = FALSE],
                                       L = as.matrix(g$species_loadings))
    }
  }

  # Check 1: each covariate from the first k gradients together.
  targets <- covariate_targets(cov)
  Y <- do.call(cbind, targets)
  check1 <- list()
  for (set in names(sets)) for (k in seq_len(r)) {
    fit <- regress_out(Y, sets[[set]]$Z[, seq_len(k), drop = FALSE])
    r2 <- 1 - fit$rss / fit$tss
    check1[[length(check1) + 1L]] <- data.frame(
      covariate = names(targets), chain = set, n_gradients = k,
      r_squared = unname(r2), adj_r_squared = unname(adjusted_share(r2, m, fit$terms)),
      stringsAsFactors = FALSE)
  }
  check1 <- do.call(rbind, check1)
  check1 <- check1[order(match(check1$covariate, names(targets))), , drop = FALSE]
  rownames(check1) <- NULL

  # Check 2: the gradients from the covariates. G = Z (L'L)^(1/2) has the sums
  # of squares of the effect the gradients hold, Z L', whatever the rotation
  # of (Z, L): Z L' = G O' with O = L (L'L)^(-1/2), whose columns are
  # orthonormal.
  weighted <- function(set) {
    e <- eigen(crossprod(set$L), symmetric = TRUE)
    root <- sqrt(pmax(e$values, 0))
    list(G = set$Z %*% e$vectors %*% diag(root, nrow = r) %*% t(e$vectors),
         inv_root = e$vectors %*% diag(ifelse(root > 0, 1 / root, 0), nrow = r) %*% t(e$vectors))
  }
  explained_by <- function(G, cols) {
    fit <- regress_out(G, F_cov[, cols, drop = FALSE])
    1 - sum(fit$rss) / sum(fit$tss)
  }
  all_cols <- seq_len(ncol(F_cov))
  overall <- do.call(rbind, lapply(names(sets), function(set) {
    share <- explained_by(weighted(sets[[set]])$G, all_cols)
    data.frame(chain = set, explained = share,
               adj_explained = adjusted_share(share, m, terms), stringsAsFactors = FALSE)
  }))

  pooled <- weighted(sets$pooled)
  G <- pooled$G
  explained <- overall$explained[1]
  by_covariate <- do.call(rbind, lapply(seq_along(cov), function(j) {
    mine <- which(owner == j)
    unique_share <- explained - explained_by(G, setdiff(all_cols, mine))
    data.frame(covariate = names(cov)[j], alone = explained_by(G, mine),
               unique = max(unique_share, 0), stringsAsFactors = FALSE)
  }))

  # Gradient by gradient, for the gradients as given: this table alone depends
  # on their rotation.
  L <- sets$pooled$L
  each <- regress_out(sweep(sets$pooled$Z, 2, sqrt(colSums(L^2)), `*`), F_cov)
  r2 <- 1 - each$rss / each$tss
  by_gradient <- data.frame(gradient = paste0("gradient", seq_len(r)),
                            share = unname(each$tss / sum(each$tss)),
                            r_squared = unname(r2),
                            adj_r_squared = unname(adjusted_share(r2, m, terms)),
                            stringsAsFactors = FALSE)

  # Check 3: the leading directions of what the covariates leave of G.
  left <- regress_out(G, F_cov)
  s <- svd(left$resid)
  held <- s$d^2 / sum(left$tss)
  keep <- which(held > 1e-10)
  u_scores <- sqrt(m) * s$u[, keep, drop = FALSE]
  u_loadings <- L %*% pooled$inv_root %*% s$v[, keep, drop = FALSE] %*%
    diag(s$d[keep], nrow = length(keep)) / sqrt(m)
  for (k in seq_along(keep)) {
    if (u_loadings[which.max(abs(u_loadings[, k])), k] < 0) {
      u_scores[, k] <- -u_scores[, k]
      u_loadings[, k] <- -u_loadings[, k]
    }
  }
  u_names <- if (length(keep) > 0L) paste0("unmeasured", seq_along(keep)) else character(0)
  dimnames(u_scores) <- list(rownames(Z)[rows], u_names)
  dimnames(u_loadings) <- list(rownames(L), u_names)
  unmeasured <- list(site_scores = u_scores, species_loadings = u_loadings,
                     share = data.frame(gradient = u_names, share = held[keep],
                                        cumulative = cumsum(held[keep]),
                                        stringsAsFactors = FALSE))

  list(covariate_from_gradients = check1,
       gradients_from_covariates = list(overall = overall, by_covariate = by_covariate,
                                        by_gradient = by_gradient),
       unmeasured = unmeasured,
       sites = matched$sites)
}
