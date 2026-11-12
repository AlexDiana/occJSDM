# The occupancy coefficient prior. b_prior = "fixed" is the original Normal(0,
# sigma_b^2) on every coefficient. b_prior = "grouped" lets a block of
# embedding covariates (embeddingVars) share a scale, sigma_emb, apart from the
# measured covariates. Naming embeddingVars is enough to switch it on; naming
# them with b_prior = "fixed" is an error. Embedding columns that are in
# occCovariates but not in embeddingVars are ordinary measured covariates.
# `columns` has one row per design column: column, variable, kind.
read_coef_prior <- function(priors, params, S, gt, columns) {
  emb_vars <- get_param(params, "embeddingVars", NULL)
  has_emb <- length(emb_vars) > 0L
  type <- get_param(priors, "b_prior", NULL)
  if (is.null(type)) {
    type <- if (has_emb) "grouped" else "fixed"
  } else if (!is.character(type) || length(type) != 1L || is.na(type) ||
             !type %in% c("fixed", "grouped")) {
    stop("b_prior must be 'fixed' or 'grouped'")
  }
  type <- unname(type)
  if (type == "fixed") {
    if (has_emb) stop("embeddingVars cannot be used with b_prior = 'fixed'")
    for (nm in c("sigma_b_scale", "sigma_emb_scale", "axis_scale"))
      if (nm %in% names(priors)) stop(nm, " requires the grouped prior")
    return(list(type = "fixed"))
  }
  if (S < 2) stop("the grouped prior needs more than one species")
  if (has_emb) {
    unknown <- setdiff(emb_vars, columns$variable)
    if (length(unknown))
      stop("embeddingVars must name occupancy covariates: ",
           paste(unknown, collapse = ", "))
    bad_kind <- unique(columns$variable[columns$variable %in% emb_vars &
                                          columns$kind != "numeric"])
    if (length(bad_kind))
      stop("embeddingVars must be numeric covariates without splines: ",
           paste(bad_kind, collapse = ", "))
  } else {
    for (nm in c("sigma_emb_scale", "axis_scale"))
      if (nm %in% names(priors)) stop(nm, " requires embeddingVars")
  }
  is_emb <- columns$variable %in% emb_vars
  p_emb <- sum(is_emb)
  scale <- function(name, default) {
    value <- get_param(priors, name, default)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be a finite positive number")
    unname(value)
  }
  sigma_b_scale <- scale("sigma_b_scale", 1)
  sigma_emb_scale <- NA_real_
  axis_scale <- NA_real_
  if (has_emb) {
    sigma_emb_scale <- scale("sigma_emb_scale", 1)
    # 0.7 / K, not 1 / sqrt(K): the half-Cauchy tails of K axis scales make
    # the embedding's prior contribution grow in proportion to K times this
    # scale (dev/simstudy/grouped-prior-checks-20261111).
    axis_scale <- scale("axis_scale", 0.7 / p_emb)
  }
  c_base <- rep(1, length(is_emb))
  if (gt > 0) c_base[is_emb] <- 1 / sqrt(gt)
  list(type = "grouped", is_emb = is_emb, sigma_b_scale = sigma_b_scale,
       sigma_emb_scale = sigma_emb_scale, axis_scale = axis_scale,
       c_base = c_base)
}

# One row per column of the occupancy design matrix X: its name, the covariate
# it came from and that covariate's kind ("numeric", "spline" or "factor").
# process_covariates() does not report which columns a covariate produced, so
# the blocks are rebuilt from the record create_covariates_matrix() keeps
# (list_matrix) and its standardised data frame X0. The blocks are in names_df
# order, as covariate_response_grid() also relies on: a numeric covariate is
# one column, or its spline basis when bs_info holds one, and a factor is one
# column per contrast. A numeric covariate that splineVars selects is a spline
# even when too few distinct values made it fall back to a single linear
# column, so that its kind does not depend on the number of sites.
occ_design_columns <- function(X, list_matrix, X0) {
  if (is.null(list_matrix)) {
    variable <- kind <- character(0)
  } else {
    blocks <- lapply(list_matrix$names_df, function(name) {
      if (!list_matrix$is_numeric[[name]])
        return(list(n = ncol(stats::contrasts(X0[[name]])), kind = "factor"))
      spline <- list_matrix$bs_info[[name]]
      list(n = if (is.null(spline)) 1L else
             length(spline$knots) + spline$degree + spline$intercept,
           kind = if (name %in% list_matrix$target_spline_vars) "spline" else "numeric")
    })
    n <- vapply(blocks, function(b) as.integer(b$n), integer(1))
    variable <- rep(list_matrix$names_df, n)
    kind <- rep(vapply(blocks, function(b) b$kind, character(1)), n)
  }
  if (length(variable) != ncol(X))
    stop("the occupancy design columns could not be matched to their covariates")
  data.frame(column = as.character(colnames(X)), variable = variable, kind = kind)
}

# Embedding covariates are standardised with the other numeric occupancy
# covariates, so one that does not vary across sites cannot be used. Checked
# before the design matrix is built, on the values it is built from: the first
# row of each site. Names that are not numeric occupancy covariates are left
# for read_coef_prior() to reject.
check_embedding_varies <- function(data_info, occCovariates, params) {
  emb_vars <- intersect(get_param(params, "embeddingVars", NULL), occCovariates)
  first <- !duplicated(data_info$Site)
  constant <- emb_vars[vapply(emb_vars, function(v) {
    x <- data_info[[v]]
    is.numeric(x) && length(unique(x[first])) < 2L
  }, logical(1))]
  if (length(constant))
    stop("embeddingVars columns must vary across sites: ",
         paste(constant, collapse = ", "))
  invisible(NULL)
}

# "measured" or "embedding" for every design column, named by column.
coef_prior_groups <- function(coef_prior, columns) {
  is_emb <- if (coef_prior$type == "grouped") coef_prior$is_emb else
    rep(FALSE, nrow(columns))
  stats::setNames(c("measured", "embedding")[is_emb + 1L], columns$column)
}

# Starting values of the embedding scale and the covariate scales. A covariate
# scale is exactly 1 on every measured column, at the start and always:
# sample_coef_scales() only ever redraws the embedding columns' scales.
coef_scales_start <- function(coef_prior, p) {
  sigma_emb <- NA_real_
  lambda <- rep(1, p)
  if (coef_prior$type == "grouped") {
    sigma_emb <- coef_prior$sigma_emb_scale
    lambda[coef_prior$is_emb] <- coef_prior$axis_scale
  }
  list(sigma_emb = sigma_emb, lambda = lambda)
}

# Starting value of C, the responses of the latent traits. `C` is the chain's
# usual start, ones on and above the diagonal, and `lambda` the starting
# covariate scales. Under the grouped prior the embedding columns are put at
# their prior scale, c_base times the starting covariate scale. The first
# scale draw treats the starting C as a draw from its prior: left at 1, the
# term sum(C[, k]^2) / c_base[k]^2 is about gt * gt and throws every axis
# scale to about gt / sqrt(S + g_b + gt) in the first sweep, whatever
# axis_scale is, and the chains take thousands of sweeps to come back.
# Measured columns, and every column under the fixed prior, are returned
# untouched. No random numbers are drawn.
coef_C_start <- function(coef_prior, C, lambda) {
  if (!identical(coef_prior$type, "grouped") || !any(coef_prior$is_emb)) return(C)
  emb <- which(coef_prior$is_emb)
  C[, emb] <- C[, emb, drop = FALSE] *
    rep(coef_prior$c_base[emb] * lambda[emb], each = nrow(C))
  C
}

# The prior standard deviations of B must be p finite positive numbers before
# they reach the joint draw and the range step: a NULL aborts R inside the C++
# draw, and the R range functions would recycle a wrong length silently.
check_coef_sds <- function(sd_b, p) {
  if (!is.numeric(sd_b) || length(sd_b) != p || !all(is.finite(sd_b)) ||
      any(sd_b <= 0))
    stop("the prior standard deviations of the occupancy coefficients must be ",
         p, " finite positive numbers")
  invisible(NULL)
}

# Per-covariate prior standard deviations from the current scale values.
# Fixed: sigma_b on every coefficient, no covariate scales. Grouped: the
# covariate scales lambda multiply the measured or embedding scale.
coef_prior_sds <- function(coef_prior, sigma_b, sigma_emb, lambda, g_b, p) {
  if (coef_prior$type == "fixed")
    return(list(b = rep(sigma_b, p), g = matrix(1, g_b, p), c = rep(1, p)))
  list(b = ifelse(coef_prior$is_emb, sigma_emb, sigma_b) * lambda,
       g = matrix(lambda, g_b, p, byrow = TRUE),
       c = coef_prior$c_base * lambda)
}

# One draw of a standard deviation with a half-Cauchy(scale) prior, given `n`
# Normal(0, sd^2) quantities whose squares sum to `ss`. Same auxiliary-variable
# step as sample_spatial_sd_half_cauchy(): refresh the auxiliary given the
# current value `sigma`, then draw the new variance.
sample_scale_half_cauchy <- function(ss, n, sigma, scale) {
  aux <- rinvgamma_cpp(1, 1 / sigma^2 + 1 / scale^2)
  sqrt(rinvgamma_cpp((n + 1) / 2, ss / 2 + 1 / aux))
}

# Update the scales of the grouped prior. R is the species by covariates
# matrix of leftovers, B minus its trait regression; G and C have one column
# per covariate. Drawn in this order: sigma_b from the measured columns,
# sigma_emb from the embedding columns with each divided by its covariate
# scale, then each embedding column's covariate scale given the new sigma_emb.
# A covariate scale multiplies the prior standard deviation of its column of
# leftovers, of G and of C, so all three inform it. With no embedding columns
# sigma_emb and lambda are returned as given.
sample_coef_scales <- function(R, G, C, coef_prior, sigma_b, sigma_emb, lambda) {
  S <- nrow(R)
  is_emb <- coef_prior$is_emb
  sigma_b <- sample_scale_half_cauchy(sum(R[, !is_emb]^2), S * sum(!is_emb),
                                      sigma_b, coef_prior$sigma_b_scale)
  if (any(is_emb)) {
    emb <- which(is_emb)
    ss_emb <- sum(sweep(R[, emb, drop = FALSE], 2, lambda[emb], "/")^2)
    sigma_emb <- sample_scale_half_cauchy(ss_emb, S * length(emb), sigma_emb,
                                          coef_prior$sigma_emb_scale)
    for (k in emb) {
      ss_k <- sum(R[, k]^2) / sigma_emb^2 + sum(G[, k]^2) +
        sum(C[, k]^2) / coef_prior$c_base[k]^2
      lambda[k] <- sample_scale_half_cauchy(ss_k, S + nrow(G) + nrow(C),
                                            lambda[k], coef_prior$axis_scale)
    }
  }
  list(sigma_b = sigma_b, sigma_emb = sigma_emb, lambda = lambda)
}

#' Grouped prior for occupancy coefficients (experimental)
#'
#' An experimental prior that lets [runOccJSDM()] use an embedding, a block of
#' tens of numeric columns per site from a satellite foundation model, as
#' occupancy covariates beside a few measured covariates.
#'
#' @section What it is for:
#' Satellite embeddings describe the habitat at a site in many numeric columns
#' that have no individual meaning. With the usual prior, which gives every
#' occupancy coefficient the same Normal prior, many columns would swamp the
#' few covariates you measured. The grouped prior gives the measured covariates
#' and the embedding columns scales of their own, learned from the data, so the
#' embedding can act as a predictor for maps, or as a control that makes sites
#' comparable in habitat, without drowning the measured covariates.
#'
#' @section How to switch it on:
#' Name the embedding columns in `listParams$embeddingVars`, a character
#' vector of names taken from `occCovariates`. That alone switches the prior
#' on. Set `listPriors$b_prior = "grouped"` only if you want the prior with no
#' embedding. Setting `b_prior = "fixed"` together with `embeddingVars` is an
#' error. Columns that are in `occCovariates` but not in `embeddingVars` are
#' treated as ordinary measured covariates.
#'
#' Each embedding variable must be numeric, must vary across sites and must
#' not be in `splineVars` (so `splineVars = TRUE` rejects any
#' `embeddingVars`). The fit needs more than one species.
#'
#' @section What else changes:
#' The prior on the measured covariates becomes a learned measured-group
#' scale, with a half-Cauchy prior. Every covariate also gets a
#' community-average response, the response of the average species.
#'
#' A species' coefficient on a covariate is a centre plus a leftover. The
#' centre is the community average, plus the effects of species traits, plus
#' the effects of latent traits through blends. A latent trait is an unmeasured
#' species trait that the model estimates, and a blend is a weighted
#' combination of the embedding columns that the latent traits respond to. The
#' leftover is what the centre does not explain. Its standard deviation is the
#' measured-group scale for a measured covariate, and the embedding leftover
#' scale times that column's axis scale for an embedding column. Each embedding
#' column's axis scale also multiplies that column's community-average, trait
#' and blend terms, so an axis that no species responds to shrinks out of the
#' model.
#'
#' @section Settings:
#' Set these in `listPriors`. Each is a half-Cauchy scale and must be a finite
#' positive number.
#' * `sigma_b_scale`: scale of the measured-group scale. Default 1.
#' * `sigma_emb_scale`: scale of the embedding leftover scale. Default 1.
#' * `axis_scale`: scale of each embedding column's axis scale. Default
#'   `0.7 / K`, where K is the number of embedding columns.
#'
#' The number of blends is `listParams$n_lattrait`, the number of latent
#' traits.
#'
#' @section What the fit holds:
#' In `results_output$jsdm_output`, `lambda_output` holds the axis scales and
#' `sigmaemb_output` the embedding leftover scale; both are present only when
#' the fit has embedding columns. `sigmab_output` holds the measured-group
#' scale. When every occupancy covariate is an embedding column there is
#' nothing for that scale to learn from, so `sigmab_output` is then a draw from
#' its prior and should be ignored. `G_output` gains a leading row, the
#' community-average response. `infos$coef_prior` records the prior and its
#' settings, and `infos$covariate_group` labels each occupancy design column
#' `"measured"` or `"embedding"`.
#'
#' Four functions summarise these:
#' * [returnCovariateScales()]: posterior summaries of the measured-group,
#'   embedding leftover and axis scales.
#' * [returnCommunityResponse()]: the community-average response of every
#'   occupancy design column.
#' * [returnEmbeddingBlends()]: the blends, as weights on the embedding
#'   columns and as scores for each site. Blends are defined only up to
#'   rotation, so do not interpret them one by one.
#' * [returnEmbeddingShare()]: for each species, the share of the variation
#'   (across sites) in the covariates' contribution to its linear predictor
#'   that comes from the embedding columns.
#'
#' @section Reading a management covariate:
#' A management covariate (for example protected-area status or time since
#' clearcut) fitted alongside the embedding has a coefficient that is a direct
#' effect: the effect of management that does not pass through habitat the
#' satellite can see. It is not the total effect, because the embedding is
#' measured after the management and partly records its results. For the total
#' effect, fit again without the embedding and without other covariates that
#' management can change, and compare the two fits on the probability scale.
#'
#' @section Status:
#' This prior is experimental and has not been through its validation study.
#' The default `axis_scale` was changed on 12 November 2026, from 1 over the
#' square root of the number of embedding columns to `0.7 / K`, after a prior
#' check found the old default far too wide. The check looks at the prior
#' alone, before any data. What it measures is the spread across sites, as a
#' standard deviation, of the embedding's contribution to a species' linear
#' predictor, on the logit scale. The design target for the prior median of
#' that spread is 1 to 2. Under the present default the prior median is about
#' 1.4 with 64 embedding columns and about 1.5 with 128, both with 8 latent
#' traits. The upper tail is wide: the 95% point is about 28 at both column
#' counts, so one prior draw in twenty gives a spread larger than that. The
#' validation study is to examine that tail.
#'
#' @examples
#' \dontrun{
#' # Simulated data with 2 measured covariates and an 8-column embedding
#' set.seed(1)
#' n <- 60; S <- 12
#' sim <- simulateOccJSDMData(
#'   list_datasettings = list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L,
#'                       K = rep(1L, n), ncov_psi = 2L, ncov_theta = 1L),
#'   list_params = list(),
#'   list_jsdmParams = list(gt = 2L, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5,
#'                     sigma_ts = 0.5, sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
#'                     tau = rep(1, S), useSpatField = FALSE,
#'                     embedding = list(n_emb = 8L)),
#'   model = "binary")
#' dat <- sim$data_list
#' measured <- c("X_psi.EnvCov.1", "X_psi.EnvCov.2")
#' embedding <- dat$embeddingVars          # "X_emb.1" to "X_emb.8"
#' fit <- runOccJSDM(dat,
#'   listParams = list(n_factors = 1L, n_lattrait = 2L,
#'                     embeddingVars = embedding),
#'   occCovariates = c(measured, embedding), spatCovariates = NULL,
#'   MCMCparams = list(nchain = 2, nburn = 500, niter = 500, nthin = 1))
#' returnCovariateScales(fit)
#' returnCommunityResponse(fit)
#' }
#'
#' @md
#' @name grouped-prior
#' @aliases groupedPrior embeddingVars b_prior
NULL
