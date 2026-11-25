# The occupancy coefficient prior. b_prior = "fixed" is the original Normal(0,
# sigma_b^2) on every coefficient. b_prior = "grouped" lets a block of
# embedding covariates (embeddingVars) share a scale, sigma_emb, apart from the
# measured covariates. Naming embeddingVars is enough to switch it on; naming
# them with b_prior = "fixed" is an error. Embedding columns that are in
# occCovariates but not in embeddingVars are ordinary measured covariates.
# `columns` has one row per design column: column, variable, kind.
read_coef_prior <- function(priors, params, S, gt, columns, model = "binary") {
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
    for (nm in c("sigma_b_scale", "sigma_emb_scale", "axis_scale", "slab_scale"))
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
    for (nm in c("sigma_emb_scale", "axis_scale", "slab_scale"))
      if (nm %in% names(priors)) stop(nm, " requires embeddingVars")
  }
  is_emb <- columns$variable %in% emb_vars
  p_emb <- sum(is_emb)
  read_scale <- function(name, default) {
    value <- get_param(priors, name, default)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be a finite positive number")
    unname(value)
  }
  sigma_b_scale <- read_scale("sigma_b_scale", 1)
  sigma_emb_scale <- NA_real_
  axis_scale <- NA_real_
  slab_scale <- NA_real_
  if (has_emb) {
    sigma_emb_scale <- read_scale("sigma_emb_scale", 1)
    # 0.7 / K, not 1 / sqrt(K): the half-Cauchy tails of K axis scales make
    # the embedding's prior contribution grow in proportion to K times this
    # scale (dev/simstudy/grouped-prior-checks-20261111).
    axis_scale <- read_scale("axis_scale", 0.7 / p_emb)
    # The cap on the embedding scales (see coef_capped()). Inf is no cap. The
    # model type is accepted for the default that a later task sets by model.
    slab_scale <- get_param(priors, "slab_scale", Inf)
    if (!is.numeric(slab_scale) || length(slab_scale) != 1L || is.na(slab_scale) ||
        slab_scale <= 0)
      stop("slab_scale must be a positive number or Inf")
    slab_scale <- unname(slab_scale)
  }
  c_base <- rep(1, length(is_emb))
  if (gt > 0) c_base[is_emb] <- 1 / sqrt(gt)
  list(type = "grouped", is_emb = is_emb, sigma_b_scale = sigma_b_scale,
       sigma_emb_scale = sigma_emb_scale, axis_scale = axis_scale,
       slab_scale = slab_scale, c_base = c_base)
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
check_embedding_varies <- function(data_info, occCovariates, params,
                                   group_by_col = "Site") {
  emb_vars <- intersect(get_param(params, "embeddingVars", NULL), occCovariates)
  # The rows process_covariates() will use: the first row of each group, or
  # every row as given when the data already have one row per site (NULL).
  first <- if (is.null(group_by_col)) rep(TRUE, nrow(data_info)) else
    !duplicated(data_info[[group_by_col]])
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
  if (is.infinite(coef_slab(coef_prior))) {
    C[, emb] <- C[, emb, drop = FALSE] *
      rep(coef_prior$c_base[emb] * lambda[emb], each = nrow(C))
  } else {
    axis <- coef_effective_scales(coef_prior, NA_real_, lambda)$axis
    C[, emb] <- C[, emb, drop = FALSE] *
      rep(coef_prior$c_base[emb] * axis[emb], each = nrow(C))
  }
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
  if (is.infinite(coef_slab(coef_prior)))
    return(list(b = ifelse(coef_prior$is_emb, sigma_emb, sigma_b) * lambda,
                g = matrix(lambda, g_b, p, byrow = TRUE),
                c = coef_prior$c_base * lambda))
  eff <- coef_effective_scales(coef_prior, sigma_emb, lambda)
  list(b = ifelse(coef_prior$is_emb, eff$leftover, sigma_b * lambda),
       g = matrix(eff$axis, g_b, p, byrow = TRUE),
       c = coef_prior$c_base * eff$axis)
}

# The cap on the embedding scales of a grouped prior record: its slab_scale,
# or Inf (no cap) when the record has none, as the hand-built records of the
# older tests do, or when it is NA (no embedding columns).
coef_slab <- function(coef_prior) {
  slab <- coef_prior$slab_scale
  if (is.null(slab) || is.na(slab)) Inf else slab
}

# A standard deviation `s` capped at `slab`: the value whose precision is
# 1 / s^2 + 1 / slab^2, so s * slab / sqrt(s^2 + slab^2). It is computed from
# the smaller of the two times a factor between 0.7 and 1, so that nothing
# overflows or underflows however extreme `s` is. With slab = Inf it is `s`
# itself.
coef_capped <- function(s, slab) {
  if (is.infinite(slab)) return(s)
  small <- pmin(s, slab)
  big <- pmax(s, slab)
  small / sqrt(1 + (small / big)^2)
}

# The scales the prior actually uses on the embedding columns, once the cap is
# applied. `axis` is the effective axis scale: lambda on measured columns, the
# capped lambda on embedding columns. `leftover` is the effective standard
# deviation of the leftovers of B: NA on measured columns, the capped
# sigma_emb * lambda on embedding columns. Both have one element per column.
coef_effective_scales <- function(coef_prior, sigma_emb, lambda) {
  slab <- coef_slab(coef_prior)
  is_emb <- coef_prior$is_emb
  axis <- lambda
  leftover <- rep(NA_real_, length(lambda))
  if (any(is_emb)) {
    axis[is_emb] <- coef_capped(lambda[is_emb], slab)
    leftover[is_emb] <- coef_capped(sigma_emb * lambda[is_emb], slab)
  }
  list(axis = axis, leftover = leftover)
}

# One draw of a standard deviation with a half-Cauchy(scale) prior, given `n`
# Normal(0, sd^2) quantities whose squares sum to `ss`. Same auxiliary-variable
# step as sample_spatial_sd_half_cauchy(): refresh the auxiliary given the
# current value `sigma`, then draw the new variance.
sample_scale_half_cauchy <- function(ss, n, sigma, scale) {
  aux <- rinvgamma_cpp(1, 1 / sigma^2 + 1 / scale^2)
  sqrt(rinvgamma_cpp((n + 1) / 2, ss / 2 + 1 / aux))
}

# One univariate slice draw (Neal 2003, "Slice sampling", Annals of
# Statistics 31: 705-767; stepping out as in his Figure 3, shrinkage as in his
# Figure 5) from the density whose log is `f`, started at `u0`, which is
# returned unchanged only in the degenerate case below. The draw leaves that
# density exactly invariant and needs no tuning. Used on the log scale of the
# embedding scales under a cap (sample_coef_scales()).
#
# The slice is every u with f(u) above a level drawn below f(u0). An interval
# of width `w` is placed at random around u0, and its ends step out by `w`
# while they are inside the slice, for at most `max_steps` steps in all, split
# at random between the two ends as Neal's Figure 3 does (the random split is
# what keeps a limited stepping out exact). A point is then drawn uniformly
# from the interval, and on each rejection the interval shrinks to that point
# on the side away from u0. A value of f that is not finite (an overflow far
# out on the log scale, or a density of zero) counts as outside the slice.
# If the level rounds to f(u0) itself, which can happen only when |f(u0)| is
# huge, no point may ever be accepted; after 200 rejections the interval has
# shrunk far below the precision of a double and u0 is returned. Random
# numbers come from R's generator, all of them uniforms: the level, the
# placement, the split, then one per proposal.
slice_sample_log <- function(f, u0, w = 1, max_steps = 50L) {
  f0 <- f(u0)
  if (!is.finite(f0))
    stop("a slice draw must start at a point where the log density is finite")
  level <- f0 + log(stats::runif(1))
  inside <- function(u) {
    fu <- f(u)
    is.finite(fu) && fu > level
  }
  left <- u0 - w * stats::runif(1)
  right <- left + w
  steps_left <- floor(max_steps * stats::runif(1))
  steps_right <- max_steps - 1 - steps_left
  while (steps_left > 0 && inside(left)) {
    left <- left - w
    steps_left <- steps_left - 1
  }
  while (steps_right > 0 && inside(right)) {
    right <- right + w
    steps_right <- steps_right - 1
  }
  for (i in seq_len(200L)) {
    u <- left + stats::runif(1) * (right - left)
    if (inside(u)) return(u)
    if (u < u0) left <- u else right <- u
  }
  u0
}

# log(exp(x) + exp(y)), elementwise, without overflow: the log of a sum of two
# precisions given their logs. With y = -Inf (no cap) it is x.
log_add_exp <- function(x, y) {
  big <- pmax(x, y)
  big + log1p(exp(pmin(x, y) - big))
}

# log(1 + exp(x)), elementwise, without overflow.
log1p_exp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# The log of the conditional density of u = log(lambda_k), an embedding
# column's axis scale, under a cap `slab`, up to a constant: f_k of the
# design's "The draw of an axis scale" (dev/superpowers/specs/
# 2026-11-25-axis-slab-design.md, section "Mathematics"):
#
#   f_k(u) = u - log(1 + e^(2u) / a^2) + (m_gc / 2) log P_g - P_g ss_gc / 2
#            + (S / 2) log P_e - P_e ss_r / 2,
#   P_g = e^(-2u) + c^-2,  P_e = e^(-2u) / sigma_emb^2 + c^-2,
#
# with a = axis_scale and c = slab. P_g is the precision of column k of G and
# (divided by c_base[k]^2) of C, and P_e that of column k of the leftovers.
# ss_gc is the sum of squares of the g_b rows of G[, k] plus sum(C[, k]^2) /
# c_base[k]^2, and m_gc = g_b + q is their number; ss_r is the sum of squares
# of the S leftovers of column k. The leading u is the Jacobian of the log.
# The parts c^-2 ss_gc / 2 and c^-2 ss_r / 2 of the two quadratic terms do
# not depend on u and are dropped, so the quadratic terms are
# e^(-2u) ss_gc / 2 and e^(-2u) ss_r / (2 sigma_emb^2); keeping them would
# only add a large constant that costs precision when a sum of squares is
# large. The log precisions and the quadratic terms are computed through
# logs, so that every term stays finite for u from about -350 to 350; beyond
# that a term overflows to -Inf, which a slice draw treats as outside the
# slice. A sum of squares of zero gives a term of exactly zero. With
# slab = Inf this is the uncapped density.
axis_scale_logdens <- function(u, ss_gc, m_gc, ss_r, S, sigma_emb, axis_scale, slab) {
  log_cap <- -2 * log(slab)
  log_pg <- log_add_exp(-2 * u, log_cap)
  log_pe <- log_add_exp(-2 * u - 2 * log(sigma_emb), log_cap)
  u - log1p_exp(2 * (u - log(axis_scale))) +
    m_gc / 2 * log_pg - exp(-2 * u + log(ss_gc)) / 2 +
    S / 2 * log_pe - exp(-2 * u - 2 * log(sigma_emb) + log(ss_r)) / 2
}

# The log of the conditional density of v = log(sigma_emb) under a cap
# `slab`, up to a constant: f_sigma of the design's "The draw of the leftover
# scale" (same file and section as above):
#
#   f_sigma(v) = v - log(1 + e^(2v) / b^2)
#                + sum over embedding columns k of
#                  [(S / 2) log P_ek - P_ek ss_rk / 2],
#   P_ek = e^(-2v) / lambda_k^2 + c^-2,
#
# with b = scale and c = slab. `ss_r` and `lambda` have one element per
# embedding column: the sum of squares of its S leftovers, and its uncapped
# axis scale. Computed as axis_scale_logdens() is, with the constant parts
# c^-2 ss_rk / 2 of the quadratic terms dropped, so that each quadratic term
# is e^(-2v) ss_rk / (2 lambda_k^2).
sigma_emb_logdens <- function(v, ss_r, S, lambda, scale, slab) {
  log_lam <- 2 * log(lambda)
  log_pe <- log_add_exp(-2 * v - log_lam, -2 * log(slab))
  v - log1p_exp(2 * (v - log(scale))) +
    sum(S / 2 * log_pe - exp(-2 * v - log_lam + log(ss_r)) / 2)
}

# Update the scales of the grouped prior. R is the species by covariates
# matrix of leftovers, B minus its trait regression; G and C have one column
# per covariate. Drawn in this order: sigma_b from the measured columns,
# sigma_emb from the embedding columns with each divided by its covariate
# scale, then each embedding column's covariate scale given the new sigma_emb.
# A covariate scale multiplies the prior standard deviation of its column of
# leftovers, of G and of C, so all three inform it. With no embedding columns
# sigma_emb and lambda are returned as given. With `hold_sigma_emb` the draw
# of sigma_emb is skipped (no random number is drawn for it) and the covariate
# scales are drawn given the value passed in: the staged burn-in of
# axis_hold_sweeps().
#
# Under a finite cap (coef_slab()) the embedding scales no longer have
# inverse-gamma conditionals. sigma_b is drawn as above, then sigma_emb and
# each axis scale, in the same order, by one slice draw each on the log scale
# (slice_sample_log()) from sigma_emb_logdens() and axis_scale_logdens(),
# started at the current value. The three sums of squares per column are the
# ones the uncapped draw forms. With an infinite cap, or a record without
# one, the uncapped draws run exactly as before and draw the same random
# numbers.
sample_coef_scales <- function(R, G, C, coef_prior, sigma_b, sigma_emb, lambda,
                               hold_sigma_emb = FALSE) {
  S <- nrow(R)
  is_emb <- coef_prior$is_emb
  sigma_b <- sample_scale_half_cauchy(sum(R[, !is_emb]^2), S * sum(!is_emb),
                                      sigma_b, coef_prior$sigma_b_scale)
  slab <- coef_slab(coef_prior)
  if (any(is_emb) && is.finite(slab)) {
    emb <- which(is_emb)
    ss_r <- colSums(R[, emb, drop = FALSE]^2)
    if (!hold_sigma_emb) {
      lambda_emb <- lambda[emb]
      f_sigma <- function(v) sigma_emb_logdens(v, ss_r, S, lambda_emb,
                                               coef_prior$sigma_emb_scale, slab)
      sigma_emb <- exp(slice_sample_log(f_sigma, log(sigma_emb)))
    }
    m_gc <- nrow(G) + nrow(C)
    for (j in seq_along(emb)) {
      k <- emb[j]
      ss_gc <- sum(G[, k]^2) + sum(C[, k]^2) / coef_prior$c_base[k]^2
      f_k <- function(u) axis_scale_logdens(u, ss_gc, m_gc, ss_r[[j]], S, sigma_emb,
                                            coef_prior$axis_scale, slab)
      lambda[k] <- exp(slice_sample_log(f_k, log(lambda[k])))
    }
  } else if (any(is_emb)) {
    emb <- which(is_emb)
    if (!hold_sigma_emb) {
      ss_emb <- sum(sweep(R[, emb, drop = FALSE], 2, lambda[emb], "/")^2)
      sigma_emb <- sample_scale_half_cauchy(ss_emb, S * length(emb), sigma_emb,
                                            coef_prior$sigma_emb_scale)
    }
    for (k in emb) {
      ss_k <- sum(R[, k]^2) / sigma_emb^2 + sum(G[, k]^2) +
        sum(C[, k]^2) / coef_prior$c_base[k]^2
      lambda[k] <- sample_scale_half_cauchy(ss_k, S + nrow(G) + nrow(C),
                                            lambda[k], coef_prior$axis_scale)
    }
  }
  list(sigma_b = sigma_b, sigma_emb = sigma_emb, lambda = lambda)
}

# A scale-expansion Metropolis move on each embedding column's axis scale.
# The axis scales mix slowly under the Gibbs draws alone, because each is
# drawn given the coefficients it scales (the centred parameterisation), and
# this move is what lets them travel. Experimental with the prior. The record
# of how it was found and tested is dev/simstudy/grouped-mixing-20261118 (fix
# attempt 1), kept under the tag archive/grouped-mixing.
#
# For embedding column k, propose log(lambda_k) + delta, delta ~ Normal(0,
# step^2), and multiply the column's leftovers, its column of G and its column
# of C by r = exp(delta). B[k, ] = Tr_b G[, k] + A C[, k] + leftovers is linear
# in those, so it is multiplied by exactly r. All three prior standard
# deviations are proportional to lambda_k (coef_prior_sds()), so the prior
# density of the scaled coefficients at the new scale is r^-m times that of
# the old ones at the old scale, m being the number of coefficients scaled,
# and the Jacobian of the scaling is r^m: they cancel. The acceptance ratio is
# what is left: the likelihood ratio, times the ratio of the half-Cauchy
# density of the axis scale expressed in log scale,
# p(lambda') lambda' / (p(lambda) lambda), p(l) = 1 / (1 + (l / axis_scale)^2).
#
# The likelihood is the one the coefficient draw uses, given the current
# Polya-Gamma weights Omega: sum(kappa * eta - Omega * eta^2 / 2). Omega is
# conditioned on, not integrated out, because the factor-score draw later in
# the same sweep still conditions on it. Scaling B[k, ] adds
# (r - 1) * X[, k] B[k, ] to eta, so the log likelihood ratio is
# (r - 1) * lin - (r - 1)^2 * quad / 2 with
# lin = sum((kappa - Omega * eta) * X[, k] B[k, ]) and
# quad = sum(Omega * (X[, k] B[k, ])^2).
#
# First version: one row per site and a binary response only (model type
# "binary"). axis_move_start() switches it on for those fits alone; fits of
# every other model type, and every fit under the fixed prior, skip it and
# draw exactly the random numbers they drew before. `step` is the proposal
# standard deviation, one value or one per embedding column. Columns are
# visited in order; all proposals are drawn first, then all uniforms. Returns
# the updated B, G, C and lambda, with `accepted` and `log_ratio` per design
# column (FALSE and NA on measured columns).
#
# Under a finite cap (coef_slab()) the prior standard deviations are no
# longer proportional to lambda_k: G and C of column k have the effective
# axis scale (lambda_k capped) and its leftovers lambda_k * sigma_emb capped
# (coef_effective_scales()). The move is then the design's "The axis scale
# move" (dev/superpowers/specs/2026-11-25-axis-slab-design.md, section
# "Mathematics"), in axis_scale_move_capped(): the same proposal on
# log(lambda_k), but each part of the column is multiplied by the ratio of
# its own prior standard deviation, new over old. Column k of G and of C is
# multiplied by rho_g, the ratio of the effective axis scales, and the
# column's leftovers by rho_e, the ratio of the leftover's capped standard
# deviations. Each part keeps its value in units of its own prior standard
# deviation, so the prior densities and the Jacobian cancel as above, and the
# acceptance ratio is again the likelihood ratio times the half-Cauchy term.
# Row k of B moves by d = (rho_g - 1) (B[k, ] - e) + (rho_e - 1) e, e being
# the column's leftovers, and the log likelihood ratio is lin - quad / 2 with
# lin = sum((kappa - Omega * eta) * X[, k] d) and
# quad = sum(Omega * (X[, k] d)^2). This needs the leftovers `R` of the `B`
# given (species by covariates) and the current `sigma_emb`, which are
# required under a finite cap and not used otherwise. With an infinite cap,
# or a record without one, rho_g = rho_e = exp(delta) and the code below runs
# unchanged, drawing the same random numbers, whether or not `R` and
# `sigma_emb` are given. The random numbers are drawn in the same order under
# a cap. Returns the same list.
AXIS_MOVE_STEP <- 0.5

# The acceptance rate a one-dimensional random walk is tuned to.
AXIS_MOVE_TARGET <- 0.44

axis_scale_move <- function(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step,
                            R = NULL, sigma_emb = NULL) {
  slab <- coef_slab(coef_prior)
  if (is.finite(slab))
    return(axis_scale_move_capped(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step,
                                  R, sigma_emb, slab))
  p <- nrow(B)
  accepted <- rep(FALSE, p)
  log_ratio <- rep(NA_real_, p)
  emb <- which(coef_prior$is_emb)
  if (length(emb) > 0L) {
    delta <- stats::rnorm(length(emb), 0, step)
    logu <- log(stats::runif(length(emb)))
    grad <- kappa - Omega * eta
    for (j in seq_along(emb)) {
      k <- emb[j]
      x <- X[, k]
      b <- B[k, ]
      r1 <- expm1(delta[j])
      lin <- sum(crossprod(x, grad) * b)
      quad <- sum(crossprod(x^2, Omega) * b^2)
      proposed <- lambda[k] * exp(delta[j])
      log_ratio[k] <- r1 * lin - r1^2 * quad / 2 + delta[j] -
        log1p((proposed / coef_prior$axis_scale)^2) +
        log1p((lambda[k] / coef_prior$axis_scale)^2)
      if (logu[j] < log_ratio[k]) {
        r <- exp(delta[j])
        grad <- grad - Omega * tcrossprod(x, r1 * b)
        B[k, ] <- b * r
        G[, k] <- G[, k] * r
        C[, k] <- C[, k] * r
        lambda[k] <- proposed
        accepted[k] <- TRUE
      }
    }
  }
  list(B = B, G = G, C = C, lambda = lambda, accepted = accepted, log_ratio = log_ratio)
}

# The log of capped(s * exp(delta), slab) / capped(s, slab), elementwise,
# given log(s): how a proposal of delta on the log scale changes a capped
# standard deviation (coef_capped()). Both capped values are written as in
# coef_capped(), through the smaller and the larger of s and the cap, so
# that the log of the cap cancels exactly when s is far above it and the
# ratio stays finite and accurate for s from far below the cap to far above
# it. Taking log(s) rather than s keeps a product of two small scales from
# underflowing.
coef_capped_log_ratio <- function(log_s, delta, slab) {
  log_c <- log(slab)
  lo0 <- pmin(log_s, log_c); hi0 <- pmax(log_s, log_c)
  lo1 <- pmin(log_s + delta, log_c); hi1 <- pmax(log_s + delta, log_c)
  (lo1 - lo0) - (log1p(exp(2 * (lo1 - hi1))) - log1p(exp(2 * (lo0 - hi0)))) / 2
}

# axis_scale_move() under a finite cap `slab`: see the comment above it. The
# factors rho_g and rho_e are formed on the log scale (coef_capped_log_ratio())
# and rho - 1 by expm1(), so that a column far above the cap, whose factors
# are within rounding of 1, moves B by a correspondingly tiny amount rather
# than by rounding error. `R` is read for column k only when column k is
# visited, and the move changes no other column's leftovers, so it need not
# be updated after an acceptance.
axis_scale_move_capped <- function(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step,
                                   R, sigma_emb, slab) {
  if (is.null(R) || is.null(sigma_emb))
    stop("under a cap the axis scale move needs the leftovers R and sigma_emb")
  p <- nrow(B)
  accepted <- rep(FALSE, p)
  log_ratio <- rep(NA_real_, p)
  emb <- which(coef_prior$is_emb)
  if (length(emb) > 0L) {
    delta <- stats::rnorm(length(emb), 0, step)
    logu <- log(stats::runif(length(emb)))
    grad <- kappa - Omega * eta
    log_sigma <- log(sigma_emb)
    for (j in seq_along(emb)) {
      k <- emb[j]
      x <- X[, k]
      b <- B[k, ]
      e <- R[, k]
      log_rho_g <- coef_capped_log_ratio(log(lambda[k]), delta[j], slab)
      log_rho_e <- coef_capped_log_ratio(log_sigma + log(lambda[k]), delta[j], slab)
      d <- expm1(log_rho_g) * (b - e) + expm1(log_rho_e) * e
      lin <- sum(crossprod(x, grad) * d)
      quad <- sum(crossprod(x^2, Omega) * d^2)
      proposed <- lambda[k] * exp(delta[j])
      log_ratio[k] <- lin - quad / 2 + delta[j] -
        log1p((proposed / coef_prior$axis_scale)^2) +
        log1p((lambda[k] / coef_prior$axis_scale)^2)
      if (logu[j] < log_ratio[k]) {
        grad <- grad - Omega * tcrossprod(x, d)
        B[k, ] <- b + d
        G[, k] <- G[, k] * exp(log_rho_g)
        C[, k] <- C[, k] * exp(log_rho_g)
        lambda[k] <- proposed
        accepted[k] <- TRUE
      }
    }
  }
  list(B = B, G = G, C = C, lambda = lambda, accepted = accepted, log_ratio = log_ratio)
}

# The same scale-expansion move for sigma_emb, the embedding leftover scale,
# which mixes slowly for the same reason: it is drawn given every leftover of
# the embedding block and they are then drawn given it. Record:
# dev/simstudy/grouped-mixing-20261118 (fix attempt 2), tag
# archive/grouped-mixing.
#
# Propose log(sigma_emb) + delta, delta ~ Normal(0, step^2), and multiply
# every leftover of the embedding block, R[, embedding], by r = exp(delta).
# G, C, A and the axis scales stay where they are, so each embedding row of B
# moves by (r - 1) times its leftovers. A leftover's prior standard deviation
# is sigma_emb times its column's axis scale, so the prior density of the
# scaled leftovers at the new sigma_emb is r^-m times that of the old ones at
# the old, m being the number of leftovers scaled, and the Jacobian of the
# scaling is r^m: they cancel. The acceptance ratio is the likelihood ratio
# times p(sigma') sigma' / (p(sigma) sigma), p the half-Cauchy density with
# scale sigma_emb_scale.
#
# The likelihood is conditioned on the Polya-Gamma weights, as in
# axis_scale_move(). The move adds (r - 1) * D to eta, with
# D = X[, embedding] %*% t(R[, embedding]), so the log likelihood ratio is
# (r - 1) * lin - (r - 1)^2 * quad / 2, lin = sum((kappa - Omega * eta) * D)
# and quad = sum(Omega * D^2).
#
# `R` is the species by covariates matrix of leftovers of the `B` given. One
# proposal: delta is drawn, then the uniform. Returns the updated B,
# sigma_emb and eta, with `accepted` and `log_ratio` (FALSE and NA, with no
# random number drawn, when there are no embedding columns). It runs in the
# same fits as axis_scale_move(), just before it.
#
# Under a finite cap (coef_slab()) a leftover's prior standard deviation is
# sigma_emb * lambda_k capped, which is not proportional to sigma_emb. The
# move is then the design's "The leftover scale move" (same file and section
# as for axis_scale_move()), in sigma_emb_move_capped(): the same proposal,
# but the leftovers of embedding column k are multiplied by their own factor
# rho_k, the ratio of the capped standard deviations at the proposed and
# current sigma_emb, so that each leftover keeps its value in units of its
# prior standard deviation and the prior and the Jacobian cancel as above.
# Each embedding row of B moves by (rho_k - 1) times its leftovers, and with
# D = X[, embedding] %*% diag(rho - 1) %*% t(R[, embedding]) the log
# likelihood ratio is lin - quad / 2, lin = sum((kappa - Omega * eta) * D)
# and quad = sum(Omega * D^2). This needs the current axis scales `lambda`
# (uncapped, one per design column), which are required under a finite cap
# and not used otherwise. With an infinite cap, or a record without one,
# every rho_k is exp(delta) and the code below runs unchanged, drawing the
# same random numbers, whether or not `lambda` is given. Same random numbers
# in the same order under a cap, and the same return.
sigma_emb_move <- function(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step,
                           lambda = NULL) {
  slab <- coef_slab(coef_prior)
  if (is.finite(slab))
    return(sigma_emb_move_capped(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step,
                                 lambda, slab))
  emb <- which(coef_prior$is_emb)
  if (length(emb) == 0L)
    return(list(B = B, sigma_emb = sigma_emb, eta = eta, accepted = FALSE,
                log_ratio = NA_real_))
  delta <- stats::rnorm(1, 0, step)
  logu <- log(stats::runif(1))
  r1 <- expm1(delta)
  D <- X[, emb, drop = FALSE] %*% t(R[, emb, drop = FALSE])
  lin <- sum((kappa - Omega * eta) * D)
  quad <- sum(Omega * D^2)
  proposed <- sigma_emb * exp(delta)
  log_ratio <- r1 * lin - r1^2 * quad / 2 + delta -
    log1p((proposed / coef_prior$sigma_emb_scale)^2) +
    log1p((sigma_emb / coef_prior$sigma_emb_scale)^2)
  accepted <- logu < log_ratio
  if (accepted) {
    B[emb, ] <- B[emb, , drop = FALSE] + r1 * t(R[, emb, drop = FALSE])
    eta <- eta + r1 * D
    sigma_emb <- proposed
  }
  list(B = B, sigma_emb = sigma_emb, eta = eta, accepted = accepted, log_ratio = log_ratio)
}

# sigma_emb_move() under a finite cap `slab`: see the comment above it. The
# factors are formed on the log scale (coef_capped_log_ratio()) and rho - 1
# by expm1(), as in axis_scale_move_capped().
sigma_emb_move_capped <- function(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step,
                                  lambda, slab) {
  emb <- which(coef_prior$is_emb)
  if (length(emb) == 0L)
    return(list(B = B, sigma_emb = sigma_emb, eta = eta, accepted = FALSE,
                log_ratio = NA_real_))
  if (is.null(lambda))
    stop("under a cap the sigma_emb move needs the axis scales lambda")
  delta <- stats::rnorm(1, 0, step)
  logu <- log(stats::runif(1))
  rho1 <- expm1(coef_capped_log_ratio(log(sigma_emb) + log(lambda[emb]), delta, slab))
  # Row j of dB is rho_j - 1 times the leftovers of embedding column j: the
  # change of that row of B.
  dB <- t(R[, emb, drop = FALSE]) * rho1
  D <- X[, emb, drop = FALSE] %*% dB
  lin <- sum((kappa - Omega * eta) * D)
  quad <- sum(Omega * D^2)
  proposed <- sigma_emb * exp(delta)
  log_ratio <- lin - quad / 2 + delta -
    log1p((proposed / coef_prior$sigma_emb_scale)^2) +
    log1p((sigma_emb / coef_prior$sigma_emb_scale)^2)
  accepted <- logu < log_ratio
  if (accepted) {
    B[emb, ] <- B[emb, , drop = FALSE] + dB
    eta <- eta + D
    sigma_emb <- proposed
  }
  list(B = B, sigma_emb = sigma_emb, eta = eta, accepted = accepted, log_ratio = log_ratio)
}

# The same scale-expansion move for sigma_b, the measured group's scale. Under
# the grouped prior sigma_b has a half-Cauchy prior with no floor, and chains
# were seen to drift to sigma_b near zero and stay there for hundreds to
# thousands of sweeps: it is drawn given the measured leftovers and they are
# then drawn given it, so near zero each draw moves it by only a few percent.
# Design: dev/superpowers/specs/2026-11-22-sigma-b-move-design.md.
#
# Propose log(sigma_b) + delta, delta ~ Normal(0, step^2), and multiply every
# leftover of the measured columns, R[, measured], by r = exp(delta). G, C, A
# and the embedding block stay where they are, so each measured row of B
# moves by (r - 1) times its leftovers. The axis scale of a measured column is
# 1, so the prior of the leftovers and the Jacobian cancel as in
# sigma_emb_move(), and the acceptance ratio is
#   (r - 1) * lin - (r - 1)^2 * quad / 2 + delta
#     - log1p((sigma_b' / sigma_b_scale)^2) + log1p((sigma_b / sigma_b_scale)^2)
# with D = X[, measured] %*% t(R[, measured]), lin = sum((kappa - Omega * eta) * D)
# and quad = sum(Omega * D^2), the likelihood conditioned on the Polya-Gamma
# weights. Far below the ordinary range D is proportional to sigma_b, so the
# ratio is close to delta and every upward proposal is accepted.
#
# `step` has one element per proposal, made in order, each from the state the
# one before it left. update_jSDMcoef() gives two: the step adapted during
# burn-in and SIGMA_B_MOVE_LARGE_STEP. All the deltas are drawn first, then
# all the uniforms. Returns B, R (the leftovers of the returned B), sigma_b
# and eta, with `accepted` and `log_ratio` of one element per proposal (FALSE
# and NA, with no random number drawn, when there are no measured columns).
sigma_b_move <- function(B, R, sigma_b, X, kappa, Omega, eta, coef_prior, step) {
  meas <- which(!coef_prior$is_emb)
  n_prop <- length(step)
  if (length(meas) == 0L)
    return(list(B = B, R = R, sigma_b = sigma_b, eta = eta,
                accepted = rep(FALSE, n_prop), log_ratio = rep(NA_real_, n_prop)))
  delta <- stats::rnorm(n_prop, 0, step)
  logu <- log(stats::runif(n_prop))
  accepted <- logical(n_prop); log_ratio <- numeric(n_prop)
  D <- X[, meas, drop = FALSE] %*% t(R[, meas, drop = FALSE])
  for (j in seq_len(n_prop)) {
    r1 <- expm1(delta[j])
    lin <- sum((kappa - Omega * eta) * D)
    quad <- sum(Omega * D^2)
    proposed <- sigma_b * exp(delta[j])
    log_ratio[j] <- r1 * lin - r1^2 * quad / 2 + delta[j] -
      log1p((proposed / coef_prior$sigma_b_scale)^2) +
      log1p((sigma_b / coef_prior$sigma_b_scale)^2)
    accepted[j] <- logu[j] < log_ratio[j]
    if (accepted[j]) {
      B[meas, ] <- B[meas, , drop = FALSE] + r1 * t(R[, meas, drop = FALSE])
      eta <- eta + r1 * D
      R[, meas] <- R[, meas, drop = FALSE] * exp(delta[j])
      D <- D * exp(delta[j])
      sigma_b <- proposed
    }
  }
  list(B = B, R = R, sigma_b = sigma_b, eta = eta, accepted = accepted, log_ratio = log_ratio)
}

# The fixed step of the second sigma_b proposal in each sweep, on the log
# scale. The adapted step settles where the ordinary state accepts about half
# of the proposals, which is small for a chain far below that state; this one
# is large enough to climb back in a few tens of sweeps and is mostly
# rejected elsewhere.
SIGMA_B_MOVE_LARGE_STEP <- 1

# Burn-in adaptation of the proposal steps, one per embedding column, from
# that sweep's log acceptance ratios: the log of each step moves by
# (acceptance probability - target) over the square root of the number of
# adapting sweeps so far, `t`, the Robbins-Monro recursion of b0_shift_step().
# A ratio that is not finite counts as a rejection. Steps stay within loose
# bounds. The chain adapts only during burn-in, so the kept sweeps use a
# fixed kernel.
axis_step_adapt <- function(step, log_ratio, t) {
  alpha <- ifelse(is.finite(log_ratio), pmin(1, exp(log_ratio)), 0)
  pmin(pmax(step * exp((alpha - AXIS_MOVE_TARGET) / sqrt(t)), 0.01), 5)
}

# A staged burn-in (record: dev/simstudy/grouped-mixing-20261118, fix attempt
# 3, tag archive/grouped-mixing). In a low-rank truth chains settle in one of
# two configurations and do not cross: the community response and blends carry
# the signal, or the leftovers do and nothing is pooled across species. Which
# one is chosen early in burn-in. The hypothesis is that the leftovers win because they can take
# up the signal in a few sweeps while the blends need many to organise. So for
# the first half of burn-in sigma_emb is held at a tenth of its prior scale:
# the chain starts there, its draw is skipped (sample_coef_scales()) and so is
# its scale move. From then on the sampler is unchanged. Burn-in is discarded,
# so the kept draws come from the same kernel as before. Same fits as the two
# scale moves. The length and the value are fixed, not tuned.
#
# The hold is applied only from AXIS_HOLD_MIN_BURN burn-in sweeps up; below
# that there is none and the sampler is that of the scale moves alone
# (Doug's decision of 19 November 2026). After the release the blends need
# time to settle. With 5,000 burn-in sweeps (2,500 after the release) every
# low-rank chain reached the blend configuration and the active coefficients
# were recovered much better in two data sets of three, at no cost in the
# third; with 1,500 (750 after the release) recovery was worse than without
# the hold in two of three. Those are the only two lengths that were tried,
# so the minimum is the shorter of the lengths known to work, not a measured
# threshold.
AXIS_HOLD_MIN_BURN <- 5000L
axis_hold_sweeps <- function(nburn)
  if (nburn >= AXIS_HOLD_MIN_BURN) as.integer(nburn %/% 2) else 0L

# The note runOccJSDM() prints for a fit that would have the staged burn-in
# but has too little burn-in for it, or NULL: a grouped fit of model type
# "binary" with embedding columns, the fits axis_move_start() gives the scale
# moves to.
axis_hold_note <- function(coef_prior, model, nburn) {
  if (is.null(axis_move_start(coef_prior, model)) || axis_hold_sweeps(nburn) > 0L)
    return(NULL)
  paste0("Grouped prior: with fewer than ",
         format(AXIS_HOLD_MIN_BURN, big.mark = ","),
         " burn-in iterations the embedding leftover scale is not held small at ",
         "the start of each chain, so chains may settle where the leftovers and ",
         "not the blends carry the species' responses. See ?\"grouped-prior\".")
}
axis_hold_value <- function(coef_prior) 0.1 * coef_prior$sigma_emb_scale

# The state of the three scale moves that a chain carries, or NULL when they
# are off: they run only in a grouped fit with embedding columns of model type
# "binary". Per embedding column: `step`, and `accepted`, the last sweep's
# outcome, of the axis scale move. `sigma_step` and `sigma_accepted` are those
# of the sigma_emb move. `sigma_b_step` is the adapted step of the sigma_b
# move and `sigma_b_accepted` the last sweep's outcomes of its two proposals
# (adapted, then large). `adapt` is set by the chain (TRUE during burn-in) and
# `t` counts the adapting sweeps; `sigma_t` counts those in which the
# sigma_emb move ran, and `sigma_b_t` those in which the sigma_b move adapted
# (all adapting sweeps, when there is a measured column). `hold` is set by
# the chain (TRUE while sigma_emb is held).
axis_move_start <- function(coef_prior, model) {
  if (!identical(model, "binary") || !identical(coef_prior$type, "grouped") ||
      !any(coef_prior$is_emb)) return(NULL)
  n_emb <- sum(coef_prior$is_emb)
  list(step = rep(AXIS_MOVE_STEP, n_emb), accepted = rep(FALSE, n_emb),
       sigma_step = AXIS_MOVE_STEP, sigma_accepted = FALSE,
       sigma_b_step = AXIS_MOVE_STEP, sigma_b_accepted = c(FALSE, FALSE),
       adapt = FALSE, t = 0L, sigma_t = 0L, sigma_b_t = 0L, hold = FALSE)
}

# infos$axis_move from the chains' records, each a matrix with one row per
# embedding column and one column per chain: the proposal step the kept sweeps
# used, and the share of sweeps in which the move was accepted, after and
# during burn-in. For the sigma_emb move, one value per chain: its step and
# its acceptance rates (the burn-in rate is over all burn-in sweeps, those in
# which sigma_emb was held included). For the sigma_b move: its adapted step,
# one per chain, its large step, and its acceptance rates, each a matrix with
# rows "adapted" and "large" and one column per chain. `hold_sweeps` and
# `hold_value` describe the staged burn-in. NULL when the moves are off.
axis_move_record <- function(parts, columns) {
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) return(NULL)
  rate <- function(phase) {
    out <- vapply(parts, function(p) p[[paste0("accepted_", phase)]] /
                    p[[paste0("sweeps_", phase)]], numeric(length(columns)))
    matrix(out, length(columns), length(parts), dimnames = list(columns, NULL))
  }
  step <- vapply(parts, function(p) p$step, numeric(length(columns)))
  sigma_rate <- function(phase) vapply(parts, function(p)
    p[[paste0("sigma_accepted_", phase)]] / p[[paste0("sweeps_", phase)]], numeric(1))
  sigma_b_rate <- function(phase) {
    out <- vapply(parts, function(p) p[[paste0("sigma_b_accepted_", phase)]] /
                    p[[paste0("sweeps_", phase)]], numeric(2))
    matrix(out, 2L, length(parts), dimnames = list(c("adapted", "large"), NULL))
  }
  list(step = matrix(step, length(columns), length(parts), dimnames = list(columns, NULL)),
       acceptance = rate("kept"), acceptance_burnin = rate("burnin"),
       sigma_emb_step = vapply(parts, function(p) p$sigma_step, numeric(1)),
       sigma_emb_acceptance = sigma_rate("kept"),
       sigma_emb_acceptance_burnin = sigma_rate("burnin"),
       sigma_b_step = vapply(parts, function(p) p$sigma_b_step, numeric(1)),
       sigma_b_large_step = SIGMA_B_MOVE_LARGE_STEP,
       sigma_b_acceptance = sigma_b_rate("kept"),
       sigma_b_acceptance_burnin = sigma_b_rate("burnin"),
       hold_sweeps = parts[[1]]$hold_sweeps, hold_value = parts[[1]]$hold_value)
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
#' nothing for that scale to learn from, so `sigmab_output` then holds draws
#' from its prior. They stay in the fit, but [returnCovariateScales()] leaves
#' the measured-group scale out and [computeDiagnostics()] does not check it.
#' `G_output` gains a leading row, the community-average response.
#' `infos$coef_prior` records the prior and its settings, and
#' `infos$covariate_group` labels each occupancy design column `"measured"` or
#' `"embedding"`. In a `binary` fit with embedding columns, `infos$axis_move`
#' records the sampler's extra steps for the scales (see Sampling and chain
#' length): their step sizes and acceptance rates, and how many burn-in
#' iterations the embedding leftover scale was held for.
#'
#' Seven functions summarise these:
#' * [returnGroupedConvergence()]: the convergence checks for a fit with
#'   embedding columns (see Checking convergence).
#' * [returnCovariateScales()]: posterior summaries of the measured-group,
#'   embedding leftover and axis scales.
#' * [returnCommunityResponse()]: the community-average response of every
#'   occupancy design column.
#' * [returnEmbeddingBlends()]: the blends, as weights on the embedding
#'   columns and as scores for each site. Blends are defined only up to
#'   rotation, so do not interpret them one by one.
#' * [returnEmbeddingGradients()]: the gradients of the species' whole
#'   response to the embedding, as scores for each site and loadings for each
#'   species (see Interpreting the gradients).
#' * [compareEmbeddingGradients()]: those gradients compared with covariates
#'   known for the same sites, and what the covariates leave unexplained.
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
#' @section Interpreting the gradients:
#' The embedding columns have no meaning one by one. What can be read is the
#' small set of gradients that the species share in their responses to the
#' embedding. Two gradients mark out a plane and three a block. Any rotation of
#' them marks out the same plane or block and fits the data equally well, so
#' the data determine the plane, not the particular directions a function
#' returns. Read the gradients as a set, never one by one.
#'
#' [returnEmbeddingGradients()] takes the gradients from each species' whole
#' response to the embedding, the shared part plus the leftover. How a response
#' is divided between those two is poorly determined, and at weak signal a fit
#' can put most of it in the leftovers although the species share gradients.
#' The leading directions of the whole response are far more stable. In
#' simulations they recovered the true gradients across sites, and chains
#' agreed on them, even in fits where the division was wrong and the scales
#' failed their convergence checks. [returnEmbeddingBlends()] works from the
#' shared part only, so rely on it only when the `state` row of
#' [returnGroupedConvergence()] is 0.
#'
#' These gradients are hidden environmental covariates of a second kind. Those
#' set by `n_factors` are free scores for each site, learned only from the
#' species found there. The embedding gradients are weighted sums of a site's
#' embedding columns, so they can be computed and mapped wherever the embedding
#' exists.
#'
#' [compareEmbeddingGradients()] relates the gradients to covariates known for
#' the same sites. It makes three comparisons.
#' * Covariate from gradients: how well do the site scores of the first few
#'   gradients, taken together, predict each covariate?
#' * Gradients from covariates: how much of the gradients do all the covariates
#'   together account for, and how much does each covariate account for alone
#'   and uniquely?
#' * Unmeasured gradients: the leading directions of what the covariates leave
#'   unexplained. This is habitat variation that the embedding records and the
#'   species respond to, but that the covariates do not describe, which
#'   includes a response to a measured covariate that is not a straight line.
#'   Each has a score at every site used, which can be mapped for those sites,
#'   and a loading per species. Read their shares, not how many there are:
#'   something is always left over, and a small total share means the
#'   covariates leave little.
#'
#' The covariates can also be the site scores of the hidden environmental
#' covariates, from this fit or from a fit without the embedding. The
#' comparison then shows whether a factor duplicates an embedding gradient, or
#' which hidden environmental covariates the satellite sees.
#' [returnOrdinationScores()] returns quantiles by sites by factors, so take
#' its second slice, the median, and add the site names:
#'
#' ```
#' scores <- returnOrdinationScores(fit)
#' factors <- data.frame(Site = dimnames(scores)[[2]],
#'                       as.data.frame(matrix(scores[2, , ], nrow = dim(scores)[2])))
#' compareEmbeddingGradients(returnEmbeddingGradients(fit), factors)
#' ```
#'
#' Factors are also defined only up to rotation, so read the overall share,
#' not single factors.
#'
#' Cautions.
#' * The number of gradients is uncertain. Check that a conclusion holds with
#'   one more and one fewer.
#' * A match is a hypothesis about a driver, not evidence of one. Anything
#'   correlated with the covariate across your sites would match as well. The
#'   function gives no test and no interval, because sites near one another
#'   are not independent.
#' * The number of gradients is not the number of drivers. One covariate with
#'   a humped response can take up two gradients, and one gradient can be a
#'   mixture of several covariates. Set `squared = TRUE` when responses may be
#'   humped.
#' * Species' positions on the gradients are less well determined than the
#'   sites' scores, most of all at weak signal.
#' * A covariate that was in the model can still match the gradients, because
#'   the embedding columns are correlated with it. The model holds it fixed
#'   when estimating the embedding coefficients, but that does not make the
#'   gradients' site scores unrelated to it. Such a match shows overlap between
#'   the covariate and the embedding, not a new driver. To ask what a new
#'   covariate adds beyond the model's covariates, give both and read `unique`.
#'
#' @section Sampling and chain length:
#' The scales of this prior are hard for the sampler. In a `binary` fit with
#' embedding columns it takes four extra steps for them. They are
#' experimental, and they do not change the model: the draws that are kept
#' come from the same posterior as without them.
#' * Every iteration, each embedding column's axis scale is proposed a new
#'   value together with that column's coefficients, rescaled by the same
#'   factor. Without this step the axis scale of a column that no species
#'   responds to moved very slowly.
#' * The embedding leftover scale gets the same kind of step, together with
#'   all the leftovers.
#' * The scale of the measured covariates' species leftovers (`sigma_b`) gets
#'   the same kind of step, together with those leftovers, twice in each
#'   iteration: once with a step tuned during burn-in and once with a large
#'   fixed step. Without it a chain could stay for thousands of iterations
#'   with this scale near zero, and disagree with the other chains on the
#'   measured coefficients meanwhile.
#' * With 5,000 or more burn-in iterations, the embedding leftover scale is
#'   held at a tenth of its prior scale for the first half of burn-in and then
#'   released.
#'
#' The reason for the last step is that the species' responses to the
#' embedding can be explained in two ways: by the shared part (the community
#' average, the traits and the blends), or by the leftovers alone, species by
#' species. Early in a chain the leftovers can take up the responses within a
#' few iterations, before the blends have formed, and the chain then stays
#' there. Such a fit reports almost no shared structure and shrinks the
#' embedding coefficients too far, and because its chains agree with each
#' other the convergence diagnostics do not show it. Holding the leftovers
#' small at first lets the blends form. In simulations in which species shared
#' three blends, every chain then found blends in all but one of 21 fits at
#' 5,000 burn-in and 5,000 kept iterations, and the embedding coefficients
#' were recovered much better in two data sets of three and equally well in
#' the third. With 1,500 burn-in iterations the blends had too little time to
#' settle after the release and the fits were worse, so below 5,000 there is
#' no hold and [runOccJSDM()] prints a note. A fit that short can still end
#' with the leftovers carrying the responses.
#'
#' So use at least 5,000 burn-in and 5,000 kept iterations, which are the
#' defaults, and at least four chains. Even then the chains of some fits
#' settle on different blends.
#'
#' The `occupancy`, `two_stage` and `continuous` model types do not have these
#' steps yet. There the scales still mix slowly, and nothing guards against
#' the leftovers taking the responses.
#'
#' @section Checking convergence:
#' The usual test, R-hat on every coefficient, misleads for these fits. Single
#' blends and single embedding coefficients are not pinned down, and their
#' draws are heavy-tailed, so it fails fits that are fine. And chains that all
#' ended with the leftovers carrying the responses agree with each other, so
#' it passes a fit that found no gradients. [returnGroupedConvergence()] asks
#' instead which state each chain is in, whether the measured coefficients and
#' the two scales shared by all columns have converged, whether the chains
#' found the same set of gradients, and how far apart their fitted effects of
#' the embedding are.
#' [runOccJSDM()] prints its result at the end of a fit made with
#' `embeddingVars`, and warns when a check is flagged. In the same report the
#' embedding coefficients are listed with their R-hat but without a warning.
#'
#' A flagged fit is marked "review", not failed, because one check can flag a
#' fit that is right: where species truly respond independently to the
#' embedding, the leftovers should carry the responses.
#'
#' Two limits of these checks. The comparisons of gradients and of fitted
#' effects use each chain's posterior mean, so they say whether the chains
#' reached the same answer, not whether each explored its own uncertainty.
#' And they are made at the fitted sites. Chains that agree there can differ
#' at sites unlike them, such as parts of a map whose habitat was not
#' sampled, because the data say little about directions of the embedding in
#' which the fitted sites barely vary. The function reports one number for
#' that, without a limit. A fit that passes is not thereby safe to
#' extrapolate.
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
#' The sampler's extra steps for `binary` fits were added in November 2026
#' and have been tried only on simulated data of one size (256 sites, 48
#' species, 64 embedding columns). In those simulations, when no species
#' responded to the embedding, every scale and every measured coefficient
#' converged at 5,000 burn-in and 5,000 kept iterations. Before the steps
#' were added the axis scales did not converge at either length that was
#' tried there (1,500 and 5,000 of each).
#'
#' Each embedding column has a scale of its own, so the prior does best when
#' the species respond to a few of the columns as they are delivered. In the
#' same simulations, turning the columns of an embedding (which changes
#' nothing in what it encodes) raised the error of the embedding's
#' contribution by about a third when the responses had been placed in 10 of
#' 64 columns, and by little or nothing when they were spread over all of
#' them. A real embedding has no such favoured columns, so simulations that
#' place the responses in a few columns overstate how well the prior will do.
#'
#' The prior was built and tested for the `binary`, `occupancy` and
#' `two_stage` model types. A `continuous` fit runs with it, but that
#' combination is outside the design and no study has tested it.
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
