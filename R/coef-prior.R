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
    axis_scale <- scale("axis_scale", 1 / sqrt(p_emb))
  }
  c_base <- rep(1, length(is_emb))
  if (gt > 0) c_base[is_emb] <- 1 / sqrt(gt)
  list(type = "grouped", is_emb = is_emb, sigma_b_scale = sigma_b_scale,
       sigma_emb_scale = sigma_emb_scale, axis_scale = axis_scale,
       c_base = c_base)
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
