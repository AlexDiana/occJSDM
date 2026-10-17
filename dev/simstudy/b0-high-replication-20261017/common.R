# Shared by choose.R, mixing_check.R and analyse.R: the per-fit files of refit_spread.R as one row per species-fit,
# the bands of true B0, and replicate-clustered means and slopes. Each script sets `here` to this directory first.
`%||%` <- function(a, b) if (is.null(a)) b else a

read_fits <- function(dirs) {
  files <- unlist(lapply(dirs, list.files, pattern = "[.]rds$", full.names = TRUE))
  lapply(files, readRDS)
}

# Bands of true B0, as in the ladder report. Extreme bands: the two outer ones in the ordinary community; the two
# lowest in the rare-dominated one, where its error is largest.
is_rare <- function(arm) grepl("rare", arm)
band_breaks <- function(arm) if (is_rare(arm)) c(-Inf, -4.25, -3.25, -2.25, -1.25, Inf) else
  c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf)
extreme_bands <- function(arm) if (is_rare(arm)) c(1L, 2L) else c(1L, 5L)

# One row per species-fit. e = J + V - D exactly (README.md); P is the spread-only term, Q the approximate pull of
# B0's Normal prior towards its centre (0 under the default prior, the posterior mean of mu_b0 under the community
# mean). Coverage is the study's own, from the stored rows.
species_table <- function(fits) {
  do.call(rbind, lapply(fits, function(x) {
    sp <- x$species
    ip <- x$intercept_prior
    community <- identical(ip$type, "community")
    centre <- if (community) x$hyper$mu_b0[["mean"]] else ip$mean %||% 0
    w <- sp$B0_sd^2 / ip$sd^2
    b0_rows <- x$rows[x$rows$block == "B0", ]
    b0_rows <- b0_rows[order(b0_rows$element), ]
    occ <- x$occupancy[order(x$occupancy$species), ]
    out <- data.frame(
      job = x$job, arm = x$arm, variant = if (nzchar(x$variant)) x$variant else "base",
      replicate = x$replicate, species = sp$species, B0_true = sp$B0_true, B0_est = sp$B0_mean, B0_sd = sp$B0_sd,
      covered = b0_rows$covered, width = b0_rows$upper - b0_rows$lower,
      e = sp$B0_mean - sp$B0_true, J = sp$aJ_mean - sp$a_true, V = sp$a_fit_mean - sp$aJ_mean,
      D = sp$m_r_mean - sp$m_true, P = sp$aS_mean - sp$a_true, Q = w / (1 - w) * (centre - sp$B0_mean),
      a_err = sp$a_fit_mean - sp$a_true,
      sd_ratio = sp$sd_r_mean / sp$sd_r_true, sd_XB_ratio = sp$sd_XB_mean / sp$sd_XB_true,
      sd_UL_ratio = sp$sd_UL_mean / sp$sd_UL_true, sd_SP_ratio = sp$sd_SP_mean / sp$sd_SP_true,
      sd_site_ratio = sp$sd_site_mean / sp$sd_site_true,
      sd_r_true = sp$sd_r_true, var_r_true = sp$sd_r_true^2,
      var_r_fit = sp$var_XB_mean + sp$var_site_mean + 2 * sp$cov_XB_site_mean,
      dvar_XB = sp$var_XB_true - sp$var_XB_mean, dvar_site = sp$var_site_true - sp$var_site_mean,
      dcov = sp$cov_XB_site_true - sp$cov_XB_site_mean,
      m_true = sp$m_true, m_UL_err = sp$m_UL_mean - sp$m_UL_true, m_SP_err = sp$m_SP_mean - sp$m_SP_true,
      cor_r = sp$cor_r_mean, psibar_true = sp$psibar_true, psibar_est = sp$psibar_mean,
      zbar_true = sp$zbar_true, occ_true = occ$occ_true, occ_est = occ$occ_est,
      B0_rhat = sp$B0_rhat, B0_ess = sp$B0_ess, sd_r_rhat = sp$sd_r_rhat, sd_r_ess = sp$sd_r_ess,
      aJ_rhat = sp$aJ_rhat, aJ_ess = sp$aJ_ess, m_r_rhat = sp$m_r_rhat, m_r_ess = sp$m_r_ess,
      mu_b0 = if (community) centre else NA_real_, sigma_b0 = ip$sd,
      stringsAsFactors = FALSE)
    out$band_index <- cut(out$B0_true, band_breaks(x$arm), labels = FALSE)
    out$band <- as.character(cut(out$B0_true, band_breaks(x$arm)))
    out$extreme <- out$band_index %in% extreme_bands(x$arm)
    out
  }))
}

# Mean over species-fits with a standard error treating each replicate as one cluster.
cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  c(mean = m, se = if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA)
}

# Least-squares slope of y on x with a standard error clustered by g, and the R-squared.
cluster_slope <- function(y, x, g) {
  xc <- x - mean(x); yc <- y - mean(y)
  b <- sum(xc * yc) / sum(xc^2)
  u <- yc - b * xc
  s <- tapply(xc * u, g, sum); G <- length(s)
  c(slope = b, se = sqrt(G / (G - 1) * sum(s^2)) / sum(xc^2), r2 = 1 - sum(u^2) / sum(yc^2))
}
