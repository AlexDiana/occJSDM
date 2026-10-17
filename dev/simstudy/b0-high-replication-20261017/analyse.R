# Scores the B0 high-replication run as README.md sets out. Run from the repository root:
#   Rscript dev/simstudy/b0-high-replication-20261017/analyse.R
# Reads the per-fit files in draws/ (diag, interventions and, if present, long) and the stored rows of the arms'
# original runs. Writes reproduction.csv, species-fits.csv, decomposition.csv, variance-shares.csv,
# regression.csv, spread-criterion.csv, interventions.csv, hyper.csv and, if the long step ran, long-chains.csv.
here <- if (file.exists("refit_spread.R")) "." else "dev/simstudy/b0-high-replication-20261017"
source(file.path(here, "common.R"))
root <- normalizePath(file.path(here, "..", "..", ".."))
arms <- c("M20", "M20_cm", "base_rare_cm_m20")
# DRAWS_DIR and OUT_DIR, if set, move where the per-fit files are read from and the tables are written to.
draws_dir <- Sys.getenv("DRAWS_DIR", file.path(here, "draws"))
out <- Sys.getenv("OUT_DIR", here)
diag_fits <- read_fits(file.path(draws_dir, "diag"))
int_fits <- read_fits(file.path(draws_dir, "interventions"))
long_fits <- read_fits(file.path(draws_dir, "long"))
fits <- c(diag_fits, int_fits)
wcsv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)

# 1. Reproduction: every diagnosis refit against its arm's stored rows, every block and element.
stored <- list(
  M20 = readRDS(file.path(root, "dev/simstudy/field-sample-ladder/results/simstudy-20261016-031826.rds")),
  cm = readRDS(file.path(root, "dev/simstudy/community-mean-ladder-20261016/results/simstudy-20261017-025051.rds")))
cols <- c("truth", "lower", "upper", "post_mean", "covered")
rep_rows <- do.call(rbind, lapply(diag_fits, function(x) {
  st <- if (x$arm == "M20") stored$M20 else stored$cm
  r0 <- if (!is.null(st$rows)) st$rows else st$params
  r0 <- r0[r0$scenario == x$arm & r0$replicate == x$replicate, ]
  m <- match(paste(x$rows$block, x$rows$element), paste(r0$block, r0$element))
  o0 <- st$occupancy[st$occupancy$scenario == x$arm & st$occupancy$replicate == x$replicate, ]
  o0 <- o0[order(o0$species), ]
  data.frame(arm = x$arm, replicate = x$replicate, rows = nrow(x$rows), stored_rows = nrow(r0),
             unmatched = sum(is.na(m)),
             identical = all(vapply(cols, function(k) identical(x$rows[[k]], r0[[k]][m]), TRUE)),
             max_abs_diff = max(vapply(cols, function(k) max(abs(x$rows[[k]] - r0[[k]][m])), 0)),
             occupancy_identical = identical(x$occupancy$occ_est, o0$occ_est) &&
               identical(x$occupancy$occ_true, o0$occ_true),
             rebuild_error = x$rebuild_error, stringsAsFactors = FALSE)
}))
reproduction <- do.call(rbind, lapply(split(rep_rows, rep_rows$arm), function(y) data.frame(
  arm = y$arm[1], replicates = nrow(y), rows_compared = sum(y$rows), unmatched = sum(y$unmatched),
  replicates_identical = sum(y$identical & y$rows == y$stored_rows & y$unmatched == 0),
  max_abs_diff = max(y$max_abs_diff), occupancy_identical = sum(y$occupancy_identical),
  max_rebuild_error = max(y$rebuild_error))))
all_rebuild <- max(vapply(c(fits, long_fits), `[[`, 0, "rebuild_error"))
print(reproduction, row.names = FALSE)
message(sprintf("largest rebuild error over all fits: %.1e", all_rebuild))
wcsv(reproduction, "reproduction.csv")

# 2. One row per species-fit.
d <- species_table(fits)
# The fitted and true standard deviations of r's parts, for ratios of means: a species with no latent-factor
# loading in the truth has a true latent-factor term of 0, so per-species ratios can be infinite.
parts <- c("r", "XB", "UL", "SP", "site")
raw <- do.call(rbind, lapply(fits, function(x) {
  sp <- x$species
  out <- data.frame(job = x$job, replicate = x$replicate, species = sp$species)
  for (p in parts) {
    out[[paste0("fit_sd_", p)]] <- sp[[paste0("sd_", p, "_mean")]]
    out[[paste0("true_sd_", p)]] <- sp[[paste0("sd_", p, "_true")]]
  }
  out
}))
d <- merge(d, raw, by = c("job", "replicate", "species"), sort = FALSE)
d$arm <- factor(d$arm, arms)
d <- d[order(d$arm, d$variant, d$replicate, d$species), ]
d$arm <- as.character(d$arm)
wcsv(d, "species-fits.csv")
variants <- unique(d$variant)

# 3. Decomposition by arm, variant and band of true B0 (and over all bands).
summ_band <- function(y) {
  m <- function(v) cl_mean(v, y$replicate)
  data.frame(
    n = nrow(y), B0_true = mean(y$B0_true),
    e = m(y$e)[1], e_se = m(y$e)[2], J = mean(y$J), V = mean(y$V), D = mean(y$D),
    P = mean(y$P), P_se = m(y$P)[2], Q = mean(y$Q), a_err = mean(y$a_err),
    coverage = mean(y$covered), coverage_se = m(y$covered)[2], width = mean(y$width),
    sd_ratio = mean(y$sd_ratio), var_ratio = mean(y$var_r_fit / y$var_r_true),
    sd_XB_ratio = sum(y$fit_sd_XB) / sum(y$true_sd_XB), sd_UL_ratio = sum(y$fit_sd_UL) / sum(y$true_sd_UL),
    sd_SP_ratio = sum(y$fit_sd_SP) / sum(y$true_sd_SP), sd_site_ratio = sum(y$fit_sd_site) / sum(y$true_sd_site),
    sd_r_fit = mean(y$fit_sd_r), sd_r_true = mean(y$sd_r_true), sd_XB_fit = mean(y$fit_sd_XB),
    sd_XB_true = mean(y$true_sd_XB), sd_site_fit = mean(y$fit_sd_site), sd_site_true = mean(y$true_sd_site),
    dvar_r = mean(y$var_r_true - y$var_r_fit), dvar_XB = mean(y$dvar_XB), dvar_site = mean(y$dvar_site),
    dcov2 = mean(2 * y$dcov), m_UL_err = mean(y$m_UL_err), m_SP_err = mean(y$m_SP_err), cor_r = mean(y$cor_r),
    psibar_err = mean(y$psibar_est - y$psibar_true), occ_err = mean(y$occ_est - y$occ_true),
    zbar_minus_psibar = mean(y$zbar_true - y$psibar_true))
}
decomp <- do.call(rbind, lapply(arms, function(k) do.call(rbind, lapply(variants, function(v) {
  x <- d[d$arm == k & d$variant == v, ]
  if (!nrow(x)) return(NULL)
  bands <- split(x, x$band_index)
  rbind(cbind(arm = k, variant = v, band = "all", band_index = 0L, extreme = NA, summ_band(x)),
        do.call(rbind, lapply(bands, function(y)
          cbind(arm = k, variant = v, band = y$band[1], band_index = y$band_index[1], extreme = y$extreme[1],
                summ_band(y)))))
}))))
rownames(decomp) <- NULL
wcsv(decomp, "decomposition.csv")

# 4. Shares of the between-species variance of e: cov(e, term) / var(e), which add to 1 for J, V and -D.
shares <- do.call(rbind, lapply(arms, function(k) do.call(rbind, lapply(variants, function(v) {
  x <- d[d$arm == k & d$variant == v, ]
  if (!nrow(x)) return(NULL)
  ve <- var(x$e)
  data.frame(arm = k, variant = v, n = nrow(x), sd_e = sqrt(ve), mean_e = mean(x$e),
             share_J = cov(x$e, x$J) / ve, share_V = cov(x$e, x$V) / ve, share_minus_D = -cov(x$e, x$D) / ve,
             share_P = cov(x$e, x$P) / ve, share_Q = cov(x$e, x$Q) / ve,
             cor_e_P = cor(x$e, x$P), cor_e_Q = cor(x$e, x$Q), cor_P_Q = cor(x$P, x$Q))
}))))
wcsv(shares, "variance-shares.csv")

# 5. Regressions of e on P, on Q, and on both, with standard errors clustered by replicate.
cluster_lm <- function(y, X, g) {
  X <- cbind(1, X)
  XtXi <- solve(crossprod(X))
  b <- XtXi %*% crossprod(X, y)
  u <- as.vector(y - X %*% b)
  meat <- Reduce(`+`, lapply(split(seq_along(y), g), function(i) {
    s <- crossprod(X[i, , drop = FALSE], u[i]); s %*% t(s)
  }))
  G <- length(unique(g))
  V <- G / (G - 1) * XtXi %*% meat %*% XtXi
  list(b = as.vector(b), se = sqrt(diag(V)), r2 = 1 - sum(u^2) / sum((y - mean(y))^2))
}
regs <- do.call(rbind, lapply(arms, function(k) do.call(rbind, lapply(variants, function(v) {
  x <- d[d$arm == k & d$variant == v, ]
  if (!nrow(x)) return(NULL)
  f <- function(model, X) {
    r <- cluster_lm(x$e, X, x$replicate)
    data.frame(arm = k, variant = v, model = model, intercept = r$b[1], intercept_se = r$se[1],
               b_P = if (model != "Q") r$b[2] else NA, b_P_se = if (model != "Q") r$se[2] else NA,
               b_Q = if (model == "Q") r$b[2] else if (model == "P+Q") r$b[3] else NA,
               b_Q_se = if (model == "Q") r$se[2] else if (model == "P+Q") r$se[3] else NA, r2 = r$r2)
  }
  rbind(f("P", x$P), f("Q", x$Q), f("P+Q", cbind(x$P, x$Q)))
}))))
wcsv(regs, "regression.csv")

# 6. The spread criterion, recomputed as choose.R computes it, on the diagnosis refits.
crit <- do.call(rbind, lapply(arms, function(k) {
  x <- d[d$arm == k & d$variant == "base", ]
  sl <- cluster_slope(x$e, x$P, x$replicate)
  a_holds <- unname(sl["slope"] >= 0.5 && sl["slope"] <= 1.5 && abs(sl["slope"]) > 1.96 * sl["se"])
  b <- lapply(extreme_bands(k), function(i) { y <- x[x$band_index == i, ]; c(e = mean(y$e), P = mean(y$P)) })
  b_holds <- all(vapply(b, function(v) sign(v[["P"]]) == sign(v[["e"]]) && abs(v[["P"]]) >= 0.5 * abs(v[["e"]]),
                        TRUE))
  data.frame(arm = k, slope = unname(sl["slope"]), slope_lo = unname(sl["slope"] - 1.96 * sl["se"]),
             slope_hi = unname(sl["slope"] + 1.96 * sl["se"]), r2 = unname(sl["r2"]),
             e_ext1 = b[[1]][["e"]], P_ext1 = b[[1]][["P"]], e_ext2 = b[[2]][["e"]], P_ext2 = b[[2]][["P"]],
             a_holds = a_holds, b_holds = b_holds, holds = a_holds && b_holds)
}))
wcsv(crit, "spread-criterion.csv")
print(crit, row.names = FALSE)

# 7. Interventions against the diagnosis refits, paired by species-fit.
rows_of <- function(fs) do.call(rbind, lapply(fs, function(x)
  cbind(job = x$job, arm = x$arm, variant = if (nzchar(x$variant)) x$variant else "base", x$rows)))
rr <- rows_of(fits)
bt <- rr$block == "beta_theta"
rr$block[bt] <- ifelse(rr$element[bt] %% 2 == 1, "beta_theta_intercept", "beta_theta_slope")
int <- do.call(rbind, lapply(arms, function(k) do.call(rbind, lapply(setdiff(variants, "base"), function(v) {
  b <- d[d$arm == k & d$variant == "base", ]
  x <- d[d$arm == k & d$variant == v, ]
  if (!nrow(x)) return(NULL)
  key <- function(y) paste(y$replicate, y$species)
  x <- x[match(key(b), key(x)), ]
  stopifnot(identical(x$B0_true, b$B0_true))
  pd <- function(f) cl_mean(f(x) - f(b), b$replicate)
  per_rep <- function(y, f) tapply(f(y), y$replicate, mean)
  comm_err <- function(y) per_rep(y, function(z) z$occ_est - z$occ_true)
  ext <- extreme_bands(k)
  band_e <- function(y, i) mean(y$e[y$band_index == i])
  band_sd <- function(y, i) mean(y$sd_ratio[y$band_index == i])
  mu_err <- function(y) { u <- tapply(y$mu_b0, y$replicate, `[`, 1) - tapply(y$B0_true, y$replicate, mean); mean(u) }
  blk <- function(vv, block, f) { z <- rr[rr$arm == k & rr$variant == vv & rr$block == block, ]; f(z) }
  out <- data.frame(
    arm = k, variant = v, species_fits = nrow(b),
    coverage_base = mean(b$covered), coverage_variant = mean(x$covered),
    coverage_diff = pd(function(y) y$covered)[1], coverage_diff_se = pd(function(y) y$covered)[2],
    width_base = mean(b$width), width_variant = mean(x$width),
    mean_e_base = mean(b$e), mean_e_variant = mean(x$e), mean_abs_e_base = mean(abs(b$e)),
    mean_abs_e_variant = mean(abs(x$e)),
    ext1_band = b$band[b$band_index == ext[1]][1], ext1_e_base = band_e(b, ext[1]), ext1_e_variant = band_e(x, ext[1]),
    ext2_band = b$band[b$band_index == ext[2]][1], ext2_e_base = band_e(b, ext[2]), ext2_e_variant = band_e(x, ext[2]),
    ext1_sd_ratio_base = band_sd(b, ext[1]), ext1_sd_ratio_variant = band_sd(x, ext[1]),
    ext2_sd_ratio_base = band_sd(b, ext[2]), ext2_sd_ratio_variant = band_sd(x, ext[2]),
    sd_ratio_base = mean(b$sd_ratio), sd_ratio_variant = mean(x$sd_ratio),
    Q_base = mean(b$Q), Q_variant = mean(x$Q), P_base = mean(b$P), P_variant = mean(x$P),
    occ_mae_base = mean(abs(b$occ_est - b$occ_true)), occ_mae_variant = mean(abs(x$occ_est - x$occ_true)),
    occ_mae_diff_se = pd(function(y) abs(y$occ_est - y$occ_true))[2],
    comm_err_base = mean(comm_err(b)), comm_err_variant = mean(comm_err(x)),
    mu_b0_err_base = if (all(is.na(b$mu_b0))) NA else mu_err(b),
    mu_b0_err_variant = if (all(is.na(x$mu_b0))) NA else mu_err(x),
    B_coverage_base = blk("base", "B", function(z) mean(z$covered)),
    B_coverage_variant = blk(v, "B", function(z) mean(z$covered)),
    theta_int_bias_base = blk("base", "beta_theta_intercept", function(z) mean(z$post_mean - z$truth)),
    theta_int_bias_variant = blk(v, "beta_theta_intercept", function(z) mean(z$post_mean - z$truth)),
    stringsAsFactors = FALSE)
  out$coverage_rises <- out$coverage_diff > 2 * out$coverage_diff_se
  out$ext_shrink_third <- abs(out$ext1_e_variant) <= 2 / 3 * abs(out$ext1_e_base) &
    abs(out$ext2_e_variant) <= 2 / 3 * abs(out$ext2_e_base)
  out$ext_halved <- abs(out$ext1_e_variant) <= 0.5 * abs(out$ext1_e_base) &
    abs(out$ext2_e_variant) <= 0.5 * abs(out$ext2_e_base)
  out$ext_fall_under_quarter <- abs(out$ext1_e_variant) > 0.75 * abs(out$ext1_e_base) &
    abs(out$ext2_e_variant) > 0.75 * abs(out$ext2_e_base)
  out$sd_towards_1 <- abs(out$ext1_sd_ratio_variant - 1) < abs(out$ext1_sd_ratio_base - 1) &
    abs(out$ext2_sd_ratio_variant - 1) < abs(out$ext2_sd_ratio_base - 1)
  out$guard_ok <- out$occ_mae_variant - out$occ_mae_base <= 0.005 &
    abs(out$comm_err_variant - out$comm_err_base) <= 0.01
  out$remedy_works <- out$coverage_rises & out$ext_shrink_third & out$sd_towards_1 & out$guard_ok
  out
}))))
if (!is.null(int)) {
  wcsv(int, "interventions.csv")
  print(int[, c("arm", "variant", "coverage_base", "coverage_variant", "ext1_e_base", "ext1_e_variant",
                "ext2_e_base", "ext2_e_variant", "remedy_works", "ext_halved")], row.names = FALSE)
}

# 8. Hyperparameters: averages over replicates of their posterior means, and their median bulk ESS.
hyper <- do.call(rbind, lapply(fits, function(x) do.call(rbind, lapply(names(x$hyper), function(h) data.frame(
  job = x$job, arm = x$arm, variant = if (nzchar(x$variant)) x$variant else "base", replicate = x$replicate,
  parameter = h, mean = x$hyper[[h]][["mean"]], lo = x$hyper[[h]][["lo"]], hi = x$hyper[[h]][["hi"]],
  rhat = x$hyper[[h]][["rhat"]], ess = x$hyper[[h]][["ess"]], stringsAsFactors = FALSE)))))
hyper_summ <- do.call(rbind, lapply(split(hyper, list(hyper$job, hyper$parameter), drop = TRUE), function(y)
  data.frame(arm = y$arm[1], variant = y$variant[1], parameter = y$parameter[1], replicates = nrow(y),
             mean = mean(y$mean), lo = mean(y$lo), hi = mean(y$hi), ess_median = median(y$ess),
             rhat_over_1.05 = mean(y$rhat > 1.05), stringsAsFactors = FALSE)))
hyper_summ <- hyper_summ[order(match(hyper_summ$arm, arms), hyper_summ$variant, hyper_summ$parameter), ]
wcsv(hyper_summ, "hyper.csv")

# 9. Long chains, if the mixing rule fired: each extreme band's mean e, P and sd ratio at both lengths.
if (length(long_fits)) {
  dl <- species_table(long_fits)
  dl$key <- paste(dl$job, dl$replicate, dl$species)
  ds <- d; ds$key <- paste(ds$job, ds$replicate, ds$species)
  ds <- ds[match(dl$key, ds$key), ]
  lc <- do.call(rbind, lapply(split(seq_len(nrow(dl)), list(dl$job, dl$band_index), drop = TRUE), function(i) {
    data.frame(job = dl$job[i[1]], band = dl$band[i[1]], extreme = dl$extreme[i[1]], n = length(i),
               e_short = mean(ds$e[i]), e_long = mean(dl$e[i]), P_short = mean(ds$P[i]), P_long = mean(dl$P[i]),
               sd_ratio_short = mean(ds$sd_ratio[i]), sd_ratio_long = mean(dl$sd_ratio[i]),
               coverage_short = mean(ds$covered[i]), coverage_long = mean(dl$covered[i]))
  }))
  wcsv(lc, "long-chains.csv")
}
