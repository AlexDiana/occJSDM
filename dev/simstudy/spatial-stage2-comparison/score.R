# Scores for stage 2: headroom share against the pilot oracle, pattern, slope, aggregation,
# paired comparisons and the development reading. Requires config.R.
PILOT_CODE <- "dev/simstudy/spatial-recoverability-pilot"
source(file.path(PILOT_CODE, "kernel.R")); source(file.path(PILOT_CODE, "predict.R"))

rmse <- function(a, b) sqrt(mean((a - b)^2))

score_arm <- function(prediction, env_prediction, oracle, truth, dataset_id, package, arm, status) {
  S <- length(truth$species)
  rows <- lapply(seq_len(S), function(s) {
    psi <- truth$psi_lattice[, s]; unif <- uninformed_reference(truth$eta_lattice[, s], truth$spatial_sd[s])
    r_env <- if (is.null(env_prediction)) NA_real_ else rmse(env_prediction$p_full[, s], psi); r_or <- rmse(oracle$probability_mean[, s], psi)
    base <- data.frame(dataset_id = dataset_id, package = package, arm = arm, key = arm_key(package, arm), species = truth$species[s],
                       occupancy = truth$occupancy[s], status = status, rmse_env = r_env, rmse_oracle = r_or, oracle_status = oracle$status,
                       stringsAsFactors = FALSE)
    if (is.null(prediction)) return(cbind(base, rmse = NA_real_, mae = NA_real_, headroom_share = NA_real_, pattern = NA_real_,
                                          field_slope = NA_real_, field_correlation = NA_real_))
    p <- prediction$p_full[, s]; r <- rmse(p, psi)
    share <- if (is.finite(r_env) && r_env - r_or > 1e-12) (r_env - r) / (r_env - r_or) else NA_real_
    true_contrast <- psi - unif; pattern <- if (sd(prediction$delta_p[, s]) > 0) cor(prediction$delta_p[, s], true_contrast) else NA_real_
    slope <- corr <- NA_real_
    if (identical(prediction$link, "logit")) {
      fc <- prediction$field_mean[, s] - mean(prediction$field_mean[, s]); tc <- truth$field_lattice[, s] - mean(truth$field_lattice[, s])
      slope <- sum(fc * tc) / sum(tc^2); corr <- if (sd(fc) > 0) cor(fc, tc) else NA_real_
    }
    cbind(base, rmse = r, mae = mean(abs(p - psi)), headroom_share = share, pattern = pattern, field_slope = slope, field_correlation = corr)
  })
  do.call(rbind, rows)
}

aggregate_arm <- function(scores, gate_occupancy) {
  s <- scores[scores$occupancy %in% gate_occupancy, ]
  key <- c("dataset_id", "package", "arm", "key")
  metrics <- c("headroom_share", "pattern", "field_slope", "field_correlation", "rmse", "rmse_env", "rmse_oracle")
  by_level <- aggregate(s[metrics], s[c(key, "occupancy")], function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_)
  out <- aggregate(by_level[metrics], by_level[key], function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_)
  skipped <- aggregate(list(skipped_species = !is.finite(s$headroom_share)), s[key], sum)
  status <- aggregate(list(status = s$status), s[key], function(x) x[1])
  out <- merge(merge(out, skipped, by = key), status, by = key)
  out$cell <- sub("-rep[0-9]+$", "", out$dataset_id); out
}

paired_comparisons <- function(agg, config) {
  rows <- list()
  for (cell in unique(agg$cell)) for (cmp in config$comparisons) {
    a <- agg[agg$cell == cell & agg$key == cmp[1], ]; b <- agg[agg$cell == cell & agg$key == cmp[2], ]
    m <- merge(a, b, by = "dataset_id", suffixes = c("_a", "_b"))
    ok <- m$status_a == "pass" & m$status_b == "pass" & is.finite(m$headroom_share_a) & is.finite(m$headroom_share_b)
    d <- m$headroom_share_b[ok] - m$headroom_share_a[ok]; k <- length(d)
    est <- if (k) mean(d) else NA_real_; se <- if (k > 1) sd(d) / sqrt(k) else NA_real_
    lo <- est - qt(0.975, max(k - 1, 1)) * se; hi <- est + qt(0.975, max(k - 1, 1)) * se
    verdict <- if (k < config$min_paired) "inconclusive" else if (est >= config$delta && lo > 0) "b clearly better" else
      if (est <= -config$delta && hi < 0) "a clearly better" else "no clear difference"
    rows[[length(rows) + 1L]] <- data.frame(cell = cell, a = cmp[1], b = cmp[2], paired = k, excluded = nrow(m) - k,
                                            estimate = est, lower = lo, upper = hi, verdict = verdict, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# `pairs` holds every comparison by cell; `max_share` is the largest headroom share any arm reached in any cell.
# Branches follow the spec's decision rule at the 900-site cells. "Any factor arm" is the operationalisation of "factor arms".
development_reading <- function(pairs, max_share) {
  gp <- pairs[pairs$a == "occJSDM/gp" & pairs$b %in% c("Hmsc/spatial2", "Hmsc/spatial6", "gllvm/spatial2"), ]
  t1 <- gp[grepl("^t1", gp$cell) & grepl("n900", gp$cell), ]; t2 <- gp[grepl("^t2", gp$cell) & grepl("n900", gp$cell), ]
  factors_better_t1 <- any(t1$verdict == "b clearly better")   # gp is arm a, a factor arm is arm b
  gp_better_t1 <- any(t1$verdict == "a clearly better"); gp_better_t2 <- any(t2$verdict == "a clearly better")
  nothing_anywhere <- all(pairs$verdict %in% c("no clear difference", "inconclusive"))
  reading <- if (factors_better_t1 && !gp_better_t2) "shared factors" else
    if (!factors_better_t1 && gp_better_t2) "per-species ranges" else
    if (factors_better_t1 && gp_better_t2) "trade-off" else
    if (nothing_anywhere && is.finite(max_share) && max_share < 0.5) "representation not the bottleneck" else
    if (nothing_anywhere) "no clear difference" else "no pre-specified branch"
  n900 <- pairs[grepl("n900", pairs$cell), ]
  list(reading = reading, truth1_900 = t1, truth2_900 = t2, max_share = max_share,
       # hmsc6 clearly better than hmsc2 under truth 2: factor dimension, not representation, drives the two-factor shortfall.
       factor_dimension = any(n900$a == "Hmsc/spatial6" & n900$b == "Hmsc/spatial2" & grepl("^t2", n900$cell) & n900$verdict == "a clearly better"),
       # hmsc2 clearly different from gllvm2 in either direction: a fitting-method effect that qualifies every factor-arm conclusion.
       fitting_method = any(n900$a == "Hmsc/spatial2" & n900$b == "gllvm/spatial2" & n900$verdict %in% c("a clearly better", "b clearly better")))
}
