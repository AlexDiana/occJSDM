# Per-species scores, aggregation, gate and map reading for the recoverability pilot.
# Requires kernel.R, generator.R, oracle.R and predict.R to be sourced first.
default_thresholds <- function() list(pattern_median = 0.6, pattern_floor = 0.4, pattern_count = 8L,
                                      slope_low = 0.5, slope_high = 1.5, headroom_pp = 1.0, communities = 10L)

score_one <- function(truth_field, pred_field, truth_psi, pred_psi, uninformed) {
  tc <- truth_field - mean(truth_field); fc <- pred_field - mean(pred_field)
  true_contrast <- truth_psi - uninformed; pred_contrast <- pred_psi - uninformed
  rmse <- function(a, b) sqrt(mean((a - b)^2))
  rmse_unif <- rmse(uninformed, truth_psi); rmse_or <- rmse(pred_psi, truth_psi)
  safe_cor <- function(a, b) if (sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_
  c(contrast_correlation = safe_cor(pred_contrast, true_contrast), contrast_rms = sqrt(mean(pred_contrast^2)),
    true_contrast_rms = sqrt(mean(true_contrast^2)),
    field_correlation = safe_cor(fc, tc), field_slope = sum(fc * tc) / sum(tc^2), field_rmse = rmse(fc, tc),
    zero_field_rmse = sqrt(mean(tc^2)),
    probability_rmse_oracle = rmse_or, probability_rmse_uninformed = rmse_unif,
    probability_mae_oracle = mean(abs(pred_psi - truth_psi)), probability_mae_uninformed = mean(abs(uninformed - truth_psi)),
    headroom_pp = 100 * (rmse_unif - rmse_or),
    headroom_fraction = if (rmse_unif > 0) (rmse_unif - rmse_or) / rmse_unif else NA_real_)
}

# `screen` is a data frame with columns species (integer index) and status ("pass" or "unresolved").
score_species <- function(cm, arrangement, prediction, screen) {
  lat <- cm$index$lattice; idx <- cm$index[[arrangement]]
  eta_lat <- training_offset(cm, lat); eta_site <- training_offset(cm, idx)
  sp_idx <- if (!is.null(prediction$species)) prediction$species else seq_len(nrow(cm$species))
  rows <- lapply(seq_along(sp_idx), function(j) {
    s <- sp_idx[j]
    held <- score_one(cm$spatial[lat, s], prediction$field_mean[, j], cm$psi[lat, s], prediction$probability_mean[, j],
                      uninformed_reference(eta_lat[, s], cm$spatial_sd[s]))
    train <- score_one(cm$spatial[idx, s], prediction$site_field_mean[, j], cm$psi[idx, s], prediction$site_probability_mean[, j],
                       uninformed_reference(eta_site[, s], cm$spatial_sd[s]))
    names(train) <- paste0("train_", names(train))
    data.frame(truth = cm$truth, cell = cm$cell, n = cm$n, range = cm$range, replicate = cm$replicate,
               arrangement = arrangement, species = cm$species$species[s], occupancy = cm$species$occupancy[s],
               class = if (cm$truth == "t2") cm$species$class[s] else NA_character_,
               status = screen$status[screen$species == s][1], as.list(held), as.list(train), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Average species within occupancy level, then pool the 25, 50 and 75% levels, per community,
# arrangement and (under truth 2) dispersal class.
aggregate_scores <- function(scores) {
  gate_levels <- c(0.25, 0.50, 0.75)
  s <- scores[scores$occupancy %in% gate_levels, ]
  s$group <- ifelse(s$truth == "t2", s$class, "all")
  key <- c("truth", "cell", "n", "range", "replicate", "arrangement", "group")
  s$range[is.na(s$range)] <- -1   # aggregate() drops NA keys; restore below
  metrics <- c("contrast_correlation", "field_slope", "headroom_pp", "headroom_fraction", "field_correlation", "field_rmse", "zero_field_rmse")
  by_level <- aggregate(s[metrics], s[c(key, "occupancy")], mean)
  out <- aggregate(by_level[metrics], by_level[key], mean)
  status <- aggregate(list(status = s$status), s[key], function(x) if (any(x == "unresolved")) "unresolved" else "pass")
  out <- merge(out, status, by = key)
  out$range[out$range == -1] <- NA_real_
  out
}

# One cell (and class under truth 2), uniform arrangement, one community per row.
gate_cells <- function(agg, thresholds) {
  stopifnot(nrow(agg) >= 2L)
  if (any(agg$status == "unresolved"))
    return(data.frame(decision = "unresolved", pattern_median = NA_real_, pattern_count = NA_integer_, slope_median = NA_real_,
                      headroom_mean = NA_real_, headroom_low = NA_real_, headroom_high = NA_real_, communities = nrow(agg)))
  pm <- median(agg$contrast_correlation); pc <- sum(agg$contrast_correlation >= thresholds$pattern_floor)
  sm <- median(agg$field_slope)
  hm <- mean(agg$headroom_pp); se <- sd(agg$headroom_pp) / sqrt(nrow(agg)); tq <- qt(0.975, nrow(agg) - 1L)
  lo <- hm - tq * se; hi <- hm + tq * se
  ok <- pm >= thresholds$pattern_median && pc >= thresholds$pattern_count &&
    sm >= thresholds$slope_low && sm <= thresholds$slope_high && hm >= thresholds$headroom_pp && lo > 0
  data.frame(decision = if (ok) "pass" else "fail", pattern_median = pm, pattern_count = pc, slope_median = sm,
             headroom_mean = hm, headroom_low = lo, headroom_high = hi, communities = nrow(agg))
}

# Spec section 3, "Reading the map". `map` has columns truth, range, class, n, decision.
read_map <- function(map) {
  t1 <- map[map$truth == "t1" & !is.na(map$range) & map$range == 0.08 & map$decision == "pass", ]
  t2 <- map[map$truth == "t2" & !is.na(map$class) & map$class == "medium" & map$decision == "pass", ]
  t1_n <- if (nrow(t1)) as.integer(min(t1$n)) else NA_integer_; t2_n <- if (nrow(t2)) as.integer(min(t2$n)) else NA_integer_
  if (!is.na(t1_n) && !is.na(t2_n)) list(branch = "both", recommended_n = max(t1_n, t2_n), truth1_n = t1_n, truth2_n = t2_n)
  else if (!is.na(t1_n)) list(branch = "truth1_only", recommended_n = t1_n, truth1_n = t1_n, truth2_n = NA_integer_)
  else if (!is.na(t2_n)) list(branch = "truth2_only", recommended_n = t2_n, truth1_n = NA_integer_, truth2_n = t2_n)
  else list(branch = "stop", recommended_n = NA_integer_, truth1_n = NA_integer_, truth2_n = NA_integer_)
}
