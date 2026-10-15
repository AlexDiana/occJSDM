# Summarise the rare-dominated arms against replicates 1 to 100 of their parent cells in the R = 300 paper study.
# Run from the repository root: Rscript dev/simstudy/rare-dominated-20261015/analyse.R
# Reads the run's results file, or its checkpoint while the run is in flight.
# Writes occupancy-by-band.csv, occupancy-by-band-parents-r300.csv, b0-by-truth.csv, community.csv, b0-shift-recovered.csv, pairing-check.csv, detection-paired.csv, blocks.csv and matched-truth.csv.
here <- "dev/simstudy/rare-dominated-20261015"
out <- Sys.getenv("OUT_DIR", here)  # where the CSVs go; the run directory by default
parent_file <- "dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds"
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
new <- readRDS(if (length(res)) res[1] else
  list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1])
new_rows <- if (!is.null(new$rows)) new$rows else new$params
par <- readRDS(parent_file)
pairs <- c(base = "base_rare", occupancy = "occupancy_rare", binary = "binary_rare")
shift <- -2.75
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
d <- rbind(par$rows[par$rows$scenario %in% names(pairs) & par$rows$replicate <= 100, cols],
           new_rows[new_rows$scenario %in% pairs & new_rows$block != "ERROR", cols])
occ_cols <- c("scenario", "replicate", "species", "occ_true", "occ_est")
o <- rbind(par$occupancy[par$occupancy$scenario %in% names(pairs) & par$occupancy$replicate <= 100, occ_cols],
           new$occupancy[new$occupancy$scenario %in% pairs, occ_cols])
# Use only replicates present in both arms of a pair, so a partial run still compares like with like.
keep <- function(x) do.call(rbind, lapply(names(pairs), function(p) {
  reps <- intersect(unique(d$replicate[d$scenario == p & d$block == "B0"]),
                    unique(d$replicate[d$scenario == pairs[[p]] & d$block == "B0"]))
  x[x$scenario %in% c(p, pairs[[p]]) & x$replicate %in% reps, ]
}))
d <- keep(d); o <- keep(o)
cells <- c(rbind(names(pairs), unname(pairs)))
for (p in names(pairs)) cat(sprintf("%s / %s: %d paired replicates\n", p, pairs[[p]],
  length(unique(d$replicate[d$scenario == pairs[[p]]]))))

# Mean of x over species-fits with a standard error that treats each replicate as one cluster,
# since the species within a replicate share data, latent factors and a fit.
cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  se <- if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA
  c(m, se)
}
ord <- function(x) x[order(match(x$cell, cells)), ]

# 1. Species-level occupancy (share of sites occupied) by true-occupancy band.
# occ_est is the posterior mean of z averaged over sites; binary has no rows, as z is observed there.
breaks <- c(0, .05, .1, .2, .3, .5, .7, .9, 1)
o$band <- cut(o$occ_true, breaks, include.lowest = TRUE)
o$err <- o$occ_est - o$occ_true
ob <- do.call(rbind, lapply(split(o, list(o$scenario, o$band), drop = TRUE), function(x) {
  b <- cl_mean(x$err, x$replicate)
  data.frame(cell = x$scenario[1], band = as.character(x$band[1]), n = nrow(x),
             true_mean = mean(x$occ_true), est_mean = mean(x$occ_est), bias = b[1], bias_se = b[2],
             est_over_true = mean(x$occ_est) / mean(x$occ_true), share_over = mean(x$err > 0))
}))
ob <- ord(ob[order(ob$cell, ob$true_mean), ]); rownames(ob) <- NULL
write.csv(ob, file.path(out, "occupancy-by-band.csv"), row.names = FALSE)
# The parents have few rare species in replicates 1 to 100, so also band them over all 300 (unpaired).
o3 <- par$occupancy[par$occupancy$scenario %in% c("base", "occupancy"), occ_cols]
o3$band <- cut(o3$occ_true, breaks, include.lowest = TRUE)
ob3 <- do.call(rbind, lapply(split(o3, list(o3$scenario, o3$band), drop = TRUE), function(x) {
  b <- cl_mean(x$occ_est - x$occ_true, x$replicate)
  data.frame(cell = x$scenario[1], band = as.character(x$band[1]), n = nrow(x),
             true_mean = mean(x$occ_true), est_mean = mean(x$occ_est), bias = b[1], bias_se = b[2])
}))
ob3 <- ob3[order(ob3$cell, ob3$true_mean), ]; rownames(ob3) <- NULL
write.csv(ob3, file.path(out, "occupancy-by-band-parents-r300.csv"), row.names = FALSE)

# 2. B0 bias and coverage by its true value.
b0 <- d[d$block == "B0", ]
b0$band <- cut(b0$truth, c(-Inf, -5, -4, -3, -2, -1, 0, 1, 2, Inf))
b0$err <- b0$post_mean - b0$truth
bb <- do.call(rbind, lapply(split(b0, list(b0$scenario, b0$band), drop = TRUE), function(x) {
  b <- cl_mean(x$err, x$replicate); cv <- cl_mean(x$covered, x$replicate)
  data.frame(cell = x$scenario[1], b0_band = as.character(x$band[1]), n = nrow(x),
             true_mean = mean(x$truth), est_mean = mean(x$post_mean), bias = b[1], bias_se = b[2],
             coverage = cv[1], coverage_se = cv[2], mean_width = mean(x$upper - x$lower),
             share_above = mean(x$truth < x$lower))
}))
bb <- ord(bb[order(bb$cell, bb$true_mean), ]); rownames(bb) <- NULL
write.csv(bb, file.path(out, "b0-by-truth.csv"), row.names = FALSE)

# 3. The community as a whole: per-replicate means over its 10 species, then their mean and SE over replicates.
per_rep <- function(x, f) sapply(split(x, x$replicate), f)
msd <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
cm <- do.call(rbind, lapply(cells, function(k) {
  x <- b0[b0$scenario == k, ]; y <- o[o$scenario == k, ]
  b0e <- msd(per_rep(x, function(z) mean(z$post_mean - z$truth)))
  cv <- msd(per_rep(x, function(z) mean(z$covered)))
  r <- data.frame(cell = k, R = length(unique(x$replicate)),
    b0_true = mean(x$truth), b0_est = mean(x$post_mean), b0_bias = b0e[1], b0_bias_se = b0e[2],
    b0_coverage = cv[1], b0_coverage_se = cv[2],
    b0_slope = unname(coef(lm(post_mean ~ truth, x))[2]),
    occ_true = NA, occ_est = NA, occ_bias = NA, occ_bias_se = NA, occ_rel_bias = NA, occ_slope = NA,
    share10_true = NA, share10_est = NA, share20_true = NA, share20_est = NA)
  if (nrow(y)) {
    ce <- msd(per_rep(y, function(z) mean(z$occ_est) - mean(z$occ_true)))
    rel <- per_rep(y, function(z) mean(z$occ_est) / mean(z$occ_true) - 1)
    r[c("occ_true", "occ_est", "occ_bias", "occ_bias_se", "occ_rel_bias", "occ_slope")] <-
      list(mean(y$occ_true), mean(y$occ_est), ce[1], ce[2], mean(rel),
           unname(coef(lm(occ_est ~ occ_true, y))[2]))
    r[c("share10_true", "share10_est", "share20_true", "share20_est")] <-
      list(mean(y$occ_true <= .1), mean(y$occ_est <= .1), mean(y$occ_true <= .2), mean(y$occ_est <= .2))
  }
  r
}))
write.csv(cm, file.path(out, "community.csv"), row.names = FALSE)

# 4. Paired: how much of the 2.75 drop in true B0 do the fits recover, species by species?
sr <- do.call(rbind, lapply(names(pairs), function(p) {
  m <- merge(b0[b0$scenario == p, ], b0[b0$scenario == pairs[[p]], ], by = c("replicate", "element"))
  dv <- per_rep(data.frame(replicate = m$replicate, v = m$post_mean.y - m$post_mean.x), function(z) mean(z$v))
  e <- msd(dv)
  data.frame(parent = p, arm = pairs[[p]], n = nrow(m), est_shift = e[1], est_shift_se = e[2],
             share_recovered = e[1] / shift, share_recovered_se = e[2] / abs(shift))
}))
write.csv(sr, file.path(out, "b0-shift-recovered.csv"), row.names = FALSE)

# 5. Pairing: true B0 is the parent's plus the shift; every other truth is the parent's.
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q")
pk <- do.call(rbind, lapply(names(pairs), function(p) do.call(rbind, lapply(blocks, function(b) {
  m <- merge(d[d$scenario == p & d$block == b, ], d[d$scenario == pairs[[p]] & d$block == b, ],
             by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  expected <- if (b == "B0") m$truth.x + shift else m$truth.x
  data.frame(arm = pairs[[p]], block = b, n = nrow(m), max_abs_diff_from_expected = max(abs(m$truth.y - expected)))
}))))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

# 6. The detection side, whose truth is identical in each pair: does the rare arm's fit move it?
# beta_theta is 2 x S, flattened by column, so odd elements are the collection intercepts.
det <- d[d$block %in% c("beta_theta", "theta0", "p", "q"), ]
det$block[det$block == "beta_theta"] <- ifelse(det$element[det$block == "beta_theta"] %% 2 == 1,
                                               "beta_theta_intercept", "beta_theta_slope")
dt <- do.call(rbind, lapply(names(pairs), function(p) do.call(rbind, lapply(unique(det$block), function(b) {
  x <- det[det$scenario == p & det$block == b, ]; y <- det[det$scenario == pairs[[p]] & det$block == b, ]
  m <- merge(x, y, by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  bx <- cl_mean(m$post_mean.x - m$truth.x, m$replicate); by <- cl_mean(m$post_mean.y - m$truth.y, m$replicate)
  dd <- cl_mean(m$post_mean.y - m$post_mean.x, m$replicate)
  data.frame(parent = p, arm = pairs[[p]], block = b, n = nrow(m),
             bias_parent = bx[1], bias_arm = by[1], bias_arm_se = by[2], paired_diff = dd[1], paired_diff_se = dd[2],
             coverage_parent = mean(m$covered.x), coverage_arm = mean(m$covered.y))
}))))
write.csv(dt, file.path(out, "detection-paired.csv"), row.names = FALSE)

# 7. Every scored block, parent against arm, for side effects of the shift beyond B0.
blk <- do.call(rbind, lapply(cells, function(k) do.call(rbind, lapply(c("B0", "B", "G", "beta_theta", "theta0", "p", "q"), function(b) {
  x <- d[d$scenario == k & d$block == b, ]; if (!nrow(x)) return(NULL)
  cv <- msd(per_rep(x, function(z) mean(z$covered)))
  data.frame(cell = k, block = b, n = nrow(x), coverage = cv[1], coverage_se = cv[2],
             bias = mean(x$post_mean - x$truth), rmse = sqrt(mean((x$post_mean - x$truth)^2)),
             mean_width = mean(x$upper - x$lower))
}))))
write.csv(blk, file.path(out, "blocks.csv"), row.names = FALSE)

# 8. Like for like: bias at the same true value, in narrow bands where parent and arm overlap.
# Parents over all 300 replicates here, unpaired, to give them enough species at low values.
o3all <- rbind(par$occupancy[par$occupancy$scenario %in% c("base", "occupancy"), occ_cols], o[o$scenario %in% pairs, occ_cols])
o3all$band <- cut(o3all$occ_true, seq(.1, .4, by = .05))
b3all <- rbind(par$rows[par$rows$block == "B0" & par$rows$scenario %in% names(pairs), cols], b0[b0$scenario %in% pairs, cols])
b3all$band <- cut(b3all$truth, seq(-3, -0.5, by = .5))
mt <- rbind(
  do.call(rbind, lapply(split(o3all, list(o3all$scenario, o3all$band), drop = TRUE), function(x) {
    b <- cl_mean(x$occ_est - x$occ_true, x$replicate)
    data.frame(quantity = "occupancy", cell = x$scenario[1], band = as.character(x$band[1]), n = nrow(x),
               true_mean = mean(x$occ_true), bias = b[1], bias_se = b[2], coverage = NA)
  })),
  do.call(rbind, lapply(split(b3all, list(b3all$scenario, b3all$band), drop = TRUE), function(x) {
    b <- cl_mean(x$post_mean - x$truth, x$replicate)
    data.frame(quantity = "B0", cell = x$scenario[1], band = as.character(x$band[1]), n = nrow(x),
               true_mean = mean(x$truth), bias = b[1], bias_se = b[2], coverage = mean(x$covered))
  })))
mt <- mt[order(mt$quantity, mt$band, match(mt$cell, cells)), ]; rownames(mt) <- NULL
write.csv(mt, file.path(out, "matched-truth.csv"), row.names = FALSE)

options(width = 160)
cat("\n== occupancy by true band\n"); print(ob, digits = 3, row.names = FALSE)
cat("\n== parents over all 300 replicates\n"); print(ob3, digits = 3, row.names = FALSE)
cat("\n== B0 by true value\n"); print(bb, digits = 3, row.names = FALSE)
cat("\n== community\n"); print(cm, digits = 3, row.names = FALSE)
cat("\n== B0 shift recovered\n"); print(sr, digits = 3, row.names = FALSE)
cat("\n== pairing check\n"); print(pk, digits = 3, row.names = FALSE)
cat("\n== detection side, paired\n"); print(dt, digits = 3, row.names = FALSE)
cat("\n== blocks\n"); print(blk, digits = 3, row.names = FALSE)
cat("\n== matched truth\n"); print(mt, digits = 3, row.names = FALSE)
