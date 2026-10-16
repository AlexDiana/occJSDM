# Compare the six community-mean arms with their default-prior counterparts on identical datasets.
# Run from the repository root: Rscript dev/simstudy/b0-community-mean-20261016/analyse.R
# Reads this run's results file (or its checkpoint while the run is in flight), the rare-dominated
# run's results and replicates 1 to 100 of the R = 300 paper study.
# Writes pairing-check.csv, community.csv, occupancy-by-band.csv, b0-by-truth.csv,
# detection-paired.csv, blocks.csv, tradeoff-by-band.csv and criteria.csv (C1 to C11 of README.md).
here <- "dev/simstudy/b0-community-mean-20261016"
out <- Sys.getenv("OUT_DIR", here)
paper_file <- "dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds"
rare_file <- "dev/simstudy/rare-dominated-20261015/results/simstudy-20261016-003851.rds"
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
new <- readRDS(if (length(res)) res[1] else
  list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1])
new_rows <- if (!is.null(new$rows)) new$rows else new$params
paper <- readRDS(paper_file); rare <- readRDS(rare_file)
pairs <- c(base_cm = "base", occupancy_cm = "occupancy", binary_cm = "binary",
           base_rare_cm = "base_rare", occupancy_rare_cm = "occupancy_rare", binary_rare_cm = "binary_rare")
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
occ_cols <- c("scenario", "replicate", "species", "occ_true", "occ_est")
d <- rbind(paper$rows[paper$rows$scenario %in% pairs & paper$rows$replicate <= 100, cols],
           rare$rows[rare$rows$scenario %in% pairs, cols],
           new_rows[new_rows$scenario %in% names(pairs) & new_rows$block != "ERROR", cols])
o <- rbind(paper$occupancy[paper$occupancy$scenario %in% pairs & paper$occupancy$replicate <= 100, occ_cols],
           rare$occupancy[rare$occupancy$scenario %in% pairs, occ_cols],
           new$occupancy[new$occupancy$scenario %in% names(pairs), occ_cols])
# Only replicates present in both arms of a pair, so a partial run still compares like with like.
keep <- function(x) do.call(rbind, lapply(names(pairs), function(a) {
  reps <- intersect(unique(d$replicate[d$scenario == a & d$block == "B0"]),
                    unique(d$replicate[d$scenario == pairs[[a]] & d$block == "B0"]))
  x[x$scenario %in% c(a, pairs[[a]]) & x$replicate %in% reps, ]
}))
d <- keep(d); o <- keep(o)
cells <- c(rbind(unname(pairs), names(pairs)))
for (a in names(pairs)) cat(sprintf("%s / %s: %d paired replicates\n", pairs[[a]], a,
  length(unique(d$replicate[d$scenario == a]))))

# Mean over species-fits with a standard error treating each replicate as one cluster.
cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  c(m, if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA)
}
per_rep <- function(x, f) sapply(split(x, x$replicate), f)
msd <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
ord <- function(x) x[order(match(x$cell, cells)), ]

# 1. Pairing: every truth, and each species' true occupancy, identical within each pair.
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q")
pk <- do.call(rbind, lapply(names(pairs), function(a) rbind(
  do.call(rbind, lapply(blocks, function(b) {
    m <- merge(d[d$scenario == pairs[[a]] & d$block == b, ], d[d$scenario == a & d$block == b, ],
               by = c("replicate", "element"))
    if (!nrow(m)) return(NULL)
    data.frame(arm = a, quantity = b, n = nrow(m), identical = all(m$truth.x == m$truth.y))
  })),
  if (any(o$scenario == a)) {
    m <- merge(o[o$scenario == pairs[[a]], ], o[o$scenario == a, ], by = c("replicate", "species"))
    data.frame(arm = a, quantity = "occ_true", n = nrow(m), identical = all(m$occ_true.x == m$occ_true.y))
  })))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

# 2. The community as a whole, per replicate first. The collection intercept is beta_theta's
# odd elements (2 x S, flattened by column); binary cells have none.
b0 <- d[d$block == "B0", ]
bt0 <- d[d$block == "beta_theta" & d$element %% 2 == 1, ]
mu <- d[d$block == "mu_b0", ]
cm <- do.call(rbind, lapply(cells, function(k) {
  x <- b0[b0$scenario == k, ]; y <- o[o$scenario == k, ]; z <- bt0[bt0$scenario == k, ]; u <- mu[mu$scenario == k, ]
  b0e <- msd(per_rep(x, function(v) mean(v$post_mean - v$truth)))
  cv <- msd(per_rep(x, function(v) mean(v$covered)))
  r <- data.frame(cell = k, R = length(unique(x$replicate)),
    b0_true = mean(x$truth), b0_est = mean(x$post_mean), b0_bias = b0e[1], b0_bias_se = b0e[2],
    b0_coverage = cv[1], b0_coverage_se = cv[2], b0_slope = unname(coef(lm(post_mean ~ truth, x))[2]),
    bt0_bias = NA, bt0_bias_se = NA, bt0_b0_cor = NA,
    mu_b0_truth = NA, mu_b0_est = NA, mu_b0_coverage = NA, mu_b0_width = NA,
    occ_true = NA, occ_est = NA, occ_bias = NA, occ_bias_se = NA, occ_slope = NA,
    share10_true = NA, share10_est = NA, share20_true = NA, share20_est = NA)
  if (nrow(z)) {
    ze <- per_rep(z, function(v) mean(v$post_mean - v$truth))
    r[c("bt0_bias", "bt0_bias_se")] <- msd(ze)
    # Across replicates: does a community whose B0 is overestimated get its collection intercepts underestimated?
    r$bt0_b0_cor <- cor(ze, per_rep(x, function(v) mean(v$post_mean - v$truth))[names(ze)])
  }
  if (nrow(u)) r[c("mu_b0_truth", "mu_b0_est", "mu_b0_coverage", "mu_b0_width")] <-
    list(mean(u$truth), mean(u$post_mean), mean(u$covered), mean(u$upper - u$lower))
  if (nrow(y)) {
    ce <- msd(per_rep(y, function(v) mean(v$occ_est) - mean(v$occ_true)))
    r[c("occ_true", "occ_est", "occ_bias", "occ_bias_se", "occ_slope")] <-
      list(mean(y$occ_true), mean(y$occ_est), ce[1], ce[2], unname(coef(lm(occ_est ~ occ_true, y))[2]))
    r[c("share10_true", "share10_est", "share20_true", "share20_est")] <-
      list(mean(y$occ_true <= .1), mean(y$occ_est <= .1), mean(y$occ_true <= .2), mean(y$occ_est <= .2))
  }
  r
}))
write.csv(cm, file.path(out, "community.csv"), row.names = FALSE)

# 3. Species occupancy by true band, each arm against its counterpart on the same species.
breaks <- c(0, .05, .1, .2, .3, .5, .7, .9, 1)
ob <- do.call(rbind, lapply(names(pairs), function(a) {
  m <- merge(o[o$scenario == pairs[[a]], ], o[o$scenario == a, ], by = c("replicate", "species"))
  if (!nrow(m)) return(NULL)
  m$band <- cut(m$occ_true.x, breaks, include.lowest = TRUE)
  do.call(rbind, lapply(split(m, m$band, drop = TRUE), function(x) {
    bd <- cl_mean(x$occ_est.x - x$occ_true.x, x$replicate); bc <- cl_mean(x$occ_est.y - x$occ_true.x, x$replicate)
    df <- cl_mean(x$occ_est.y - x$occ_est.x, x$replicate)
    data.frame(arm = a, band = as.character(x$band[1]), n = nrow(x), true_mean = mean(x$occ_true.x),
               est_default = mean(x$occ_est.x), est_community = mean(x$occ_est.y),
               bias_default = bd[1], bias_default_se = bd[2], bias_community = bc[1], bias_community_se = bc[2],
               paired_diff = df[1], paired_diff_se = df[2])
  }))
}))
rownames(ob) <- NULL
write.csv(ob, file.path(out, "occupancy-by-band.csv"), row.names = FALSE)

# 4. B0 by its true value.
b0$band <- cut(b0$truth, c(-Inf, -5, -4, -3, -2, -1, 0, 1, 2, Inf))
bb <- do.call(rbind, lapply(split(b0, list(b0$scenario, b0$band), drop = TRUE), function(x) {
  b <- cl_mean(x$post_mean - x$truth, x$replicate); cv <- cl_mean(x$covered, x$replicate)
  data.frame(cell = x$scenario[1], b0_band = as.character(x$band[1]), n = nrow(x),
             true_mean = mean(x$truth), est_mean = mean(x$post_mean), bias = b[1], bias_se = b[2],
             coverage = cv[1], coverage_se = cv[2], mean_width = mean(x$upper - x$lower))
}))
bb <- ord(bb[order(bb$cell, bb$true_mean), ]); rownames(bb) <- NULL
write.csv(bb, file.path(out, "b0-by-truth.csv"), row.names = FALSE)

# 5. The detection side and B0, each arm against its counterpart, element by element.
det <- d[d$block %in% c("B0", "beta_theta", "theta0", "p", "q"), ]
det$block[det$block == "beta_theta"] <- ifelse(det$element[det$block == "beta_theta"] %% 2 == 1,
                                               "beta_theta_intercept", "beta_theta_slope")
dt <- do.call(rbind, lapply(names(pairs), function(a) do.call(rbind, lapply(unique(det$block), function(b) {
  m <- merge(det[det$scenario == pairs[[a]] & det$block == b, ], det[det$scenario == a & det$block == b, ],
             by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  bx <- cl_mean(m$post_mean.x - m$truth.x, m$replicate); by <- cl_mean(m$post_mean.y - m$truth.y, m$replicate)
  dd <- cl_mean(m$post_mean.y - m$post_mean.x, m$replicate); dc <- cl_mean(m$covered.y - m$covered.x, m$replicate)
  data.frame(arm = a, block = b, n = nrow(m), bias_default = bx[1], bias_community = by[1],
             bias_community_se = by[2], paired_diff = dd[1], paired_diff_se = dd[2],
             coverage_default = mean(m$covered.x), coverage_community = mean(m$covered.y),
             coverage_diff = dc[1], coverage_diff_se = dc[2],
             width_default = mean(m$upper.x - m$lower.x), width_community = mean(m$upper.y - m$lower.y))
}))))
write.csv(dt, file.path(out, "detection-paired.csv"), row.names = FALSE)

# 6. Every scored block.
blk <- do.call(rbind, lapply(cells, function(k) do.call(rbind, lapply(
  c("B0", "mu_b0", "B", "G", "sigma_b", "beta_theta", "theta0", "p", "q", "resid_cor"), function(b) {
  x <- d[d$scenario == k & d$block == b, ]
  if (b == "resid_cor") x <- {
    y <- rbind(paper$rows, rare$rows, new_rows)
    y[y$scenario == k & y$block == b & y$replicate %in% unique(d$replicate[d$scenario == k]), cols]
  }
  if (!nrow(x)) return(NULL)
  cv <- msd(per_rep(x, function(v) mean(v$covered)))
  data.frame(cell = k, block = b, n = nrow(x), coverage = cv[1], coverage_se = cv[2],
             bias = mean(x$post_mean - x$truth), rmse = sqrt(mean((x$post_mean - x$truth)^2)),
             mean_width = mean(x$upper - x$lower))
}))))
write.csv(blk, file.path(out, "blocks.csv"), row.names = FALSE)

# 8. The trade-off species by species: B0 and collection-intercept errors by true occupancy band,
# in the eDNA arms. Collection intercept element 2s - 1 belongs to species s.
tb <- merge(merge(
  setNames(b0[, c("scenario", "replicate", "element", "post_mean", "truth")],
           c("scenario", "replicate", "species", "b0_est", "b0_true")),
  transform(setNames(bt0[, c("scenario", "replicate", "element", "post_mean", "truth")],
                     c("scenario", "replicate", "element", "bt0_est", "bt0_true")),
            species = (element + 1) / 2)[, c("scenario", "replicate", "species", "bt0_est", "bt0_true")]),
  o, by = c("scenario", "replicate", "species"))
tb$band <- cut(tb$occ_true, c(0, .1, .2, .3, .5, .7, 1), include.lowest = TRUE)
to <- do.call(rbind, lapply(split(tb, list(tb$scenario, tb$band), drop = TRUE), function(x) {
  if (nrow(x) < 10) return(NULL)
  b <- cl_mean(x$b0_est - x$b0_true, x$replicate); t <- cl_mean(x$bt0_est - x$bt0_true, x$replicate)
  data.frame(cell = x$scenario[1], band = as.character(x$band[1]), n = nrow(x),
             b0_bias = b[1], b0_bias_se = b[2], bt0_bias = t[1], bt0_bias_se = t[2],
             error_cor = cor(x$b0_est - x$b0_true, x$bt0_est - x$bt0_true))
}))
to <- ord(to[order(to$cell, to$band), ]); rownames(to) <- NULL
write.csv(to, file.path(out, "tradeoff-by-band.csv"), row.names = FALSE)

# 7. Success criteria C1 to C11 (README.md). C12 and C13 are in mixing.R.
g <- function(df, ...) { f <- list(...); for (n in names(f)) df <- df[df[[n]] == f[[n]], ]; df }
crit <- list()
add <- function(id, arm, quantity, default, community, pass, rule) crit[[length(crit) + 1]] <<-
  data.frame(criterion = id, arm = arm, quantity = quantity, default = default, community = community,
             pass = pass, rule = rule)
for (a in c("base_rare_cm", "occupancy_rare_cm")) {
  cd <- g(cm, cell = pairs[[a]]); cc <- g(cm, cell = a)
  add("C1", a, "community occupancy error", cd$occ_bias, cc$occ_bias, abs(cc$occ_bias) <= .05, "within 0.05")
  add("C1 halved", a, "community occupancy error", cd$occ_bias, cc$occ_bias, abs(cc$occ_bias) <= abs(cd$occ_bias) / 2, "at most half the default's")
  add("C2", a, "community B0 error", cd$b0_bias, cc$b0_bias, abs(cc$b0_bias) <= .5, "within 0.5")
  add("C2 halved", a, "community B0 error", cd$b0_bias, cc$b0_bias, abs(cc$b0_bias) <= abs(cd$b0_bias) / 2, "at most half the default's")
  for (bd in c("[0,0.05]", "(0.05,0.1]", "(0.1,0.2]")) {
    x <- g(ob, arm = a, band = bd)
    add("C3", a, paste("occupancy bias, band", bd), x$bias_default, x$bias_community, x$bias_community < .10, "under +0.10")
  }
  add("C4", a, "B0 coverage", cd$b0_coverage, cc$b0_coverage, cc$b0_coverage >= .85, "at least 0.85")
  x <- g(dt, arm = a, block = "beta_theta_intercept")
  add("C5", a, "collection intercept bias", x$bias_default, x$bias_community, abs(x$bias_community) <= .15, "within 0.15")
  add("C5 halved", a, "collection intercept bias", x$bias_default, x$bias_community, abs(x$bias_community) <= abs(x$bias_default) / 2, "at most half the default's")
  add("C6", a, "mu_b0 coverage of -2.75", NA, cc$mu_b0_coverage, cc$mu_b0_coverage >= .85, "at least 0.85")
}
cd <- g(cm, cell = "binary_rare"); cc <- g(cm, cell = "binary_rare_cm")
add("C7", "binary_rare_cm", "B0 coverage", cd$b0_coverage, cc$b0_coverage, cc$b0_coverage >= .85, "at least 0.85")
add("C7", "binary_rare_cm", "community B0 error", cd$b0_bias, cc$b0_bias, abs(cc$b0_bias) <= .3, "within 0.3")
add("C6", "binary_rare_cm", "mu_b0 coverage of -2.75", NA, cc$mu_b0_coverage, cc$mu_b0_coverage >= .85, "at least 0.85")
for (a in c("base_cm", "occupancy_cm", "binary_cm")) {
  cd <- g(cm, cell = pairs[[a]]); cc <- g(cm, cell = a)
  add("C8", a, "B0 coverage", cd$b0_coverage, cc$b0_coverage, cc$b0_coverage >= cd$b0_coverage - .02, "no more than 0.02 lower")
  if (!is.na(cd$occ_bias)) {
    add("C9", a, "community occupancy error", cd$occ_bias, cc$occ_bias, abs(cc$occ_bias) <= abs(cd$occ_bias) + .02, "no more than 0.02 larger in absolute value")
    x <- g(ob, arm = a); x <- x[x$n >= 10, ]
    worst <- which.max(abs(x$bias_community) - abs(x$bias_default))
    add("C10", a, paste("occupancy bias, worst band", x$band[worst]), x$bias_default[worst], x$bias_community[worst],
        all(abs(x$bias_community) <= abs(x$bias_default) + .03), "no band more than 0.03 larger in absolute value")
  }
  for (b in c("beta_theta_intercept", "p", "theta0")) {
    x <- g(dt, arm = a, block = b); if (!nrow(x)) next
    add("C11", a, paste(b, "bias"), x$bias_default, x$bias_community, abs(x$bias_community - x$bias_default) <= .05, "changed by at most 0.05")
    add("C11", a, paste(b, "coverage"), x$coverage_default, x$coverage_community, x$coverage_community >= x$coverage_default - .02, "no more than 0.02 lower")
  }
  bd <- g(blk, cell = pairs[[a]], block = "B"); bc <- g(blk, cell = a, block = "B")
  add("C11", a, "B coverage", bd$coverage, bc$coverage, bc$coverage >= bd$coverage - .02, "no more than 0.02 lower")
}
crit <- do.call(rbind, crit)
write.csv(crit, file.path(out, "criteria.csv"), row.names = FALSE)

options(width = 180)
cat("\n== pairing\n"); print(pk, row.names = FALSE)
cat("\n== community\n"); print(cm, digits = 3, row.names = FALSE)
cat("\n== occupancy by true band\n"); print(ob, digits = 3, row.names = FALSE)
cat("\n== B0 by true value\n"); print(bb, digits = 3, row.names = FALSE)
cat("\n== detection side and B0, paired\n"); print(dt, digits = 3, row.names = FALSE)
cat("\n== blocks\n"); print(blk, digits = 3, row.names = FALSE)
cat("\n== trade-off by band\n"); print(to, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
