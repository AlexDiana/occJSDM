# Compare the six M = 3 rare-dominated arms with each other (the three priors fit identical datasets)
# and with their M = 2 counterparts (paired on B0, eta and the true occupancy states).
# Run from the repository root: Rscript dev/simstudy/rare-m3-20261016/analyse.R
# Reads this run's results file (or its checkpoint while the run is in flight) and the results of the
# rare-dominated run (default prior), the community-mean run (B0 only) and the both-means run.
# Writes pairing-check.csv, community.csv, occupancy-by-band.csv, b0-by-truth.csv, detection-paired.csv,
# blocks.csv, tradeoff-by-band.csv, m-effect.csv and criteria.csv (C1 to C6, C14, C17 and C18 of README.md).
here <- "dev/simstudy/rare-m3-20261016"
out <- Sys.getenv("OUT_DIR", here)
rare_file <- "dev/simstudy/rare-dominated-20261015/results/simstudy-20261016-003851.rds"
cm_file <- "dev/simstudy/b0-community-mean-20261016/results/simstudy-20261016-110650.rds"
both_file <- "dev/simstudy/theta-community-mean-20261016/results/simstudy-20261016-132957.rds"
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
new <- readRDS(if (length(res)) res[1] else
  list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1])
new_rows <- if (!is.null(new$rows)) new$rows else new$params
old <- lapply(c(rare_file, cm_file, both_file), readRDS)

stems <- c("base_rare", "occupancy_rare")
priors <- c(default = "", b0_only = "_cm", both = "_cm_both")
m3 <- setNames(paste0(rep(stems, each = 3), priors, "_m3"), paste0(rep(stems, each = 3), priors, "_m3"))
m2_of <- setNames(sub("_m3$", "", names(m3)), names(m3))
prior_of <- setNames(rep(names(priors), 2), names(m3))
stem_of <- setNames(rep(stems, each = 3), names(m3))
cells <- c(rbind(unname(m2_of), names(m3)))
pairs <- rbind(
  # Within M = 3, identical datasets.
  data.frame(arm = paste0(stems, "_cm_m3"), ref = paste0(stems, "_m3"), kind = "default_m3"),
  data.frame(arm = paste0(stems, "_cm_both_m3"), ref = paste0(stems, "_m3"), kind = "default_m3"),
  data.frame(arm = paste0(stems, "_cm_both_m3"), ref = paste0(stems, "_cm_m3"), kind = "b0only_m3"),
  # Across M, the same prior; paired on B0, eta and the true occupancy states.
  data.frame(arm = names(m3), ref = unname(m2_of), kind = "m2"))

cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
occ_cols <- c("scenario", "replicate", "species", "occ_true", "occ_est")
d <- do.call(rbind, c(lapply(old, function(x) x$rows[x$rows$scenario %in% cells, cols]),
                      list(new_rows[new_rows$scenario %in% cells & new_rows$block != "ERROR", cols])))
o <- do.call(rbind, c(lapply(old, function(x) x$occupancy[x$occupancy$scenario %in% cells, occ_cols]),
                      list(new$occupancy[new$occupancy$scenario %in% cells, occ_cols])))
# Only replicates present in every cell of a stem, so a partial run still compares like with like.
keep <- function(x) do.call(rbind, lapply(stems, function(st) {
  f <- cells[grepl(paste0("^", st), cells)]
  reps <- Reduce(intersect, lapply(f, function(k) unique(d$replicate[d$scenario == k & d$block == "B0"])))
  x[x$scenario %in% f & x$replicate %in% reps, ]
}))
d <- keep(d); o <- keep(o)
for (st in stems) cat(sprintf("%s: %d replicates in every cell\n", st,
                              length(unique(d$replicate[d$scenario == paste0(st, "_m3")]))))

cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  c(m, if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA)
}
per_rep <- function(x, f) sapply(split(x, x$replicate), f)
msd <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
ord <- function(x, col = "cell") x[order(match(x[[col]], cells)), ]

# 1. Pairing: within M = 3 every truth agrees; across M, B0 and the true occupancy agree.
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q")
pk <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  rbind(
    do.call(rbind, lapply(blocks, function(b) {
      m <- merge(d[d$scenario == r & d$block == b, ], d[d$scenario == a & d$block == b, ], by = c("replicate", "element"))
      if (!nrow(m)) return(NULL)
      data.frame(arm = a, ref = r, kind = pairs$kind[i], quantity = b, n = nrow(m),
                 identical = all(m$truth.x == m$truth.y), max_abs_diff = max(abs(m$truth.x - m$truth.y)))
    })),
    {
      m <- merge(o[o$scenario == r, ], o[o$scenario == a, ], by = c("replicate", "species"))
      data.frame(arm = a, ref = r, kind = pairs$kind[i], quantity = "occ_true", n = nrow(m),
                 identical = all(m$occ_true.x == m$occ_true.y), max_abs_diff = max(abs(m$occ_true.x - m$occ_true.y)))
    })
}))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

# 2. The community as a whole, per replicate first. The collection intercept is beta_theta's odd elements.
b0 <- d[d$block == "B0", ]
bt0 <- d[d$block == "beta_theta" & d$element %% 2 == 1, ]
mub <- d[d$block == "mu_b0", ]
mut <- d[d$block == "mu_theta", ]
cm <- do.call(rbind, lapply(cells, function(k) {
  x <- b0[b0$scenario == k, ]; y <- o[o$scenario == k, ]; z <- bt0[bt0$scenario == k, ]
  u <- mub[mub$scenario == k, ]; v <- mut[mut$scenario == k, ]
  b0e <- per_rep(x, function(w) mean(w$post_mean - w$truth))
  ze <- per_rep(z, function(w) mean(w$post_mean - w$truth))
  oe <- per_rep(y, function(w) mean(w$occ_est) - mean(w$occ_true))
  r <- data.frame(cell = k, M = if (grepl("_m3$", k)) 3L else 2L, R = length(unique(x$replicate)),
    b0_true = mean(x$truth), b0_est = mean(x$post_mean), b0_bias = mean(b0e), b0_bias_se = msd(b0e)[2],
    b0_coverage = mean(x$covered), b0_slope = unname(coef(lm(post_mean ~ truth, x))[2]),
    bt0_true = mean(z$truth), bt0_est = mean(z$post_mean), bt0_bias = mean(ze), bt0_bias_se = msd(ze)[2],
    bt0_coverage = mean(z$covered),
    occ_true = mean(y$occ_true), occ_est = mean(y$occ_est), occ_bias = mean(oe), occ_bias_se = msd(oe)[2],
    occ_slope = unname(coef(lm(occ_est ~ occ_true, y))[2]),
    share10_true = mean(y$occ_true <= .1), share10_est = mean(y$occ_est <= .1),
    share20_true = mean(y$occ_true <= .2), share20_est = mean(y$occ_est <= .2),
    mu_b0_est = NA, mu_b0_coverage = NA, mu_b0_width = NA,
    mu_theta_truth = NA, mu_theta_est = NA, mu_theta_bias = NA, mu_theta_coverage = NA, mu_theta_width = NA)
  if (nrow(u)) r[c("mu_b0_est", "mu_b0_coverage", "mu_b0_width")] <-
    list(mean(u$post_mean), mean(u$covered), mean(u$upper - u$lower))
  if (nrow(v)) r[c("mu_theta_truth", "mu_theta_est", "mu_theta_bias", "mu_theta_coverage", "mu_theta_width")] <-
    list(mean(v$truth), mean(v$post_mean), mean(v$post_mean - v$truth), mean(v$covered), mean(v$upper - v$lower))
  r
}))
write.csv(cm, file.path(out, "community.csv"), row.names = FALSE)

# 3. Species occupancy by true band, each arm against its reference on the same species.
breaks <- c(0, .05, .1, .2, .3, .5, .7, .9, 1)
ob <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  m <- merge(o[o$scenario == r, ], o[o$scenario == a, ], by = c("replicate", "species"))
  m$band <- cut(m$occ_true.x, breaks, include.lowest = TRUE)
  do.call(rbind, lapply(split(m, m$band, drop = TRUE), function(x) {
    br <- cl_mean(x$occ_est.x - x$occ_true.x, x$replicate); ba <- cl_mean(x$occ_est.y - x$occ_true.y, x$replicate)
    df <- cl_mean(x$occ_est.y - x$occ_est.x, x$replicate)
    data.frame(arm = a, ref = r, kind = pairs$kind[i], band = as.character(x$band[1]), n = nrow(x),
               true_mean = mean(x$occ_true.x), bias_ref = br[1], bias_ref_se = br[2], bias_arm = ba[1],
               bias_arm_se = ba[2], paired_diff = df[1], paired_diff_se = df[2])
  }))
}))
rownames(ob) <- NULL
write.csv(ob, file.path(out, "occupancy-by-band.csv"), row.names = FALSE)

# 4. B0 by its true value.
b0$band <- cut(b0$truth, c(-Inf, -5, -4, -3, -2, -1, 0, 1, 2, Inf))
bb <- do.call(rbind, lapply(split(b0, list(b0$scenario, b0$band), drop = TRUE), function(x) {
  b <- cl_mean(x$post_mean - x$truth, x$replicate); cv <- cl_mean(x$covered, x$replicate)
  data.frame(cell = x$scenario[1], b0_band = as.character(x$band[1]), n = nrow(x), true_mean = mean(x$truth),
             bias = b[1], bias_se = b[2], coverage = cv[1], coverage_se = cv[2], mean_width = mean(x$upper - x$lower))
}))
bb <- ord(bb[order(bb$cell, bb$true_mean), ]); rownames(bb) <- NULL
write.csv(bb, file.path(out, "b0-by-truth.csv"), row.names = FALSE)

# 5. The detection side and B0, each arm against its reference, element by element. Across M the
# collection truths differ slightly (the slopes and the covariate's centring), so each side is scored
# against its own truth.
det <- d[d$block %in% c("B0", "beta_theta", "theta0", "p", "q"), ]
det$block[det$block == "beta_theta"] <- ifelse(det$element[det$block == "beta_theta"] %% 2 == 1,
                                               "beta_theta_intercept", "beta_theta_slope")
dt <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) do.call(rbind, lapply(unique(det$block), function(b) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  m <- merge(det[det$scenario == r & det$block == b, ], det[det$scenario == a & det$block == b, ], by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  bx <- cl_mean(m$post_mean.x - m$truth.x, m$replicate); by <- cl_mean(m$post_mean.y - m$truth.y, m$replicate)
  dd <- cl_mean((m$post_mean.y - m$truth.y) - (m$post_mean.x - m$truth.x), m$replicate)
  data.frame(arm = a, ref = r, kind = pairs$kind[i], block = b, n = nrow(m), bias_ref = bx[1], bias_arm = by[1],
             bias_arm_se = by[2], error_diff = dd[1], error_diff_se = dd[2],
             coverage_ref = mean(m$covered.x), coverage_arm = mean(m$covered.y),
             width_ref = mean(m$upper.x - m$lower.x), width_arm = mean(m$upper.y - m$lower.y))
}))))
write.csv(dt, file.path(out, "detection-paired.csv"), row.names = FALSE)

# 6. Every scored block.
blk <- do.call(rbind, lapply(cells, function(k) do.call(rbind, lapply(
  c("B0", "mu_b0", "mu_theta", "B", "G", "beta_theta", "theta0", "p", "q"), function(b) {
  x <- d[d$scenario == k & d$block == b, ]
  if (!nrow(x)) return(NULL)
  cv <- msd(per_rep(x, function(v) mean(v$covered)))
  data.frame(cell = k, block = b, n = nrow(x), coverage = cv[1], coverage_se = cv[2],
             bias = mean(x$post_mean - x$truth), rmse = sqrt(mean((x$post_mean - x$truth)^2)),
             mean_width = mean(x$upper - x$lower))
}))))
write.csv(blk, file.path(out, "blocks.csv"), row.names = FALSE)

# 7. The trade-off species by species: B0 and collection-intercept errors by true occupancy band.
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
             b0_bias = b[1], b0_bias_se = b[2], bt0_bias = t[1], bt0_bias_se = t[2])
}))
to <- ord(to[order(to$cell, to$band), ]); rownames(to) <- NULL
write.csv(to, file.path(out, "tradeoff-by-band.csv"), row.names = FALSE)

# 8. The third field sample, prior by prior (C17, C18): community errors per replicate, M = 3 minus M = 2.
rep_err <- function(k) {
  y <- o[o$scenario == k, ]; z <- bt0[bt0$scenario == k, ]; x <- b0[b0$scenario == k, ]
  data.frame(replicate = sort(unique(x$replicate)),
             occ = per_rep(y, function(w) mean(w$occ_est) - mean(w$occ_true)),
             bt0 = per_rep(z, function(w) mean(w$post_mean - w$truth)),
             b0 = per_rep(x, function(w) mean(w$post_mean - w$truth)),
             b0_cov = per_rep(x, function(w) mean(w$covered)))
}
me <- do.call(rbind, lapply(names(m3), function(a) {
  m <- merge(rep_err(m2_of[[a]]), rep_err(a), by = "replicate", suffixes = c("_m2", "_m3"))
  do.call(rbind, lapply(c("occ", "bt0", "b0", "b0_cov"), function(q) {
    dd <- m[[paste0(q, "_m3")]] - m[[paste0(q, "_m2")]]
    data.frame(arm = a, prior = prior_of[[a]], quantity = q, m2 = mean(m[[paste0(q, "_m2")]]),
               m3 = mean(m[[paste0(q, "_m3")]]), paired_diff = mean(dd), paired_diff_se = sd(dd) / sqrt(nrow(m)))
  }))
}))
write.csv(me, file.path(out, "m-effect.csv"), row.names = FALSE)

# 9. Success criteria (README.md). C12 and C13 are in mixing.R.
g <- function(df, ...) { f <- list(...); for (n in names(f)) df <- df[df[[n]] == f[[n]], ]; df }
crit <- list()
# A band with no species (possible in a partial run) gives NA rather than an error.
na1 <- function(x) if (length(x)) x else NA
add <- function(id, arm, quantity, m2, value, pass, rule) crit[[length(crit) + 1]] <<-
  data.frame(criterion = id, arm = arm, quantity = quantity, m2 = na1(m2), m3 = na1(value), pass = na1(pass),
             rule = rule)
for (a in c(paste0(stems, "_cm_both_m3"), paste0(stems, "_cm_m3"))) {
  c2 <- g(cm, cell = m2_of[[a]]); c3 <- g(cm, cell = a)
  add("C1", a, "community occupancy error", c2$occ_bias, c3$occ_bias, abs(c3$occ_bias) <= .05, "within 0.05")
  add("C2", a, "community B0 error", c2$b0_bias, c3$b0_bias, abs(c3$b0_bias) <= .5, "within 0.5")
  for (bd in c("[0,0.05]", "(0.05,0.1]", "(0.1,0.2]")) {
    x <- g(ob, arm = a, kind = "m2", band = bd)
    add("C3", a, paste("occupancy bias, band", bd), x$bias_ref, x$bias_arm, x$bias_arm < .10, "under +0.10")
  }
  add("C4", a, "B0 coverage", c2$b0_coverage, c3$b0_coverage, c3$b0_coverage >= .85, "at least 0.85")
  add("C5", a, "collection intercept bias", c2$bt0_bias, c3$bt0_bias, abs(c3$bt0_bias) <= .15, "within 0.15")
  add("C6", a, "mu_b0 coverage of -2.75", c2$mu_b0_coverage, c3$mu_b0_coverage, c3$mu_b0_coverage >= .85, "at least 0.85")
  if (grepl("both", a)) {
    add("C14", a, "mu_theta coverage of the realised mean", c2$mu_theta_coverage, c3$mu_theta_coverage,
        c3$mu_theta_coverage >= .85, "at least 0.85")
    add("C17", a, "collection intercept bias, M = 3 against M = 2", c2$bt0_bias, c3$bt0_bias,
        abs(c3$bt0_bias) <= abs(c2$bt0_bias) / 2, "at most half the M = 2 value in absolute value")
  }
}
for (a in names(m3)) {
  x <- g(me, arm = a, quantity = "occ")
  add("C18", a, "community occupancy error", x$m2, x$m3, NA, "descriptive")
  x <- g(ob, arm = a, kind = "m2", band = "[0,0.05]")
  add("C18", a, "occupancy bias, band [0,0.05]", x$bias_ref, x$bias_arm, NA, "descriptive")
  x <- g(me, arm = a, quantity = "b0_cov")
  add("C18", a, "B0 coverage", x$m2, x$m3, NA, "descriptive")
}
crit <- do.call(rbind, crit)
write.csv(crit, file.path(out, "criteria.csv"), row.names = FALSE)

options(width = 200)
cat("\n== pairing\n"); print(pk, digits = 3, row.names = FALSE)
cat("\n== community\n"); print(cm, digits = 3, row.names = FALSE)
cat("\n== occupancy by true band\n"); print(ob, digits = 3, row.names = FALSE)
cat("\n== detection side and B0, paired\n"); print(dt, digits = 3, row.names = FALSE)
cat("\n== trade-off by band\n"); print(to, digits = 3, row.names = FALSE)
cat("\n== M = 3 against M = 2\n"); print(me, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
