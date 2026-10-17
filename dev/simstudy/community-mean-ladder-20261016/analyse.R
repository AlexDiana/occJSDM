# The field-sample ladder with B0's community mean: every rung (M = 2, 3, 5, 10, 20) of each community (ordinary,
# rare-dominated) under each prior (default, community mean), scored as README.md sets out.
# Run from the repository root: Rscript dev/simstudy/community-mean-ladder-20261016/analyse.R
# Reads this run's results file (or its checkpoint while the run is in flight) and replicates 1 to 100 of the
# reused arms in six earlier runs. Writes pairing-check.csv, by-m.csv, occupancy-by-band.csv, squeeze.csv,
# squeeze-answer.csv, prior-paired.csv, b0-by-truth.csv and mu-b0.csv.
here <- "dev/simstudy/community-mean-ladder-20261016"
out <- Sys.getenv("OUT_DIR", here)
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
new_file <- if (length(res)) res[1] else
  list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1]
if (is.na(new_file)) message("no results from this run yet: scoring the reused arms only")

# Every rung, where it comes from and which truths it carries. `truths` is the seed label: rungs sharing one share
# B0, eta and the true occupancy states replicate by replicate. `base` is the default prior's M = 2 on the base
# truths, the counterpart of base_cm; the ordinary default ladder's own M = 2 rung is the ladder's M2.
cells <- data.frame(
  cell = c("M2", "base", "base_m3", "M5", "M10", "M20",
           "base_cm", "base_cm_m3", "M5_cm", "M10_cm", "M20_cm",
           "base_rare", "base_rare_m3", "base_rare_m5", "base_rare_m10", "base_rare_m20",
           "base_rare_cm", "base_rare_cm_m3", "base_rare_cm_m5", "base_rare_cm_m10", "base_rare_cm_m20"),
  community = rep(c("ordinary", "rare"), c(11, 10)),
  prior = rep(c("default", "community", "default", "community"), c(6, 5, 5, 5)),
  M = c(2, 2, 3, 5, 10, 20, 2, 3, 5, 10, 20, 2, 3, 5, 10, 20, 2, 3, 5, 10, 20),
  truths = c("mladder", "base", "base", rep("mladder", 3), "base", "base", rep("mladder", 3), rep("base", 10)),
  file = c("ladder", "paper", "om3", "ladder", "ladder", "ladder", "b0cm", "om3", "new", "new", "new",
           "rare", "rm3", "new", "new", "new", "b0cm", "rm3", "new", "new", "new"),
  stringsAsFactors = FALSE)
cells$rung <- cells$cell != "base"
files <- c(
  ladder = "dev/simstudy/field-sample-ladder/results/simstudy-20261016-031826.rds",
  paper = "dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds",
  b0cm = "dev/simstudy/b0-community-mean-20261016/results/simstudy-20261016-110650.rds",
  rare = "dev/simstudy/rare-dominated-20261015/results/simstudy-20261016-003851.rds",
  om3 = "dev/simstudy/ordinary-m3-20261016/results/simstudy-20261016-183451.rds",
  rm3 = "dev/simstudy/rare-m3-20261016/results/simstudy-20261016-154536.rds",
  new = new_file)
files <- files[!is.na(files)]
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
occ_cols <- c("scenario", "replicate", "species", "occ_true", "occ_est")
d <- list(); o <- list()
for (k in names(files)) {
  x <- readRDS(files[[k]]); rows <- if (!is.null(x$rows)) x$rows else x$params
  want <- cells$cell[cells$file == k]
  d[[k]] <- rows[rows$scenario %in% want & rows$replicate <= 100 & rows$block != "ERROR", cols]
  o[[k]] <- x$occupancy[x$occupancy$scenario %in% want & x$occupancy$replicate <= 100, occ_cols]
}
d <- do.call(rbind, d); o <- do.call(rbind, o); rownames(d) <- rownames(o) <- NULL
# The collection intercepts are beta_theta's odd elements, the slopes its even ones.
bt <- d$block == "beta_theta"
d$block[bt] <- ifelse(d$element[bt] %% 2 == 1, "beta_theta_intercept", "beta_theta_slope")
d$element[bt] <- (d$element[bt] + 1) %/% 2
reps_of <- function(k) sort(unique(d$replicate[d$scenario == k & d$block == "B0"]))
cells$R <- vapply(cells$cell, function(k) length(reps_of(k)), integer(1))
print(cells[, c("cell", "community", "prior", "M", "truths", "R")], row.names = FALSE)
cells <- cells[cells$R > 0, ]

# Mean over species-fits with a standard error treating each replicate as one cluster.
cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  c(m, if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA)
}
per_rep <- function(x, f) sapply(split(x, x$replicate), f)
msd <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
# Least-squares slope of estimated on true occupancy, with a jackknife standard error over replicates.
slope <- function(y) { v <- var(y$occ_true); if (v > 0) cov(y$occ_true, y$occ_est) / v else NA }
slope_se <- function(y) {
  r <- unique(y$replicate); s <- vapply(r, function(i) slope(y[y$replicate != i, ]), 0)
  sqrt((length(r) - 1) / length(r) * sum((s - mean(s))^2))
}

# 1. Pairing. Within an M the two priors fit identical datasets, so every truth agrees. Across M within a seed
# label, everything drawn before the field samples agrees: the JSDM layer, the true occupancy states, p, q and
# theta0. The collection coefficients' stored truths do not: the slopes are redrawn with the collection covariate,
# and the intercepts are scored at that covariate's mean, which moves with it (the raw intercepts agree; see
# tests/testthat/test-simstudy-cm-ladder.R).
prior_pairs <- data.frame(
  arm = c("base_cm", "base_cm_m3", "M5_cm", "M10_cm", "M20_cm",
          "base_rare_cm", "base_rare_cm_m3", "base_rare_cm_m5", "base_rare_cm_m10", "base_rare_cm_m20"),
  ref = c("base", "base_m3", "M5", "M10", "M20",
          "base_rare", "base_rare_m3", "base_rare_m5", "base_rare_m10", "base_rare_m20"))
m_pairs <- data.frame(
  arm = c("M5", "M10", "M20", "M5_cm", "M10_cm", "M20_cm", "base_m3", "base_cm_m3",
          "base_rare_m3", "base_rare_m5", "base_rare_m10", "base_rare_m20",
          "base_rare_cm_m3", "base_rare_cm_m5", "base_rare_cm_m10", "base_rare_cm_m20"),
  ref = c("M2", "M2", "M2", "M2", "M2", "M2", "base", "base",
          rep("base_rare", 4), rep("base_rare", 4)))
avail <- function(p) p[p$arm %in% cells$cell & p$ref %in% cells$cell, ]
prior_pairs <- avail(prior_pairs); m_pairs <- avail(m_pairs)
pairs <- rbind(cbind(prior_pairs, kind = "same M, other prior"), cbind(m_pairs, kind = "other M, same truths"))
blocks <- c("B0", "B", "G", "beta_theta_intercept", "beta_theta_slope", "theta0", "p", "q")
pk <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$arm[i]; r <- pairs$ref[i]; same_m <- pairs$kind[i] == "same M, other prior"
  rbind(
    do.call(rbind, lapply(blocks, function(b) {
      m <- merge(d[d$scenario == r & d$block == b, ], d[d$scenario == a & d$block == b, ],
                 by = c("replicate", "element"))
      if (!nrow(m)) return(NULL)
      data.frame(arm = a, ref = r, kind = pairs$kind[i], quantity = b, n = nrow(m),
                 identical = all(m$truth.x == m$truth.y), max_abs_diff = max(abs(m$truth.x - m$truth.y)),
                 expected_identical = same_m || !b %in% c("beta_theta_intercept", "beta_theta_slope"))
    })),
    {
      m <- merge(o[o$scenario == r, ], o[o$scenario == a, ], by = c("replicate", "species"))
      data.frame(arm = a, ref = r, kind = pairs$kind[i], quantity = "occ_true", n = nrow(m),
                 identical = all(m$occ_true.x == m$occ_true.y),
                 max_abs_diff = max(abs(m$occ_true.x - m$occ_true.y)), expected_identical = TRUE)
    })
}))
pk$as_expected <- pk$identical == pk$expected_identical
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

# 2. Every rung's summary.
b0 <- d[d$block == "B0", ]; bt0 <- d[d$block == "beta_theta_intercept", ]; pp <- d[d$block == "p", ]
mub <- d[d$block == "mu_b0", ]; bb_ <- d[d$block == "B", ]; th <- d[d$block == "theta0", ]
rep_err <- function(k) {
  y <- o[o$scenario == k, ]; x <- b0[b0$scenario == k, ]; z <- bt0[bt0$scenario == k, ]; w <- pp[pp$scenario == k, ]
  data.frame(replicate = sort(unique(x$replicate)),
             occ = per_rep(y, function(v) mean(v$occ_est) - mean(v$occ_true)),
             mae = per_rep(y, function(v) mean(abs(v$occ_est - v$occ_true))),
             b0_cov = per_rep(x, function(v) mean(v$covered)),
             b0_width = per_rep(x, function(v) mean(v$upper - v$lower)),
             b0 = per_rep(x, function(v) mean(v$post_mean - v$truth)),
             bt0 = per_rep(z, function(v) mean(v$post_mean - v$truth)),
             p = per_rep(w, function(v) mean(v$post_mean - v$truth)))
}
bm <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  k <- cells$cell[i]; y <- o[o$scenario == k, ]; x <- b0[b0$scenario == k, ]; u <- mub[mub$scenario == k, ]
  e <- rep_err(k)
  tb <- merge(setNames(x[, c("replicate", "element", "post_mean", "truth")], c("replicate", "sp", "b0e", "b0t")),
              setNames(bt0[bt0$scenario == k, c("replicate", "element", "post_mean", "truth")],
                       c("replicate", "sp", "bte", "btt")))
  r <- data.frame(cells[i, c("cell", "community", "prior", "M", "truths", "rung", "R")],
    occ_true = mean(y$occ_true), occ_est = mean(y$occ_est),
    occ_slope = slope(y), occ_slope_se = slope_se(y),
    occ_mae = mean(e$mae), occ_mae_se = msd(e$mae)[2],
    occ_comm_err = mean(e$occ), occ_comm_err_se = msd(e$occ)[2],
    share20_true = mean(y$occ_true < .2), share20_est = mean(y$occ_est < .2),
    b0_coverage = mean(e$b0_cov), b0_coverage_se = msd(e$b0_cov)[2], b0_width = mean(e$b0_width),
    b0_bias = mean(e$b0), b0_bias_se = msd(e$b0)[2],
    b0_slope = unname(coef(lm(post_mean ~ truth, x))[2]),
    bt0_bias = mean(e$bt0), bt0_bias_se = msd(e$bt0)[2], bt0_coverage = mean(bt0$covered[bt0$scenario == k]),
    p_bias = mean(e$p), p_bias_se = msd(e$p)[2], p_coverage = mean(pp$covered[pp$scenario == k]),
    tradeoff_cor = cor(tb$b0e - tb$b0t, tb$bte - tb$btt),
    B_coverage = mean(bb_$covered[bb_$scenario == k]), theta0_coverage = mean(th$covered[th$scenario == k]),
    mu_b0_est = NA, mu_b0_coverage = NA, mu_b0_width = NA)
  if (nrow(u)) r[c("mu_b0_est", "mu_b0_coverage", "mu_b0_width")] <-
    list(mean(u$post_mean), mean(u$covered), mean(u$upper - u$lower))
  r
}))
write.csv(bm, file.path(out, "by-m.csv"), row.names = FALSE)

# 3. Species occupancy by band of true occupancy.
breaks <- c(0, .05, .1, .2, .3, .5, .7, .9, 1)
o$band <- cut(o$occ_true, breaks, include.lowest = TRUE)
ob <- do.call(rbind, lapply(split(o, list(o$scenario, o$band), drop = TRUE), function(x) {
  b <- cl_mean(x$occ_est - x$occ_true, x$replicate)
  data.frame(cell = x$scenario[1], band = as.character(x$band[1]), n = nrow(x), true_mean = mean(x$occ_true),
             est_mean = mean(x$occ_est), bias = b[1], bias_se = b[2], mae = mean(abs(x$occ_est - x$occ_true)))
}))
ob <- merge(cells[, c("cell", "community", "prior", "M", "truths")], ob, by = "cell")
ob <- ob[order(match(ob$cell, cells$cell), match(ob$band, levels(o$band))), ]; rownames(ob) <- NULL
write.csv(ob, file.path(out, "occupancy-by-band.csv"), row.names = FALSE)

# 4. The rules of README.md: mostly removed (slope >= 0.80, every band with >= 10 species-fits within 0.10) and
# removed (slope >= 0.90, every such band within 0.05). Judged on point estimates; a failure within two standard
# errors of its threshold is "narrow".
sq <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  k <- cells$cell[i]; s <- bm[bm$cell == k, ]; x <- ob[ob$cell == k & ob$n >= 10, ]
  w <- which.max(abs(x$bias))
  judge <- function(slope_min, band_max) {
    gap_s <- slope_min - s$occ_slope; gap_b <- abs(x$bias) - band_max   # positive where failing
    pass <- gap_s <= 1e-12 && all(gap_b <= 1e-12)
    c(pass = pass, narrow = !pass && gap_s <= 2 * s$occ_slope_se && all(gap_b <= 2 * x$bias_se))
  }
  mr <- judge(.8, .1); rm <- judge(.9, .05)
  data.frame(cells[i, c("cell", "community", "prior", "M", "truths", "rung")],
             occ_slope = s$occ_slope, occ_slope_se = s$occ_slope_se, bands_judged = nrow(x),
             worst_band = x$band[w], worst_band_n = x$n[w], worst_band_bias = x$bias[w],
             worst_band_bias_se = x$bias_se[w],
             mostly_removed = mr[["pass"]], mostly_removed_narrow_fail = mr[["narrow"]],
             removed = rm[["pass"]], removed_narrow_fail = rm[["narrow"]])
}))
write.csv(sq, file.path(out, "squeeze.csv"), row.names = FALSE)
first_m <- function(x, col) { m <- x$M[x[[col]]]; if (length(m)) as.character(min(m)) else "more than 20" }
sa <- do.call(rbind, lapply(split(sq[sq$rung, ], list(sq$community[sq$rung], sq$prior[sq$rung])), function(x) {
  x <- x[order(x$M), ]
  data.frame(community = x$community[1], prior = x$prior[1], rungs = paste(x$M, collapse = ", "),
             mostly_removed_at = first_m(x, "mostly_removed"), removed_at = first_m(x, "removed"),
             # A rule met at one rung and failed at a higher one would make "the smallest rung" misleading.
             mostly_removed_monotone = all(diff(x$mostly_removed) >= 0),
             removed_monotone = all(diff(x$removed) >= 0))
}))
sa <- sa[order(sa$community, -xtfrm(sa$prior)), ]; rownames(sa) <- NULL
write.csv(sa, file.path(out, "squeeze-answer.csv"), row.names = FALSE)

# 5. The community mean against the default prior on identical datasets, paired by replicate.
qs <- c("occ", "mae", "b0_cov", "b0_width", "bt0", "p")
pr <- do.call(rbind, lapply(seq_len(nrow(prior_pairs)), function(i) {
  a <- prior_pairs$arm[i]; r <- prior_pairs$ref[i]
  x <- rep_err(a); y <- rep_err(r); m <- merge(y, x, by = "replicate", suffixes = c("_default", "_community"))
  ya <- o[o$scenario == a, ]; yr <- o[o$scenario == r, ]
  sd_ <- slope(ya) - slope(yr)
  jr <- vapply(m$replicate, function(j) slope(ya[ya$replicate != j, ]) - slope(yr[yr$replicate != j, ]), 0)
  rbind(
    do.call(rbind, lapply(qs, function(q) {
      dd <- m[[paste0(q, "_community")]] - m[[paste0(q, "_default")]]
      data.frame(arm = a, ref = r, community = cells$community[cells$cell == a], M = cells$M[cells$cell == a],
                 quantity = q, R = nrow(m), default = mean(m[[paste0(q, "_default")]]),
                 community_mean = mean(m[[paste0(q, "_community")]]), diff = mean(dd), diff_se = sd(dd) / sqrt(nrow(m)))
    })),
    data.frame(arm = a, ref = r, community = cells$community[cells$cell == a], M = cells$M[cells$cell == a],
               quantity = "occ_slope", R = nrow(m), default = slope(yr), community_mean = slope(ya), diff = sd_,
               diff_se = sqrt((length(jr) - 1) / length(jr) * sum((jr - mean(jr))^2))))
}))
write.csv(pr, file.path(out, "prior-paired.csv"), row.names = FALSE)

# 6. B0 by its true value, on bands placed around each community's centre.
b0 <- merge(b0, cells[, c("cell", "community")], by.x = "scenario", by.y = "cell")
b0$band <- ifelse(b0$community == "rare", as.character(cut(b0$truth, c(-Inf, -4.25, -3.25, -2.25, -1.25, Inf))),
                  as.character(cut(b0$truth, c(-Inf, -1.5, -.5, .5, 1.5, Inf))))
bb <- do.call(rbind, lapply(split(b0, list(b0$scenario, b0$band), drop = TRUE), function(x) {
  b <- cl_mean(x$post_mean - x$truth, x$replicate); cv <- cl_mean(x$covered, x$replicate)
  data.frame(cell = x$scenario[1], b0_band = x$band[1], n = nrow(x), true_mean = mean(x$truth),
             bias = b[1], bias_se = b[2], coverage = cv[1], coverage_se = cv[2], mean_width = mean(x$upper - x$lower))
}))
bb <- bb[order(match(bb$cell, cells$cell), bb$true_mean), ]; rownames(bb) <- NULL
write.csv(bb, file.path(out, "b0-by-truth.csv"), row.names = FALSE)

# 7. mu_b0 against the generating mean (0, or -2.75 in the rare community) and the realised mean of the true B0.
mu <- do.call(rbind, lapply(cells$cell[cells$prior == "community"], function(k) {
  u <- mub[mub$scenario == k, ]; x <- b0[b0$scenario == k, ]
  rm <- per_rep(x, function(w) mean(w$truth)); em <- per_rep(x, function(w) mean(w$post_mean))
  m <- merge(u, data.frame(replicate = as.integer(names(rm)), realised = rm, b0_est_mean = em), by = "replicate")
  data.frame(cell = k, community = cells$community[cells$cell == k], M = cells$M[cells$cell == k], R = nrow(m),
             generating_mean = m$truth[1], coverage_of_generating = mean(m$covered), estimate = mean(m$post_mean),
             realised_mean = mean(m$realised), error_vs_realised = mean(m$post_mean - m$realised),
             covers_realised = mean(m$lower <= m$realised & m$realised <= m$upper),
             mean_b0_estimate = mean(m$b0_est_mean), width = mean(m$upper - m$lower))
}))
write.csv(mu, file.path(out, "mu-b0.csv"), row.names = FALSE)

options(width = 220)
cat("\n== pairing: rows not as expected\n"); print(pk[!pk$as_expected, ], digits = 3, row.names = FALSE)
cat("pairing rows:", nrow(pk), " as expected:", sum(pk$as_expected), "\n")
cat("\n== by M\n")
print(bm[, c("cell", "community", "prior", "M", "R", "occ_slope", "occ_slope_se", "occ_mae", "occ_comm_err",
             "share20_true", "share20_est", "b0_coverage", "b0_coverage_se", "b0_width", "bt0_bias", "p_bias",
             "tradeoff_cor", "mu_b0_est", "mu_b0_coverage")], digits = 3, row.names = FALSE)
cat("\n== occupancy by true band\n")
print(ob[, c("cell", "band", "n", "true_mean", "est_mean", "bias", "bias_se")], digits = 3, row.names = FALSE)
cat("\n== squeeze rules\n"); print(sq, digits = 3, row.names = FALSE)
cat("\n== answer\n"); print(sa, row.names = FALSE)
cat("\n== community minus default, paired\n"); print(pr, digits = 3, row.names = FALSE)
cat("\n== B0 by truth\n"); print(bb, digits = 3, row.names = FALSE)
cat("\n== mu_b0\n"); print(mu, digits = 3, row.names = FALSE)
