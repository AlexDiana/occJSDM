# Compare the six both-means arms with the default-prior and B0-only arms fitted to identical datasets,
# and the two _wide arms with their sd-2 twins.
# Run from the repository root: Rscript dev/simstudy/theta-community-mean-20261016/analyse.R
# Reads this run's results file (or its checkpoint while the run is in flight), the community-mean run's
# results, the rare-dominated run's results and replicates 1 to 100 of the R = 300 paper study.
# Writes pairing-check.csv, community.csv, occupancy-by-band.csv, b0-by-truth.csv, detection-paired.csv,
# blocks.csv, tradeoff-by-band.csv, wide-paired.csv and criteria.csv (C1 to C11, C14 and C15 of README.md).
here <- "dev/simstudy/theta-community-mean-20261016"
out <- Sys.getenv("OUT_DIR", here)
paper_file <- "dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds"
rare_file <- "dev/simstudy/rare-dominated-20261015/results/simstudy-20261016-003851.rds"
cm_file <- "dev/simstudy/b0-community-mean-20261016/results/simstudy-20261016-110650.rds"
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
new <- readRDS(if (length(res)) res[1] else
  list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1])
new_rows <- if (!is.null(new$rows)) new$rows else new$params
paper <- readRDS(paper_file); rare <- readRDS(rare_file); cmr <- readRDS(cm_file)

# Each arm with its default-prior and B0-only counterparts; the wide arms' reference is the sd-2 arm.
both <- c("base_cm_both", "occupancy_cm_both", "base_rare_cm_both", "occupancy_rare_cm_both")
default_of <- c(base_cm_both = "base", occupancy_cm_both = "occupancy",
                base_rare_cm_both = "base_rare", occupancy_rare_cm_both = "occupancy_rare")
b0only_of <- setNames(paste0(default_of, "_cm"), both)
wide_of <- c(base_rare_cm_both_wide = "base_rare_cm_both",
             occupancy_rare_cm_both_wide = "occupancy_rare_cm_both")
pairs <- rbind(
  data.frame(arm = both, ref = unname(default_of[both]), kind = "default"),
  data.frame(arm = both, ref = unname(b0only_of[both]), kind = "b0_only"),
  data.frame(arm = names(wide_of), ref = unname(wide_of), kind = "sd2"))
cells <- c("base", "base_cm", "base_cm_both", "occupancy", "occupancy_cm", "occupancy_cm_both",
           "base_rare", "base_rare_cm", "base_rare_cm_both", "base_rare_cm_both_wide",
           "occupancy_rare", "occupancy_rare_cm", "occupancy_rare_cm_both", "occupancy_rare_cm_both_wide")

cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
occ_cols <- c("scenario", "replicate", "species", "occ_true", "occ_est")
d <- rbind(paper$rows[paper$rows$scenario %in% cells & paper$rows$replicate <= 100, cols],
           rare$rows[rare$rows$scenario %in% cells, cols],
           cmr$rows[cmr$rows$scenario %in% cells, cols],
           new_rows[new_rows$scenario %in% cells & new_rows$block != "ERROR", cols])
o <- rbind(paper$occupancy[paper$occupancy$scenario %in% cells & paper$occupancy$replicate <= 100, occ_cols],
           rare$occupancy[rare$occupancy$scenario %in% cells, occ_cols],
           cmr$occupancy[cmr$occupancy$scenario %in% cells, occ_cols],
           new$occupancy[new$occupancy$scenario %in% cells, occ_cols])
# Only replicates present in every cell of a family, so a partial run still compares like with like.
family <- list(base = c("base", "base_cm", "base_cm_both"),
               occupancy = c("occupancy", "occupancy_cm", "occupancy_cm_both"),
               base_rare = c("base_rare", "base_rare_cm", "base_rare_cm_both", "base_rare_cm_both_wide"),
               occupancy_rare = c("occupancy_rare", "occupancy_rare_cm", "occupancy_rare_cm_both",
                                  "occupancy_rare_cm_both_wide"))
keep <- function(x) do.call(rbind, lapply(family, function(f) {
  reps <- Reduce(intersect, lapply(f, function(k) unique(d$replicate[d$scenario == k & d$block == "B0"])))
  x[x$scenario %in% f & x$replicate %in% reps, ]
}))
d <- keep(d); o <- keep(o)
for (f in family) cat(sprintf("%s: %d replicates in every cell\n", paste(f, collapse = ", "),
                              length(unique(d$replicate[d$scenario == f[1]]))))

# Mean over species-fits with a standard error treating each replicate as one cluster.
cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  c(m, if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA)
}
per_rep <- function(x, f) sapply(split(x, x$replicate), f)
msd <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
ord <- function(x, col = "cell") x[order(match(x[[col]], cells)), ]

# 1. Pairing: every truth, and each species' true occupancy, identical within each pair.
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q")
pk <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  rbind(
    do.call(rbind, lapply(blocks, function(b) {
      m <- merge(d[d$scenario == r & d$block == b, ], d[d$scenario == a & d$block == b, ],
                 by = c("replicate", "element"))
      if (!nrow(m)) return(NULL)
      data.frame(arm = a, ref = r, quantity = b, n = nrow(m), identical = all(m$truth.x == m$truth.y))
    })),
    {
      m <- merge(o[o$scenario == r, ], o[o$scenario == a, ], by = c("replicate", "species"))
      data.frame(arm = a, ref = r, quantity = "occ_true", n = nrow(m), identical = all(m$occ_true.x == m$occ_true.y))
    })
}))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

# 2. The community as a whole, per replicate first. The collection intercept is beta_theta's
# odd elements (2 x S, flattened by column).
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
  r <- data.frame(cell = k, R = length(unique(x$replicate)),
    b0_true = mean(x$truth), b0_est = mean(x$post_mean), b0_bias = mean(b0e), b0_bias_se = msd(b0e)[2],
    b0_coverage = mean(x$covered), b0_slope = unname(coef(lm(post_mean ~ truth, x))[2]),
    bt0_true = mean(z$truth), bt0_est = mean(z$post_mean), bt0_bias = mean(ze), bt0_bias_se = msd(ze)[2],
    bt0_coverage = mean(z$covered), bt0_b0_cor = cor(ze, b0e[names(ze)]),
    occ_true = mean(y$occ_true), occ_est = mean(y$occ_est), occ_bias = mean(oe), occ_bias_se = msd(oe)[2],
    occ_slope = unname(coef(lm(occ_est ~ occ_true, y))[2]),
    share10_true = mean(y$occ_true <= .1), share10_est = mean(y$occ_est <= .1),
    share20_true = mean(y$occ_true <= .2), share20_est = mean(y$occ_est <= .2),
    mu_b0_truth = NA, mu_b0_est = NA, mu_b0_coverage = NA, mu_b0_width = NA,
    mu_theta_truth = NA, mu_theta_est = NA, mu_theta_bias = NA, mu_theta_bias_se = NA,
    mu_theta_coverage = NA, mu_theta_width = NA, mu_errors_cor = NA)
  if (nrow(u)) r[c("mu_b0_truth", "mu_b0_est", "mu_b0_coverage", "mu_b0_width")] <-
    list(mean(u$truth), mean(u$post_mean), mean(u$covered), mean(u$upper - u$lower))
  if (nrow(v)) {
    ve <- v$post_mean - v$truth
    r[c("mu_theta_truth", "mu_theta_est", "mu_theta_bias", "mu_theta_bias_se", "mu_theta_coverage",
        "mu_theta_width")] <- list(mean(v$truth), mean(v$post_mean), mean(ve), msd(ve)[2], mean(v$covered),
                                   mean(v$upper - v$lower))
    # Across replicates: do the two community means' errors move against each other?
    if (nrow(u)) {
      m <- merge(u, v, by = "replicate")
      r$mu_errors_cor <- cor(m$post_mean.x - m$truth.x, m$post_mean.y - m$truth.y)
    }
  }
  r
}))
write.csv(cm, file.path(out, "community.csv"), row.names = FALSE)

# 3. Species occupancy by true band, each arm against its reference on the same species.
breaks <- c(0, .05, .1, .2, .3, .5, .7, .9, 1)
ob <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  m <- merge(o[o$scenario == r, ], o[o$scenario == a, ], by = c("replicate", "species"))
  if (!nrow(m)) return(NULL)
  m$band <- cut(m$occ_true.x, breaks, include.lowest = TRUE)
  do.call(rbind, lapply(split(m, m$band, drop = TRUE), function(x) {
    br <- cl_mean(x$occ_est.x - x$occ_true.x, x$replicate); ba <- cl_mean(x$occ_est.y - x$occ_true.x, x$replicate)
    df <- cl_mean(x$occ_est.y - x$occ_est.x, x$replicate)
    data.frame(arm = a, ref = r, kind = pairs$kind[i], band = as.character(x$band[1]), n = nrow(x),
               true_mean = mean(x$occ_true.x), est_ref = mean(x$occ_est.x), est_arm = mean(x$occ_est.y),
               bias_ref = br[1], bias_ref_se = br[2], bias_arm = ba[1], bias_arm_se = ba[2],
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

# 5. The detection side and B0, each arm against its reference, element by element.
det <- d[d$block %in% c("B0", "beta_theta", "theta0", "p", "q"), ]
det$block[det$block == "beta_theta"] <- ifelse(det$element[det$block == "beta_theta"] %% 2 == 1,
                                               "beta_theta_intercept", "beta_theta_slope")
dt <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) do.call(rbind, lapply(unique(det$block), function(b) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  m <- merge(det[det$scenario == r & det$block == b, ], det[det$scenario == a & det$block == b, ],
             by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  bx <- cl_mean(m$post_mean.x - m$truth.x, m$replicate); by <- cl_mean(m$post_mean.y - m$truth.y, m$replicate)
  dd <- cl_mean(m$post_mean.y - m$post_mean.x, m$replicate); dc <- cl_mean(m$covered.y - m$covered.x, m$replicate)
  data.frame(arm = a, ref = r, kind = pairs$kind[i], block = b, n = nrow(m), bias_ref = bx[1],
             bias_arm = by[1], bias_arm_se = by[2], paired_diff = dd[1], paired_diff_se = dd[2],
             coverage_ref = mean(m$covered.x), coverage_arm = mean(m$covered.y),
             coverage_diff = dc[1], coverage_diff_se = dc[2],
             width_ref = mean(m$upper.x - m$lower.x), width_arm = mean(m$upper.y - m$lower.y))
}))))
write.csv(dt, file.path(out, "detection-paired.csv"), row.names = FALSE)

# 6. Every scored block.
all_rows <- rbind(paper$rows[, cols], rare$rows[, cols], cmr$rows[, cols], new_rows[new_rows$block != "ERROR", cols])
blk <- do.call(rbind, lapply(cells, function(k) do.call(rbind, lapply(
  c("B0", "mu_b0", "mu_theta", "B", "G", "sigma_b", "beta_theta", "theta0", "p", "q", "resid_cor"), function(b) {
  x <- if (b == "resid_cor")
    all_rows[all_rows$scenario == k & all_rows$block == b & all_rows$replicate %in% unique(d$replicate[d$scenario == k]), ]
  else d[d$scenario == k & d$block == b, ]
  if (!nrow(x)) return(NULL)
  cv <- msd(per_rep(x, function(v) mean(v$covered)))
  data.frame(cell = k, block = b, n = nrow(x), coverage = cv[1], coverage_se = cv[2],
             bias = mean(x$post_mean - x$truth), rmse = sqrt(mean((x$post_mean - x$truth)^2)),
             mean_width = mean(x$upper - x$lower))
}))))
write.csv(blk, file.path(out, "blocks.csv"), row.names = FALSE)

# 7. The trade-off species by species: B0 and collection-intercept errors by true occupancy band.
# Collection intercept element 2s - 1 belongs to species s.
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

# 8. The wide hyperpriors against sd 2, replicate by replicate (C15).
rep_summary <- function(k) {
  y <- o[o$scenario == k, ]; z <- bt0[bt0$scenario == k, ]; x <- b0[b0$scenario == k, ]
  u <- mub[mub$scenario == k, ]; v <- mut[mut$scenario == k, ]
  reps <- sort(unique(x$replicate))
  data.frame(replicate = reps,
    occ_err = per_rep(y, function(w) mean(w$occ_est) - mean(w$occ_true))[as.character(reps)],
    bt0_err = per_rep(z, function(w) mean(w$post_mean - w$truth))[as.character(reps)],
    b0_err = per_rep(x, function(w) mean(w$post_mean - w$truth))[as.character(reps)],
    mu_b0 = u$post_mean[match(reps, u$replicate)], mu_b0_width = (u$upper - u$lower)[match(reps, u$replicate)],
    mu_theta = v$post_mean[match(reps, v$replicate)], mu_theta_width = (v$upper - v$lower)[match(reps, v$replicate)])
}
wp <- do.call(rbind, lapply(names(wide_of), function(a) {
  w <- rep_summary(a); s <- rep_summary(wide_of[[a]])
  m <- merge(s, w, by = "replicate", suffixes = c("_sd2", "_sd5"))
  q <- c("occ_err", "bt0_err", "b0_err", "mu_b0", "mu_theta")
  rbind(
    do.call(rbind, lapply(q, function(k) {
      dd <- m[[paste0(k, "_sd5")]] - m[[paste0(k, "_sd2")]]
      data.frame(arm = a, quantity = k, mean_sd2 = mean(m[[paste0(k, "_sd2")]]), mean_sd5 = mean(m[[paste0(k, "_sd5")]]),
                 paired_diff = mean(dd), paired_diff_se = sd(dd) / sqrt(nrow(m)), mean_abs_diff = mean(abs(dd)))
    })),
    do.call(rbind, lapply(c("mu_b0_width", "mu_theta_width"), function(k) {
      ratio <- m[[paste0(k, "_sd5")]] / m[[paste0(k, "_sd2")]]
      data.frame(arm = a, quantity = k, mean_sd2 = mean(m[[paste0(k, "_sd2")]]), mean_sd5 = mean(m[[paste0(k, "_sd5")]]),
                 paired_diff = mean(ratio) - 1, paired_diff_se = sd(ratio) / sqrt(nrow(m)),
                 mean_abs_diff = mean(m[[paste0(k, "_sd5")]]) / mean(m[[paste0(k, "_sd2")]]) - 1)
    })))
}))
write.csv(wp, file.path(out, "wide-paired.csv"), row.names = FALSE)

# 9. Success criteria C1 to C11, C14 and C15 (README.md). C12, C13 and C16 are in mixing.R.
g <- function(df, ...) { f <- list(...); for (n in names(f)) df <- df[df[[n]] == f[[n]], ]; df }
crit <- list()
add <- function(id, arm, quantity, default, b0_only, value, pass, rule) crit[[length(crit) + 1]] <<-
  data.frame(criterion = id, arm = arm, quantity = quantity, default = default, b0_only = b0_only,
             both = value, pass = pass, rule = rule)
for (a in c("base_rare_cm_both", "occupancy_rare_cm_both")) {
  cd <- g(cm, cell = default_of[[a]]); c0 <- g(cm, cell = b0only_of[[a]]); cc <- g(cm, cell = a)
  add("C1", a, "community occupancy error", cd$occ_bias, c0$occ_bias, cc$occ_bias, abs(cc$occ_bias) <= .05, "within 0.05")
  add("C2", a, "community B0 error", cd$b0_bias, c0$b0_bias, cc$b0_bias, abs(cc$b0_bias) <= .5, "within 0.5")
  for (bd in c("[0,0.05]", "(0.05,0.1]", "(0.1,0.2]")) {
    x <- g(ob, arm = a, kind = "default", band = bd); x0 <- g(ob, arm = a, kind = "b0_only", band = bd)
    add("C3", a, paste("occupancy bias, band", bd), x$bias_ref, x0$bias_ref, x$bias_arm, x$bias_arm < .10, "under +0.10")
  }
  add("C4", a, "B0 coverage", cd$b0_coverage, c0$b0_coverage, cc$b0_coverage, cc$b0_coverage >= .85, "at least 0.85")
  x <- g(dt, arm = a, kind = "default", block = "beta_theta_intercept")
  x0 <- g(dt, arm = a, kind = "b0_only", block = "beta_theta_intercept")
  add("C5", a, "collection intercept bias", x$bias_ref, x0$bias_ref, x$bias_arm, abs(x$bias_arm) <= .15, "within 0.15")
  add("C5 halved", a, "collection intercept bias", x$bias_ref, x0$bias_ref, x$bias_arm,
      abs(x$bias_arm) <= abs(x0$bias_ref) / 2, "at most half the B0-only arm's in absolute value")
  add("C6", a, "mu_b0 coverage of -2.75", NA, c0$mu_b0_coverage, cc$mu_b0_coverage, cc$mu_b0_coverage >= .85, "at least 0.85")
}
for (a in both) {
  cc <- g(cm, cell = a)
  add("C14", a, "mu_theta coverage of the realised mean", NA, NA, cc$mu_theta_coverage,
      cc$mu_theta_coverage >= .85, "at least 0.85")
}
for (a in c("base_cm_both", "occupancy_cm_both")) {
  cd <- g(cm, cell = default_of[[a]]); c0 <- g(cm, cell = b0only_of[[a]]); cc <- g(cm, cell = a)
  add("C8", a, "B0 coverage", cd$b0_coverage, c0$b0_coverage, cc$b0_coverage, cc$b0_coverage >= cd$b0_coverage - .02, "no more than 0.02 lower")
  add("C9", a, "community occupancy error", cd$occ_bias, c0$occ_bias, cc$occ_bias,
      abs(cc$occ_bias) <= abs(cd$occ_bias) + .02, "no more than 0.02 larger in absolute value")
  x <- g(ob, arm = a, kind = "default"); x <- x[x$n >= 10, ]
  worst <- which.max(abs(x$bias_arm) - abs(x$bias_ref))
  x0 <- g(ob, arm = a, kind = "b0_only", band = x$band[worst])
  add("C10", a, paste("occupancy bias, worst band", x$band[worst]), x$bias_ref[worst], x0$bias_ref, x$bias_arm[worst],
      all(abs(x$bias_arm) <= abs(x$bias_ref) + .03), "no band more than 0.03 larger in absolute value")
  for (b in c("beta_theta_intercept", "p", "theta0")) {
    x <- g(dt, arm = a, kind = "default", block = b); if (!nrow(x)) next
    x0 <- g(dt, arm = a, kind = "b0_only", block = b)
    add("C11", a, paste(b, "bias"), x$bias_ref, x0$bias_ref, x$bias_arm, abs(x$bias_arm - x$bias_ref) <= .05, "changed by at most 0.05")
    add("C11", a, paste(b, "coverage"), x$coverage_ref, x0$coverage_ref, x$coverage_arm,
        x$coverage_arm >= x$coverage_ref - .02, "no more than 0.02 lower")
  }
  bd <- g(blk, cell = default_of[[a]], block = "B"); b0b <- g(blk, cell = b0only_of[[a]], block = "B")
  bc <- g(blk, cell = a, block = "B")
  add("C11", a, "B coverage", bd$coverage, b0b$coverage, bc$coverage, bc$coverage >= bd$coverage - .02, "no more than 0.02 lower")
}
for (a in names(wide_of)) {
  w <- g(wp, arm = a)
  add("C15", a, "change in community occupancy error, sd 5 against sd 2", NA, NA,
      g(w, quantity = "occ_err")$paired_diff, abs(g(w, quantity = "occ_err")$paired_diff) < .03, "under 0.03")
  add("C15", a, "change in collection intercept error, sd 5 against sd 2", NA, NA,
      g(w, quantity = "bt0_err")$paired_diff, abs(g(w, quantity = "bt0_err")$paired_diff) < .15, "under 0.15")
  add("C15", a, "mu_theta interval widening, sd 5 against sd 2", NA, NA,
      g(w, quantity = "mu_theta_width")$mean_abs_diff, g(w, quantity = "mu_theta_width")$mean_abs_diff < .25,
      "under 25% (ratio of mean widths)")
}
crit <- do.call(rbind, crit)
write.csv(crit, file.path(out, "criteria.csv"), row.names = FALSE)

options(width = 200)
cat("\n== pairing: all identical?", all(pk$identical), "\n"); print(pk[!pk$identical, ], row.names = FALSE)
cat("\n== community\n"); print(cm, digits = 3, row.names = FALSE)
cat("\n== occupancy by true band\n"); print(ob, digits = 3, row.names = FALSE)
cat("\n== B0 by true value\n"); print(bb, digits = 3, row.names = FALSE)
cat("\n== detection side and B0, paired\n"); print(dt, digits = 3, row.names = FALSE)
cat("\n== blocks\n"); print(blk, digits = 3, row.names = FALSE)
cat("\n== trade-off by band\n"); print(to, digits = 3, row.names = FALSE)
cat("\n== wide against sd 2\n"); print(wp, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
