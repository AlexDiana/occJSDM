# Compare the four M = 3 ordinary-community arms: B0's community mean against the default prior on identical
# datasets (C8 to C11 of README.md), and each arm against its M = 2 counterpart, paired on B0, eta and the true
# occupancy states (D1, D2).
# Run from the repository root: Rscript dev/simstudy/ordinary-m3-20261016/analyse.R
# Reads this run's results file (or its checkpoint while the run is in flight), replicates 1 to 100 of the
# R = 300 paper study (base, occupancy) and the community-mean run (base_cm, occupancy_cm).
# Writes pairing-check.csv, community.csv, occupancy-by-band.csv, b0-by-truth.csv, detection-paired.csv,
# blocks.csv, tradeoff-by-band.csv, m-effect.csv, mu-b0.csv and criteria.csv.
here <- "dev/simstudy/ordinary-m3-20261016"
out <- Sys.getenv("OUT_DIR", here)
paper_file <- "dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds"
cm_file <- "dev/simstudy/b0-community-mean-20261016/results/simstudy-20261016-110650.rds"
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
new <- readRDS(if (length(res)) res[1] else
  list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1])
new_rows <- if (!is.null(new$rows)) new$rows else new$params
old <- lapply(c(paper_file, cm_file), readRDS)

stems <- c("base", "occupancy")
# For each stem: M = 2 default, M = 2 community, M = 3 default, M = 3 community.
cells <- c(rbind(stems, paste0(stems, "_cm"), paste0(stems, "_m3"), paste0(stems, "_cm_m3")))
m2_of <- setNames(sub("_m3$", "", cells[grepl("_m3$", cells)]), cells[grepl("_m3$", cells)])
pairs <- rbind(
  # The test: within M = 3, identical datasets.
  data.frame(arm = paste0(stems, "_cm_m3"), ref = paste0(stems, "_m3"), kind = "default_m3"),
  # The same comparison at M = 2, recomputed on the replicates used here.
  data.frame(arm = paste0(stems, "_cm"), ref = stems, kind = "default_m2"),
  # Across M, the same prior; paired on B0, eta and the true occupancy states.
  data.frame(arm = names(m2_of), ref = unname(m2_of), kind = "m2"))

cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
occ_cols <- c("scenario", "replicate", "species", "occ_true", "occ_est")
d <- do.call(rbind, c(lapply(old, function(x) x$rows[x$rows$scenario %in% cells & x$rows$replicate <= 100, cols]),
                      list(new_rows[new_rows$scenario %in% cells & new_rows$block != "ERROR", cols])))
o <- do.call(rbind, c(lapply(old, function(x) x$occupancy[x$occupancy$scenario %in% cells &
                                                            x$occupancy$replicate <= 100, occ_cols]),
                      list(new$occupancy[new$occupancy$scenario %in% cells, occ_cols])))
# Only replicates present in every cell of a stem, so a partial run still compares like with like.
keep <- function(x) do.call(rbind, lapply(stems, function(st) {
  f <- cells[sub("_cm", "", sub("_m3$", "", cells)) == st]
  reps <- Reduce(intersect, lapply(f, function(k) unique(d$replicate[d$scenario == k & d$block == "B0"])))
  x[x$scenario %in% f & x$replicate %in% reps, ]
}))
d <- keep(d); o <- keep(o)
for (st in stems) cat(sprintf("%s: %d replicates in every cell\n", st,
                              length(unique(d$replicate[d$scenario == paste0(st, "_m3")]))))

# Mean over species-fits with a standard error treating each replicate as one cluster.
cl_mean <- function(x, rep) {
  s <- tapply(x, rep, sum); n <- tapply(x, rep, length)
  m <- sum(s) / sum(n); e <- s - m * n; R <- length(s)
  c(m, if (R > 1) sqrt(R / (R - 1) * sum(e^2)) / sum(n) else NA)
}
per_rep <- function(x, f) sapply(split(x, x$replicate), f)
msd <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
ord <- function(x, col = "cell") x[order(match(x[[col]], cells)), ]

# 1. Pairing: within M = 3 and within M = 2 every truth agrees; across M, B0 and the true occupancy agree.
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
cm <- do.call(rbind, lapply(cells, function(k) {
  x <- b0[b0$scenario == k, ]; y <- o[o$scenario == k, ]; z <- bt0[bt0$scenario == k, ]; u <- mub[mub$scenario == k, ]
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
    mu_b0_est = NA, mu_b0_coverage = NA, mu_b0_width = NA)
  if (nrow(u)) r[c("mu_b0_est", "mu_b0_coverage", "mu_b0_width")] <-
    list(mean(u$post_mean), mean(u$covered), mean(u$upper - u$lower))
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
    # The change in absolute bias, which C10 judges, with its own clustered standard error.
    ab <- cl_mean(abs(x$occ_est.y - x$occ_true.y) - abs(x$occ_est.x - x$occ_true.x), x$replicate)
    data.frame(arm = a, ref = r, kind = pairs$kind[i], band = as.character(x$band[1]), n = nrow(x),
               true_mean = mean(x$occ_true.x), bias_ref = br[1], bias_ref_se = br[2], bias_arm = ba[1],
               bias_arm_se = ba[2], paired_diff = df[1], paired_diff_se = df[2],
               abs_bias_change = abs(ba[1]) - abs(br[1]), species_abs_error_change = ab[1],
               species_abs_error_change_se = ab[2])
  }))
}))
rownames(ob) <- NULL
write.csv(ob, file.path(out, "occupancy-by-band.csv"), row.names = FALSE)

# 4. B0 by its true value.
b0$band <- cut(b0$truth, c(-Inf, -2, -1, 0, 1, 2, Inf))
bb <- do.call(rbind, lapply(split(b0, list(b0$scenario, b0$band), drop = TRUE), function(x) {
  b <- cl_mean(x$post_mean - x$truth, x$replicate); cv <- cl_mean(x$covered, x$replicate)
  data.frame(cell = x$scenario[1], b0_band = as.character(x$band[1]), n = nrow(x), true_mean = mean(x$truth),
             bias = b[1], bias_se = b[2], coverage = cv[1], coverage_se = cv[2], mean_width = mean(x$upper - x$lower))
}))
bb <- ord(bb[order(bb$cell, bb$true_mean), ]); rownames(bb) <- NULL
write.csv(bb, file.path(out, "b0-by-truth.csv"), row.names = FALSE)

# 5. The detection side, B0 and B, each arm against its reference, element by element. Across M the collection
# truths differ slightly (the slopes and the covariate's centring), so each side is scored against its own truth.
det <- d[d$block %in% c("B0", "B", "beta_theta", "theta0", "p", "q"), ]
det$block[det$block == "beta_theta"] <- ifelse(det$element[det$block == "beta_theta"] %% 2 == 1,
                                               "beta_theta_intercept", "beta_theta_slope")
dt <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) do.call(rbind, lapply(unique(det$block), function(b) {
  a <- pairs$arm[i]; r <- pairs$ref[i]
  m <- merge(det[det$scenario == r & det$block == b, ], det[det$scenario == a & det$block == b, ], by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  bx <- cl_mean(m$post_mean.x - m$truth.x, m$replicate); by <- cl_mean(m$post_mean.y - m$truth.y, m$replicate)
  dd <- cl_mean((m$post_mean.y - m$truth.y) - (m$post_mean.x - m$truth.x), m$replicate)
  dc <- cl_mean(m$covered.y - m$covered.x, m$replicate)
  data.frame(arm = a, ref = r, kind = pairs$kind[i], block = b, n = nrow(m), bias_ref = bx[1], bias_arm = by[1],
             bias_arm_se = by[2], error_diff = dd[1], error_diff_se = dd[2],
             coverage_ref = mean(m$covered.x), coverage_arm = mean(m$covered.y),
             coverage_diff = dc[1], coverage_diff_se = dc[2],
             width_ref = mean(m$upper.x - m$lower.x), width_arm = mean(m$upper.y - m$lower.y))
}))))
write.csv(dt, file.path(out, "detection-paired.csv"), row.names = FALSE)

# 6. Every scored block.
blk <- do.call(rbind, lapply(cells, function(k) do.call(rbind, lapply(
  c("B0", "mu_b0", "B", "G", "beta_theta", "theta0", "p", "q"), function(b) {
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

# 8. D1: the third field sample. Per replicate, each arm's community errors and the option's cost (community
# minus default), at M = 2 and M = 3; M = 3 minus M = 2 paired by replicate.
rep_err <- function(k) {
  y <- o[o$scenario == k, ]; z <- bt0[bt0$scenario == k, ]; x <- b0[b0$scenario == k, ]
  data.frame(replicate = sort(unique(x$replicate)),
             occ = per_rep(y, function(w) mean(w$occ_est) - mean(w$occ_true)),
             bt0 = per_rep(z, function(w) mean(w$post_mean - w$truth)),
             b0 = per_rep(x, function(w) mean(w$post_mean - w$truth)),
             b0_cov = per_rep(x, function(w) mean(w$covered)))
}
qs <- c("occ", "bt0", "b0", "b0_cov")
pair_rows <- function(label, prior, x2, x3) do.call(rbind, lapply(qs, function(q) {
  m <- merge(x2, x3, by = "replicate", suffixes = c("_m2", "_m3"))
  dd <- m[[paste0(q, "_m3")]] - m[[paste0(q, "_m2")]]
  data.frame(arm = label, prior = prior, quantity = q, m2 = mean(m[[paste0(q, "_m2")]]),
             m3 = mean(m[[paste0(q, "_m3")]]), paired_diff = mean(dd), paired_diff_se = sd(dd) / sqrt(nrow(m)))
}))
cost <- function(a, r) { x <- rep_err(a); y <- rep_err(r); stopifnot(identical(x$replicate, y$replicate))
  cbind(replicate = x$replicate, x[qs] - y[qs]) }
me <- do.call(rbind, c(
  lapply(names(m2_of), function(a) pair_rows(a, if (grepl("_cm", a)) "community" else "default",
                                             rep_err(m2_of[[a]]), rep_err(a))),
  lapply(stems, function(st) pair_rows(paste0(st, ": community minus default"), "cost",
                                       cost(paste0(st, "_cm"), st), cost(paste0(st, "_cm_m3"), paste0(st, "_m3"))))))
write.csv(me, file.path(out, "m-effect.csv"), row.names = FALSE)

# 9. D2: mu_b0 against the generating mean, 0, and the realised mean of the true B0.
mu <- do.call(rbind, lapply(cells[grepl("_cm", cells)], function(k) {
  u <- mub[mub$scenario == k, ]; x <- b0[b0$scenario == k, ]
  rm <- per_rep(x, function(w) mean(w$truth))
  m <- merge(u, data.frame(replicate = as.integer(names(rm)), realised = rm), by = "replicate")
  data.frame(cell = k, R = nrow(m), coverage_of_0 = mean(m$covered), estimate = mean(m$post_mean),
             realised_mean = mean(m$realised), error_vs_realised = mean(m$post_mean - m$realised),
             covers_realised = mean(m$lower <= m$realised & m$realised <= m$upper), width = mean(m$upper - m$lower))
}))
write.csv(mu, file.path(out, "mu-b0.csv"), row.names = FALSE)

# 10. Success criteria C8 to C11 (README.md); C12 and C13 are in mixing.R. "worsening" is how much worse the
# community arm is on the criterion's own scale, so a criterion passes when worsening <= threshold; "se" is the
# paired standard error of the difference, and a failure within two of them of the threshold is "narrow".
g <- function(df, ...) { f <- list(...); for (n in names(f)) df <- df[df[[n]] == f[[n]], ]; df }
crit <- list()
add <- function(id, arm, quantity, default, community, worsening, threshold, se, rule) {
  pass <- worsening <= threshold + 1e-12
  crit[[length(crit) + 1]] <<- data.frame(criterion = id, arm = arm, quantity = quantity, default = default,
    community = community, worsening = worsening, threshold = threshold, se = se, pass = pass,
    narrow = !pass & worsening - threshold <= 2 * se, rule = rule)
}
for (st in stems) for (a in c(paste0(st, "_cm_m3"), paste0(st, "_cm"))) {
  r <- if (grepl("_m3$", a)) paste0(st, "_m3") else st
  id <- function(x) if (grepl("_m3$", a)) x else paste(x, "at M = 2")
  cd <- g(cm, cell = r); cc <- g(cm, cell = a)
  ea <- rep_err(a); er <- rep_err(r)
  x <- g(dt, arm = a, ref = r, block = "B0")
  add(id("C8"), a, "B0 coverage", cd$b0_coverage, cc$b0_coverage, cd$b0_coverage - cc$b0_coverage, .02,
      x$coverage_diff_se, "no more than 0.02 lower")
  add(id("C9"), a, "community occupancy error", cd$occ_bias, cc$occ_bias, abs(cc$occ_bias) - abs(cd$occ_bias), .02,
      sd(ea$occ - er$occ) / sqrt(nrow(ea)), "no more than 0.02 larger in absolute value")
  x <- g(ob, arm = a, ref = r); x <- x[x$n >= 10, ]
  for (j in seq_len(nrow(x)))
    add(id("C10"), a, paste("occupancy bias, band", x$band[j]), x$bias_ref[j], x$bias_arm[j], x$abs_bias_change[j],
        .03, x$paired_diff_se[j], "no band more than 0.03 larger in absolute value")
  for (b in c("beta_theta_intercept", "p", "theta0")) {
    x <- g(dt, arm = a, ref = r, block = b); if (!nrow(x)) next
    add(id("C11"), a, paste(b, "bias"), x$bias_ref, x$bias_arm, abs(x$error_diff), .05, x$error_diff_se,
        "changed by at most 0.05")
    add(id("C11"), a, paste(b, "coverage"), x$coverage_ref, x$coverage_arm, -x$coverage_diff, .02,
        x$coverage_diff_se, "no more than 0.02 lower")
  }
  x <- g(dt, arm = a, ref = r, block = "B")
  add(id("C11"), a, "B coverage", x$coverage_ref, x$coverage_arm, -x$coverage_diff, .02, x$coverage_diff_se,
      "no more than 0.02 lower")
}
crit <- do.call(rbind, crit)
write.csv(crit, file.path(out, "criteria.csv"), row.names = FALSE)

options(width = 220)
cat("\n== pairing (anything not identical)\n"); print(pk[!pk$identical, ], digits = 3, row.names = FALSE)
cat("\n== community\n"); print(cm, digits = 3, row.names = FALSE)
cat("\n== occupancy by true band\n"); print(ob, digits = 3, row.names = FALSE)
cat("\n== detection side, B0 and B, paired\n"); print(dt, digits = 3, row.names = FALSE)
cat("\n== trade-off by band\n"); print(to, digits = 3, row.names = FALSE)
cat("\n== M = 3 against M = 2 (D1)\n"); print(me, digits = 3, row.names = FALSE)
cat("\n== mu_b0 (D2)\n"); print(mu, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
