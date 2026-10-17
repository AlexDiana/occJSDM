# Summarise the five sigma_b prior arms against no_spatial, no_spatial_long,
# no_spatial_gt4 and no_spatial_gt4_long, which the two earlier B runs fitted
# to the same datasets at revision c40f385 (package code).
# Run from the repository root: Rscript dev/simstudy/sigma-b-prior-20261017/analyse.R
# Reads this run's results file, or its checkpoint while the run is in flight,
# and the two earlier results files for the comparison arms.
# Writes pairing-check.csv, coverage-by-arm.csv, sigma-b-posterior.csv,
# paired-contrasts.csv, b-by-effect-size.csv, b-by-effect-size-contrasts.csv,
# chain-length.csv and criteria.csv.
here <- "dev/simstudy/sigma-b-prior-20261017"
out <- Sys.getenv("OUT_DIR", here)  # where the CSVs go; the run directory by default
cand_file <- "dev/simstudy/b-nonspatial-candidates-20261015/results/simstudy-20261016-004116.rds"
gt4_file <- "dev/simstudy/b-gt4-long-20261016/results/simstudy-20261016-074749.rds"
res_dir <- Sys.getenv("RESULTS_DIR", file.path(here, "results"))
res <- list.files(res_dir, pattern = "^simstudy-.*[.]rds$", full.names = TRUE)
run <- if (length(res)) readRDS(res[1]) else NULL
new <- if (length(res)) run$rows else
  readRDS(list.files(res_dir, pattern = "^checkpoint", full.names = TRUE)[1])$params
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
new_arms <- c("no_spatial_sbhc", "no_spatial_sbig2", "no_spatial_sbhc_long",
              "no_spatial_sbig2_long", "no_spatial_sbhc_gt4_long")
old_arms <- c("no_spatial", "no_spatial_long", "no_spatial_gt4", "no_spatial_gt4_long")
arms <- c(old_arms, new_arms)
cand <- readRDS(cand_file); gt4 <- readRDS(gt4_file)
d <- rbind(cand$rows[cand$rows$scenario %in% old_arms, cols],
           gt4$rows[gt4$rows$scenario == "no_spatial_gt4_long", cols],
           new[new$scenario %in% new_arms & new$block != "ERROR", cols])
d$width <- d$upper - d$lower; d$err <- d$post_mean - d$truth
by_arm <- vapply(arms, function(a) length(unique(d$replicate[d$scenario == a])), integer(1))
cat("replicates per arm:\n"); print(by_arm)
reps <- Reduce(intersect, lapply(arms, function(a) unique(d$replicate[d$scenario == a])))
d <- d[d$replicate %in% reps, ]
cat(sprintf("replicates present in all %d arms: %d\n", length(arms), length(reps)))
if (length(res)) {
  pv <- run$provenance
  cat(sprintf("failed replicates: %d; revision %s, dirty at start %s, at write %s, moved %s\n",
              run$failed_replicates, substr(pv$git_sha, 1, 7), pv$git_dirty,
              pv$git_dirty_at_write, pv$git_moved_during_run))
  cat("cell MCMC settings recorded:", deparse(pv$cell_MCMCparams), "\n")
  cat(sprintf("elapsed %.1f min\n", run$elapsed_min))
}
blocks <- c("B", "B0", "G", "sigma_b", "beta_theta", "theta0", "p", "q")
zq <- 2 * qnorm(0.975)
per_rep <- function(x, v) aggregate(x[[v]], list(replicate = x$replicate), mean)

# 1. Pairing, checked first: are the truths element for element the same in
# every arm, with the same rows present? And does no new arm share a
# posterior with the arm it is compared with (it should not)?
pk <- do.call(rbind, lapply(new_arms, function(a) do.call(rbind, lapply(old_arms, function(o)
  do.call(rbind, lapply(blocks, function(b) {
    x <- d[d$scenario == a & d$block == b, ]
    y <- d[d$scenario == o & d$block == b, ]
    if (!nrow(x) || !nrow(y)) return(NULL)
    m <- merge(x, y, by = c("replicate", "element"))
    data.frame(arm = a, against = o, block = b, n_new = nrow(x), n_old = nrow(y), n_matched = nrow(m),
               max_abs_truth_diff = if (nrow(m)) max(abs(m$truth.x - m$truth.y)) else NA,
               identical_post_means = if (nrow(m)) mean(m$post_mean.x == m$post_mean.y) else NA)
  }))))))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)
cat(sprintf("truths identical element for element in every block and arm: %s (%d rows compared); any shared posterior means: %s\n",
            all(pk$max_abs_truth_diff == 0 & pk$n_new == pk$n_matched & pk$n_old == pk$n_matched),
            sum(pk$n_matched), any(pk$identical_post_means > 0)))

# 2. Coverage and shape of the intervals by arm and block, with a normal 95%
# interval for coverage from the replicate-based SE.
cov_rows <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(blocks, function(b) {
  x <- d[d$scenario == a & d$block == b, ]; if (!nrow(x)) return(NULL)
  c_r <- per_rep(x, "covered")$x; se <- sd(c_r) / sqrt(length(c_r))
  data.frame(arm = a, block = b, R = length(c_r), coverage = mean(c_r), se = se,
             ci_lower = mean(c_r) - qnorm(0.975) * se, ci_upper = mean(c_r) + qnorm(0.975) * se,
             mean_width = mean(x$width), rmse = sqrt(mean(x$err^2)), bias = mean(x$err),
             slope = if (isTRUE(sd(x$truth) > 0)) unname(coef(lm(post_mean ~ truth, x))[2]) else NA,
             truth_on_mean = if (isTRUE(sd(x$truth) > 0)) unname(coef(lm(truth ~ post_mean, x))[2]) else NA,
             err_to_width = sqrt(mean((x$err / (x$width / zq))^2)))
}))))
write.csv(cov_rows, file.path(out, "coverage-by-arm.csv"), row.names = FALSE)

# 3. Where sigma_b's posterior settles, by arm. The truth is 0.35 in every
# replicate; the default prior's mean is 0.33.
sb <- d[d$block == "sigma_b", ]
sbt <- do.call(rbind, lapply(arms, function(a) {
  x <- sb[sb$scenario == a, ]
  data.frame(arm = a, R = nrow(x), truth = x$truth[1],
             mean_post_mean = mean(x$post_mean), sd_post_mean = sd(x$post_mean),
             min_post_mean = min(x$post_mean), max_post_mean = max(x$post_mean),
             mean_lower = mean(x$lower), mean_upper = mean(x$upper), mean_width = mean(x$width),
             above_truth = mean(x$post_mean > x$truth),
             coverage = mean(x$covered), coverage_se = sd(x$covered) / sqrt(nrow(x)))
}))
write.csv(sbt, file.path(out, "sigma-b-posterior.csv"), row.names = FALSE)

# 4. Paired contrasts on per-replicate coverage, width, squared error and
# posterior mean, all on identical datasets, for B, B0, G and sigma_b. The
# posterior-mean difference matters for sigma_b (does it move?); for the
# other blocks it is a signed average and reported only for completeness.
contrasts <- list(
  c("half-Cauchy, short", "no_spatial_sbhc", "no_spatial"),
  c("weak InvGamma, short", "no_spatial_sbig2", "no_spatial"),
  c("half-Cauchy against weak InvGamma, short", "no_spatial_sbhc", "no_spatial_sbig2"),
  c("half-Cauchy, long (primary)", "no_spatial_sbhc_long", "no_spatial_long"),
  c("weak InvGamma, long", "no_spatial_sbig2_long", "no_spatial_long"),
  c("half-Cauchy against weak InvGamma, long", "no_spatial_sbhc_long", "no_spatial_sbig2_long"),
  c("long chains, given half-Cauchy", "no_spatial_sbhc_long", "no_spatial_sbhc"),
  c("long chains, given weak InvGamma", "no_spatial_sbig2_long", "no_spatial_sbig2"),
  c("long chains, given the default prior", "no_spatial_long", "no_spatial"),
  c("4 latent traits, given half-Cauchy and long", "no_spatial_sbhc_gt4_long", "no_spatial_sbhc_long"),
  c("half-Cauchy, given 4 latent traits and long", "no_spatial_sbhc_gt4_long", "no_spatial_gt4_long"),
  c("half-Cauchy long against the latent-trait probe", "no_spatial_sbhc_long", "no_spatial_gt4_long"),
  c("4 latent traits, given the default prior and long", "no_spatial_gt4_long", "no_spatial_long"))
interactions <- list(
  list("interaction: prior x chain length, (sbhc_long - sbhc) - (long - no_spatial)",
       c("no_spatial_sbhc_long", "no_spatial_sbhc", "no_spatial_long", "no_spatial")),
  list("interaction: prior x chain length, weak InvGamma",
       c("no_spatial_sbig2_long", "no_spatial_sbig2", "no_spatial_long", "no_spatial")),
  list("interaction: prior x latent traits at long, (sbhc_gt4_long - sbhc_long) - (gt4_long - long)",
       c("no_spatial_sbhc_gt4_long", "no_spatial_sbhc_long", "no_spatial_gt4_long", "no_spatial_long")))
rep_stat <- function(a, b, v) {
  x <- d[d$scenario == a & d$block == b, ]; if (v == "se2") x$se2 <- x$err^2
  if (!nrow(x)) return(stats::setNames(rep(NA_real_, length(reps)), reps))  # block absent in this arm
  s <- per_rep(x, v); stats::setNames(s$x, s$replicate)[as.character(reps)]
}
contrast_blocks <- c("B", "B0", "G", "sigma_b")
pc <- do.call(rbind, lapply(contrast_blocks, function(b) {
  rows <- do.call(rbind, lapply(contrasts, function(k) {
    a1 <- d[d$scenario == k[2] & d$block == b, ]; a0 <- d[d$scenario == k[3] & d$block == b, ]
    dd <- rep_stat(k[2], b, "covered") - rep_stat(k[3], b, "covered")
    de <- rep_stat(k[2], b, "se2") - rep_stat(k[3], b, "se2")
    dm <- rep_stat(k[2], b, "post_mean") - rep_stat(k[3], b, "post_mean")
    sl <- function(x) if (isTRUE(sd(x$truth) > 0)) unname(coef(lm(post_mean ~ truth, x))[2]) else NA
    data.frame(block = b, contrast = k[1], minuend = k[2], subtrahend = k[3], R = length(dd),
               coverage_diff = mean(dd), coverage_diff_se = sd(dd) / sqrt(length(dd)),
               width_ratio = mean(a1$width) / mean(a0$width),
               mse_diff = mean(de), mse_diff_se = sd(de) / sqrt(length(de)),
               post_mean_diff = mean(dm), post_mean_diff_se = sd(dm) / sqrt(length(dm)),
               slope_minuend = sl(a1), slope_subtrahend = sl(a0))
  }))
  w <- vapply(arms, function(a) mean(d$width[d$scenario == a & d$block == b]), numeric(1))
  irows <- do.call(rbind, lapply(interactions, function(k) {
    a <- k[[2]]
    ic <- rep_stat(a[1], b, "covered") - rep_stat(a[2], b, "covered") -
          rep_stat(a[3], b, "covered") + rep_stat(a[4], b, "covered")
    ie <- rep_stat(a[1], b, "se2") - rep_stat(a[2], b, "se2") -
          rep_stat(a[3], b, "se2") + rep_stat(a[4], b, "se2")
    im <- rep_stat(a[1], b, "post_mean") - rep_stat(a[2], b, "post_mean") -
          rep_stat(a[3], b, "post_mean") + rep_stat(a[4], b, "post_mean")
    data.frame(block = b, contrast = k[[1]], minuend = "", subtrahend = "", R = length(ic),
               coverage_diff = mean(ic), coverage_diff_se = sd(ic) / sqrt(length(ic)),
               width_ratio = (w[[a[1]]] / w[[a[2]]]) / (w[[a[3]]] / w[[a[4]]]),
               mse_diff = mean(ie), mse_diff_se = sd(ie) / sqrt(length(ie)),
               post_mean_diff = mean(im), post_mean_diff_se = sd(im) / sqrt(length(im)),
               slope_minuend = NA, slope_subtrahend = NA)
  }))
  rbind(rows, irows)
}))
write.csv(pc, file.path(out, "paired-contrasts.csv"), row.names = FALSE)

# 5. B coverage by size of the true effect, by arm, and paired element by
# element for the main contrasts (the B truths are identical in every arm).
# Standard errors treat replicates as clusters.
bins <- c(0, 1, 2, 3, Inf); lab <- c("under 1", "1 to 2", "2 to 3", "over 3")
bb <- d[d$block == "B", ]; bb$bin <- cut(abs(bb$truth), bins, labels = lab, right = FALSE)
eb <- aggregate(cbind(covered, pull = abs(truth) - abs(post_mean), width) ~ scenario + bin, bb, mean)
eb$n <- aggregate(covered ~ scenario + bin, bb, length)$covered
names(eb)[names(eb) == "covered"] <- "coverage"
cluster_se <- function(v, r) {  # SE of a mean over elements, replicates as clusters
  s <- tapply(v, r, sum); n_r <- tapply(v, r, length)
  est <- sum(s) / sum(n_r); R <- length(s)
  sqrt(sum((s - est * n_r)^2) * R / (R - 1)) / sum(n_r)
}
eb$coverage_se <- mapply(function(a, l) { x <- bb[bb$scenario == a & bb$bin == l, ]
  cluster_se(x$covered, x$replicate) }, eb$scenario, as.character(eb$bin))
eb$scenario <- factor(eb$scenario, arms)
eb <- eb[order(eb$scenario, eb$bin), ]
write.csv(eb, file.path(out, "b-by-effect-size.csv"), row.names = FALSE)

flag <- function(a) { x <- bb[bb$scenario == a, ]; x <- x[order(x$replicate, x$element), ]; x }
fl <- lapply(stats::setNames(arms, arms), flag)
same_rows <- function(x, y) identical(x$replicate, y$replicate) &&
  identical(x$element, y$element) && identical(x$truth, y$truth)
stopifnot(all(vapply(fl, same_rows, logical(1), y = fl$no_spatial)))
key <- fl$no_spatial[c("replicate", "bin")]
cv <- function(a) fl[[a]]$covered
diffs <- list(
  "sbhc minus no_spatial" = cv("no_spatial_sbhc") - cv("no_spatial"),
  "sbig2 minus no_spatial" = cv("no_spatial_sbig2") - cv("no_spatial"),
  "sbhc_long minus long" = cv("no_spatial_sbhc_long") - cv("no_spatial_long"),
  "sbig2_long minus long" = cv("no_spatial_sbig2_long") - cv("no_spatial_long"),
  "sbhc_long minus sbhc" = cv("no_spatial_sbhc_long") - cv("no_spatial_sbhc"),
  "sbhc_gt4_long minus gt4_long" = cv("no_spatial_sbhc_gt4_long") - cv("no_spatial_gt4_long"),
  "sbhc_gt4_long minus sbhc_long" = cv("no_spatial_sbhc_gt4_long") - cv("no_spatial_sbhc_long"),
  "sbhc_long minus gt4_long" = cv("no_spatial_sbhc_long") - cv("no_spatial_gt4_long"),
  "sbhc_long minus no_spatial" = cv("no_spatial_sbhc_long") - cv("no_spatial"))
ebc <- do.call(rbind, lapply(names(diffs), function(nm) do.call(rbind, lapply(lab, function(l) {
  i <- key$bin == l; dc <- diffs[[nm]][i]; r <- key$replicate[i]
  s <- tapply(dc, r, sum); n_r <- tapply(dc, r, length)
  est <- sum(s) / sum(n_r); R <- length(s)
  se <- sqrt(sum((s - est * n_r)^2) * R / (R - 1)) / sum(n_r)
  data.frame(contrast = nm, bin = l, n = sum(i), coverage_diff = est, coverage_diff_se = se)
}))))
write.csv(ebc, file.path(out, "b-by-effect-size-contrasts.csv"), row.names = FALSE)

# 6. Chain length on identical datasets under each prior: how far the
# interval endpoints move relative to the interval's width, and how often the
# covered flag changes.
cl <- do.call(rbind, lapply(list(c("no_spatial_sbhc_long", "no_spatial_sbhc"),
                                 c("no_spatial_sbig2_long", "no_spatial_sbig2"),
                                 c("no_spatial_long", "no_spatial")), function(k) {
  do.call(rbind, lapply(c("B", "B0", "sigma_b"), function(b) {
    m <- merge(d[d$scenario == k[1] & d$block == b, ], d[d$scenario == k[2] & d$block == b, ],
               by = c("replicate", "element"))
    data.frame(long = k[1], short = k[2], block = b, n = nrow(m),
               lower_shift_over_width = mean(abs(m$lower.x - m$lower.y) / m$width.y),
               upper_shift_over_width = mean(abs(m$upper.x - m$upper.y) / m$width.y),
               mean_shift_over_width = mean(abs(m$post_mean.x - m$post_mean.y) / m$width.y),
               width_ratio_median = median(m$width.x / m$width.y),
               covered_flag_changed = mean(m$covered.x != m$covered.y),
               gained = mean(m$covered.x > m$covered.y), lost = mean(m$covered.x < m$covered.y))
  }))
}))
write.csv(cl, file.path(out, "chain-length.csv"), row.names = FALSE)

# 7. The pre-registered criteria (README.md), evaluated mechanically.
cov_of <- function(a, b) cov_rows[cov_rows$arm == a & cov_rows$block == b, ]
pc_of <- function(nm, b) pc[pc$contrast == nm & pc$block == b, ]
over3 <- function(a) eb[eb$scenario == a & eb$bin == "over 3", ]
verdict_rows <- list()
add <- function(...) verdict_rows[[length(verdict_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
for (k in list(c("half-Cauchy, short", "no_spatial_sbhc", "no_spatial"),
               c("weak InvGamma, short", "no_spatial_sbig2", "no_spatial"),
               c("half-Cauchy, long (primary)", "no_spatial_sbhc_long", "no_spatial_long"),
               c("weak InvGamma, long", "no_spatial_sbig2_long", "no_spatial_long"),
               c("half-Cauchy, given 4 latent traits and long", "no_spatial_sbhc_gt4_long", "no_spatial_gt4_long"))) {
  cB <- cov_of(k[2], "B"); g <- pc_of(k[1], "B"); b0 <- pc_of(k[1], "B0"); cG <- cov_of(k[2], "G")
  contains <- cB$ci_lower <= 0.95 && cB$ci_upper >= 0.95
  gain <- g$coverage_diff > 2 * g$coverage_diff_se
  g_ok <- cG$coverage <= 0.97
  b0_ok <- b0$coverage_diff >= -2 * b0$coverage_diff_se
  v <- if (gain && contains && g_ok && b0_ok) "restores" else
       if (gain) "partly restores" else
       if (abs(g$coverage_diff) <= 2 * g$coverage_diff_se) "no effect" else "worse"
  add(criterion = paste("B coverage:", k[1]), arm = k[2], against = k[3],
      estimate = cB$coverage, se = cB$se, paired_diff = g$coverage_diff, paired_se = g$coverage_diff_se,
      detail = sprintf("CI %.3f to %.3f contains 0.95: %s; gain > 2 SE: %s; G %.3f <= 0.97: %s; B0 diff %+.3f (SE %.3f) not below -2 SE: %s; over 3: %.3f",
                       cB$ci_lower, cB$ci_upper, contains, gain, cG$coverage, g_ok,
                       b0$coverage_diff, b0$coverage_diff_se, b0_ok, over3(k[2])$coverage),
      verdict = v)
}
two_se <- function(nm, b = "B", label = nm) {
  x <- pc_of(nm, b)
  add(criterion = label, arm = x$minuend, against = x$subtrahend, estimate = NA, se = NA,
      paired_diff = x$coverage_diff, paired_se = x$coverage_diff_se,
      detail = sprintf("%s coverage difference %+.3f (SE %.3f)", b, x$coverage_diff, x$coverage_diff_se),
      verdict = if (x$coverage_diff > 2 * x$coverage_diff_se) "minuend higher by > 2 SE" else
                if (x$coverage_diff < -2 * x$coverage_diff_se) "minuend lower by > 2 SE" else "within 2 SE")
}
two_se("half-Cauchy against weak InvGamma, short", label = "form against weight, short")
two_se("half-Cauchy against weak InvGamma, long", label = "form against weight, long")
two_se("long chains, given half-Cauchy", label = "chain length still needed under half-Cauchy")
two_se("long chains, given weak InvGamma", label = "chain length still needed under weak InvGamma")
two_se("4 latent traits, given half-Cauchy and long", label = "the probe is replaced (sbhc_gt4_long against sbhc_long)")
two_se("half-Cauchy long against the latent-trait probe", label = "half-Cauchy long against gt4_long")
for (k in list(c("half-Cauchy, short", "no_spatial_sbhc"), c("weak InvGamma, short", "no_spatial_sbig2"),
               c("half-Cauchy, long (primary)", "no_spatial_sbhc_long"),
               c("weak InvGamma, long", "no_spatial_sbig2_long"),
               c("half-Cauchy, given 4 latent traits and long", "no_spatial_sbhc_gt4_long"))) {
  x <- pc_of(k[1], "sigma_b"); s <- sbt[sbt$arm == k[2], ]
  add(criterion = paste("sigma_b moves:", k[1]), arm = k[2], against = x$subtrahend,
      estimate = s$mean_post_mean, se = s$sd_post_mean / sqrt(s$R),
      paired_diff = x$post_mean_diff, paired_se = x$post_mean_diff_se,
      detail = sprintf("posterior mean %.3f (SD across replicates %.3f), interval %.3f to %.3f, covers 0.35 at %.2f; paired shift %+.3f (SE %.3f)",
                       s$mean_post_mean, s$sd_post_mean, s$mean_lower, s$mean_upper, s$coverage,
                       x$post_mean_diff, x$post_mean_diff_se),
      verdict = if (abs(x$post_mean_diff) > 2 * x$post_mean_diff_se) "moves" else "does not move")
}
criteria <- do.call(rbind, verdict_rows)
write.csv(criteria, file.path(out, "criteria.csv"), row.names = FALSE)

options(width = 200)
cat("\n== pairing check (summary)\n")
print(aggregate(cbind(n_matched, max_abs_truth_diff, identical_post_means) ~ arm + against, pk,
                function(v) c(sum = sum(v), max = max(v))[if (length(v)) 1 else 1]), digits = 3, row.names = FALSE)
cat("\n== coverage by arm (B, B0, G, sigma_b)\n")
print(cov_rows[cov_rows$block %in% c("B", "B0", "G", "sigma_b"), ], digits = 3, row.names = FALSE)
cat("\n== sigma_b posterior by arm\n"); print(sbt, digits = 3, row.names = FALSE)
cat("\n== paired contrasts\n"); print(pc, digits = 3, row.names = FALSE)
cat("\n== B by effect size\n"); print(eb[order(eb$bin, eb$scenario), ], digits = 3, row.names = FALSE)
cat("\n== B by effect size, paired\n"); print(ebc, digits = 3, row.names = FALSE)
cat("\n== chain length\n"); print(cl, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(criteria[c("criterion", "verdict", "detail")], row.names = FALSE, right = FALSE)
