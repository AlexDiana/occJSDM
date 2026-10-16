# Summarise no_spatial_gt4_long against no_spatial, no_spatial_gt4 and
# no_spatial_long, which the non-spatial candidates run fitted to the same
# datasets at revision c40f385.
# Run from the repository root: Rscript dev/simstudy/b-gt4-long-20261016/analyse.R
# Reads this run's results file, or its checkpoint while the run is in flight,
# and the candidates run's results for the three paired arms.
# Writes pairing-check.csv, coverage-by-arm.csv, paired-contrasts.csv,
# b-by-effect-size.csv, b-by-effect-size-contrasts.csv and chain-length.csv.
here <- "dev/simstudy/b-gt4-long-20261016"
out <- Sys.getenv("OUT_DIR", here)  # where the CSVs go; the run directory by default
prev_file <- "dev/simstudy/b-nonspatial-candidates-20261015/results/simstudy-20261016-004116.rds"
res_dir <- Sys.getenv("RESULTS_DIR", file.path(here, "results"))  # override to test on another run
res <- list.files(res_dir, pattern = "^simstudy-.*[.]rds$", full.names = TRUE)
run <- if (length(res)) readRDS(res[1]) else NULL
new <- if (length(res)) run$rows else
  readRDS(list.files(res_dir, pattern = "^checkpoint", full.names = TRUE)[1])$params
prev <- readRDS(prev_file)
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
arms <- c("no_spatial", "no_spatial_long", "no_spatial_gt4", "no_spatial_gt4_long")
old_arms <- arms[1:3]
d <- rbind(prev$rows[prev$rows$scenario %in% old_arms, cols],
           new[new$scenario == "no_spatial_gt4_long" & new$block != "ERROR", cols])
d$width <- d$upper - d$lower; d$err <- d$post_mean - d$truth
reps <- Reduce(intersect, lapply(arms, function(a) unique(d$replicate[d$scenario == a])))
d <- d[d$replicate %in% reps, ]
cat(sprintf("replicates present in all four arms: %d\n", length(reps)))
if (length(res)) {
  pv <- run$provenance
  cat(sprintf("failed replicates: %d; revision %s, dirty at start %s, at write %s, moved %s\n",
              run$failed_replicates, substr(pv$git_sha, 1, 7), pv$git_dirty,
              pv$git_dirty_at_write, pv$git_moved_during_run))
  cat("cell MCMC settings recorded:", deparse(pv$cell_MCMCparams), "\n")
}
blocks <- c("B", "B0", "G", "sigma_b", "beta_theta", "theta0", "p", "q")
zq <- 2 * qnorm(0.975)
per_rep <- function(x, v) aggregate(x[[v]], list(replicate = x$replicate), mean)

# 1. Pairing, checked first: are the truths element for element the same in
# all four arms, with the same rows present? And does the new arm share no
# posterior with the old ones (it should not)?
pk <- do.call(rbind, lapply(old_arms, function(a) do.call(rbind, lapply(unique(d$block), function(b) {
  x <- d[d$scenario == "no_spatial_gt4_long" & d$block == b, ]
  y <- d[d$scenario == a & d$block == b, ]
  m <- merge(x, y, by = c("replicate", "element"))
  data.frame(arm = a, block = b, n_new = nrow(x), n_old = nrow(y), n_matched = nrow(m),
             max_abs_truth_diff = if (nrow(m)) max(abs(m$truth.x - m$truth.y)) else NA,
             identical_post_means = if (nrow(m)) mean(m$post_mean.x == m$post_mean.y) else NA)
}))))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)
cat(sprintf("truths identical element for element in every block and arm: %s (%d rows compared)\n",
            all(pk$max_abs_truth_diff == 0 & pk$n_new == pk$n_matched & pk$n_old == pk$n_matched),
            sum(pk$n_matched)))

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

# 3. Paired contrasts on per-replicate coverage, width and squared error, all
# on identical datasets. The last row of each block is the interaction,
# (gt4_long - gt4) - (long - no_spatial): zero when the two effects add.
contrasts <- list(
  c("both against neither", "no_spatial_gt4_long", "no_spatial"),
  c("long chains, given 4 latent traits", "no_spatial_gt4_long", "no_spatial_gt4"),
  c("4 latent traits, given long chains", "no_spatial_gt4_long", "no_spatial_long"),
  c("4 latent traits alone", "no_spatial_gt4", "no_spatial"),
  c("long chains alone", "no_spatial_long", "no_spatial"))
rep_stat <- function(a, b, v) {
  x <- d[d$scenario == a & d$block == b, ]; if (v == "se2") x$se2 <- x$err^2
  s <- per_rep(x, v); stats::setNames(s$x, s$replicate)[as.character(reps)]
}
pc <- do.call(rbind, lapply(c("B", "B0"), function(b) {
  rows <- do.call(rbind, lapply(contrasts, function(k) {
    a1 <- d[d$scenario == k[2] & d$block == b, ]; a0 <- d[d$scenario == k[3] & d$block == b, ]
    dd <- rep_stat(k[2], b, "covered") - rep_stat(k[3], b, "covered")
    de <- rep_stat(k[2], b, "se2") - rep_stat(k[3], b, "se2")
    data.frame(block = b, contrast = k[1], minuend = k[2], subtrahend = k[3], R = length(dd),
               coverage_diff = mean(dd), coverage_diff_se = sd(dd) / sqrt(length(dd)),
               width_ratio = mean(a1$width) / mean(a0$width),
               mse_diff = mean(de), mse_diff_se = sd(de) / sqrt(length(de)),
               slope_minuend = unname(coef(lm(post_mean ~ truth, a1))[2]),
               slope_subtrahend = unname(coef(lm(post_mean ~ truth, a0))[2]))
  }))
  cv <- lapply(stats::setNames(arms, arms), function(a) rep_stat(a, b, "covered"))
  sq <- lapply(stats::setNames(arms, arms), function(a) rep_stat(a, b, "se2"))
  ic <- cv$no_spatial_gt4_long - cv$no_spatial_gt4 - cv$no_spatial_long + cv$no_spatial
  ie <- sq$no_spatial_gt4_long - sq$no_spatial_gt4 - sq$no_spatial_long + sq$no_spatial
  w <- vapply(arms, function(a) mean(d$width[d$scenario == a & d$block == b]), numeric(1))
  rbind(rows, data.frame(block = b, contrast = "interaction: (gt4_long - gt4) - (long - no_spatial)",
                         minuend = "", subtrahend = "", R = length(ic),
                         coverage_diff = mean(ic), coverage_diff_se = sd(ic) / sqrt(length(ic)),
                         width_ratio = (w[["no_spatial_gt4_long"]] / w[["no_spatial_gt4"]]) /
                                       (w[["no_spatial_long"]] / w[["no_spatial"]]),
                         mse_diff = mean(ie), mse_diff_se = sd(ie) / sqrt(length(ie)),
                         slope_minuend = NA, slope_subtrahend = NA))
}))
write.csv(pc, file.path(out, "paired-contrasts.csv"), row.names = FALSE)

# 4. B coverage by size of the true effect, by arm, and paired element by
# element against each of the other three (the B truths are identical in every
# arm), with the interaction per bin. Standard errors treat replicates as
# clusters.
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
write.csv(eb[order(eb$scenario, eb$bin), ], file.path(out, "b-by-effect-size.csv"), row.names = FALSE)

flag <- function(a) { x <- bb[bb$scenario == a, ]; x <- x[order(x$replicate, x$element), ]; x }
fl <- lapply(stats::setNames(arms, arms), flag)
same_rows <- function(x, y) identical(x$replicate, y$replicate) &&
  identical(x$element, y$element) && identical(x$truth, y$truth)
stopifnot(all(vapply(fl, same_rows, logical(1), y = fl$no_spatial)))
key <- fl$no_spatial[c("replicate", "bin")]
diffs <- list(
  "gt4_long minus no_spatial" = fl$no_spatial_gt4_long$covered - fl$no_spatial$covered,
  "gt4_long minus gt4" = fl$no_spatial_gt4_long$covered - fl$no_spatial_gt4$covered,
  "gt4_long minus long" = fl$no_spatial_gt4_long$covered - fl$no_spatial_long$covered,
  "interaction" = fl$no_spatial_gt4_long$covered - fl$no_spatial_gt4$covered -
                  fl$no_spatial_long$covered + fl$no_spatial$covered)
ebc <- do.call(rbind, lapply(names(diffs), function(nm) do.call(rbind, lapply(lab, function(l) {
  i <- key$bin == l; dc <- diffs[[nm]][i]; r <- key$replicate[i]
  s <- tapply(dc, r, sum); n_r <- tapply(dc, r, length)
  est <- sum(s) / sum(n_r); R <- length(s)
  se <- sqrt(sum((s - est * n_r)^2) * R / (R - 1)) / sum(n_r)
  data.frame(contrast = nm, bin = l, n = sum(i), coverage_diff = est, coverage_diff_se = se)
}))))
write.csv(ebc, file.path(out, "b-by-effect-size-contrasts.csv"), row.names = FALSE)

# 5. Chain length on identical datasets, at 4 latent traits and, for
# reference, at 2: how far the interval endpoints move relative to the
# interval's width, and how often the covered flag changes.
cl <- do.call(rbind, lapply(list(c("no_spatial_gt4_long", "no_spatial_gt4"),
                                 c("no_spatial_long", "no_spatial")), function(k) {
  do.call(rbind, lapply(c("B", "B0"), function(b) {
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

options(width = 160)
cat("\n== pairing check\n"); print(pk, digits = 3, row.names = FALSE)
cat("\n== coverage by arm (B, B0, G, sigma_b)\n")
print(cov_rows[cov_rows$block %in% c("B", "B0", "G", "sigma_b"), ], digits = 3, row.names = FALSE)
cat("\n== paired contrasts\n"); print(pc, digits = 3, row.names = FALSE)
cat("\n== B by effect size\n"); print(eb[order(eb$bin, eb$scenario), ], digits = 3, row.names = FALSE)
cat("\n== B by effect size, paired\n"); print(ebc, digits = 3, row.names = FALSE)
cat("\n== chain length\n"); print(cl, digits = 3, row.names = FALSE)
