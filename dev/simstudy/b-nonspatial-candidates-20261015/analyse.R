# Summarise the non-spatial B-candidate arms against no_spatial.
# Run from the repository root: Rscript dev/simstudy/b-nonspatial-candidates-20261015/analyse.R
# Reads the run's results file, or its checkpoint while the run is in flight, and
# the spatial 2x2's results for the check that no_spatial reproduces.
# Writes coverage-by-arm.csv, paired-contrasts.csv, b-by-effect-size.csv,
# b-by-effect-size-contrasts.csv, chain-length.csv and pairing-check.csv.
here <- "dev/simstudy/b-nonspatial-candidates-20261015"
out <- Sys.getenv("OUT_DIR", here)  # where the CSVs go; the run directory by default
twobytwo_file <- "dev/simstudy/spatial-field-2x2-20261015/results/simstudy-20261015-222234.rds"
res_dir <- Sys.getenv("RESULTS_DIR", file.path(here, "results"))  # override to test on another run
res <- list.files(res_dir, pattern = "^simstudy-.*[.]rds$", full.names = TRUE)
run <- if (length(res)) readRDS(res[1]) else NULL
new <- if (length(res)) run$rows else
  readRDS(list.files(res_dir, pattern = "^checkpoint", full.names = TRUE)[1])$params
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
d <- new[new$block != "ERROR", cols]
d$width <- d$upper - d$lower; d$err <- d$post_mean - d$truth
arms <- c("no_spatial", "no_spatial_long", "no_factors", "no_factors_fit2",
          "no_spatial_notraits", "no_spatial_gt4")
d <- d[d$scenario %in% arms, ]
reps <- Reduce(intersect, lapply(arms, function(a) unique(d$replicate[d$scenario == a])))
d <- d[d$replicate %in% reps, ]
cat(sprintf("replicates present in all six arms: %d\n", length(reps)))
if (length(res)) cat(sprintf("failed replicates: %d\n", run$failed_replicates))
blocks <- c("B", "B0", "G", "sigma_b", "beta_theta", "theta0", "p", "q")
zq <- 2 * qnorm(0.975)

per_rep <- function(x, v) aggregate(x[[v]], list(replicate = x$replicate), mean)

# 1. Coverage and shape of the intervals by arm and block. slope is the
# posterior mean regressed on the truth; truth_on_mean is the reverse, which is
# near 1 when the posterior means are calibrated for these truths even though
# shrinkage keeps slope below 1.
cov_rows <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(blocks, function(b) {
  x <- d[d$scenario == a & d$block == b, ]; if (!nrow(x)) return(NULL)
  c_r <- per_rep(x, "covered")$x
  data.frame(arm = a, block = b, R = length(c_r), coverage = mean(c_r), se = sd(c_r) / sqrt(length(c_r)),
             mean_width = mean(x$width), rmse = sqrt(mean(x$err^2)), bias = mean(x$err),
             slope = if (isTRUE(sd(x$truth) > 0)) unname(coef(lm(post_mean ~ truth, x))[2]) else NA,
             truth_on_mean = if (isTRUE(sd(x$truth) > 0)) unname(coef(lm(truth ~ post_mean, x))[2]) else NA,
             err_to_width = sqrt(mean((x$err / (x$width / zq))^2)))
}))))
write.csv(cov_rows, file.path(out, "coverage-by-arm.csv"), row.names = FALSE)

# 2. Paired contrasts on per-replicate coverage, width and error. same_data
# says whether the two arms fit identical datasets.
contrasts <- list(
  c("longer chains", "no_spatial_long", "no_spatial", TRUE),
  c("no factors in truth or fit", "no_factors", "no_spatial", FALSE),
  c("fitting 2 factors that are not there", "no_factors_fit2", "no_factors", TRUE),
  c("factors in the truth, 2 fitted", "no_spatial", "no_factors_fit2", FALSE),
  c("traits left out of the fit", "no_spatial_notraits", "no_spatial", TRUE),
  c("4 latent traits fitted instead of 2", "no_spatial_gt4", "no_spatial", TRUE))
pc <- do.call(rbind, lapply(c("B", "B0"), function(b) do.call(rbind, lapply(contrasts, function(k) {
  a1 <- d[d$scenario == k[2] & d$block == b, ]; a0 <- d[d$scenario == k[3] & d$block == b, ]
  m <- merge(per_rep(a1, "covered"), per_rep(a0, "covered"), by = "replicate")
  w <- merge(per_rep(a1, "width"), per_rep(a0, "width"), by = "replicate")
  a1$se2 <- a1$err^2; a0$se2 <- a0$err^2
  e <- merge(per_rep(a1, "se2"), per_rep(a0, "se2"), by = "replicate")
  dd <- m$x.x - m$x.y; de <- e$x.x - e$x.y
  data.frame(block = b, contrast = k[1], minuend = k[2], subtrahend = k[3], same_data = as.logical(k[4]),
             R = nrow(m), coverage_diff = mean(dd), coverage_diff_se = sd(dd) / sqrt(nrow(m)),
             width_ratio = mean(w$x.x) / mean(w$x.y),
             mse_diff = mean(de), mse_diff_se = sd(de) / sqrt(nrow(e)),
             rmse_minuend = sqrt(mean(a1$err^2)), rmse_subtrahend = sqrt(mean(a0$err^2)),
             slope_minuend = unname(coef(lm(post_mean ~ truth, a1))[2]),
             slope_subtrahend = unname(coef(lm(post_mean ~ truth, a0))[2]))
}))))
write.csv(pc, file.path(out, "paired-contrasts.csv"), row.names = FALSE)

# 3. B coverage by size of the true effect, by arm, and paired against
# no_spatial element by element (the B truths are identical in every arm). The
# standard error of a bin's difference treats replicates as clusters.
bins <- c(0, 1, 2, 3, Inf); lab <- c("under 1", "1 to 2", "2 to 3", "over 3")
bb <- d[d$block == "B", ]; bb$bin <- cut(abs(bb$truth), bins, labels = lab, right = FALSE)
eb <- aggregate(cbind(covered, pull = abs(truth) - abs(post_mean), width) ~ scenario + bin, bb, mean)
eb$n <- aggregate(covered ~ scenario + bin, bb, length)$covered
names(eb)[names(eb) == "covered"] <- "coverage"
write.csv(eb[order(eb$scenario, eb$bin), ], file.path(out, "b-by-effect-size.csv"), row.names = FALSE)

ref <- bb[bb$scenario == "no_spatial", c("replicate", "element", "bin", "covered")]
ebc <- do.call(rbind, lapply(setdiff(arms, "no_spatial"), function(a) {
  m <- merge(bb[bb$scenario == a, c("replicate", "element", "covered")], ref, by = c("replicate", "element"))
  m$dc <- m$covered.x - m$covered.y
  do.call(rbind, lapply(lab, function(l) {
    x <- m[m$bin == l, ]
    s <- tapply(x$dc, x$replicate, sum); n_r <- tapply(x$dc, x$replicate, length)
    est <- sum(s) / sum(n_r); R <- length(s)
    se <- sqrt(sum((s - est * n_r)^2) * R / (R - 1)) / sum(n_r)
    data.frame(arm = a, bin = l, n = nrow(x), coverage_diff = est, coverage_diff_se = se)
  }))
}))
write.csv(ebc, file.path(out, "b-by-effect-size-contrasts.csv"), row.names = FALSE)

# 4. Chain length on identical datasets: how far the interval endpoints move,
# relative to the interval's width, and how often the covered flag changes.
cl <- do.call(rbind, lapply(c("B", "B0"), function(b) {
  m <- merge(d[d$scenario == "no_spatial_long" & d$block == b, ],
             d[d$scenario == "no_spatial" & d$block == b, ], by = c("replicate", "element"))
  data.frame(block = b, n = nrow(m),
             lower_shift_over_width = mean(abs(m$lower.x - m$lower.y) / m$width.y),
             upper_shift_over_width = mean(abs(m$upper.x - m$upper.y) / m$width.y),
             mean_shift_over_width = mean(abs(m$post_mean.x - m$post_mean.y) / m$width.y),
             width_ratio_median = median(m$width.x / m$width.y),
             covered_flag_changed = mean(m$covered.x != m$covered.y),
             gained = mean(m$covered.x > m$covered.y), lost = mean(m$covered.x < m$covered.y))
}))
write.csv(cl, file.path(out, "chain-length.csv"), row.names = FALSE)

# 5. Pairing: are the truths element-for-element no_spatial's? Do the
# fit-only arms share no_spatial's posterior for nothing (they should not)?
# And does this run's no_spatial reproduce the spatial 2x2's?
pk <- do.call(rbind, lapply(setdiff(arms, "no_spatial"), function(a) do.call(rbind, lapply(blocks, function(b) {
  x <- d[d$scenario == "no_spatial" & d$block == b, ]; y <- d[d$scenario == a & d$block == b, ]
  m <- merge(x, y, by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  data.frame(arm = a, block = b, n = nrow(m), max_abs_truth_diff = max(abs(m$truth.x - m$truth.y)),
             identical_post_means = mean(m$post_mean.x == m$post_mean.y))
}))))
old <- readRDS(twobytwo_file)$rows
old <- old[old$scenario == "no_spatial" & old$replicate %in% reps, cols]
mo <- merge(old, new[new$scenario == "no_spatial", cols], by = c("replicate", "block", "element"))
repro <- data.frame(arm = "no_spatial vs the 2x2's no_spatial", block = "all", n = nrow(mo),
                    max_abs_truth_diff = max(abs(mo$truth.x - mo$truth.y)),
                    identical_post_means = mean(mo$post_mean.x == mo$post_mean.y &
                                                mo$lower.x == mo$lower.y & mo$upper.x == mo$upper.y))
pk <- rbind(pk, repro)
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

# How much of the linear predictor the factors carry in no_spatial's truth.
if (length(res)) {
  occ <- run$occupancy[run$occupancy$scenario %in% c("no_spatial", "no_factors"), ]
  cat("\n== true variance shares of the linear predictor, mean over species\n")
  print(aggregate(cbind(share_env, share_lat) ~ scenario, occ, mean), digits = 3, row.names = FALSE)
}

options(width = 160)
cat("\n== coverage by arm (B, B0, G, sigma_b)\n")
print(cov_rows[cov_rows$block %in% c("B", "B0", "G", "sigma_b"), ], digits = 3, row.names = FALSE)
cat("\n== paired contrasts\n"); print(pc, digits = 3, row.names = FALSE)
cat("\n== B by effect size\n"); print(eb[order(eb$bin, eb$scenario), ], digits = 3, row.names = FALSE)
cat("\n== B by effect size, paired against no_spatial\n"); print(ebc, digits = 3, row.names = FALSE)
cat("\n== chain length\n"); print(cl, digits = 3, row.names = FALSE)
cat("\n== pairing check\n"); print(pk, digits = 3, row.names = FALSE)
