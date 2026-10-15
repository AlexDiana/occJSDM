# Compare the free-amplitude (half-Cauchy) cells with the same datasets
# fitted under the default amplitude prior.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/amplitude-prior-20261015/compare.R
#
# Default-prior rows come from replicates 1 to 50 of the R = 300 production
# run (`base`, `binary`) and from the range-inside-grid run (`base_ls075`,
# `binary_ls075`). Every pair shares its datasets. Conventions as in
# dev/simstudy/paper-study-20261015/calibration/calibration.R: coverage is the
# mean of per-replicate coverage fractions and its standard error is their SD
# over sqrt(R).

run_dir <- "dev/simstudy/amplitude-prior-20261015"
one_rds <- function(dir) {
  f <- list.files(dir, "^simstudy-.*\\.rds$", full.names = TRUE)
  stopifnot(length(f) == 1)
  readRDS(f)
}
prod <- readRDS(file.path("dev/simstudy/paper-study-20261015/results",
                          "simstudy-20261015-172834.rds"))$rows
rng <- one_rds("dev/simstudy/range-inside-grid-20261015/results")
hc <- one_rds(file.path(run_dir, "results"))
stopifnot(rng$failed_replicates == 0, hc$failed_replicates == 0)
R <- 50L

cols <- names(prod)
d <- rbind(prod[prod$scenario %in% c("base", "binary") & prod$replicate <= R, ],
           rng$rows[rng$rows$scenario %in% c("base_ls075", "binary_ls075"), cols],
           hc$rows[cols])
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q")
d <- d[d$block %in% blocks, ]
d$width <- d$upper - d$lower
d$err <- d$post_mean - d$truth
cells <- c("base", "base_hc", "base_ls075", "base_ls075_hc",
           "binary", "binary_hc", "binary_ls075", "binary_ls075_hc")
for (sc in cells) stopifnot(setequal(unique(d$replicate[d$scenario == sc]), seq_len(R)))

# Pairs: the prior's effect at each range, and prior plus range together.
pairs <- data.frame(
  scenario = c("base_hc", "base_ls075_hc", "base_ls075_hc",
               "binary_hc", "binary_ls075_hc", "binary_ls075_hc"),
  reference = c("base", "base_ls075", "base",
                "binary", "binary_ls075", "binary"),
  contrast = rep(c("prior at production range", "prior at 0.075 raw",
                   "prior and range together"), 2),
  stringsAsFactors = FALSE)

key <- function(x) paste(x$replicate, x$block, x$element, sep = "|")
pairing <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  do.call(rbind, lapply(c("B0", "B", "G"), function(bl) {
    a <- d[d$scenario == pairs$scenario[i] & d$block == bl, ]
    b <- d[d$scenario == pairs$reference[i] & d$block == bl, ]
    m <- match(key(a), key(b))
    stopifnot(!anyNA(m), nrow(a) == nrow(b))
    data.frame(scenario = pairs$scenario[i], reference = pairs$reference[i],
               block = bl, elements = nrow(a),
               max_abs_truth_diff = max(abs(a$truth - b$truth[m])))
  }))
}))
stopifnot(all(pairing$max_abs_truth_diff == 0))
write.csv(pairing, file.path(run_dir, "pairing-check.csv"), row.names = FALSE)

zsd <- function(err, w) sqrt(mean((err / (w / (2 * qnorm(0.975))))^2))
slope <- function(x) if (stats::var(x$truth) > 0)
  unname(coef(lm(post_mean ~ truth, x))[2]) else NA_real_

summ <- do.call(rbind, lapply(cells, function(sc) do.call(rbind, lapply(blocks, function(bl) {
  sub <- d[d$scenario == sc & d$block == bl, ]
  if (!nrow(sub)) return(NULL)
  c_r <- tapply(sub$covered, sub$replicate, mean)
  w_r <- tapply(sub$width, sub$replicate, mean)
  est <- mean(c_r); se <- sd(c_r) / sqrt(length(c_r))
  data.frame(scenario = sc, block = bl, R = length(c_r),
             coverage = est, se = se,
             z_vs_nominal = if (se > 0) (est - 0.95) / se else NA_real_,
             mean_width = mean(w_r), bias = mean(sub$err),
             rmse = sqrt(mean(sub$err^2)), err_to_width = zsd(sub$err, sub$width),
             slope = slope(sub))
}))))
rownames(summ) <- NULL
write.csv(summ, file.path(run_dir, "coverage-r50.csv"), row.names = FALSE)

paired <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) do.call(rbind, lapply(blocks, function(bl) {
  a <- d[d$scenario == pairs$scenario[i] & d$block == bl, ]
  b <- d[d$scenario == pairs$reference[i] & d$block == bl, ]
  if (!nrow(a)) return(NULL)
  ca <- tapply(a$covered, a$replicate, mean)
  cb <- tapply(b$covered, b$replicate, mean)[names(ca)]
  wa <- tapply(a$width, a$replicate, mean)
  wb <- tapply(b$width, b$replicate, mean)[names(wa)]
  dd <- ca - cb
  data.frame(scenario = pairs$scenario[i], reference = pairs$reference[i],
             contrast = pairs$contrast[i], block = bl,
             reference_coverage = mean(cb), coverage = mean(ca),
             diff = mean(dd), diff_se = sd(dd) / sqrt(length(dd)),
             width_ratio = mean(wa) / mean(wb),
             width_ratio_se = sd(wa / wb) / sqrt(length(wa)),
             reference_rmse = sqrt(mean(b$err^2)), rmse = sqrt(mean(a$err^2)),
             reference_slope = slope(b), slope = slope(a),
             reference_err_to_width = zsd(b$err, b$width),
             err_to_width = zsd(a$err, a$width))
}))))
rownames(paired) <- NULL
write.csv(paired, file.path(run_dir, "paired-vs-default-prior.csv"), row.names = FALSE)

band <- function(x) cut(abs(x), c(0, 1, 2, 3, Inf), right = FALSE,
                        labels = c("<1", "1-2", "2-3", ">=3"))
by_size <- do.call(rbind, lapply(c("B", "B0"), function(bl) do.call(rbind, lapply(cells, function(sc) {
  sub <- d[d$scenario == sc & d$block == bl, ]
  sub$band <- band(sub$truth)
  do.call(rbind, lapply(levels(sub$band), function(lv) {
    s <- sub[sub$band == lv, ]
    if (!nrow(s)) return(NULL)
    data.frame(block = bl, scenario = sc, abs_truth = lv, elements = nrow(s),
               coverage = mean(s$covered), mean_width = mean(s$width),
               mean_shrink = mean(sign(s$truth) * (s$truth - s$post_mean)))
  }))
}))))
write.csv(by_size, file.path(run_dir, "b-b0-by-effect-size.csv"), row.names = FALSE)

show <- summ[summ$block %in% c("B0", "B"), c("scenario", "block", "coverage", "se",
                                              "mean_width", "slope", "err_to_width", "rmse")]
show[-(1:2)] <- lapply(show[-(1:2)], round, 3)
print(show[order(show$block, match(show$scenario, cells)), ], row.names = FALSE)
pd <- paired[paired$block %in% c("B0", "B"),
             c("scenario", "reference", "block", "diff", "diff_se", "width_ratio",
               "reference_slope", "slope")]
pd[-(1:3)] <- lapply(pd[-(1:3)], round, 3)
print(pd, row.names = FALSE)
op <- paired[!paired$block %in% c("B0", "B"), c("scenario", "reference", "block", "diff", "diff_se", "width_ratio")]
op[-(1:3)] <- lapply(op[-(1:3)], round, 3)
print(op, row.names = FALSE)
cat("written to", run_dir, "\n")
