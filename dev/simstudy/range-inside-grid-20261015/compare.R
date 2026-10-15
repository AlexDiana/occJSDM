# Compare the range-inside-grid cells with production base and binary.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/range-inside-grid-20261015/compare.R [results.rds]
#
# Without an argument it reads the single simstudy-*.rds in results/. The
# production rows are replicates 1 to 50 of `base` and `binary` from the
# R = 300 run of 15 October 2026. Summaries follow the conventions of
# dev/simstudy/paper-study-20261015/calibration/calibration.R: coverage is
# the mean of per-replicate coverage fractions and its standard error is
# their SD over sqrt(R), never a pooled-element binomial SE.

args <- commandArgs(trailingOnly = TRUE)
run_dir <- "dev/simstudy/range-inside-grid-20261015"
res_path <- if (length(args) >= 1) args[1] else
  list.files(file.path(run_dir, "results"), "^simstudy-.*\\.rds$", full.names = TRUE)
stopifnot(length(res_path) == 1)

prod_all <- readRDS(file.path("dev/simstudy/paper-study-20261015/results",
                              "simstudy-20261015-172834.rds"))$rows
diag <- readRDS(res_path)
new <- diag$rows
R <- 50L
stopifnot(diag$failed_replicates == 0)

parent_of <- c(base_ls05 = "base", base_ls075 = "base",
               binary_ls05 = "binary", binary_ls075 = "binary")
l_s_of <- c(base = 0.15, binary = 0.15, base_ls05 = 0.05, binary_ls05 = 0.05,
            base_ls075 = 0.075, binary_ls075 = 0.075)
cells <- c("base", "base_ls075", "base_ls05", "binary", "binary_ls075", "binary_ls05")
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q")

prod <- prod_all[prod_all$scenario %in% c("base", "binary") &
                 prod_all$replicate <= R, ]
d <- rbind(prod, new[names(prod)])
d <- d[d$block %in% blocks, ]
d$width <- d$upper - d$lower
d$err <- d$post_mean - d$truth
for (sc in cells) stopifnot(setequal(unique(d$replicate[d$scenario == sc]), seq_len(R)))

# --- Pairing: the new cells share their parent's coefficients ---------------

key <- function(x) paste(x$replicate, x$block, x$element, sep = "|")
pairing <- do.call(rbind, lapply(names(parent_of), function(sc) {
  do.call(rbind, lapply(c("B0", "B", "G"), function(bl) {
    a <- d[d$scenario == sc & d$block == bl, ]
    b <- d[d$scenario == parent_of[[sc]] & d$block == bl, ]
    m <- match(key(a), key(b))
    stopifnot(!anyNA(m), nrow(a) == nrow(b))
    data.frame(scenario = sc, parent = parent_of[[sc]], block = bl,
               elements = nrow(a), max_abs_truth_diff = max(abs(a$truth - b$truth[m])))
  }))
}))
stopifnot(all(pairing$max_abs_truth_diff == 0))
write.csv(pairing, file.path(run_dir, "pairing-check.csv"), row.names = FALSE)

# --- Replicate-based summaries ---------------------------------------------

zsd <- function(err, w) sqrt(mean((err / (w / (2 * qnorm(0.975))))^2))

summarise_cell_block <- function(sub) {
  c_r <- tapply(sub$covered, sub$replicate, mean)
  w_r <- tapply(sub$width, sub$replicate, mean)
  e_r <- tapply(sub$err, sub$replicate, mean)
  est <- mean(c_r); se <- sd(c_r) / sqrt(length(c_r))
  data.frame(
    scenario = sub$scenario[1], block = sub$block[1],
    l_s_raw = l_s_of[[sub$scenario[1]]], R = length(c_r),
    elements_per_rep = nrow(sub) / length(c_r),
    coverage = est, se = se,
    z_vs_nominal = if (se > 0) (est - 0.95) / se else NA_real_,
    mean_width = mean(w_r), width_se = sd(w_r) / sqrt(length(w_r)),
    bias = mean(sub$err), bias_se = sd(e_r) / sqrt(length(e_r)),
    rmse = sqrt(mean(sub$err^2)),
    err_to_width = zsd(sub$err, sub$width),
    slope = if (stats::var(sub$truth) > 0)
      unname(coef(lm(post_mean ~ truth, sub))[2]) else NA_real_)
}

summ <- do.call(rbind, lapply(cells, function(sc) do.call(rbind, lapply(blocks, function(bl) {
  sub <- d[d$scenario == sc & d$block == bl, ]
  if (!nrow(sub)) return(NULL)
  summarise_cell_block(sub)
}))))
rownames(summ) <- NULL
write.csv(summ, file.path(run_dir, "coverage-r50.csv"), row.names = FALSE)

# --- Paired differences against the parent ----------------------------------

paired <- do.call(rbind, lapply(names(parent_of), function(sc) do.call(rbind, lapply(blocks, function(bl) {
  a <- d[d$scenario == sc & d$block == bl, ]
  b <- d[d$scenario == parent_of[[sc]] & d$block == bl, ]
  if (!nrow(a)) return(NULL)
  ca <- tapply(a$covered, a$replicate, mean)
  cb <- tapply(b$covered, b$replicate, mean)[names(ca)]
  wa <- tapply(a$width, a$replicate, mean)
  wb <- tapply(b$width, b$replicate, mean)[names(wa)]
  dd <- ca - cb
  data.frame(
    scenario = sc, parent = parent_of[[sc]], block = bl,
    parent_coverage = mean(cb), coverage = mean(ca),
    diff = mean(dd), diff_se = sd(dd) / sqrt(length(dd)),
    parent_width = mean(wb), width = mean(wa),
    width_ratio = mean(wa) / mean(wb),
    width_ratio_se = sd(wa / wb) / sqrt(length(wa)),
    parent_rmse = sqrt(mean(b$err^2)), rmse = sqrt(mean(a$err^2)),
    parent_slope = unname(coef(lm(post_mean ~ truth, b))[2]),
    slope = unname(coef(lm(post_mean ~ truth, a))[2]),
    parent_err_to_width = zsd(b$err, b$width),
    err_to_width = zsd(a$err, a$width))
}))))
rownames(paired) <- NULL
write.csv(paired, file.path(run_dir, "paired-vs-production.csv"), row.names = FALSE)

# --- B and B0 coverage by size of the true effect ---------------------------

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

# --- Console digest -----------------------------------------------------------

show <- summ[summ$block %in% c("B0", "B"), c("scenario", "block", "coverage", "se",
                                              "mean_width", "slope", "err_to_width", "rmse")]
show[-(1:2)] <- lapply(show[-(1:2)], round, 3)
print(show[order(show$block, match(show$scenario, cells)), ], row.names = FALSE)
pd <- paired[paired$block %in% c("B0", "B"), c("scenario", "block", "diff", "diff_se",
                                                "width_ratio", "parent_slope", "slope")]
pd[-(1:2)] <- lapply(pd[-(1:2)], round, 3)
print(pd, row.names = FALSE)
cat("written to", run_dir, "\n")
