# Summarise amp-diag.R's output into amp-diag.csv.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/amplitude-prior-20261015/amp-diag-summary.R \
#     dev/simstudy/amplitude-prior-20261015/results/amp-diag.rds
#
# Also checks that each refit is the same fit as the matching run replicate:
# B and B0 posterior means must equal the rows of the R = 300 production run
# (base, binary), the range-inside-grid run (base_ls075, binary_ls075) or this
# run (the half-Cauchy cells).
a <- commandArgs(TRUE)
run_dir <- "dev/simstudy/amplitude-prior-20261015"
d <- readRDS(a[1])
one_rds <- function(dir) {
  f <- list.files(dir, "^simstudy-.*\\.rds$", full.names = TRUE)
  stopifnot(length(f) == 1)
  readRDS(f)$rows
}
prod <- readRDS(file.path("dev/simstudy/paper-study-20261015/results",
                          "simstudy-20261015-172834.rds"))$rows
cols <- names(prod)
ref <- rbind(prod[prod$scenario %in% c("base", "binary") & prod$replicate <= 10, ],
             one_rds("dev/simstudy/range-inside-grid-20261015/results")[cols],
             one_rds(file.path(run_dir, "results"))[cols])
ref <- ref[ref$replicate <= 10 & ref$block %in% c("B0", "B"), ]
key <- function(x) paste(x$scenario, x$replicate, x$block, x$element, sep = "|")
# Draws are stored chain after chain (niter x nchain, flattened), 2 chains.
rhat <- function(x) posterior::rhat(matrix(x, ncol = 2))
ess <- function(x) posterior::ess_bulk(matrix(x, ncol = 2))

out <- do.call(rbind, lapply(d, function(x) {
  rows <- x$rows[x$rows$block %in% c("B0", "B"), ]
  m <- match(key(rows), key(ref))
  stopifnot(!anyNA(m))
  grid <- x$l_s_grid
  data.frame(scenario = x$scenario, prior = x$prior, replicate = x$replicate,
             l_s_raw = x$l_s_raw, l_s_std = mean(x$l_s_std),
             fitted_l_s_mean = mean(x$l_s_draws),
             share_at_grid_max = mean(x$l_s_draws == max(grid)),
             share_at_grid_min = mean(x$l_s_draws == min(grid)),
             sigma_bs_mean = mean(x$sigma_bs_draws),
             sigma_bs_q025 = unname(quantile(x$sigma_bs_draws, 0.025)),
             sigma_bs_q975 = unname(quantile(x$sigma_bs_draws, 0.975)),
             sigma_bs_rhat = rhat(x$sigma_bs_draws), sigma_bs_ess = ess(x$sigma_bs_draws),
             l_s_rhat = rhat(x$l_s_draws), l_s_ess = ess(x$l_s_draws),
             spatial_share_true = mean(x$spatial_share_true),
             spatial_share_fit = mean(x$spatial_share_fit),
             field_sd_true = x$field_sd_true,
             max_abs_post_mean_diff_vs_run = max(abs(rows$post_mean - ref$post_mean[m])))
}))
order_cells <- c("base", "base_hc", "base_ls075", "base_ls075_hc",
                 "binary", "binary_hc", "binary_ls075", "binary_ls075_hc")
out <- out[order(match(out$scenario, order_cells), out$replicate), ]
rownames(out) <- NULL
write.csv(out, file.path(run_dir, "amp-diag.csv"), row.names = FALSE)

agg <- do.call(rbind, lapply(order_cells, function(sc) {
  x <- out[out$scenario == sc, ]
  data.frame(scenario = sc, true_std = round(mean(x$l_s_std), 2),
             sigma_bs = round(mean(x$sigma_bs_mean), 3),
             sigma_bs_rep_min = round(min(x$sigma_bs_mean), 3),
             sigma_bs_rep_max = round(max(x$sigma_bs_mean), 3),
             fitted_l_s = round(mean(x$fitted_l_s_mean), 3),
             at_max = round(mean(x$share_at_grid_max), 3),
             at_min = round(mean(x$share_at_grid_min), 3),
             sigma_bs_rhat_median = round(median(x$sigma_bs_rhat), 2),
             sigma_bs_ess_median = round(median(x$sigma_bs_ess), 1),
             l_s_rhat_median = round(median(x$l_s_rhat), 2),
             l_s_ess_median = round(median(x$l_s_ess), 1),
             share_true = round(mean(x$spatial_share_true), 3),
             share_fit = round(mean(x$spatial_share_fit), 3),
             max_diff_vs_run = signif(max(x$max_abs_post_mean_diff_vs_run), 2))
}))
print(agg, row.names = FALSE)
