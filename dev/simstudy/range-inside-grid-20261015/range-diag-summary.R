# Summarise range-diag.R's output into range-diag.csv.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/range-inside-grid-20261015/range-diag-summary.R \
#     dev/simstudy/range-inside-grid-20261015/results/range-diag.rds
#
# Also checks that each refit is the same fit as the run's replicate: B and B0
# posterior means must match the R = 300 production rows (parents) or this
# run's rows (diagnostic cells).
a <- commandArgs(TRUE)
run_dir <- "dev/simstudy/range-inside-grid-20261015"
d <- readRDS(a[1])
prod <- readRDS(file.path("dev/simstudy/paper-study-20261015/results",
                          "simstudy-20261015-172834.rds"))$rows
runr <- readRDS(list.files(file.path(run_dir, "results"), "^simstudy-.*\\.rds$",
                           full.names = TRUE))$rows
ref <- rbind(prod[prod$scenario %in% c("base", "binary") & prod$replicate <= 10, ],
             runr[runr$replicate <= 10, names(prod)])
key <- function(x) paste(x$scenario, x$replicate, x$block, x$element, sep = "|")
# Draws are stored chain after chain (niter x nchain, flattened), 2 chains.
rhat <- function(x) posterior::rhat(matrix(x, ncol = 2))
ess <- function(x) posterior::ess_bulk(matrix(x, ncol = 2))

out <- do.call(rbind, lapply(d, function(x) {
  rows <- x$rows[x$rows$block %in% c("B0", "B"), ]
  m <- match(key(rows), key(ref))
  stopifnot(!anyNA(m))
  grid <- x$l_s_grid
  data.frame(scenario = x$scenario, replicate = x$replicate,
             l_s_raw = x$l_s_raw, l_s_std = mean(x$l_s_std),
             fitted_l_s_mean = mean(x$l_s_draws),
             fitted_l_s_median = median(x$l_s_draws),
             share_at_grid_max = mean(x$l_s_draws == max(grid)),
             share_at_grid_min = mean(x$l_s_draws == min(grid)),
             share_within_one_step = mean(abs(x$l_s_draws - mean(x$l_s_std)) <=
                                            diff(grid)[1] + 1e-12),
             l_s_rhat = rhat(x$l_s_draws), l_s_ess = ess(x$l_s_draws),
             sigma_bs_mean = mean(x$sigma_bs_draws),
             sigma_bs_rhat = rhat(x$sigma_bs_draws), sigma_bs_ess = ess(x$sigma_bs_draws),
             max_abs_post_mean_diff_vs_run = max(abs(rows$post_mean - ref$post_mean[m])))
}))
out <- out[order(match(out$scenario, c("base", "base_ls075", "base_ls05", "binary",
                                       "binary_ls075", "binary_ls05")), out$replicate), ]
rownames(out) <- NULL
write.csv(out, file.path(run_dir, "range-diag.csv"), row.names = FALSE)

agg <- do.call(rbind, lapply(split(out, out$scenario), function(x) data.frame(
  scenario = x$scenario[1], true_std = round(mean(x$l_s_std), 3),
  fitted_mean = round(mean(x$fitted_l_s_mean), 3),
  fitted_rep_min = round(min(x$fitted_l_s_mean), 3),
  fitted_rep_max = round(max(x$fitted_l_s_mean), 3),
  at_max = round(mean(x$share_at_grid_max), 3),
  at_min = round(mean(x$share_at_grid_min), 3),
  within_step = round(mean(x$share_within_one_step), 3),
  l_s_rhat_median = round(median(x$l_s_rhat), 2),
  l_s_ess_median = round(median(x$l_s_ess), 1),
  sigma_bs = round(mean(x$sigma_bs_mean), 3),
  sigma_bs_rhat_median = round(median(x$sigma_bs_rhat), 2),
  max_diff_vs_run = signif(max(x$max_abs_post_mean_diff_vs_run), 2))))
print(agg[order(match(agg$scenario, unique(out$scenario))), ], row.names = FALSE)
