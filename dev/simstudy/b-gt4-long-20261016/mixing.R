# Mixing of B and B0 in no_spatial_gt4_long fits, against the candidates
# run's mixing figures for no_spatial, no_spatial_long and no_spatial_gt4.
# Run from the repository root:
#   Rscript dev/simstudy/b-gt4-long-20261016/mixing.R
# Refits replicates 1 to 5 of no_spatial_gt4_long in this process, keeping the
# draws, and checks each refit reproduces the run's posterior means exactly
# before reporting its bulk ESS, tail ESS and R-hat. Also re-simulates all 100
# replicates to check the cell fits exactly no_spatial's datasets.
# Writes mixing.csv and dataset-pairing.csv.
here <- "dev/simstudy/b-gt4-long-20261016"
prev_mix <- "dev/simstudy/b-nonspatial-candidates-20261015/mixing.csv"
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-simstudy.R")
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1"); RcppParallel::setThreadOptions(numThreads = 1)
# The runner's workers run L'Ecuyer streams (clusterSetRNGStream); under that
# generator set.seed() in draw_truth() gives the same datasets and chain seeds.
RNGkind("L'Ecuyer-CMRG")
run <- readRDS(list.files(file.path(here, "results"), "^simstudy-.*[.]rds$", full.names = TRUE)[1])$rows
sc <- simstudy_scenarios(); names(sc) <- vapply(sc, `[[`, "", "label")
reps <- as.integer(Sys.getenv("MIX_REPS", "5"))

one <- function(lab, r) {
  s <- sc[[lab]]
  mc <- s$MCMCparams %||% list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1)
  truth <- draw_truth(s, simstudy_seed_for(s, r)); sim <- simstudy_simulate(truth)
  fit <- simstudy_fit(sim, s, mc)
  rows <- simstudy_rows(fit, sim, truth, s, r)
  ref <- run[run$scenario == lab & run$replicate == r, ]
  m <- merge(rows, ref, by = c("block", "element"))
  same <- nrow(m) == nrow(ref) && all(m$post_mean.x == m$post_mean.y)
  jo <- fit$results_output$jsdm_output
  do.call(rbind, lapply(c("B", "B0"), function(b) {
    arr <- if (b == "B") jo$B_output else jo$B0_output
    dm <- posterior_draw_matrix(arr)   # elements x (niter * nchain), chains in blocks
    nit <- dim(arr)[length(dim(arr)) - 1]; nch <- dim(arr)[length(dim(arr))]
    st <- t(apply(dm, 1, function(v) {
      x <- matrix(v, nit, nch)
      c(ess_bulk = posterior::ess_bulk(x), ess_tail = posterior::ess_tail(x),
        rhat = posterior::rhat(x))
    }))
    data.frame(arm = lab, replicate = r, block = b, element = seq_len(nrow(dm)),
               draws = nit * nch, reproduces_run = same, st)
  }))
}
# Datasets, all 100 replicates, simulation only: the cell should hand the
# fitter exactly no_spatial's data.
simdat <- function(lab, r) {
  s <- sc[[lab]]; simstudy_fit_args(simstudy_simulate(draw_truth(s, simstudy_seed_for(s, r))), s)$data
}
dp <- data.frame(replicate = 1:100, gt4_long = vapply(1:100, function(r)
  identical(simdat("no_spatial_gt4_long", r), simdat("no_spatial", r)), logical(1)))
cat("no_spatial_gt4_long fits no_spatial's datasets at all 100 replicates:", all(dp$gt4_long), "\n")
write.csv(dp, file.path(here, "dataset-pairing.csv"), row.names = FALSE)

mx <- do.call(rbind, lapply(seq_len(reps), function(r) one("no_spatial_gt4_long", r)))
write.csv(mx, file.path(here, "mixing.csv"), row.names = FALSE)
cat("refits reproduce the run:", all(mx$reproduces_run), "\n")

# The candidates run's figures, at the same replicates and at all it refitted.
old <- read.csv(prev_mix)
q <- function(v) sprintf("%.0f / %.0f / %.0f", quantile(v, .1), median(v), quantile(v, .9))
show <- function(x, tag) for (a in unique(x$arm)) for (b in c("B", "B0")) {
  y <- x[x$arm == a & x$block == b, ]
  cat(sprintf("%-20s %-3s %-12s draws %5d  bulk ESS p10/med/p90 %s  tail ESS %s  R-hat > 1.01: %.2f  max %.3f\n",
              a, b, tag, y$draws[1], q(y$ess_bulk), q(y$ess_tail), mean(y$rhat > 1.01), max(y$rhat)))
}
show(mx, sprintf("reps 1-%d", reps))
show(old[old$replicate <= reps, ], sprintf("reps 1-%d", reps))
show(old, "all refitted")
