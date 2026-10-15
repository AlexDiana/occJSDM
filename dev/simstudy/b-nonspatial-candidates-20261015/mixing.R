# Mixing of B and B0 in no_spatial fits at the study's settings and at 5,000 + 5,000,
# and in no_spatial_gt4 fits at the study's settings.
# Run from the repository root:
#   Rscript dev/simstudy/b-nonspatial-candidates-20261015/mixing.R
# Refits replicates 1 to 10 of no_spatial and no_spatial_gt4, and 1 to 5 of
# no_spatial_long, in this process, keeping the draws, and checks each refit reproduces the run's
# posterior means exactly before reporting its bulk ESS, tail ESS and R-hat.
# Also re-simulates all 100 replicates of each arm to check which share datasets.
# Writes mixing.csv and dataset-pairing.csv.
here <- "dev/simstudy/b-nonspatial-candidates-20261015"
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-simstudy.R")
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1"); RcppParallel::setThreadOptions(numThreads = 1)
# The runner's workers run L'Ecuyer streams (clusterSetRNGStream); under that
# generator set.seed() in draw_truth() gives the same datasets and chain seeds.
RNGkind("L'Ecuyer-CMRG")
run <- readRDS(list.files(file.path(here, "results"), "^simstudy-.*[.]rds$", full.names = TRUE)[1])$rows
sc <- simstudy_scenarios(); names(sc) <- vapply(sc, `[[`, "", "label")
reps <- as.integer(Sys.getenv("MIX_REPS", "10"))

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
# Datasets, all 100 replicates, simulation only: the fit-only arms should hand
# the fitter no_spatial's data, the d = 0 arms each other's.
simdat <- function(lab, r) {
  s <- sc[[lab]]; simstudy_fit_args(simstudy_simulate(draw_truth(s, simstudy_seed_for(s, r))), s)$data
}
dp <- do.call(rbind, lapply(1:100, function(r) {
  ns <- simdat("no_spatial", r); d0 <- simdat("no_factors", r)
  data.frame(replicate = r,
             long = identical(simdat("no_spatial_long", r), ns),
             gt4 = identical(simdat("no_spatial_gt4", r), ns),
             notraits = identical(simdat("no_spatial_notraits", r), ns[setdiff(names(ns), "traits")]),
             fit2 = identical(simdat("no_factors_fit2", r), d0),
             d0_differs = !identical(d0[c("OTU", "info")], ns[c("OTU", "info")]))
}))
cat("datasets as documented at all 100 replicates:",
    vapply(dp[-1], all, logical(1)), "\n")
write.csv(dp, file.path(here, "dataset-pairing.csv"), row.names = FALSE)

mx <- do.call(rbind, c(lapply(seq_len(reps), function(r) one("no_spatial", r)),
                       lapply(seq_len(reps), function(r) one("no_spatial_gt4", r)),
                       lapply(seq_len(min(reps, 5L)), function(r) one("no_spatial_long", r))))
write.csv(mx, file.path(here, "mixing.csv"), row.names = FALSE)
cat("refits reproduce the run:", all(mx$reproduces_run), "\n")
q <- function(v) sprintf("%.0f / %.0f / %.0f", quantile(v, .1), median(v), quantile(v, .9))
for (a in unique(mx$arm)) for (b in c("B", "B0")) {
  x <- mx[mx$arm == a & mx$block == b, ]
  cat(sprintf("%-16s %-3s draws %5d  bulk ESS p10/med/p90 %s  tail ESS %s  R-hat > 1.01: %.2f  max %.3f\n",
              a, b, x$draws[1], q(x$ess_bulk), q(x$ess_tail), mean(x$rhat > 1.01), max(x$rhat)))
}
