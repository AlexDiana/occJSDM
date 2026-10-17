# Mixing of B, B0 and sigma_b in the five sigma_b prior arms and, for sigma_b
# under the default prior, in no_spatial and no_spatial_long.
# Run from the repository root, after the study run has finished (it uses 5
# worker processes, the laptop's whole fitting budget):
#   Rscript dev/simstudy/sigma-b-prior-20261017/mixing.R
# Refits replicates 1 to 5 of each arm with the draws kept, checks each refit
# reproduces the run's posterior means exactly before reporting its bulk ESS,
# tail ESS and R-hat, and re-simulates all 100 replicates to check every new
# cell fits exactly no_spatial's datasets.
# Writes mixing.csv and dataset-pairing.csv.
here <- "dev/simstudy/sigma-b-prior-20261017"
cand_file <- "dev/simstudy/b-nonspatial-candidates-20261015/results/simstudy-20261016-004116.rds"
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-simstudy.R")
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1"); RcppParallel::setThreadOptions(numThreads = 1)
# The runner's workers run L'Ecuyer streams (clusterSetRNGStream); under that
# generator set.seed() in draw_truth() gives the same datasets and chain seeds.
RNGkind("L'Ecuyer-CMRG")
run <- readRDS(list.files(file.path(here, "results"), "^simstudy-.*[.]rds$", full.names = TRUE)[1])$rows
cand <- readRDS(cand_file)$rows
ref_rows <- rbind(run, cand[cand$scenario %in% c("no_spatial", "no_spatial_long"), ])
sc <- simstudy_scenarios(); names(sc) <- vapply(sc, `[[`, "", "label")
reps <- as.integer(Sys.getenv("MIX_REPS", "5"))
new_arms <- c("no_spatial_sbhc", "no_spatial_sbig2", "no_spatial_sbhc_long",
              "no_spatial_sbig2_long", "no_spatial_sbhc_gt4_long")
mix_arms <- c(new_arms, "no_spatial", "no_spatial_long")

# Datasets, all 100 replicates, simulation only: every new cell should hand
# the fitter exactly no_spatial's data.
simdat <- function(lab, r) {
  s <- sc[[lab]]; simstudy_fit_args(simstudy_simulate(draw_truth(s, simstudy_seed_for(s, r))), s)$data
}
dp <- data.frame(replicate = 1:100)
for (a in new_arms) dp[[a]] <- vapply(1:100, function(r)
  identical(simdat(a, r), simdat("no_spatial", r)), logical(1))
cat("every new cell fits no_spatial's datasets at all 100 replicates:", all(as.matrix(dp[new_arms])), "\n")
write.csv(dp, file.path(here, "dataset-pairing.csv"), row.names = FALSE)

one <- function(lab, r) {
  s <- sc[[lab]]
  mc <- s$MCMCparams %||% list(nchain = 2, nburn = 1000, niter = 1000, nthin = 1)
  truth <- draw_truth(s, simstudy_seed_for(s, r)); sim <- simstudy_simulate(truth)
  fit <- simstudy_fit(sim, s, mc)
  rows <- simstudy_rows(fit, sim, truth, s, r)
  ref <- ref_rows[ref_rows$scenario == lab & ref_rows$replicate == r, ]
  m <- merge(rows, ref, by = c("block", "element"))
  same <- nrow(m) == nrow(ref) && all(m$post_mean.x == m$post_mean.y)
  jo <- fit$results_output$jsdm_output
  stat <- function(x) c(ess_bulk = posterior::ess_bulk(x), ess_tail = posterior::ess_tail(x),
                        rhat = posterior::rhat(x))
  out <- do.call(rbind, lapply(c("B", "B0"), function(b) {
    arr <- if (b == "B") jo$B_output else jo$B0_output
    dm <- posterior_draw_matrix(arr)   # elements x (niter * nchain), chains in blocks
    nit <- dim(arr)[length(dim(arr)) - 1]; nch <- dim(arr)[length(dim(arr))]
    st <- t(apply(dm, 1, function(v) stat(matrix(v, nit, nch))))
    data.frame(arm = lab, replicate = r, block = b, element = seq_len(nrow(dm)),
               draws = nit * nch, reproduces_run = same, st)
  }))
  out$post_mean <- NA_real_; out$post_median <- NA_real_; out$q025 <- NA_real_; out$q975 <- NA_real_
  sbo <- jo$sigmab_output   # niter x nchain; its posterior is heavy-tailed under
  # the half-Cauchy, so the median and quantiles are kept beside the mean
  rbind(out, data.frame(arm = lab, replicate = r, block = "sigma_b", element = 1L,
                        draws = length(sbo), reproduces_run = same,
                        t(stat(matrix(sbo, nrow(sbo), ncol(sbo)))),
                        post_mean = mean(sbo), post_median = median(sbo),
                        q025 = unname(quantile(sbo, .025)), q975 = unname(quantile(sbo, .975))))
}

jobs <- expand.grid(replicate = seq_len(reps), arm = mix_arms, stringsAsFactors = FALSE)
n_workers <- min(5L, nrow(jobs))
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("one", "sc", "ref_rows", "simdat"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1"); RcppParallel::setThreadOptions(numThreads = 1)
  suppressMessages(devtools::load_all(".", quiet = TRUE))
  source("tests/testthat/helper-simstudy.R")
  RNGkind("L'Ecuyer-CMRG")
  NULL
}))
t0 <- Sys.time()
res <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(k, jobs) one(jobs$arm[k], jobs$replicate[k]),
                             jobs = jobs)
mx <- do.call(rbind, res)
cat(sprintf("%d refits in %.1f min on %d workers\n", nrow(jobs),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_workers))
write.csv(mx, file.path(here, "mixing.csv"), row.names = FALSE)
cat("refits reproduce the run:", all(mx$reproduces_run), "\n")
if (!all(mx$reproduces_run)) print(unique(mx[!mx$reproduces_run, c("arm", "replicate")]))

q <- function(v) sprintf("%.0f / %.0f / %.0f", quantile(v, .1), median(v), quantile(v, .9))
for (a in mix_arms) for (b in c("B", "B0", "sigma_b")) {
  y <- mx[mx$arm == a & mx$block == b, ]
  if (b == "sigma_b") {
    cat(sprintf("%-26s %-7s draws %5d  bulk ESS per fit %s  tail ESS %s  R-hat max %.3f\n",
                a, b, y$draws[1], paste(round(y$ess_bulk), collapse = " "),
                paste(round(y$ess_tail), collapse = " "), max(y$rhat)))
    cat(sprintf("%-26s         posterior mean per fit %s  median %s  2.5%% %s  97.5%% %s\n", a,
                paste(sprintf("%.2f", y$post_mean), collapse = " "), paste(sprintf("%.2f", y$post_median), collapse = " "),
                paste(sprintf("%.2f", y$q025), collapse = " "), paste(sprintf("%.2f", y$q975), collapse = " ")))
  } else {
    cat(sprintf("%-26s %-7s draws %5d  bulk ESS p10/med/p90 %s  tail ESS %s  R-hat > 1.01: %.2f  max %.3f\n",
                a, b, y$draws[1], q(y$ess_bulk), q(y$ess_tail), mean(y$rhat > 1.01), max(y$rhat)))
  }
}
