# Re-check of the stage 2 gp arm's per-iteration cost under Apple's Accelerate BLAS, 14 October 2026 (SHELVED.md section 9).
# Not a study attempt. Runs under system R 4.6.1 with Accelerate (not run-r, whose pinned R 4.5.0 loads the reference BLAS)
# and the installed occJSDM. Usage:  Rscript probe-accelerate-gp.R <support_points> [sites]   (sites defaults to 900)
# Mirrors dev/profiling/2026-10-14-cpu-profile/make-binary900.R's community (12 species, binary, range 0.08) and the stage 2
# adapter's call shape (n_factors 0, n_lattrait 0, one chain) at the ladder's support-point settings. Three schedules per setting
# separate the fixed cost from the per-iteration cost by least squares, as results/probe-gp-slope.txt did with two.
suppressMessages(library(occJSDM))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n"); cat("occJSDM:", as.character(packageVersion("occJSDM")), "\n")
set.seed(5162)
n <- { a <- commandArgs(trailingOnly = TRUE); if (length(a) >= 2) as.integer(a[2]) else 900L }; S <- 12; g <- 3; gt <- 2; d <- 2; ds <- 2
M <- rep(1, n); P <- 1; N <- sum(M); K <- rep(1, P * N)
sim <- simulateOccJSDMData(list(n = n, S = S, g = g, M = M, P = P, K = K, ncov_psi = 2, ncov_theta = 2),
  list(p = matrix(0.5, P, S), q = matrix(0.03, P, S), theta0 = rep(0.05, S), theta_baseline = rep(0.3, S)),
  list(gt = gt, d = d, ds = ds, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5, sigma_h = 1, sigma_s = 1, l_s = 0.08, useSpatField = TRUE),
  model = "binary")
data <- sim$data_list
sp <- as.integer(commandArgs(trailingOnly = TRUE)[1])
scheds <- list(c(nburn = 5L, niter = 10L), c(nburn = 10L, niter = 20L), c(nburn = 20L, niter = 40L))
secs <- c(); its <- c()
for (sched in scheds) {
  set.seed(20261014)
  t <- system.time(fit <- suppressMessages(runOccJSDM(data = data,
    listParams = list(n_factors = 0L, n_lattrait = 0L, n_supportpoints = sp),
    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"), spatCovariates = c("Xs.1", "Xs.2"),
    MCMCparams = list(nchain = 1, nburn = sched[["nburn"]], niter = sched[["niter"]], nthin = 1, ncores = 1),
    summarisedLatentPresences = TRUE)))[["elapsed"]]
  stopifnot(nrow(fit$infos$list_Xs$X_tilde) == sp)
  secs <- c(secs, t); its <- c(its, sum(sched))
  cat(sprintf("support points %d, iterations %d: %.1f s whole fit\n", sp, sum(sched), t))
}
fit_ls <- lm(secs ~ its)
cat(sprintf("support points %d: per-iteration cost %.3f s (least squares over three schedules), fixed cost about %.1f s\n",
  sp, coef(fit_ls)[["its"]], coef(fit_ls)[["(Intercept)"]]))
