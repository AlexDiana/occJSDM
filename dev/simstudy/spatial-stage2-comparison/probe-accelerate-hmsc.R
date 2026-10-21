# Re-check of the stage 2 Hmsc arms' per-iteration cost under Apple's Accelerate BLAS, 14 October 2026 (SHELVED.md section 9).
# Not a study attempt. Runs under system R 4.6.1 with Accelerate, Hmsc 3.3-7 from the R 4.6.1 library of the Hmsc longer-runs
# check (laptop path below). Full method, probit, one covariate, nf spatial factors fixed, one chain, as upstream/adapter-hmsc.R
# calls it. Usage:  Rscript probe-accelerate-hmsc.R <sites> <factors> <iterations...>   (a third of each schedule is warm-up)
a <- commandArgs(trailingOnly = TRUE); n <- as.integer(a[1]); nf <- as.integer(a[2]); iters <- as.integer(a[-(1:2)])
.libPaths(c("/Users/douglasyu/src/occJSDM/dev/simstudy/results/hmsc-longer-runs-20261013/library", .libPaths()))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
suppressMessages({ library(occJSDM); library(Hmsc) })
cat("BLAS:", extSoftVersion()[["BLAS"]], " Hmsc:", as.character(packageVersion("Hmsc")), "from", dirname(system.file(package = "Hmsc")), "\n")
set.seed(5162); S <- 12; M <- rep(1, n); P <- 1; K <- rep(1, n)
sim <- simulateOccJSDMData(list(n = n, S = S, g = 3, M = M, P = P, K = K, ncov_psi = 2, ncov_theta = 2),
  list(p = matrix(0.5, P, S), q = matrix(0.03, P, S), theta0 = rep(0.05, S), theta_baseline = rep(0.3, S)),
  list(gt = 2, d = 2, ds = 2, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5, sigma_h = 1, sigma_s = 1, l_s = 0.08, useSpatField = TRUE),
  model = "binary")
data <- sim$data_list; Y <- as.matrix(data$OTU); storage.mode(Y) <- "numeric"
ids <- sprintf("s%04d", seq_len(n)); rownames(Y) <- ids
xy <- cbind(data$info$Xs.1, data$info$Xs.2); rownames(xy) <- ids
X <- data.frame(environment = as.numeric(scale(data$info$X_psi.EnvCov.1)))
level <- setPriors(HmscRandomLevel(sData = xy, sMethod = "Full"), nfMin = nf, nfMax = nf)
m <- Hmsc(Y = Y, XData = X, XFormula = ~., XScale = FALSE, distr = "probit",
          studyDesign = data.frame(site = factor(ids, levels = ids)), ranLevels = list(site = level))
secs <- c()
for (it in iters) {
  set.seed(1)
  t <- system.time(sampleMcmc(m, samples = it - it %/% 3, transient = it %/% 3, thin = 1L, nChains = 1L, nParallel = 1L, verbose = 0))[["elapsed"]]
  secs <- c(secs, t); cat(sprintf("Hmsc Full, %d sites, %d factors, iterations %d: %.1f s\n", n, nf, it, t))
}
f <- lm(secs ~ iters)
cat(sprintf("Hmsc Full, %d sites, %d factors: per-iteration cost %.3f s, fixed cost about %.1f s\n", n, nf, coef(f)[["iters"]], coef(f)[["(Intercept)"]]))
