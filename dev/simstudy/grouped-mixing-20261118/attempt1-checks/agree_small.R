# A same-answer check with power: a small low-rank problem where both arms can converge.
# One simulated dataset, 4 chains of 3,000 + 12,000 per arm. usage: Rscript agree_small.R <lib> <out.rds>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
set.seed(11); n <- 150L; S <- 16L; K <- 8L
ds <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n), ncov_psi = 2L, ncov_theta = 1L)
jp <- list(gt = 2L, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5, sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
           tau = rep(1, S), useSpatField = FALSE, embedding = list(n_emb = K, n_blends = 2L, n_active = 3L))
sim <- simulateOccJSDMData(ds, list(), jp, model = "binary"); dat <- sim$data_list
occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars), dat$embeddingVars)
set.seed(if (length(a) > 2) as.integer(a[3]) else 12L)
t0 <- proc.time()
fit <- suppressMessages(suppressWarnings(runOccJSDM(dat, listParams = list(n_factors = 1L, n_lattrait = 2L, embeddingVars = dat$embeddingVars),
  occCovariates = occ, MCMCparams = list(nchain = 4, nburn = 3000, niter = 12000, nthin = 1, ncores = 1))))
el <- (proc.time() - t0)
jo <- fit$results_output$jsdm_output
emb <- unname(which(fit$infos$covariate_group == "embedding")); B <- jo$B_output; X <- fit$X_psi
summ <- function(x) { ess <- as.vector(occJSDM:::computeESSparams(x)); rh <- as.vector(occJSDM:::computeRhat(x))
  d <- dim(x); m <- matrix(x, prod(d[seq_len(length(d) - 2)]), d[length(d) - 1] * d[length(d)])
  data.frame(mean = rowMeans(m), sd = apply(m, 1, stats::sd), ess = ess, rhat = rh,
             q10 = apply(m, 1, quantile, .1), q50 = apply(m, 1, quantile, .5), q90 = apply(m, 1, quantile, .9)) }
spread <- array(NA_real_, c(S, dim(B)[3], dim(B)[4]))
for (ch in seq_len(dim(B)[4])) for (it in seq_len(dim(B)[3])) spread[, it, ch] <- apply(X[, emb, drop = FALSE] %*% B[emb, , it, ch], 2, stats::sd)
saveRDS(list(commit = readLines(file.path(a[1], "COMMIT"), 1), seconds = el[["elapsed"]], cpu = el[["user.self"]] + el[["sys.self"]],
  B_meas = summ(B[-emb, , , , drop = FALSE]), B_emb = summ(B[emb, , , , drop = FALSE]), spread = summ(spread),
  log_lambda = summ(log(jo$lambda_output)), log_sigma_emb = summ(log(jo$sigmaemb_output)), log_sigma_b = summ(log(jo$sigmab_output)),
  B0 = summ(jo$B0_output), axis_move = fit$infos$axis_move), a[2])
