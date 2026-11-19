#!/usr/bin/env Rscript
# One full-length fit, and the pieces the grouped prior's convergence check is made from, kept so
# that the check can be applied to them afterwards. The fits of the rotation and orientation tests
# kept only each chain's mean effect, so the state and R-hat checks were never run on a full-length
# fit; this is for that. Data as in repro.R and prior_scale_experiment.R (256 sites, 48 species,
# 64 embedding columns). All data are simulated here; no dataset content.
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/checks/convergence_check_on_new_fits.R \
#       <library> <null|lowrank|unstructured> <seed> [nchain=4] [nburn=5000] [niter=5000]
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
truth <- a[2]; seed <- as.integer(a[3])
nchain <- if (length(a) > 3) as.integer(a[4]) else 4L
nburn <- if (length(a) > 4) as.integer(a[5]) else 5000L; niter <- if (length(a) > 5) as.integer(a[6]) else 5000L
n <- 256L; S <- 48L; K <- 64L; gt <- 4L
out_dir <- "dev/simstudy/grouped-mixing-20261118/results/convergence-check"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
embedding <- list(n_emb = K)
if (truth == "null") embedding$B <- matrix(0, K, S)
if (truth == "lowrank") embedding <- c(embedding, list(n_blends = 3L, n_active = 10L))
if (truth == "unstructured") { set.seed(1000L + seed); Bu <- matrix(0, K, S); Bu[sample(K, 10L), ] <- rnorm(10L * S, sd = 0.9); embedding$B <- Bu }
set.seed(seed)
ds <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n), ncov_psi = 5L, ncov_theta = 1L)
jp <- list(gt = gt, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5, sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
           tau = rep(1, S), useSpatField = FALSE, embedding = embedding)
sim <- simulateOccJSDMData(ds, list(), jp, model = "binary"); dat <- sim$data_list
occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars), dat$embeddingVars)
set.seed(100L + seed)
t0 <- Sys.time()
fit <- suppressMessages(suppressWarnings(runOccJSDM(dat, listParams = list(n_factors = 2L, n_lattrait = gt, embeddingVars = dat$embeddingVars),
  occCovariates = occ, MCMCparams = list(nchain = nchain, nburn = nburn, niter = niter, nthin = 1, ncores = 1))))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
jo <- fit$results_output$jsdm_output; is_emb <- fit$infos$coef_prior$is_emb; B <- jo$B_output
rhat_of <- function(arr, n_params) { b <- occJSDM:::blockRhatESS(arr, n_params, niter, nchain); ifelse(b$variable, b$rhat, NA_real_) }
true_B <- sim$true_params$jsdmParams_true$B
saveRDS(list(truth = truth, seed = seed, nchain = nchain, nburn = nburn, niter = niter, seconds = secs,
  commit = tryCatch(readLines(file.path(a[1], "COMMIT"), 1), error = function(e) NA_character_),
  sigma_emb_median = apply(matrix(jo$sigmaemb_output, niter, nchain), 2, stats::median),
  rhat_measured = rhat_of(B[!is_emb, , , , drop = FALSE], sum(!is_emb) * S),
  rhat_log_sigma_emb = rhat_of(array(log(jo$sigmaemb_output), c(1L, niter, nchain)), 1L),
  rhat_log_sigma_b = rhat_of(array(log(jo$sigmab_output), c(1L, niter, nchain)), 1L),
  rhat_log_axis = rhat_of(log(jo$lambda_output), K),
  B_emb_mean = vapply(seq_len(nchain), function(ch) rowMeans(B[is_emb, , , ch], dims = 2), matrix(0, K, S)),
  B_meas_mean = vapply(seq_len(nchain), function(ch) rowMeans(B[!is_emb, , , ch], dims = 2), matrix(0, sum(!is_emb), S)),
  X_emb = unname(fit$X_psi[, is_emb]), true_B_emb = true_B[is_emb, ], true_B_meas = true_B[!is_emb, ],
  as_built = returnGroupedConvergence(fit)),
  file.path(out_dir, sprintf("%s_s%d_b%d.rds", truth, seed, nburn)))
cat(sprintf("%s seed %d done in %.0f s\n", truth, seed, secs))
