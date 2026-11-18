# Is the share of embedding B entries over R-hat 1.1 a property of the sampler or of the diagnostic?
# Refits one repro.R configuration and compares three R-hats on the same draws: coda's
# gelman.diag point estimate (what the package and repro.R use, which carries a degrees-of-freedom
# correction that grows when the draws are heavy-tailed), the plain sqrt(V / W) without that
# correction, and the rank-normalised split R-hat (the larger of bulk and folded) of Vehtari et al. 2021.
#   Rscript rhat_check.R <lib> <truth> <seed> <nburn> <niter> <out.rds>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(occJSDM))
truth <- a[2]; seed <- as.integer(a[3]); nburn <- as.integer(a[4]); niter <- as.integer(a[5])
n <- 256L; S <- 48L; K <- 64L; gt <- 4L
set.seed(seed)
ds <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n), ncov_psi = 5L, ncov_theta = 1L)
embedding <- list(n_emb = K)
if (truth == "null") embedding$B <- matrix(0, K, S) else embedding <- c(embedding, list(n_blends = 3L, n_active = 10L))
jp <- list(gt = gt, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5, sigma_h = 1, sigma_s = 0.5, l_s = 0.3,
           tau = rep(1, S), useSpatField = FALSE, embedding = embedding)
sim <- simulateOccJSDMData(ds, list(), jp, model = "binary"); dat <- sim$data_list
occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars), dat$embeddingVars)
fit <- suppressMessages(suppressWarnings(runOccJSDM(dat, listParams = list(n_factors = 2L, n_lattrait = gt, embeddingVars = dat$embeddingVars),
  occCovariates = occ, MCMCparams = list(nchain = 2, nburn = nburn, niter = niter, nthin = 1, ncores = 1))))
jo <- fit$results_output$jsdm_output
emb <- unname(which(fit$infos$covariate_group == "embedding"))
plain <- function(x) { nn <- nrow(x); W <- mean(apply(x, 2, var)); B <- nn * var(colMeans(x)); sqrt(((nn - 1) / nn * W + B / nn) / W) }
split2 <- function(x) { h <- floor(nrow(x) / 2); cbind(x[seq_len(h), , drop = FALSE], x[h + seq_len(h), , drop = FALSE]) }
znorm <- function(x) matrix(qnorm((rank(x) - 3 / 8) / (length(x) + 1 / 4)), nrow(x))
rank_rhat <- function(x) { x <- split2(x); max(plain(znorm(x)), plain(znorm(abs(x - median(x))))) }
coda_rhat <- function(x) coda::gelman.diag(coda::mcmc.list(lapply(seq_len(ncol(x)), function(i) coda::mcmc(x[, i]))), autoburnin = FALSE)$psrf[1, 1]
three <- function(x) c(coda = coda_rhat(x), plain = plain(x), rank = rank_rhat(x))
kurt <- function(v) mean((v - mean(v))^4) / var(v)^2
B <- jo$B_output
rB <- t(sapply(seq_len(length(emb) * S), function(i) { k <- emb[(i - 1) %% length(emb) + 1]; s <- (i - 1) %/% length(emb) + 1; x <- B[k, s, , ]; c(three(x), kurtosis = kurt(as.vector(x))) }))
meas <- setdiff(seq_len(dim(B)[1]), emb)
rM <- t(sapply(seq_len(length(meas) * S), function(i) { k <- meas[(i - 1) %% length(meas) + 1]; s <- (i - 1) %/% length(meas) + 1; three(B[k, s, , ]) }))
rL <- t(apply(log(jo$lambda_output), 1, three))
spread <- array(NA_real_, c(S, dim(B)[3], dim(B)[4])); X <- fit$X_psi
for (ch in 1:2) for (it in seq_len(dim(B)[3])) spread[, it, ch] <- apply(X[, emb, drop = FALSE] %*% B[emb, , it, ch], 2, stats::sd)
rS <- t(apply(spread, 1, three))
out <- list(truth = truth, seed = seed, nburn = nburn, commit = readLines(file.path(a[1], "COMMIT"), 1),
            B_emb = rB, B_meas = rM, log_lambda = rL, spread = rS,
            log_sigma_emb = three(log(jo$sigmaemb_output)), log_sigma_b = three(log(jo$sigmab_output)))
saveRDS(out, a[6])
sh <- function(m) round(colMeans(m[, 1:3, drop = FALSE] > 1.1), 3)
cat(sprintf("%s seed %d %d+%d (%s)\n", truth, seed, nburn, niter, out$commit)); cat("share over 1.1 (coda, plain, rank-normalised split):\n")
print(rbind(B_embedding = sh(rB), B_measured = sh(rM), log_axis_scales = sh(rL), embedding_spread = sh(rS)))
cat("median kurtosis of an embedding B entry's draws (3 = Normal):", round(median(rB[, "kurtosis"]), 1), "\n")
cat("log sigma_emb:", round(out$log_sigma_emb, 3), " log sigma_b:", round(out$log_sigma_b, 3), "\n")
