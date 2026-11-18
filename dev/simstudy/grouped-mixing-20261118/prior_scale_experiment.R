#!/usr/bin/env Rscript
# A settings-only experiment on the prior of sigma_emb. No package code is changed: it sets
# listPriors$sigma_emb_scale, which is also the chain's starting value of sigma_emb, so a result
# here does not separate a better prior from a better start.
#
# Question: in the low-rank truth chains settle in one of two configurations (README, "The two
# low-rank configurations") and none crosses. Does a smaller half-Cauchy scale on sigma_emb put
# every chain in the blend configuration, and what does it cost in the null truth and in a truth
# whose species respond independently?
#
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/prior_scale_experiment.R \
#       --truth=null|lowrank|unstructured --seed=1 --scale=1 [--nchain=4] [--nburn=1500] [--niter=1500] \
#       --lib=<library> [--out=<dir>]
#
# Data are simulated exactly as in repro.R (same seeds give the same null and low-rank data sets):
# 256 sites, 48 species, 64 embedding columns, 5 measured covariates, 4 latent traits. The
# "unstructured" truth has 10 active embedding columns whose coefficients are independent
# Normal(0, 0.9^2) across species (no community average, no blends), and zero elsewhere.
# All data are simulated here; no dataset content.
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
truth <- getopt("truth", "lowrank"); seed <- as.integer(getopt("seed", "1"))
scale <- as.numeric(getopt("scale", "1")); nchain <- as.integer(getopt("nchain", "4"))
nburn <- as.integer(getopt("nburn", "1500")); niter <- as.integer(getopt("niter", "1500"))
n <- as.integer(getopt("n", "256")); S <- as.integer(getopt("S", "48")); K <- as.integer(getopt("K", "64"))
lib <- getopt("lib", NA_character_)
out_dir <- getopt("out", "dev/simstudy/grouped-mixing-20261118/results/prior-scale")
stopifnot(truth %in% c("null", "lowrank", "unstructured"))
if (!is.na(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
suppressMessages(library(occJSDM))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gt <- 4L; p_meas <- 5L
tag <- sprintf("%s_s%d_scale%s_b%d", truth, seed, format(scale), nburn)

embedding <- list(n_emb = K)
if (truth == "null") embedding$B <- matrix(0, K, S)
if (truth == "lowrank") embedding <- c(embedding, list(n_blends = 3L, n_active = min(10L, K)))
if (truth == "unstructured") {
  set.seed(1000L + seed)
  B_un <- matrix(0, K, S); B_un[sample(K, min(10L, K)), ] <- rnorm(min(10L, K) * S, sd = 0.9)
  embedding$B <- B_un
}
set.seed(seed)
datasettings <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n),
                     ncov_psi = p_meas, ncov_theta = 1L)
jsdmParams <- list(gt = gt, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                   sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, S),
                   useSpatField = FALSE, embedding = embedding)
sim <- simulateOccJSDMData(datasettings, list(), jsdmParams, model = "binary")
dat <- sim$data_list
occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars), dat$embeddingVars)

t0 <- Sys.time()
fit <- suppressMessages(suppressWarnings(runOccJSDM(
  dat, listParams = list(n_factors = 2L, n_lattrait = gt, embeddingVars = dat$embeddingVars),
  listPriors = if (scale == 1) list() else list(sigma_emb_scale = scale),
  occCovariates = occ,
  MCMCparams = list(nchain = nchain, nburn = nburn, niter = niter, nthin = 1, ncores = 1))))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
jo <- fit$results_output$jsdm_output
stopifnot(isTRUE(all.equal(fit$infos$coef_prior$sigma_emb_scale, scale)))

# ---- diagnostics: coda's R-hat (the package's) and the rank-normalised split R-hat ----------
plain <- function(x) { nn <- nrow(x); W <- mean(apply(x, 2, var)); Bv <- nn * var(colMeans(x)); sqrt(((nn - 1) / nn * W + Bv / nn) / W) }
split2 <- function(x) { h <- floor(nrow(x) / 2); cbind(x[seq_len(h), , drop = FALSE], x[h + seq_len(h), , drop = FALSE]) }
znorm <- function(x) matrix(qnorm((rank(x) - 3 / 8) / (length(x) + 1 / 4)), nrow(x))
rank_rhat <- function(x) { x <- split2(x); max(plain(znorm(x)), plain(znorm(abs(x - median(x))))) }
coda_rhat <- function(x) coda::gelman.diag(coda::mcmc.list(lapply(seq_len(ncol(x)), function(i) coda::mcmc(x[, i]))), autoburnin = FALSE)$psrf[[1, 1]]
both <- function(x) c(coda = coda_rhat(x), rank = rank_rhat(x))
over <- function(m) colMeans(m > 1.1)

B <- jo$B_output                                   # p x S x iter x chain
emb <- unname(which(fit$infos$covariate_group == "embedding")); meas <- setdiff(seq_len(dim(B)[1]), emb)
stopifnot(length(emb) == K)
true_B <- sim$true_params$jsdmParams_true$B
active <- emb[rowSums(true_B[emb, , drop = FALSE]^2) > 0]; inactive <- setdiff(emb, active)
X <- fit$X_psi
entries <- function(rows) t(sapply(seq_len(length(rows) * S), function(i) both(B[rows[(i - 1) %% length(rows) + 1], (i - 1) %/% length(rows) + 1, , ])))
r_emb <- entries(emb); r_meas <- entries(meas)
r_lam <- t(apply(log(jo$lambda_output), 1, both))
spread <- array(NA_real_, c(S, dim(B)[3], dim(B)[4]))
for (ch in seq_len(nchain)) for (it in seq_len(dim(B)[3])) spread[, it, ch] <- apply(X[, emb, drop = FALSE] %*% B[emb, , it, ch], 2, stats::sd)
r_spread <- t(apply(spread, 1, both))
log_sig <- log(jo$sigmaemb_output)                 # iter x chain
r_sig <- both(log_sig)

rms <- function(x) sqrt(mean(x^2))
slope <- function(est, tru) if (all(tru == 0)) NA_real_ else sum(est * tru) / sum(tru^2)
per_chain <- do.call(rbind, lapply(seq_len(nchain), function(ch) {
  Bc <- apply(B[, , , ch], c(1, 2), mean)
  contribution <- X[, emb, drop = FALSE] %*% (Bc[emb, ] - true_B[emb, ])
  data.frame(tag = tag, truth = truth, seed = seed, scale = scale, chain = ch,
             log_sigma_emb_mean = mean(log_sig[, ch]), log_sigma_emb_min = min(log_sig[, ch]),
             log_sigma_emb_max = max(log_sig[, ch]),
             mean_log_axis_active = if (length(active)) mean(log(jo$lambda_output[match(active, emb), , ch])) else NA_real_,
             mean_log_axis_inactive = mean(log(jo$lambda_output[match(inactive, emb), , ch])),
             err_active = if (length(active)) rms(Bc[active, ] - true_B[active, ]) else NA_real_,
             slope_active = if (length(active)) slope(Bc[active, ], true_B[active, ]) else NA_real_,
             err_inactive = rms(Bc[inactive, ] - true_B[inactive, ]),
             err_measured = rms(Bc[meas, ] - true_B[meas, ]),
             err_contribution = rms(contribution))
}))
Bm <- apply(B, c(1, 2), mean)
fit_row <- data.frame(tag = tag, truth = truth, seed = seed, scale = scale, nchain = nchain, nburn = nburn, niter = niter,
  seconds = round(secs), package_commit = tryCatch(readLines(file.path(dirname(getNamespaceInfo("occJSDM", "path")), "COMMIT"), 1), error = function(e) NA_character_),
  logsig_rhat_coda = r_sig[["coda"]], logsig_rhat_rank = r_sig[["rank"]],
  logsig_chain_gap = diff(range(colMeans(log_sig))),
  axis_over_coda = sum(r_lam[, "coda"] > 1.1), axis_over_rank = sum(r_lam[, "rank"] > 1.1), axis_rhat_max_rank = max(r_lam[, "rank"]),
  B_emb_over_coda = over(r_emb)[["coda"]], B_emb_over_rank = over(r_emb)[["rank"]],
  B_meas_over_coda = over(r_meas)[["coda"]], B_meas_over_rank = over(r_meas)[["rank"]], B_meas_rhat_max_rank = max(r_meas[, "rank"]),
  spread_over_rank = over(r_spread)[["rank"]], spread_rhat_max_rank = max(r_spread[, "rank"]),
  err_active = if (length(active)) rms(Bm[active, ] - true_B[active, ]) else NA_real_,
  slope_active = if (length(active)) slope(Bm[active, ], true_B[active, ]) else NA_real_,
  err_inactive = rms(Bm[inactive, ] - true_B[inactive, ]), err_measured = rms(Bm[meas, ] - true_B[meas, ]),
  err_contribution = rms(X[, emb, drop = FALSE] %*% (Bm[emb, ] - true_B[emb, ])),
  truth_rms_active = if (length(active)) rms(true_B[active, ]) else NA_real_,
  truth_rms_contribution = rms(X[, emb, drop = FALSE] %*% true_B[emb, ]),
  sigma_emb_median = stats::median(jo$sigmaemb_output),
  sigma_move_acceptance = paste(round(fit$infos$axis_move$sigma_emb_acceptance, 2), collapse = " "),
  sigma_move_step = paste(round(fit$infos$axis_move$sigma_emb_step, 2), collapse = " "))
utils::write.csv(per_chain, file.path(out_dir, paste0(tag, "_chains.csv")), row.names = FALSE)
utils::write.csv(fit_row, file.path(out_dir, paste0(tag, ".csv")), row.names = FALSE)
num <- vapply(per_chain, is.numeric, logical(1)); per_chain[num] <- lapply(per_chain[num], round, 3)
print(per_chain[, -(1:4)], row.names = FALSE); print(t(fit_row), quote = FALSE)
