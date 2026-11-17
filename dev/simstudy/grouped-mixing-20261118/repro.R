#!/usr/bin/env Rscript
# Synthetic reproduction of the grouped prior's slow mixing, for occJSDM's own debugging.
# Contains no dataset content: sites, covariates and species are all simulated here.
#
# Why this exists. The AlphaEco simulation run (17 Nov 2026, package 32847ce) found every grouped fit
# failing the convergence rule at 5,000 + 5,000: the log axis scales, the embedding leftover scale
# and 25 to 65 percent of the B entries over R-hat 1.1 in the null and low-rank truths. Two readings
# are open (see the README in this folder). This script measures the evidence for each.
#
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/repro.R \
#       [--truth=null|lowrank] [--n=256] [--S=48] [--K=64] [--p_meas=5] [--gt=4] \
#       [--nburn=1500] [--niter=1500] [--seed=1] [--lib=<study library>] [--out=<dir>]
#
# Run with the INSTALLED package (never devtools::load_all(), which compiles without optimisation),
# one worker per run, and only while the fitting-worker budget has room (five workers in total).
# --lib points at a library holding the package under test; the commit is recorded with the result.
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
truth   <- getopt("truth", "null")
n       <- as.integer(getopt("n", "256"));    S      <- as.integer(getopt("S", "48"))
K       <- as.integer(getopt("K", "64"));     p_meas <- as.integer(getopt("p_meas", "5"))
gt      <- as.integer(getopt("gt", "4"));     seed   <- as.integer(getopt("seed", "1"))
nburn   <- as.integer(getopt("nburn", "1500")); niter <- as.integer(getopt("niter", "1500"))
lib     <- getopt("lib", NA_character_)
out_dir <- getopt("out", "dev/simstudy/grouped-mixing-20261118/results")
stopifnot(truth %in% c("null", "lowrank"))
if (!is.na(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
suppressMessages(library(occJSDM))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag <- sprintf("%s_n%d_S%d_K%d_gt%d_b%d_s%d", truth, n, S, K, gt, nburn, seed)

# ---- simulated data: measured covariates, a correlated embedding block, species ------------
set.seed(seed)
datasettings <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n),
                     ncov_psi = p_meas, ncov_theta = 1L)
embedding <- list(n_emb = K)
if (truth == "null") embedding$B <- matrix(0, K, S) else embedding <- c(embedding, list(n_blends = 3L, n_active = min(10L, K)))
jsdmParams <- list(gt = gt, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                   sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, S),
                   useSpatField = FALSE, embedding = embedding)
sim <- simulateOccJSDMData(datasettings, list(), jsdmParams, model = "binary")
dat <- sim$data_list
occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars), dat$embeddingVars)

# ---- fit: grouped prior, package defaults, one worker --------------------------------------
t0 <- Sys.time()
fit <- suppressMessages(suppressWarnings(runOccJSDM(
  dat, listParams = list(n_factors = 2L, n_lattrait = gt, embeddingVars = dat$embeddingVars),
  occCovariates = occ,
  MCMCparams = list(nchain = 2, nburn = nburn, niter = niter, nthin = 1, ncores = 1))))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
jo <- fit$results_output$jsdm_output

# ---- diagnostics ---------------------------------------------------------------------------
as_iter_chain <- function(x) { d <- dim(x); matrix(x, nrow = prod(d[-length(d)]), ncol = d[length(d)]) }
rhat <- function(x) as.vector(occJSDM:::computeRhat(x))
ess  <- function(x) as.vector(occJSDM:::computeESSparams(x))
lam  <- jo$lambda_output                       # K x iter x chain
log_lam <- log(lam)
sig  <- as_iter_chain(jo$sigmaemb_output)      # iter x chain
nit  <- dim(lam)[2]; nch <- dim(lam)[3]
true_B <- sim$true_params$jsdmParams_true$B
active <- if (truth == "null") rep(FALSE, K) else {
  emb_rows <- which(sim$true_params$jsdmParams_true$covariate_group == "embedding")
  rowSums(true_B[emb_rows, , drop = FALSE]^2) > 0
}

# Autocorrelation time = kept draws per effective draw (both chains pooled).
tau <- function(e) (nit * nch) / e
rh_ll <- rhat(log_lam); es_ll <- ess(log_lam)
# H3 check: the product sigma_emb * lambda_k is the prior sd of the leftovers, and what the
# leftover entries identify. If its ESS is much higher than that of lambda_k alone, the two scales
# trade off along a ridge.
prod_arr <- array(NA_real_, dim(lam))
for (k in seq_len(K)) prod_arr[k, , ] <- lam[k, , ] * sig
rh_pr <- rhat(log(prod_arr)); es_pr <- ess(log(prod_arr))
mean_log_lam <- apply(log_lam, c(2, 3), mean)  # iter x chain
cor_sig_lam <- cor(as.vector(log(sig)), as.vector(mean_log_lam))
B <- jo$B_output                               # p x S x iter x chain
emb_rows <- seq(p_meas + 1L, dim(B)[1])   # measured columns come first, then the embedding
rh_B <- rhat(B[emb_rows, , , , drop = FALSE]); es_B <- ess(B[emb_rows, , , , drop = FALSE])

row <- data.frame(tag = tag, truth = truth, n = n, S = S, K = K, gt = gt, nburn = nburn, niter = niter,
                  seed = seed, seconds = round(secs), package_commit = tryCatch(
                    readLines(file.path(dirname(getNamespaceInfo("occJSDM", "path")), "COMMIT"), 1),
                    error = function(e) NA_character_),
                  n_active = sum(active),
                  loglam_rhat_gt1.1 = sum(rh_ll > 1.1, na.rm = TRUE), loglam_rhat_max = max(rh_ll, na.rm = TRUE),
                  loglam_tau_med_inactive = median(tau(es_ll[!active])),
                  loglam_tau_med_active = if (any(active)) median(tau(es_ll[active])) else NA_real_,
                  prod_tau_med = median(tau(es_pr)), prod_rhat_gt1.1 = sum(rh_pr > 1.1, na.rm = TRUE),
                  logsig_rhat = rhat(log(sig)), logsig_tau = tau(ess(log(sig))),
                  cor_logsig_meanloglam = cor_sig_lam,
                  B_emb_share_rhat_gt1.1 = mean(rh_B > 1.1, na.rm = TRUE), B_emb_tau_med = median(tau(es_B)))
utils::write.csv(row, file.path(out_dir, paste0(tag, ".csv")), row.names = FALSE)
cols <- data.frame(tag = tag, k = seq_len(K), active = active, rhat_loglam = rh_ll, tau_loglam = tau(es_ll),
                   rhat_logprod = rh_pr, tau_logprod = tau(es_pr),
                   median_lambda = apply(lam, 1, median))
utils::write.csv(cols, file.path(out_dir, paste0(tag, "_columns.csv")), row.names = FALSE)
# Draws of a few axis scales and sigma_emb, thinned, for looking at traces and autocorrelation.
keep <- seq(1L, nit, by = 5L)
saveRDS(list(log_lambda = log_lam[sample(K, min(K, 8L)), keep, , drop = FALSE], log_sigma_emb = log(sig[keep, , drop = FALSE]),
             mean_log_lambda = mean_log_lam[keep, , drop = FALSE]), file.path(out_dir, paste0(tag, "_traces.rds")))
print(t(row), quote = FALSE)
