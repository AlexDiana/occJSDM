#!/usr/bin/env Rscript
# Orientation test of the grouped prior: does turning the embedding to its principal components
# before fitting recover what a random rotation loses? A copy of rotation_test.R (kept apart so
# that fits already running from that file are not disturbed) with one more option.
#
#   Rscript --vanilla dev/simstudy/grouped-mixing-20261118/orientation_test.R \
#       --truth=lowrank --seed=1 --orient=pca [--rotate=0] [--active=10] [--nchain=4] \
#       [--nburn=5000] [--niter=5000] --lib=<library> [--out=<dir>]
#
# --orient=pca replaces the embedding columns handed to the fit by their principal component
# scores: the columns are centred (not rescaled one by one, which would depend on how they are
# rotated) and multiplied by the right singular vectors. All components are kept, so this is again
# only a rotation and the model can represent the same truth exactly. The scores are the same,
# up to the sign of each, whatever rotation the embedding arrives in; the script checks that
# against rotation 1 of rotation_test.R before fitting. The fit then standardises each score
# column as it does any column, so a component with little variance across sites becomes a
# unit-variance column whose true coefficient is small. --orient=asis with --rotate behaves as
# rotation_test.R. Scoring is that of rotation_test.R: only quantities a rotation does not change.
# All data are simulated here; no dataset content.
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
truth <- getopt("truth", "lowrank"); seed <- as.integer(getopt("seed", "1")); rotate <- as.integer(getopt("rotate", "0"))
nchain <- as.integer(getopt("nchain", "4"))
nburn <- as.integer(getopt("nburn", "5000")); niter <- as.integer(getopt("niter", "5000"))
n <- as.integer(getopt("n", "256")); S <- as.integer(getopt("S", "48")); K <- as.integer(getopt("K", "64"))
lib <- getopt("lib", NA_character_); n_active <- as.integer(getopt("active", "10"))
orient <- getopt("orient", "pca"); stopifnot(orient %in% c("asis", "pca"))
out_dir <- getopt("out", "dev/simstudy/grouped-mixing-20261118/results/orientation")
stopifnot(truth %in% c("null", "lowrank"))
if (!is.na(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
suppressMessages(library(occJSDM))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gt <- 4L; p_meas <- 5L
tag <- sprintf("%s%s_s%d_rot%d_%s_b%d", truth, if (n_active == 10L) "" else paste0("-active", n_active), seed, rotate, orient, nburn)

embedding <- list(n_emb = K)
if (truth == "null") embedding$B <- matrix(0, K, S) else embedding <- c(embedding, list(n_blends = 3L, n_active = min(n_active, K)))
set.seed(seed)
datasettings <- list(n = n, S = S, g = 2L, M = rep(1L, n), P = 1L, K = rep(1L, n),
                     ncov_psi = p_meas, ncov_theta = 1L)
jsdmParams <- list(gt = gt, d = 1L, ds = 0, sigma_b = 0.5, sigma_bs = 0.5, sigma_ts = 0.5,
                   sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, S),
                   useSpatField = FALSE, embedding = embedding)
sim <- simulateOccJSDMData(datasettings, list(), jsdmParams, model = "binary")
dat <- sim$data_list
emb_vars <- dat$embeddingVars
true_B <- sim$true_params$jsdmParams_true$B
emb_true <- which(sim$true_params$jsdmParams_true$covariate_group == "embedding")
# The truth's contribution, in the columns as simulated (they are already centred and scaled).
X_sim <- as.matrix(dat$info[, emb_vars])
true_contribution <- unname(X_sim %*% true_B[emb_true, ])
if (rotate > 0L) {
  set.seed(5000L + rotate)
  qrq <- qr(matrix(rnorm(K * K), K, K)); Q <- qr.Q(qrq) %*% diag(sign(diag(qr.R(qrq))))
  stopifnot(isTRUE(all.equal(crossprod(Q), diag(K))))
  X_rot <- X_sim %*% Q; colnames(X_rot) <- emb_vars
  dat$info[, emb_vars] <- X_rot
}
pc_scores <- function(x) { x <- scale(x, center = TRUE, scale = FALSE); sv <- svd(x); sv$u %*% diag(sv$d) }
if (orient == "pca") {
  X_pc <- pc_scores(as.matrix(dat$info[, emb_vars])); colnames(X_pc) <- emb_vars
  # The same scores, up to sign, from the embedding in another rotation (rotation 1).
  set.seed(5001L); qr1 <- qr(matrix(rnorm(K * K), K, K)); Q1 <- qr.Q(qr1) %*% diag(sign(diag(qr.R(qr1))))
  other <- pc_scores(X_sim %*% Q1)
  stopifnot(isTRUE(all.equal(abs(unname(X_pc)), abs(other), tolerance = 1e-6)))
  component_share <- colSums(X_pc^2) / sum(X_pc^2)
  dat$info[, emb_vars] <- X_pc
}
occ <- c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), emb_vars), emb_vars)

set.seed(100L + seed)      # the chain seeds, the same for every rotation of a data set
t0 <- Sys.time()
fit <- suppressMessages(suppressWarnings(runOccJSDM(
  dat, listParams = list(n_factors = 2L, n_lattrait = gt, embeddingVars = emb_vars),
  occCovariates = occ,
  MCMCparams = list(nchain = nchain, nburn = nburn, niter = niter, nthin = 1, ncores = 1))))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
jo <- fit$results_output$jsdm_output

plain <- function(x) { nn <- nrow(x); W <- mean(apply(x, 2, var)); Bv <- nn * var(colMeans(x)); sqrt(((nn - 1) / nn * W + Bv / nn) / W) }
split2 <- function(x) { h <- floor(nrow(x) / 2); cbind(x[seq_len(h), , drop = FALSE], x[h + seq_len(h), , drop = FALSE]) }
znorm <- function(x) matrix(qnorm((rank(x) - 3 / 8) / (length(x) + 1 / 4)), nrow(x))
rank_rhat <- function(x) { x <- split2(x); max(plain(znorm(x)), plain(znorm(abs(x - median(x))))) }
rms <- function(x) sqrt(mean(x^2))

B <- jo$B_output                                   # p x S x iter x chain
emb <- unname(which(fit$infos$covariate_group == "embedding")); meas <- setdiff(seq_len(dim(B)[1]), emb)
stopifnot(length(emb) == K)
X <- fit$X_psi; Xe <- unname(X[, emb, drop = FALSE])
if (rotate == 0L && orient == "asis") stopifnot(isTRUE(all.equal(Xe, unname(X_sim), tolerance = 1e-8)))
# Rotation-free quantities, by draw: the contribution at a fixed sample of 300 site-by-species
# cells (for R-hat), and the spread across sites per species.
set.seed(77L); cells <- cbind(sample(n, 300L, replace = TRUE), sample(S, 300L, replace = TRUE))
nit <- dim(B)[3]
cell_draws <- array(NA_real_, c(300L, nit, nchain)); spread <- array(NA_real_, c(S, nit, nchain))
contribution_mean <- array(0, c(n, S, nchain))
for (ch in seq_len(nchain)) for (it in seq_len(nit)) {
  contrib <- Xe %*% B[emb, , it, ch]
  cell_draws[, it, ch] <- contrib[cells]
  spread[, it, ch] <- apply(contrib, 2, stats::sd)
  contribution_mean[, , ch] <- contribution_mean[, , ch] + contrib / nit
}
r_cells <- apply(cell_draws, 1, rank_rhat); r_spread <- apply(spread, 1, rank_rhat)
r_meas <- sapply(seq_len(length(meas) * S), function(i) rank_rhat(B[meas[(i - 1) %% length(meas) + 1], (i - 1) %/% length(meas) + 1, , ]))
r_lam <- apply(log(jo$lambda_output), 1, rank_rhat)
log_sig <- log(jo$sigmaemb_output)
slope <- function(est) sum(est * true_contribution) / sum(true_contribution^2)

per_chain <- do.call(rbind, lapply(seq_len(nchain), function(ch) {
  Bc <- apply(B[meas, , , ch, drop = FALSE], c(1, 2), mean)
  data.frame(tag = tag, truth = truth, active = n_active, seed = seed, rotate = rotate, orient = orient, chain = ch,
             log_sigma_emb_mean = mean(log_sig[, ch]), log_sigma_emb_max = max(log_sig[, ch]),
             mean_log_axis = mean(log(jo$lambda_output[, , ch])), sd_log_axis_across_columns = stats::sd(rowMeans(log(jo$lambda_output[, , ch]))),
             err_contribution = rms(contribution_mean[, , ch] - true_contribution),
             slope_contribution = if (truth == "null") NA_real_ else slope(contribution_mean[, , ch]),
             spread_mean = mean(spread[, , ch]), err_measured = rms(Bc - true_B[meas, ]))
}))
pooled <- apply(contribution_mean, c(1, 2), mean)
between <- if (nchain > 1) mean(sapply(seq_len(nchain), function(ch) rms(contribution_mean[, , ch] - pooled))) else NA_real_
fit_row <- data.frame(tag = tag, truth = truth, active = n_active, seed = seed, rotate = rotate, orient = orient, nchain = nchain, nburn = nburn, niter = niter,
  seconds = round(secs), package_commit = tryCatch(readLines(file.path(dirname(getNamespaceInfo("occJSDM", "path")), "COMMIT"), 1), error = function(e) NA_character_),
  err_contribution = rms(pooled - true_contribution), slope_contribution = if (truth == "null") NA_real_ else slope(pooled),
  truth_rms_contribution = rms(true_contribution), chain_rms_from_pooled = between,
  err_measured = rms(apply(B[meas, , , , drop = FALSE], c(1, 2), mean) - true_B[meas, ]),
  cells_over = mean(r_cells > 1.1), cells_rhat_max = max(r_cells), spread_over = mean(r_spread > 1.1), spread_rhat_max = max(r_spread),
  meas_over = mean(r_meas > 1.1), meas_rhat_max = max(r_meas), axis_over = sum(r_lam > 1.1), axis_rhat_max = max(r_lam),
  logsig_rhat = rank_rhat(log_sig), sigma_emb_median = stats::median(jo$sigmaemb_output),
  chains_in_blend = sum(colMeans(log_sig) < 0.5))
utils::write.csv(per_chain, file.path(out_dir, paste0(tag, "_chains.csv")), row.names = FALSE)
utils::write.csv(fit_row, file.path(out_dir, paste0(tag, ".csv")), row.names = FALSE)
saveRDS(list(tag = tag, contribution_mean = contribution_mean, true_contribution = true_contribution,
             spread_chain_mean = apply(spread, c(1, 3), mean),
             B_measured_chain_mean = apply(B[meas, , , , drop = FALSE], c(1, 2, 4), mean), true_B_measured = true_B[meas, ],
             component_share = if (orient == "pca") component_share else NULL,
             log_axis_mean = apply(log(jo$lambda_output), 1, mean)),
        file.path(out_dir, paste0(tag, "_summary.rds")))
num <- vapply(per_chain, is.numeric, logical(1)); per_chain[num] <- lapply(per_chain[num], round, 3)
print(per_chain[, -(1:6)], row.names = FALSE); print(t(fit_row), quote = FALSE)
