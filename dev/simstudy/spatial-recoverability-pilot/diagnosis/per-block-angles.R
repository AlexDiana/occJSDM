# Per-block mixing diagnostic for truth 1: the accepted ellipse angle of the one-angle (joint)
# sampler against the per-block angles of the blockwise study sampler, with per-block ESS for
# both, on one job at a given budget. Run from the worktree root:
#   run-r -f dev/simstudy/spatial-recoverability-pilot/diagnosis/per-block-angles.R --args JOB NBURN NITER
args <- commandArgs(TRUE); job <- args[1]; nburn <- as.integer(args[2]); niter <- as.integer(args[3])
PILOT <- "dev/simstudy/spatial-recoverability-pilot"
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R")) source(file.path(PILOT, f))
suppressWarnings(library(posterior))
compile_pilot_samplers(".")
Rcpp::sourceCpp(file.path(PILOT, "diagnosis/joint-angles.cpp"), cacheDir = file.path(tempdir(), "recoverability-pilot-cpp"))
S0 <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"   # the study directory, as in compare.R
cm <- readRDS(file.path(S0, "inputs", paste0(sub("-(uniform|clustered)$", "", job), ".rds"))); arr <- sub(".*-(uniform|clustered)$", "\\1", job)
idx <- cm$index[[arr]]; xy <- cm$points[idx, , drop = FALSE]; n <- length(idx); Sp <- nrow(cm$species); blocks <- 2L + Sp
f <- stable_factor(sq_exp(xy, xy, cm$range)); offset <- training_offset(cm, idx); w <- joint_weights(cm); ci <- sqrt(1 - cm$shared_fraction)
truth <- cm$spatial[idx, , drop = FALSE]; tc <- sweep(truth, 2, colMeans(truth), "-")
set.seed(1); initial <- cbind(rep(0, blocks * n), do.call(rbind, replicate(blocks, f$lower %*% matrix(rnorm(n * 3), n, 3), simplify = FALSE)))
Y <- cm$z[idx, , drop = FALSE]
set.seed(2); t0 <- Sys.time(); a <- joint_angles_cpp(f$lower, Y, offset, w, 2L, ci, tc, initial, nburn, niter, 10L); ta <- as.numeric(Sys.time() - t0, units = "secs")
set.seed(2); t0 <- Sys.time(); b <- block_ellipse_draws_cpp(f$lower, Y, offset, w, 2L, ci, tc, initial, nburn, niter, 10L); tb <- as.numeric(Sys.time() - t0, units = "secs")
fold <- function(th) th - 2 * pi * round(th / (2 * pi))
cat("\n==", job, "nburn", nburn, "niter", niter, " joint secs", round(ta), " block secs", round(tb), "\n")
th <- fold(a$angles)
cat("JOINT accepted angle, folded into (-pi, pi]: median |theta|", round(median(abs(th)), 4), " mean |theta|", round(mean(abs(th)), 4),
    " mean cos", round(mean(cos(th)), 5), " implied autocorrelation time (1+rho)/(1-rho) for a prior-like block:", round((1 + mean(cos(th))) / (1 - mean(cos(th)))), "\n")
cat("JOINT proposals per step", round(mean(a$mean_proposals), 2), "\n")
bn <- c("shared1", "shared2", paste0("indep", seq_len(Sp)))
tab <- cbind(proposals = rowMeans(b$block_proposals), mean_abs_angle = rowMeans(b$block_mean_abs_angle)); rownames(tab) <- bn
cat("BLOCK per-block proposals per block step and mean |accepted angle|:\n"); print(round(tab, 3))
nk <- dim(a$draws)[2]
per_block <- function(x) t(sapply(seq_len(blocks), function(bb) {
  d <- matrix(x$draws[(bb - 1) * n + 1:n, , ], n, nk * 4); mm <- array(colMeans(d), c(nk, 4)); ss <- array(apply(d, 2, sd), c(nk, 4))
  p <- array(colSums(d * rowMeans(d)) / sum(rowMeans(d)^2), c(nk, 4))
  c(ess_mean = ess_bulk(mm), ess_sd = ess_bulk(ss), ess_proj = ess_bulk(p), rhat_proj = rhat(p)) }))
pa <- per_block(a); pb <- per_block(b); rownames(pa) <- rownames(pb) <- bn
cat("per-block ESS of", nk * 4, "thinned draws, JOINT:\n"); print(round(pa, 1))
cat("per-block ESS of", nk * 4, "thinned draws, BLOCK:\n"); print(round(pb, 1))
summ <- function(tr) { sc <- screen_traces(tr); c(ess_med = median(sc$ess_bulk), ess_min = min(sc$ess_bulk), tail_min = min(sc$ess_tail), rhat_med = median(sc$rhat), rhat_max = max(sc$rhat), frac_pass = mean(sc$pass)) }
cat("36-trace study screen, JOINT:\n"); print(round(summ(a$traces), 3)); cat("36-trace study screen, BLOCK:\n"); print(round(summ(b$traces), 3))
cat("median posterior SD at sites (block run): shared1", round(median(apply(matrix(b$draws[1:n, , ], n, nk * 4), 1, sd)), 3),
    " indep1", round(median(apply(matrix(b$draws[2 * n + 1:n, , ], n, nk * 4), 1, sd)), 3), " (prior SD 1 for every block)\n")
