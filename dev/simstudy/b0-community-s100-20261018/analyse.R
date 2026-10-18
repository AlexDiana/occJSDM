# Mixing of B0's community mean at 100 species against 10 (README.md).
# Run from the repository root: Rscript dev/simstudy/b0-community-s100-20261018/analyse.R
# Reads the kept draws of the main step (draws/main, 100 species) and the s10
# step (draws/s10, 10 species), and the timing step's fits (draws/timing).
# Writes:
#   reproduction.csv  the s10 refits against the M = 3 runs' long-chains.csv
#   fits.csv          one row per fit: mu_b0, B0, collection intercepts, the
#                     slow direction, occupancy error, elapsed time
#   summary.csv       the same averaged by arm, with medians and ranges
#   rule.csv          README.md's rule, strict (every replicate) and median
#   timing.csv        elapsed time of the 200 + 200 timing fits
#   traces.rds        per-draw mu_b0, mean B0 and mean collection intercept
#   traces-mu-b0.png  trace plots of mu_b0, both chains, every fit
here <- Sys.getenv("RUN_DIR", "dev/simstudy/b0-community-s100-20261018")
out <- Sys.getenv("OUT_DIR", here)
rhat <- posterior::rhat; ess_bulk <- posterior::ess_bulk; ess_tail <- posterior::ess_tail

arms <- data.frame(
  arm = c("base_cm_s100", "base_rare_cm_s100", "base_cm_m3", "base_rare_cm_m3"),
  species = c(100L, 100L, 10L, 10L),
  community = c("ordinary", "rare", "ordinary", "rare"),
  step = c("main", "main", "s10", "s10"),
  stringsAsFactors = FALSE)
load_step <- function(step) {
  f <- sort(list.files(file.path(here, "draws", step), "[.]rds$", full.names = TRUE))
  setNames(lapply(f, readRDS), sub("[.]rds$", "", basename(f)))
}
fits <- c(load_step("s10"), load_step("main"))
stopifnot(length(fits) > 0)

# --- Reproduction: the s10 refits against the stored long-chain rows -------
# The same statistics as the M = 3 runs' long.R, compared with their "long"
# column for the same arm and replicate.
long_stat <- function(x) {
  mb <- x$draws$mu_b0
  c(occ_error = mean(x$occupancy$occ_est) - mean(x$occupancy$occ_true),
    mu_b0 = mean(mb), mu_b0_sd = sd(as.vector(mb)), mu_b0_rhat = rhat(mb), mu_b0_ess = ess_bulk(mb),
    b0_rhat_max = max(apply(x$draws$B0, 1, rhat)),
    b0_ess_median = median(apply(x$draws$B0, 1, ess_bulk)),
    bt0_rhat_max = max(apply(x$draws$beta_theta0, 1, rhat)),
    bt0_ess_median = median(apply(x$draws$beta_theta0, 1, ess_bulk)))
}
stored <- rbind(
  read.csv(file.path("dev", "simstudy", "ordinary-m3-20261016", "long-chains.csv")),
  read.csv(file.path("dev", "simstudy", "rare-m3-20261016", "long-chains.csv")))
s10 <- Filter(function(x) x$arm %in% c("base_cm_m3", "base_rare_cm_m3"), fits)
reproduction <- do.call(rbind, lapply(s10, function(x) {
  s <- long_stat(x)
  st <- stored[stored$arm == x$arm & stored$replicate == x$replicate, ]
  stored_val <- st$long[match(names(s), st$quantity)]
  data.frame(arm = x$arm, replicate = x$replicate, quantity = names(s),
             refit = unname(s), stored = stored_val, difference = unname(s) - stored_val)
}))
if (!is.null(reproduction)) {
  write.csv(reproduction, file.path(out, "reproduction.csv"), row.names = FALSE)
  cat(sprintf("reproduction: %d values, %d missing, max |difference| %.3g\n",
              nrow(reproduction), sum(is.na(reproduction$stored)),
              max(abs(reproduction$difference), na.rm = TRUE)))
}

# --- One row per fit --------------------------------------------------------
per_element <- function(a, prefix) {
  rh <- apply(a, 1, rhat); eb <- apply(a, 1, ess_bulk)
  setNames(c(median(rh), max(rh), mean(rh > 1.05), median(eb), min(eb)),
           paste0(prefix, c("_rhat_median", "_rhat_max", "_share_rhat_over_1.05",
                            "_ess_median", "_ess_min")))
}
chain_cor <- function(a, b) mean(vapply(seq_len(ncol(a)), function(k) cor(a[, k], b[, k]), 0))
traces <- list()
fit_row <- function(key) {
  x <- fits[[key]]
  mb <- x$draws$mu_b0                          # draws x chains
  B0 <- x$draws$B0                             # species x draws x chains
  bt <- x$draws$beta_theta0                    # species x draws x chains
  mB0 <- apply(B0, c(2, 3), mean); mbt <- apply(bt, c(2, 3), mean)
  traces[[key]] <<- lapply(list(mu_b0 = mb, mean_B0 = mB0, mean_bt0 = mbt), signif, 7)
  b0 <- x$rows[x$rows$block == "B0", ]
  bti <- x$rows[x$rows$block == "beta_theta" & x$rows$element %% 2 == 1, ]
  q <- quantile(mb, c(0.025, 0.975))
  a <- arms[arms$arm == x$arm, ]
  cbind(a, data.frame(
    replicate = x$replicate, draws = length(mb),
    mu_b0_rhat = rhat(mb), mu_b0_ess_bulk = ess_bulk(mb), mu_b0_ess_tail = ess_tail(mb),
    mu_b0_mean = mean(mb), mu_b0_sd = sd(as.vector(mb)), mu_b0_q2.5 = q[[1]], mu_b0_q97.5 = q[[2]],
    true_b0_mean = mean(b0$truth),
    mu_b0_covers_true_mean = q[[1]] <= mean(b0$truth) && mean(b0$truth) <= q[[2]],
    t(per_element(B0, "b0")), t(per_element(bt, "bt0")),
    cor_mu_b0_mean_b0 = chain_cor(mb, mB0), cor_mean_b0_mean_bt0 = chain_cor(mB0, mbt),
    mean_b0_rhat = rhat(mB0), mean_b0_ess_bulk = ess_bulk(mB0),
    mean_bt0_rhat = rhat(mbt), mean_bt0_ess_bulk = ess_bulk(mbt),
    b0_error = mean(b0$post_mean - b0$truth), b0_coverage = mean(b0$covered),
    bt0_error = mean(bti$post_mean - bti$truth), bt0_coverage = mean(bti$covered),
    occ_true = mean(x$occupancy$occ_true),
    occ_error = mean(x$occupancy$occ_est) - mean(x$occupancy$occ_true),
    elapsed_minutes = x$elapsed / 60, check.names = FALSE))
}
fr <- do.call(rbind, lapply(names(fits), fit_row))
fr <- fr[order(match(fr$arm, arms$arm), fr$replicate), ]; rownames(fr) <- NULL
write.csv(fr, file.path(out, "fits.csv"), row.names = FALSE)
saveRDS(traces[paste0(fr$arm, "-", sprintf("%02d", fr$replicate))], file.path(out, "traces.rds"),
        compress = "xz")

# --- By arm -----------------------------------------------------------------
num <- names(fr)[vapply(fr, is.numeric, TRUE) & !names(fr) %in% c("species", "replicate", "draws")]
summ <- do.call(rbind, lapply(split(fr, factor(fr$arm, arms$arm)), function(d) {
  if (!nrow(d)) return(NULL)
  do.call(rbind, lapply(num, function(v) data.frame(
    arm = d$arm[1], species = d$species[1], community = d$community[1], replicates = nrow(d),
    quantity = v, mean = mean(d[[v]]), median = median(d[[v]]), min = min(d[[v]]), max = max(d[[v]]))))
}))
rownames(summ) <- NULL
write.csv(summ, file.path(out, "summary.csv"), row.names = FALSE)

# --- The rule ---------------------------------------------------------------
rule <- do.call(rbind, lapply(split(fr, factor(fr$arm, arms$arm)), function(d) {
  if (!nrow(d)) return(NULL)
  data.frame(arm = d$arm[1], species = d$species[1], community = d$community[1],
             replicates = nrow(d),
             rhat_under_1.05 = sum(d$mu_b0_rhat < 1.05), rhat_max = max(d$mu_b0_rhat),
             ess_at_least_100 = sum(d$mu_b0_ess_bulk >= 100), ess_min = min(d$mu_b0_ess_bulk),
             ess_median = median(d$mu_b0_ess_bulk),
             pass_strict = all(d$mu_b0_rhat < 1.05) && all(d$mu_b0_ess_bulk >= 100),
             pass_median = all(d$mu_b0_rhat < 1.05) && median(d$mu_b0_ess_bulk) >= 100)
}))
rownames(rule) <- NULL
write.csv(rule, file.path(out, "rule.csv"), row.names = FALSE)

# --- Timing fits ------------------------------------------------------------
tf <- load_step("timing")
if (length(tf)) {
  tm <- do.call(rbind, lapply(tf, function(x) data.frame(
    arm = x$arm, replicate = x$replicate, nburn = x$mcmc$nburn, niter = x$mcmc$niter,
    elapsed_minutes = x$elapsed / 60,
    projected_minutes_5000 = x$elapsed / 60 * (5000 + 5000) / (x$mcmc$nburn + x$mcmc$niter))))
  rownames(tm) <- NULL
  write.csv(tm, file.path(out, "timing.csv"), row.names = FALSE)
}

# --- Trace plots of mu_b0 ---------------------------------------------------
# One panel per fit, both chains (categorical slots 1 and 2 of the dataviz
# reference palette), with the realised mean of the true B0 as a dashed line.
ink <- "#0b0b0b"; ink2 <- "#52514e"; rule_col <- "#b5b3ad"; chain_col <- c("#2a78d6", "#eb6834")
rows_shown <- arms[arms$arm %in% fr$arm, ]
png(file.path(out, "traces-mu-b0.png"), width = 2000, height = 380 * nrow(rows_shown) + 120, res = 150)
op <- par(mfrow = c(nrow(rows_shown), 5), mar = c(3.4, 3.2, 2.4, 0.6), oma = c(0, 0, 3, 0),
          mgp = c(1.7, 0.5, 0), tcl = -0.25, col.axis = ink2, col.lab = ink2, fg = rule_col,
          cex.axis = 0.8, cex.main = 0.85)
for (i in seq_len(nrow(rows_shown))) for (r in 1:5) {
  key <- sprintf("%s-%02d", rows_shown$arm[i], r)
  if (is.null(traces[[key]])) { plot.new(); next }
  mb <- traces[[key]]$mu_b0; f <- fr[fr$arm == rows_shown$arm[i] & fr$replicate == r, ]
  plot(NA, xlim = c(1, nrow(mb)), ylim = range(mb, f$true_b0_mean), xlab = "draw",
       ylab = if (r == 1) "mu_b0" else "", bty = "l", las = 1)
  for (k in seq_len(ncol(mb))) lines(seq_len(nrow(mb)), mb[, k], col = chain_col[k], lwd = 0.5)
  abline(h = f$true_b0_mean, col = ink2, lty = 2, lwd = 1)
  title(main = sprintf("%d species, %s, replicate %d\nR-hat %.3f, bulk ESS %.0f",
                       rows_shown$species[i], rows_shown$community[i], r, f$mu_b0_rhat,
                       f$mu_b0_ess_bulk), col.main = ink, font.main = 1)
}
mtext("Trace plots of mu_b0 at 5,000 + 5,000 (chain 1 blue, chain 2 orange; dashed: mean of the true B0)",
      outer = TRUE, line = 1, col = ink, cex = 0.95)
par(op); invisible(dev.off())

options(width = 200)
cat("\n== rule\n"); print(rule, digits = 4, row.names = FALSE)
show <- c("arm", "replicate", "mu_b0_rhat", "mu_b0_ess_bulk", "mu_b0_ess_tail", "mu_b0_sd",
          "mu_b0_mean", "true_b0_mean", "b0_rhat_median", "b0_share_rhat_over_1.05", "b0_ess_median",
          "bt0_share_rhat_over_1.05", "bt0_ess_median", "cor_mu_b0_mean_b0", "cor_mean_b0_mean_bt0",
          "mean_b0_ess_bulk", "occ_error", "elapsed_minutes")
cat("\n== fits\n"); print(fr[, show], digits = 3, row.names = FALSE)
