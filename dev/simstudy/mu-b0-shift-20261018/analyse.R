# Mixing of B0's community mean with the collapsed shift step on (README.md).
# Run from the repository root:
#   Rscript dev/simstudy/mu-b0-shift-20261018/analyse.R
# Reads draws/drift (refit-drift.R, the three diagnostic fits) and draws/main
# (refit-s100.R, the 100-species check), whichever exist, and the committed
# tables of the diagnostic (../mu-b0-drift-20261018) and of the check
# (../b0-community-s100-20261018) for the "before".
# Writes:
#   drift.csv          per diagnostic fit and quantity: chain 2 minus chain 1
#                      of the mean over draws, the posterior sd, R-hat and
#                      bulk ESS, before (the diagnostic) and after (the step on)
#   drift-species.csv  per diagnostic fit: the across-species sd of the B0
#                      shift before and after, the step's acceptance rates and
#                      final proposal sd per chain, and the elapsed time
#   fits.csv           one row per fit of the main step: the statistics of
#                      b0-community-s100-20261018/analyse.R plus the step's
#                      acceptance rates and proposal sd
#   before-after.csv   per fit of the main step: mu_b0's R-hat, bulk and tail
#                      ESS and posterior sd, the mean B0's bulk ESS and the
#                      elapsed time, before and after
#   rule.csv           the check's rule by arm, strict and median
#   traces.rds         per-draw mu_b0, mean B0 and mean collection intercept
#                      of every fit of both steps
#   traces-drift.png   mu_b0 and the mean B0 in the diagnostic fits, before
#                      and after
#   traces-mu-b0.png   mu_b0 in the main step's fits, before and after
here <- Sys.getenv("RUN_DIR", "dev/simstudy/mu-b0-shift-20261018")
out <- Sys.getenv("OUT_DIR", here)
diag_dir <- file.path(dirname(here), "mu-b0-drift-20261018")
check_dir <- file.path(dirname(here), "b0-community-s100-20261018")
rhat <- posterior::rhat; ess_bulk <- posterior::ess_bulk; ess_tail <- posterior::ess_tail

load_step <- function(step) {
  f <- sort(list.files(file.path(here, "draws", step), "[.]rds$", full.names = TRUE))
  setNames(lapply(f, readRDS), sub("[.]rds$", "", basename(f)))
}
drift <- load_step("drift")
main <- load_step("main")
stopifnot(length(drift) > 0 || length(main) > 0)
traces <- list()
ink <- "#0b0b0b"; ink2 <- "#52514e"; rule_col <- "#b5b3ad"; chain_col <- c("#2a78d6", "#eb6834")
shift_fields <- function(x) {
  b <- x$b0_shift
  data.frame(acceptance_chain1 = b$acceptance[1], acceptance_chain2 = b$acceptance[2],
             acceptance_burnin_chain1 = b$acceptance_burnin[1], acceptance_burnin_chain2 = b$acceptance_burnin[2],
             sd_chain1 = b$sd[1], sd_chain2 = b$sd[2], proposals = b$proposals)
}

# --- The diagnostic fits -----------------------------------------------------
if (length(drift)) {
  before_diff <- read.csv(file.path(diag_dir, "chain-difference.csv"))
  before_part <- read.csv(file.path(diag_dir, "partners.csv"))
  before_stat <- function(fit, stat, term) {
    v <- before_part$value[before_part$fit == fit & before_part$statistic == stat & before_part$term == term]
    if (length(v)) v[1] else NA_real_
  }
  per_draw <- function(x) {
    d <- x$draws
    S <- dim(d$B0)[1]; Tn <- dim(d$B0)[2]; nch <- dim(d$B0)[3]
    q <- list(mu_b0 = d$mu_b0, mean_B0 = apply(d$B0, c(2, 3), mean),
              mean_psi = apply(d$psi, c(3, 4), mean), occupied = apply(d$z, c(3, 4), sum),
              mean_bt0 = apply(d$beta_theta0, c(2, 3), mean),
              mean_logit_theta0 = apply(qlogis(d$theta0), c(2, 3), mean))
    lapply(q, function(m) matrix(m, Tn, nch))
  }
  drift_rows <- list(); species_rows <- list()
  for (key in names(drift)) {
    x <- drift[[key]]; q <- per_draw(x)
    traces[[paste0("drift:", key)]] <- lapply(q[c("mu_b0", "mean_B0", "mean_bt0")], signif, 7)
    drift_rows[[key]] <- do.call(rbind, lapply(names(q), function(nm) {
      m <- q[[nm]]
      b <- before_diff[before_diff$fit == key & before_diff$quantity == nm, ]
      data.frame(fit = key, quantity = nm,
                 chain1_after = mean(m[, 1]), chain2_after = mean(m[, 2]),
                 difference_before = if (nrow(b)) b$difference[1] else NA_real_,
                 difference_after = mean(m[, 2]) - mean(m[, 1]),
                 posterior_sd_before = if (nrow(b)) b$posterior_sd[1] else NA_real_,
                 posterior_sd_after = sd(as.vector(m)),
                 rhat_before = before_stat(key, "rhat", nm), rhat_after = rhat(m),
                 ess_bulk_before = before_stat(key, "ess_bulk", nm), ess_bulk_after = ess_bulk(m),
                 ess_tail_after = ess_tail(m))
    }))
    d_B0 <- apply(x$draws$B0[, , 2], 1, mean) - apply(x$draws$B0[, , 1], 1, mean)
    species_rows[[key]] <- cbind(
      data.frame(fit = key, sd_across_species_of_d_B0_before = before_stat(key, "across_species_sd_of_d", "B0"),
                 sd_across_species_of_d_B0_after = sd(d_B0), mean_d_B0_after = mean(d_B0),
                 max_abs_d_B0_after = max(abs(d_B0))),
      shift_fields(x), data.frame(elapsed_minutes = x$elapsed / 60))
  }
  drift_tab <- do.call(rbind, drift_rows); rownames(drift_tab) <- NULL
  drift_tab$difference_in_sds_before <- drift_tab$difference_before / drift_tab$posterior_sd_before
  drift_tab$difference_in_sds_after <- drift_tab$difference_after / drift_tab$posterior_sd_after
  write.csv(drift_tab, file.path(out, "drift.csv"), row.names = FALSE)
  species_tab <- do.call(rbind, species_rows); rownames(species_tab) <- NULL
  write.csv(species_tab, file.path(out, "drift-species.csv"), row.names = FALSE)

  # Trace plots: before from the check's traces.rds (the same chains, every
  # draw; thinned by 10 here to match), after from this run's stored draws.
  before_tr <- readRDS(file.path(check_dir, "traces.rds"))
  png(file.path(out, "traces-drift.png"), width = 2000, height = 1500, res = 150)
  op <- par(mfcol = c(4, length(drift)), mar = c(3.4, 3.4, 2.4, 0.6), oma = c(0, 0, 3, 0),
            mgp = c(1.7, 0.5, 0), tcl = -0.25, col.axis = ink2, col.lab = ink2, fg = rule_col,
            cex.axis = 0.8, cex.main = 0.85)
  panel <- function(m, label, stat) {
    plot(NA, xlim = c(1, nrow(m)), ylim = range(m), xlab = "stored draw", ylab = label, bty = "l", las = 1)
    for (k in seq_len(ncol(m))) lines(seq_len(nrow(m)), m[, k], col = adjustcolor(chain_col[k], 0.6), lwd = 0.6)
    title(main = stat, col.main = ink, font.main = 1)
  }
  for (key in names(drift)) {
    thin <- function(m) m[seq(10, nrow(m), by = 10), , drop = FALSE]
    bm <- thin(before_tr[[key]]$mu_b0); bb <- thin(before_tr[[key]]$mean_B0)
    am <- traces[[paste0("drift:", key)]]$mu_b0; ab <- traces[[paste0("drift:", key)]]$mean_B0
    r <- drift_tab[drift_tab$fit == key, ]
    panel(bm, "mu_b0", sprintf("%s, before\nR-hat %.3f, bulk ESS %.0f", key,
                               r$rhat_before[r$quantity == "mu_b0"], r$ess_bulk_before[r$quantity == "mu_b0"]))
    panel(am, "mu_b0", sprintf("after\nR-hat %.3f, bulk ESS %.0f",
                               r$rhat_after[r$quantity == "mu_b0"], r$ess_bulk_after[r$quantity == "mu_b0"]))
    panel(bb, "mean B0", sprintf("before, chain difference %+.3f", r$difference_before[r$quantity == "mean_B0"]))
    panel(ab, "mean B0", sprintf("after, chain difference %+.3f", r$difference_after[r$quantity == "mean_B0"]))
  }
  mtext("The diagnostic's three fits without (before) and with (after) the shift step; chain 1 blue, chain 2 orange; every 10th draw",
        outer = TRUE, line = 1, col = ink, cex = 0.95)
  par(op); invisible(dev.off())
  options(width = 200)
  cat("\n== drift: chain 2 minus chain 1, before and after\n")
  print(drift_tab[, c("fit", "quantity", "difference_before", "difference_after", "difference_in_sds_before",
                      "difference_in_sds_after", "rhat_before", "rhat_after", "ess_bulk_before", "ess_bulk_after")],
        digits = 3, row.names = FALSE)
  cat("\n== drift: species and the step\n"); print(species_tab, digits = 3, row.names = FALSE)
}

# --- The 100-species check ---------------------------------------------------
if (length(main)) {
  arms <- data.frame(arm = c("base_cm_s100", "base_rare_cm_s100"), species = 100L,
                     community = c("ordinary", "rare"), stringsAsFactors = FALSE)
  before_fits <- read.csv(file.path(check_dir, "fits.csv"))
  per_element <- function(a, prefix) {
    rh <- apply(a, 1, rhat); eb <- apply(a, 1, ess_bulk)
    setNames(c(median(rh), max(rh), mean(rh > 1.05), median(eb), min(eb)),
             paste0(prefix, c("_rhat_median", "_rhat_max", "_share_rhat_over_1.05", "_ess_median", "_ess_min")))
  }
  batch_sd <- function(m, b = 500) sd(as.vector(apply(m, 2, function(v) colMeans(matrix(v, nrow = b)))))
  chain_cor <- function(a, b) mean(vapply(seq_len(ncol(a)), function(k) cor(a[, k], b[, k]), 0))
  fit_row <- function(key) {
    x <- main[[key]]
    mb <- x$draws$mu_b0; B0 <- x$draws$B0; bt <- x$draws$beta_theta0
    mB0 <- apply(B0, c(2, 3), mean); mbt <- apply(bt, c(2, 3), mean)
    traces[[paste0("main:", key)]] <<- lapply(list(mu_b0 = mb, mean_B0 = mB0, mean_bt0 = mbt), signif, 7)
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
      mu_b0_batch_sd = batch_sd(mb), mu_b0_batch_sd_ratio = batch_sd(mb) / sd(as.vector(mb)),
      b0_error = mean(b0$post_mean - b0$truth), b0_coverage = mean(b0$covered),
      bt0_error = mean(bti$post_mean - bti$truth), bt0_coverage = mean(bti$covered),
      occ_true = mean(x$occupancy$occ_true),
      occ_error = mean(x$occupancy$occ_est) - mean(x$occupancy$occ_true),
      elapsed_minutes = x$elapsed / 60, check.names = FALSE), shift_fields(x))
  }
  fr <- do.call(rbind, lapply(names(main), fit_row))
  fr <- fr[order(match(fr$arm, arms$arm), fr$replicate), ]; rownames(fr) <- NULL
  write.csv(fr, file.path(out, "fits.csv"), row.names = FALSE)

  compare <- c("mu_b0_rhat", "mu_b0_ess_bulk", "mu_b0_ess_tail", "mu_b0_sd", "mu_b0_mean",
               "mean_b0_rhat", "mean_b0_ess_bulk", "b0_ess_median", "b0_share_rhat_over_1.05",
               "bt0_ess_median", "cor_mu_b0_mean_b0", "mu_b0_batch_sd_ratio", "occ_error",
               "b0_coverage", "elapsed_minutes")
  ba <- do.call(rbind, lapply(seq_len(nrow(fr)), function(i) {
    b <- before_fits[before_fits$arm == fr$arm[i] & before_fits$replicate == fr$replicate[i], ]
    data.frame(arm = fr$arm[i], community = fr$community[i], replicate = fr$replicate[i], quantity = compare,
               before = if (nrow(b)) unlist(b[1, compare]) else NA_real_, after = unlist(fr[i, compare]))
  }))
  rownames(ba) <- NULL
  write.csv(ba, file.path(out, "before-after.csv"), row.names = FALSE)

  rule <- do.call(rbind, lapply(split(fr, factor(fr$arm, arms$arm)), function(d) {
    if (!nrow(d)) return(NULL)
    data.frame(arm = d$arm[1], species = d$species[1], community = d$community[1], replicates = nrow(d),
               rhat_under_1.05 = sum(d$mu_b0_rhat < 1.05), rhat_max = max(d$mu_b0_rhat),
               ess_at_least_100 = sum(d$mu_b0_ess_bulk >= 100), ess_min = min(d$mu_b0_ess_bulk),
               ess_median = median(d$mu_b0_ess_bulk),
               pass_strict = all(d$mu_b0_rhat < 1.05) && all(d$mu_b0_ess_bulk >= 100),
               pass_median = all(d$mu_b0_rhat < 1.05) && median(d$mu_b0_ess_bulk) >= 100)
  }))
  rownames(rule) <- NULL
  write.csv(rule, file.path(out, "rule.csv"), row.names = FALSE)

  before_tr <- readRDS(file.path(check_dir, "traces.rds"))
  reps <- sort(unique(fr$replicate))
  png(file.path(out, "traces-mu-b0.png"), width = 2000, height = 380 * 2 * nrow(arms) + 120, res = 150)
  op <- par(mfrow = c(2 * nrow(arms), length(reps)), mar = c(3.4, 3.2, 2.4, 0.6), oma = c(0, 0, 3, 0),
            mgp = c(1.7, 0.5, 0), tcl = -0.25, col.axis = ink2, col.lab = ink2, fg = rule_col,
            cex.axis = 0.8, cex.main = 0.85)
  for (i in seq_len(nrow(arms))) for (when in c("before", "after")) for (r in reps) {
    key <- sprintf("%s-%02d", arms$arm[i], r)
    f <- fr[fr$arm == arms$arm[i] & fr$replicate == r, ]
    b <- before_fits[before_fits$arm == arms$arm[i] & before_fits$replicate == r, ]
    mb <- if (when == "before") before_tr[[key]]$mu_b0 else traces[[paste0("main:", key)]]$mu_b0
    if (is.null(mb) || !nrow(f)) { plot.new(); next }
    plot(NA, xlim = c(1, nrow(mb)), ylim = range(mb, f$true_b0_mean), xlab = "draw",
         ylab = if (r == reps[1]) "mu_b0" else "", bty = "l", las = 1)
    for (k in seq_len(ncol(mb))) lines(seq_len(nrow(mb)), mb[, k], col = adjustcolor(chain_col[k], 0.55), lwd = 0.5)
    abline(h = f$true_b0_mean, col = ink2, lty = 2, lwd = 1)
    st <- if (when == "before") c(b$mu_b0_rhat[1], b$mu_b0_ess_bulk[1]) else c(f$mu_b0_rhat, f$mu_b0_ess_bulk)
    title(main = sprintf("%s, replicate %d, %s\nR-hat %.3f, bulk ESS %.0f", arms$community[i], r, when, st[1], st[2]),
          col.main = ink, font.main = 1)
  }
  mtext("mu_b0 at 100 species and 5,000 + 5,000, without (before) and with (after) the shift step; chain 1 blue, chain 2 orange; dashed: mean of the true B0",
        outer = TRUE, line = 1, col = ink, cex = 0.95)
  par(op); invisible(dev.off())

  options(width = 200)
  cat("\n== rule\n"); print(rule, digits = 4, row.names = FALSE)
  show <- c("arm", "replicate", "mu_b0_rhat", "mu_b0_ess_bulk", "mu_b0_ess_tail", "mu_b0_sd", "mu_b0_mean",
            "true_b0_mean", "b0_share_rhat_over_1.05", "b0_ess_median", "bt0_ess_median", "cor_mu_b0_mean_b0",
            "mean_b0_ess_bulk", "mu_b0_batch_sd_ratio", "occ_error", "b0_coverage",
            "acceptance_chain1", "acceptance_chain2", "sd_chain1", "sd_chain2", "elapsed_minutes")
  cat("\n== fits\n"); print(fr[, show], digits = 3, row.names = FALSE)
  cat("\n== before and after\n")
  wide <- reshape(ba, direction = "wide", idvar = c("arm", "community", "replicate"), timevar = "quantity",
                  v.names = c("before", "after"))
  print(wide[, c("arm", "replicate", "before.mu_b0_rhat", "after.mu_b0_rhat", "before.mu_b0_ess_bulk",
                 "after.mu_b0_ess_bulk", "before.mean_b0_ess_bulk", "after.mean_b0_ess_bulk",
                 "before.elapsed_minutes", "after.elapsed_minutes")], digits = 3, row.names = FALSE)
}
saveRDS(traces, file.path(out, "traces.rds"), compress = "xz")
