# Interval-calibration assessment of the 15 October 2026 R = 300 paper study.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/paper-study-20261015/calibration/calibration.R \
#     [path/to/simstudy-20260729-143756.rds] [path/to/simstudy-20260810-121954.rds]
#
# The July and August baselines are gitignored and live only in the main
# checkout's dev/simstudy/results/. Without them the script still writes the
# R = 300 summaries and skips the comparisons.
#
# Coverage uncertainty is computed across replicates, never across pooled
# elements: the elements of one block in one replicate share a dataset and an
# MCMC run, so they are not independent checks (PLAN.md section 9).

args <- commandArgs(trailingOnly = TRUE)
run_dir <- "dev/simstudy/paper-study-20261015"
out_dir <- file.path(run_dir, "calibration")
jul_path <- if (length(args) >= 1) args[1] else NA
aug_path <- if (length(args) >= 2) args[2] else NA

new <- readRDS(file.path(run_dir, "results", "simstudy-20261015-172834.rds"))$rows
new$width <- new$upper - new$lower
new$err <- new$post_mean - new$truth

cells <- c("base", "spatial_isolated", "traits_isolated", "primers_3",
           "d_underfit", "d_overfit", "species_20", "occupancy", "binary",
           "low_information")
blocks <- c("B0", "B", "G", "beta_theta", "theta0", "p", "q", "sigma_b",
            "resid_cor")
spatial_cells <- setdiff(cells, "traits_isolated")

key <- function(d) paste(d$scenario, d$replicate, d$block, d$element, sep = "|")

# Per-replicate means of a column, one row per (cell, block, replicate).
per_rep <- function(d, col) {
  a <- aggregate(d[[col]], d[c("scenario", "block", "replicate")], mean)
  names(a)[4] <- col
  a
}

# --- 1. Replicate-based coverage, width, bias and RMSE at R = 300 ----------

summarise_cell_block <- function(d) {
  rc <- per_rep(d, "covered")
  rw <- per_rep(d, "width")
  re <- per_rep(d, "err")
  out <- lapply(split(seq_len(nrow(rc)), paste(rc$scenario, rc$block)), function(ix) {
    sc <- rc$scenario[ix[1]]; bl <- rc$block[ix[1]]
    sub <- d[d$scenario == sc & d$block == bl, ]
    c_r <- rc$covered[ix]
    R <- length(c_r)
    est <- mean(c_r)
    se <- sd(c_r) / sqrt(R)
    pooled <- mean(sub$covered)
    n_el <- nrow(sub)
    se_naive <- sqrt(pooled * (1 - pooled) / n_el)
    w_r <- rw$width[match(paste(sc, bl, rc$replicate[ix]),
                          paste(rw$scenario, rw$block, rw$replicate))]
    e_r <- re$err[match(paste(sc, bl, rc$replicate[ix]),
                        paste(re$scenario, re$block, re$replicate))]
    # Standardised error: error over the posterior SD implied by the interval.
    # Near 1 for a calibrated, roughly normal posterior; above 1 means the
    # intervals are narrow for the size of the errors they make.
    zsd <- sqrt(mean((sub$err / (sub$width / (2 * qnorm(0.975))))^2))
    data.frame(
      scenario = sc, block = bl, R = R,
      elements_per_rep = n_el / R,
      coverage = est, se = se,
      lo95 = est - 1.96 * se, hi95 = est + 1.96 * se,
      z_vs_nominal = if (se > 0) (est - 0.95) / se else NA_real_,
      pooled_coverage = pooled, naive_se = se_naive,
      design_effect = if (se_naive > 0) (se / se_naive)^2 else NA_real_,
      mean_width = mean(w_r), width_se = sd(w_r) / sqrt(R),
      median_width = median(sub$width),
      bias = mean(sub$err), bias_se = sd(e_r) / sqrt(R),
      rmse = sqrt(mean(sub$err^2)),
      err_to_width = zsd,
      # Slope of posterior mean on truth: 1 means no shrinkage of the point
      # estimates, 0 means they ignore the truth entirely.
      slope = if (stats::var(sub$truth) > 0)
        unname(coef(lm(post_mean ~ truth, sub))[2]) else NA_real_,
      stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  res <- res[order(match(res$scenario, cells), match(res$block, blocks)), ]
  rownames(res) <- NULL
  res
}

cov300 <- summarise_cell_block(new)
write.csv(cov300, file.path(out_dir, "coverage-r300.csv"), row.names = FALSE)

# Nine-cell pooled view per block, excluding low_information. Cells are
# independent, so the SE combines the per-cell replicate SEs.
pool9 <- do.call(rbind, lapply(blocks, function(bl) {
  x <- cov300[cov300$block == bl & cov300$scenario != "low_information", ]
  if (!nrow(x)) return(NULL)
  data.frame(block = bl, cells = nrow(x), coverage = mean(x$coverage),
             se = sqrt(sum(x$se^2)) / nrow(x),
             min_cell = x$scenario[which.min(x$coverage)],
             min_coverage = min(x$coverage),
             max_cell = x$scenario[which.max(x$coverage)],
             max_coverage = max(x$coverage),
             mean_width = mean(x$mean_width), stringsAsFactors = FALSE)
}))
write.csv(pool9, file.path(out_dir, "coverage-r300-nine-cells.csv"), row.names = FALSE)

# Within-block detail for beta_theta: intercept row and slope row separately.
# Element order is (row, species) column-major on a 2 x S matrix, so odd
# elements are intercepts and even elements are slopes.
bt <- new[new$block == "beta_theta", ]
bt$row <- ifelse(bt$element %% 2 == 1, "intercept", "slope")
bt_rows <- do.call(rbind, lapply(split(bt, list(bt$scenario, bt$row), drop = TRUE), function(d) {
  c_r <- tapply(d$covered, d$replicate, mean)
  data.frame(scenario = d$scenario[1], row = d$row[1], coverage = mean(c_r),
             se = sd(c_r) / sqrt(length(c_r)), mean_width = mean(d$width),
             bias = mean(d$err), stringsAsFactors = FALSE)
}))
bt_rows <- bt_rows[order(match(bt_rows$scenario, cells), bt_rows$row), ]
write.csv(bt_rows, file.path(out_dir, "beta-theta-rows-r300.csv"), row.names = FALSE)

# --- 2. Paired comparison with the 29 July 2026 R = 100 baseline -----------

if (!is.na(jul_path) && file.exists(jul_path)) {
  jul <- readRDS(jul_path)$rows
  jul$width <- jul$upper - jul$lower
  n100 <- new[new$replicate <= 100, ]

  # Pairing check: same elements, same truth where the scoring did not change.
  m <- match(key(n100), key(jul))
  stopifnot(!anyNA(m), nrow(jul) == nrow(n100))
  jt <- jul$truth[m]
  pairing <- do.call(rbind, lapply(blocks, function(bl) {
    ix <- n100$block == bl
    data.frame(block = bl, elements = sum(ix),
               max_abs_truth_diff = max(abs(n100$truth[ix] - jt[ix])),
               share_identical = mean(abs(n100$truth[ix] - jt[ix]) < 1e-12),
               stringsAsFactors = FALSE)
  }))
  write.csv(pairing, file.path(out_dir, "pairing-check.csv"), row.names = FALSE)

  # July intervals rescored against the October truth for beta_theta, p and
  # q, whose truth changed purely in scoring units while the fitted scale did
  # not: the transformation depends only on the simulated data, which are the
  # same in both runs. G is left on July's own truth, because July's fit used
  # raw traits, so its G intervals are in raw units and July's raw truth is
  # the right comparison for them. Every other block's truth is identical.
  rescore <- c("beta_theta", "p", "q")
  jul_r <- jul[m, ]
  jul_r$truth_new <- ifelse(jul_r$block %in% rescore, n100$truth, jul_r$truth)
  jul_r$covered_rescored <- as.integer(jul_r$lower <= jul_r$truth_new &
                                       jul_r$truth_new <= jul_r$upper)

  paired <- do.call(rbind, lapply(split(seq_len(nrow(n100)),
                                        paste(n100$scenario, n100$block)), function(ix) {
    a <- n100[ix, ]; b <- jul_r[ix, ]
    ca <- tapply(a$covered, a$replicate, mean)
    cb <- tapply(b$covered, b$replicate, mean)
    cbr <- tapply(b$covered_rescored, b$replicate, mean)
    dd <- ca - cb[names(ca)]
    ddr <- ca - cbr[names(ca)]
    data.frame(
      scenario = a$scenario[1], block = a$block[1],
      jul_coverage = mean(cb), jul_coverage_rescored = mean(cbr),
      oct_coverage = mean(ca),
      diff = mean(dd), diff_se = sd(dd) / sqrt(length(dd)),
      diff_vs_rescored = mean(ddr), diff_vs_rescored_se = sd(ddr) / sqrt(length(ddr)),
      jul_width = mean(b$width), oct_width = mean(a$width),
      width_ratio = mean(a$width) / mean(b$width),
      jul_bias = mean(b$post_mean - b$truth_new), oct_bias = mean(a$err),
      jul_rmse = sqrt(mean((b$post_mean - b$truth_new)^2)),
      oct_rmse = sqrt(mean(a$err^2)),
      flips_lost = sum(b$covered_rescored == 1 & a$covered == 0),
      flips_gained = sum(b$covered_rescored == 0 & a$covered == 1),
      cor_post_mean = cor(a$post_mean, b$post_mean),
      stringsAsFactors = FALSE)
  }))
  paired <- paired[order(match(paired$block, blocks), match(paired$scenario, cells)), ]
  rownames(paired) <- NULL
  write.csv(paired, file.path(out_dir, "paired-july-r100.csv"), row.names = FALSE)

  # Shrinkage diagnostics for B and B0: regression slope of posterior mean on
  # truth (1 = no attenuation), and the standardised error, per cell and run.
  shrink <- do.call(rbind, lapply(c("B", "B0"), function(bl) {
    do.call(rbind, lapply(cells, function(sc) {
      a <- n100[n100$block == bl & n100$scenario == sc, ]
      b <- jul_r[jul_r$block == bl & jul_r$scenario == sc, ]
      if (!nrow(a)) return(NULL)
      zs <- function(err, w) sqrt(mean((err / (w / (2 * qnorm(0.975))))^2))
      data.frame(block = bl, scenario = sc,
                 jul_slope = unname(coef(lm(post_mean ~ truth, b))[2]),
                 oct_slope = unname(coef(lm(post_mean ~ truth, a))[2]),
                 jul_err_to_width = zs(b$post_mean - b$truth, b$width),
                 oct_err_to_width = zs(a$err, a$width),
                 stringsAsFactors = FALSE)
    }))
  }))
  write.csv(shrink, file.path(out_dir, "shrinkage-b-b0.csv"), row.names = FALSE)
}

# --- 3. Where the change happened: the 10 August 2026 run ------------------

if (!is.na(aug_path) && file.exists(aug_path)) {
  aug <- readRDS(aug_path)$rows
  aug$width <- aug$upper - aug$lower
  n100 <- new[new$replicate <= 100, ]
  tl <- do.call(rbind, lapply(c("B", "B0"), function(bl) {
    do.call(rbind, lapply(cells, function(sc) {
      a <- aug[aug$block == bl & aug$scenario == sc, ]
      b <- n100[n100$block == bl & n100$scenario == sc, ]
      data.frame(block = bl, scenario = sc,
                 aug10_coverage = mean(a$covered), oct_r100_coverage = mean(b$covered),
                 aug10_width = mean(a$width), oct_r100_width = mean(b$width),
                 aug10_bias = mean(a$post_mean - a$truth),
                 oct_r100_bias = mean(b$post_mean - b$truth),
                 stringsAsFactors = FALSE)
    }))
  }))
  write.csv(tl, file.path(out_dir, "august-10-b-b0.csv"), row.names = FALSE)
}

cat("written to", out_dir, "\n")
