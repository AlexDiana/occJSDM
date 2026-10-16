# Score the field-sample ladder: does more replication stop the prior deciding occupancy?
# Run from the repository root: Rscript dev/simstudy/field-sample-ladder/analyse.R
# Environment overrides (all optional):
#   RESULTS  the run's results file (default: newest results/simstudy-*.rds in the run directory)
#   JULY     the 29 July 2026 ladder results for the paired B0 comparison (default: the main
#            checkout's dev/simstudy/results/simstudy-20260729-212525.rds; skipped if absent)
#   ARMS     comma-separated cells to score (default M2,M5,M10,M20,K30)
#   OUT_DIR  where the CSVs go (default: the run directory)
here <- "dev/simstudy/field-sample-ladder"
env <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
res_file <- env("RESULTS", { f <- list.files(file.path(here, "results"), "^simstudy-.*[.]rds$", full.names = TRUE); if (length(f)) tail(sort(f), 1) else NA })
july_file <- env("JULY", "/Users/douglasyu/src/occJSDM/dev/simstudy/results/simstudy-20260729-212525.rds")
arms <- strsplit(env("ARMS", "M2,M5,M10,M20,K30"), ",")[[1]]
out <- env("OUT_DIR", here)
stopifnot(!is.na(res_file), file.exists(res_file))
r <- readRDS(res_file); x <- r$rows; x <- x[x$block != "ERROR" & x$scenario %in% arms, ]
x$err <- x$post_mean - x$truth; x$width <- x$upper - x$lower
occ <- r$occupancy; if (!is.null(occ)) occ <- occ[occ$scenario %in% arms, ]
cat("results:", res_file, "\narms present:", paste(intersect(arms, unique(x$scenario)), collapse = ", "), "\n")
slope <- function(d, y = "post_mean", t = "truth") if (nrow(d) > 2 && var(d[[t]]) > 0) unname(coef(lm(d[[y]] ~ d[[t]]))[2]) else NA
rep_cov <- function(d) { c_r <- tapply(d$covered, d$replicate, mean); c(cov = mean(c_r), se = sd(c_r) / sqrt(length(c_r)), R = length(c_r)) }

# species-level B0 and collection-intercept errors, for the trade-off
species_err <- function(a) {
  b0 <- x[x$scenario == a & x$block == "B0", ]; if (!nrow(b0)) return(NULL)
  b0 <- data.frame(replicate = b0$replicate, species = b0$element, e_b0 = b0$err)
  bt <- x[x$scenario == a & x$block == "beta_theta", ]
  if (nrow(bt)) {  # (ncov_theta + 1) rows per species, column-major: the first row of each species is its intercept
    k <- nrow(bt) / nrow(b0); bt <- bt[(bt$element - 1) %% k == 0, ]
    b0 <- merge(b0, data.frame(replicate = bt$replicate, species = (bt$element - 1) %/% k + 1, e_th = bt$err), by = c("replicate", "species"))
  }
  b0
}

# 1. one row per arm
summ <- do.call(rbind, lapply(intersect(arms, unique(x$scenario)), function(a) {
  g <- function(b) x[x$scenario == a & x$block == b, ]
  b0 <- g("B0"); cb0 <- rep_cov(b0); cb <- rep_cov(g("B"))
  o <- if (!is.null(occ)) occ[occ$scenario == a, ] else NULL
  se <- species_err(a)
  row <- data.frame(arm = a, R = cb0[["R"]],
    B0_coverage = cb0[["cov"]], B0_se = cb0[["se"]], B0_mean_error = mean(b0$err), B0_slope = slope(b0), B0_width = mean(b0$width),
    occ_slope = if (!is.null(o) && nrow(o)) slope(o, "occ_est", "occ_true") else NA,
    occ_mean_error = if (!is.null(o) && nrow(o)) mean(o$occ_est - o$occ_true) else NA,
    occ_mean_abs_error = if (!is.null(o) && nrow(o)) mean(abs(o$occ_est - o$occ_true)) else NA,
    B0_vs_collection_cor = if (!is.null(se) && "e_th" %in% names(se)) cor(se$e_b0, se$e_th) else NA,
    B_coverage = cb[["cov"]], B_slope = slope(g("B")))
  for (b in c("theta0", "p", "q")) { d <- g(b); row[[paste0(b, "_coverage")]] <- if (nrow(d)) rep_cov(d)[["cov"]] else NA; row[[paste0(b, "_mean_error")]] <- if (nrow(d)) mean(d$err) else NA }
  row
}))
write.csv(summ, file.path(out, "ladder-by-arm.csv"), row.names = FALSE)

# 2. species occupancy error by true occupancy band
ob <- NULL
if (!is.null(occ) && nrow(occ)) {
  occ$band <- cut(occ$occ_true, c(0, .2, .3, .5, .7, .8, 1), include.lowest = TRUE)
  ob <- aggregate(cbind(occ_true, occ_est, error = occ_est - occ_true) ~ scenario + band, occ, mean)
  ob$n <- aggregate(occ_true ~ scenario + band, occ, length)$occ_true
  ob <- ob[order(match(ob$scenario, arms), ob$band), ]
  write.csv(ob, file.path(out, "occupancy-by-band.csv"), row.names = FALSE)
}

# 3. B0 error by true B0 band
b0 <- x[x$block == "B0", ]; b0$band <- cut(b0$truth, c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf))
bb <- aggregate(cbind(error = err, covered) ~ scenario + band, b0, mean)
bb$n <- aggregate(err ~ scenario + band, b0, length)$err
bb <- bb[order(match(bb$scenario, arms), bb$band), ]
write.csv(bb, file.path(out, "b0-by-band.csv"), row.names = FALSE)

# 4. truths identical across arms at each replicate (B0 and B, against the first arm)
pk <- do.call(rbind, lapply(setdiff(intersect(arms, unique(x$scenario)), arms[1]), function(a) do.call(rbind, lapply(c("B0", "B"), function(b) {
  m <- merge(x[x$scenario == arms[1] & x$block == b, c("replicate", "element", "truth")], x[x$scenario == a & x$block == b, c("replicate", "element", "truth")], by = c("replicate", "element"))
  data.frame(arm = a, block = b, n = nrow(m), max_abs_truth_diff = if (nrow(m)) max(abs(m$truth.x - m$truth.y)) else NA)
}))))
if (!is.null(pk)) write.csv(pk, file.path(out, "pairing-across-arms.csv"), row.names = FALSE)

# 5. paired B0 against the 29 July ladder (same seeds, older code)
jp <- NULL
if (file.exists(july_file)) {
  j <- readRDS(july_file)$rows; j <- j[j$block == "B0" & j$scenario %in% arms, ]
  jp <- do.call(rbind, lapply(intersect(unique(j$scenario), unique(x$scenario)), function(a) {
    n <- x[x$scenario == a & x$block == "B0", ]; o <- j[j$scenario == a, ]
    m <- merge(o[, c("replicate", "element", "truth", "lower", "upper", "post_mean", "covered")], n[, c("replicate", "element", "truth", "lower", "upper", "post_mean", "covered")], by = c("replicate", "element"), suffixes = c(".jul", ".now"))
    if (!nrow(m)) return(NULL)
    d <- tapply(m$covered.now - m$covered.jul, m$replicate, mean)
    data.frame(arm = a, R = length(d), max_abs_truth_diff = max(abs(m$truth.now - m$truth.jul)),
      july_coverage = mean(tapply(m$covered.jul, m$replicate, mean)), now_coverage = mean(tapply(m$covered.now, m$replicate, mean)),
      coverage_diff = mean(d), coverage_diff_se = sd(d) / sqrt(length(d)),
      width_ratio = mean(m$upper.now - m$lower.now) / mean(m$upper.jul - m$lower.jul),
      july_slope = unname(coef(lm(post_mean.jul ~ truth.jul, m))[2]), now_slope = unname(coef(lm(post_mean.now ~ truth.now, m))[2]),
      july_mean_error = mean(m$post_mean.jul - m$truth.jul), now_mean_error = mean(m$post_mean.now - m$truth.now))
  }))
  if (!is.null(jp)) write.csv(jp, file.path(out, "b0-paired-july.csv"), row.names = FALSE)
} else cat("July ladder file not found; skipping the paired comparison\n")

options(width = 160)
cat("\n== by arm\n"); print(summ, digits = 3, row.names = FALSE)
if (!is.null(ob)) { cat("\n== species occupancy error by true occupancy\n"); print(ob, digits = 2, row.names = FALSE) }
cat("\n== B0 error by true B0\n"); print(bb, digits = 2, row.names = FALSE)
if (!is.null(pk)) { cat("\n== truths across arms\n"); print(pk, row.names = FALSE) }
if (!is.null(jp)) { cat("\n== B0 paired with the 29 July ladder\n"); print(jp, digits = 3, row.names = FALSE) }
