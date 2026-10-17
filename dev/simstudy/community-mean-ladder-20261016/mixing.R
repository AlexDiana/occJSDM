# Mixing of mu_b0, B0 and the collection intercepts in the refits of replicates 1 to 10 of the six community-mean
# arms (launch.sh draws), set beside the M = 2 and M = 3 arms' values from the earlier runs' mixing tables; and,
# if the long-chain step ran, how far each refitted replicate moves between 1,000 + 1,000 and 5,000 + 5,000.
# Run from the repository root: Rscript dev/simstudy/community-mean-ladder-20261016/mixing.R
# First checks that every refit reproduces its arm's stored rows exactly.
# Writes reproduction.csv, mixing.csv (one row per element), mixing-summary.csv, mu-b0-by-m.csv and, if the long
# chains ran, long-chains.csv and long-chains-summary.csv.
here <- "dev/simstudy/community-mean-ladder-20261016"
out <- Sys.getenv("OUT_DIR", here)
fits <- lapply(list.files(file.path(here, "draws", "short"), "[.]rds$", full.names = TRUE), readRDS)
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
ck <- list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)
main <- if (length(res)) readRDS(res[1]) else if (length(ck)) readRDS(ck[1]) else NULL
stored <- if (is.null(main)) data.frame(scenario = character(0), replicate = integer(0)) else
  if (!is.null(main$rows)) main$rows else main$params
cols <- c("truth", "lower", "upper", "post_mean", "covered")

rp <- do.call(rbind, lapply(fits, function(f) {
  ref <- stored[stored$scenario == f$arm & stored$replicate == f$replicate, ]
  if (!nrow(ref)) return(data.frame(arm = f$arm, replicate = f$replicate, rows = nrow(f$rows),
                                    stored_rows = 0L, identical = NA))
  m <- merge(f$rows, ref, by = c("block", "element"))
  same <- nrow(m) == nrow(ref) && nrow(m) == nrow(f$rows) &&
    all(vapply(cols, function(k) identical(m[[paste0(k, ".x")]], m[[paste0(k, ".y")]]), logical(1)))
  data.frame(arm = f$arm, replicate = f$replicate, rows = nrow(f$rows), stored_rows = nrow(ref), identical = same)
}))
rp <- rp[order(rp$arm, rp$replicate), ]
write.csv(rp, file.path(out, "reproduction.csv"), row.names = FALSE)

# Rank-normalised R-hat and bulk and tail ESS per element. B0 is S x niter x nchain, mu_b0 niter x nchain.
diag3 <- function(arr) t(apply(arr, 1, function(x) c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
                                                    ess_tail = posterior::ess_tail(x))))
mx_of <- function(fits) do.call(rbind, lapply(fits, function(f) {
  dr <- f$draws
  parts <- list(mu_b0 = if (!is.null(dr$mu_b0)) array(dr$mu_b0, c(1, dim(dr$mu_b0))),
                B0 = dr$B0, beta_theta0 = dr$beta_theta0)
  do.call(rbind, lapply(names(parts), function(b) {
    a <- parts[[b]]; if (is.null(a)) return(NULL)
    st <- diag3(a)
    data.frame(arm = f$arm, replicate = f$replicate, block = b, element = seq_len(nrow(st)),
               draws = prod(dim(a)[2:3]), st)
  }))
}))
mx <- mx_of(fits)
write.csv(mx, file.path(out, "mixing.csv"), row.names = FALSE)

# The M = 2 and M = 3 community-mean arms, from the earlier runs' refits at the same chain length.
earlier <- list(
  b0cm = "dev/simstudy/b0-community-mean-20261016/mixing.csv",
  om3 = "dev/simstudy/ordinary-m3-20261016/mixing.csv",
  rm3 = "dev/simstudy/rare-m3-20261016/mixing.csv")
want <- c("base_cm", "base_rare_cm", "base_cm_m3", "base_rare_cm_m3")
old <- do.call(rbind, lapply(earlier, function(f) {
  x <- read.csv(f); x[x$arm %in% want & x$block %in% c("mu_b0", "B0", "beta_theta0"), names(mx)]
}))
arms <- data.frame(
  arm = c("base_cm", "base_cm_m3", "M5_cm", "M10_cm", "M20_cm",
          "base_rare_cm", "base_rare_cm_m3", "base_rare_cm_m5", "base_rare_cm_m10", "base_rare_cm_m20"),
  community = rep(c("ordinary", "rare"), each = 5), M = rep(c(2, 3, 5, 10, 20), 2),
  source = rep(c("b0-community-mean-20261016", "ordinary-m3-20261016", rep("this run", 3),
                 "b0-community-mean-20261016", "rare-m3-20261016", rep("this run", 3)), 1))
summ <- function(x) do.call(rbind, lapply(split(x, list(x$arm, x$block), drop = TRUE), function(y) {
  data.frame(arm = y$arm[1], block = y$block[1], replicates = length(unique(y$replicate)), elements = nrow(y),
             rhat_median = median(y$rhat), rhat_max = max(y$rhat), share_rhat_over_1.05 = mean(y$rhat > 1.05),
             replicates_rhat_under_1.05 = if (nrow(y) == length(unique(y$replicate))) sum(y$rhat < 1.05) else NA,
             ess_bulk_p10 = unname(quantile(y$ess_bulk, .1)), ess_bulk_median = median(y$ess_bulk),
             ess_tail_median = median(y$ess_tail))
}))
sm <- merge(arms, summ(rbind(old, mx)), by = "arm")
sm <- sm[order(match(sm$block, c("mu_b0", "B0", "beta_theta0")), sm$community, sm$M), ]; rownames(sm) <- NULL
write.csv(sm, file.path(out, "mixing-summary.csv"), row.names = FALSE)
mu <- sm[sm$block == "mu_b0", ]
mu$converged <- mu$replicates_rhat_under_1.05 >= 9 & mu$ess_bulk_median >= 100
write.csv(mu[, c("community", "M", "arm", "source", "replicates", "replicates_rhat_under_1.05", "rhat_median",
                 "rhat_max", "ess_bulk_median", "ess_bulk_p10", "ess_tail_median", "converged")],
          file.path(out, "mu-b0-by-m.csv"), row.names = FALSE)

# The long chains: each refitted replicate's errors at 5,000 + 5,000 against the main run's at 1,000 + 1,000.
long_files <- list.files(file.path(here, "draws", "long"), "[.]rds$", full.names = TRUE)
if (length(long_files) && !is.null(main)) {
  lf <- lapply(long_files, readRDS)
  occ_main <- main$occupancy
  errs <- function(rows, occ) {
    b0 <- rows[rows$block == "B0", ]; bt <- rows[rows$block == "beta_theta" & rows$element %% 2 == 1, ]
    u <- rows[rows$block == "mu_b0", ]
    data.frame(occ = mean(occ$occ_est) - mean(occ$occ_true), mae = mean(abs(occ$occ_est - occ$occ_true)),
               slope = cov(occ$occ_true, occ$occ_est) / var(occ$occ_true), b0_cov = mean(b0$covered),
               b0_width = mean(b0$upper - b0$lower), bt0 = mean(bt$post_mean - bt$truth),
               mu_b0 = if (nrow(u)) u$post_mean else NA)
  }
  lc <- do.call(rbind, lapply(lf, function(f) {
    s <- errs(stored[stored$scenario == f$arm & stored$replicate == f$replicate, ],
              occ_main[occ_main$scenario == f$arm & occ_main$replicate == f$replicate, ])
    l <- errs(f$rows, f$occupancy)
    a <- f$draws$mu_b0
    data.frame(arm = f$arm, replicate = f$replicate, setNames(s, paste0(names(s), "_short")),
               setNames(l, paste0(names(l), "_long")), mu_b0_rhat_long = posterior::rhat(a),
               mu_b0_ess_bulk_long = posterior::ess_bulk(a))
  }))
  lc <- lc[order(lc$arm, lc$replicate), ]
  write.csv(lc, file.path(out, "long-chains.csv"), row.names = FALSE)
  qn <- c("occ", "mae", "slope", "b0_cov", "b0_width", "bt0", "mu_b0")
  lsum <- do.call(rbind, lapply(split(lc, lc$arm), function(x) do.call(rbind, lapply(qn, function(q) {
    dd <- x[[paste0(q, "_long")]] - x[[paste0(q, "_short")]]
    data.frame(arm = x$arm[1], quantity = q, replicates = nrow(x), short = mean(x[[paste0(q, "_short")]]),
               long = mean(x[[paste0(q, "_long")]]), diff = mean(dd), diff_se = sd(dd) / sqrt(nrow(x)),
               max_abs_diff = max(abs(dd)))
  }))))
  lsum <- rbind(lsum, do.call(rbind, lapply(split(lc, lc$arm), function(x) data.frame(
    arm = x$arm[1], quantity = c("mu_b0 R-hat (long)", "mu_b0 bulk ESS (long)"), replicates = nrow(x),
    short = NA, long = c(max(x$mu_b0_rhat_long), median(x$mu_b0_ess_bulk_long)), diff = NA, diff_se = NA,
    max_abs_diff = NA))))
  lsum <- lsum[order(match(lsum$arm, arms$arm)), ]; rownames(lsum) <- NULL
  write.csv(lsum, file.path(out, "long-chains-summary.csv"), row.names = FALSE)
}

options(width = 200)
cat("refits:", nrow(rp), " reproducing their stored rows exactly:", sum(rp$identical, na.rm = TRUE),
    " without stored rows yet:", sum(is.na(rp$identical)), "\n")
print(rp[!is.na(rp$identical) & !rp$identical, ], row.names = FALSE)
cat("\n== mu_b0 by M\n"); print(mu, digits = 3, row.names = FALSE)
cat("\n== mixing summary\n"); print(sm, digits = 3, row.names = FALSE)
if (exists("lsum")) { cat("\n== long chains\n"); print(lsum, digits = 3, row.names = FALSE) }
