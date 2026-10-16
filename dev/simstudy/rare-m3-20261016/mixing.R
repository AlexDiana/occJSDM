# Mixing of mu_b0, mu_theta, B0 and the collection intercepts in the refits of replicates 1 to 10 of the
# six M = 3 arms (launch.sh draws), the community-mean arms against the default prior on the same datasets.
# Run from the repository root: Rscript dev/simstudy/rare-m3-20261016/mixing.R
# First checks that every refit reproduces its arm's stored rows exactly.
# Writes reproduction.csv, mixing.csv (one row per element), mixing-summary.csv, identification.csv,
# collection-rate.csv and criteria-mixing.csv (C12 and C13 of README.md).
here <- "dev/simstudy/rare-m3-20261016"
out <- Sys.getenv("OUT_DIR", here)
fits <- lapply(list.files(file.path(here, "draws", "short"), "[.]rds$", full.names = TRUE), readRDS)
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
ck <- list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)
stored <- if (length(res)) readRDS(res[1])$rows else if (length(ck)) readRDS(ck[1])$rows else
  data.frame(scenario = character(0), replicate = integer(0))
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
write.csv(rp, file.path(out, "reproduction.csv"), row.names = FALSE)

diag3 <- function(arr) t(apply(arr, 1, function(x) c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
                                                    ess_tail = posterior::ess_tail(x))))
mx <- do.call(rbind, lapply(fits, function(f) {
  dr <- f$draws
  parts <- list(mu_b0 = if (!is.null(dr$mu_b0)) array(dr$mu_b0, c(1, dim(dr$mu_b0))),
                mu_theta = if (!is.null(dr$mu_theta)) array(dr$mu_theta, c(1, dim(dr$mu_theta))),
                B0 = dr$B0, beta_theta0 = dr$beta_theta0, theta0 = dr$theta0)
  do.call(rbind, lapply(names(parts), function(b) {
    a <- parts[[b]]; if (is.null(a)) return(NULL)
    st <- diag3(a)
    data.frame(arm = f$arm, replicate = f$replicate, block = b, element = seq_len(nrow(st)),
               draws = prod(dim(a)[2:3]), st)
  }))
}))
write.csv(mx, file.path(out, "mixing.csv"), row.names = FALSE)

arms <- c("base_rare_m3", "base_rare_cm_m3", "base_rare_cm_both_m3",
          "occupancy_rare_m3", "occupancy_rare_cm_m3", "occupancy_rare_cm_both_m3")
sm <- do.call(rbind, lapply(split(mx, list(mx$arm, mx$block), drop = TRUE), function(x) {
  data.frame(arm = x$arm[1], block = x$block[1], replicates = length(unique(x$replicate)), elements = nrow(x),
             rhat_median = median(x$rhat), rhat_max = max(x$rhat), share_rhat_over_1.05 = mean(x$rhat > 1.05),
             replicates_rhat_under_1.05 = if (nrow(x) == length(unique(x$replicate))) sum(x$rhat < 1.05) else NA,
             ess_bulk_p10 = unname(quantile(x$ess_bulk, .1)), ess_bulk_median = median(x$ess_bulk),
             ess_tail_median = median(x$ess_tail))
}))
sm <- sm[order(match(sm$block, c("mu_b0", "mu_theta", "B0", "beta_theta0", "theta0")), match(sm$arm, arms)), ]
rownames(sm) <- NULL
write.csv(sm, file.path(out, "mixing-summary.csv"), row.names = FALSE)

# How far the data pin the community means, and their within-fit correlation.
idf <- do.call(rbind, lapply(fits, function(f) {
  if (is.null(f$draws$mu_b0)) return(NULL)
  a <- as.vector(f$draws$mu_b0); b <- if (is.null(f$draws$mu_theta)) NULL else as.vector(f$draws$mu_theta)
  data.frame(arm = f$arm, replicate = f$replicate, mu_b0_sd = sd(a),
             mu_theta_sd = if (is.null(b)) NA else sd(b), cor_mu_b0_mu_theta = if (is.null(b)) NA else cor(a, b))
}))
write.csv(idf, file.path(out, "identification.csv"), row.names = FALSE)

# The collection rate on the probability scale, at the collection covariate's mean.
cr <- do.call(rbind, lapply(fits, function(f) {
  tr <- f$rows[f$rows$block == "beta_theta" & f$rows$element %% 2 == 1, ]
  tr <- tr[order(tr$element), ]
  est <- apply(plogis(f$draws$beta_theta0), 1, mean)
  data.frame(arm = f$arm, replicate = f$replicate, rate_true = mean(plogis(tr$truth)), rate_est = mean(est),
             rate_err = mean(est - plogis(tr$truth)), logit_mean_err = mean(tr$post_mean - tr$truth))
}))
crs <- do.call(rbind, lapply(split(cr, cr$arm), function(x) data.frame(
  arm = x$arm[1], replicates = nrow(x), rate_true = mean(x$rate_true), rate_est = mean(x$rate_est),
  rate_err = mean(x$rate_err), rate_err_se = sd(x$rate_err) / sqrt(nrow(x)), logit_mean_err = mean(x$logit_mean_err))))
crs <- crs[order(match(crs$arm, arms)), ]; rownames(crs) <- NULL
write.csv(crs, file.path(out, "collection-rate.csv"), row.names = FALSE)

default_of <- c(base_rare_cm_m3 = "base_rare_m3", base_rare_cm_both_m3 = "base_rare_m3",
                occupancy_rare_cm_m3 = "occupancy_rare_m3", occupancy_rare_cm_both_m3 = "occupancy_rare_m3")
crit <- do.call(rbind, lapply(names(default_of), function(a) {
  one <- function(id, quantity, default, value, pass, rule)
    data.frame(criterion = id, arm = a, quantity = quantity, default = default, value = value, pass = pass, rule = rule)
  rows <- list()
  for (b in intersect(c("mu_b0", "mu_theta"), unique(mx$block[mx$arm == a]))) {
    u <- mx[mx$arm == a & mx$block == b, ]
    rows[[length(rows) + 1]] <- one("C12", paste("replicates with", b, "R-hat under 1.05"), NA,
                                    sum(u$rhat < 1.05), sum(u$rhat < 1.05) >= 9, "at least 9 of 10")
    rows[[length(rows) + 1]] <- one("C12", paste(b, "median bulk ESS"), NA, median(u$ess_bulk),
                                    median(u$ess_bulk) >= 100, "at least 100")
  }
  for (b in c("B0", "beta_theta0")) {
    cb <- sm[sm$arm == a & sm$block == b, ]; db <- sm[sm$arm == default_of[[a]] & sm$block == b, ]
    rows[[length(rows) + 1]] <- one("C13", paste(b, "median bulk ESS"), db$ess_bulk_median, cb$ess_bulk_median,
                                    cb$ess_bulk_median >= .75 * db$ess_bulk_median, "at least 75% of the default's")
    rows[[length(rows) + 1]] <- one("C13", paste("share of", b, "with R-hat over 1.05"), db$share_rhat_over_1.05,
                                    cb$share_rhat_over_1.05, cb$share_rhat_over_1.05 <= db$share_rhat_over_1.05 + .05,
                                    "no more than 0.05 higher")
  }
  do.call(rbind, rows)
}))
write.csv(crit, file.path(out, "criteria-mixing.csv"), row.names = FALSE)

options(width = 200)
cat("refits:", nrow(rp), " reproducing their stored rows exactly:", sum(rp$identical, na.rm = TRUE),
    " without stored rows yet:", sum(is.na(rp$identical)), "\n")
print(rp[!is.na(rp$identical) & !rp$identical, ], row.names = FALSE)
cat("\n== mixing summary\n"); print(sm, digits = 3, row.names = FALSE)
cat("\n== identification\n")
print(aggregate(cbind(mu_b0_sd, mu_theta_sd, cor_mu_b0_mu_theta) ~ arm, idf, median, na.action = na.pass), digits = 3)
cat("\n== collection rate\n"); print(crs, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
