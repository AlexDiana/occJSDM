# Mixing of mu_b0, mu_theta, B0 and the collection intercepts in the refits of replicates 1 to 10
# (refit.R, launch.sh draws), against the default and B0-only arms on the same datasets, and how far
# the data pin the two community means (C16).
# Run from the repository root: Rscript dev/simstudy/theta-community-mean-20261016/mixing.R
# First checks that every refit reproduces its arm's stored rows exactly: the both-means arms this
# run's, the comparison arms the paper study's, the rare-dominated run's and the community-mean run's.
# A comparison refit that matches is this branch giving the earlier code's numbers.
# The default and B0-only arms' R-hat and ESS are read from the community-mean run's mixing.csv. Their
# draws are not committed; if OLD_DRAWS names the folder that run's refit_draws.R wrote, the
# probability-scale collection rates are computed from them too.
# Writes reproduction.csv, mixing.csv (one row per element), mixing-summary.csv, identification.csv,
# collection-rate.csv and criteria-mixing.csv (C12, C13 and C16 of README.md).
here <- "dev/simstudy/theta-community-mean-20261016"
out <- Sys.getenv("OUT_DIR", here)
old_draws <- Sys.getenv("OLD_DRAWS", "")
files <- list.files(file.path(here, "draws", "short"), "[.]rds$", full.names = TRUE)
fits <- lapply(files, readRDS)
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
stored <- rbind(
  readRDS("dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds")$rows,
  readRDS("dev/simstudy/rare-dominated-20261015/results/simstudy-20261016-003851.rds")$rows,
  readRDS("dev/simstudy/b0-community-mean-20261016/results/simstudy-20261016-110650.rds")$rows,
  if (length(res)) readRDS(res[1])$rows)
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

both <- c("base_cm_both", "occupancy_cm_both", "base_rare_cm_both", "occupancy_rare_cm_both",
          "base_rare_cm_both_wide", "occupancy_rare_cm_both_wide")
fits <- fits[vapply(fits, function(f) f$arm %in% both, logical(1))]

diag3 <- function(arr) {
  # arr is [element, niter, nchain]; returns one row per element.
  t(apply(arr, 1, function(x) c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
                                ess_tail = posterior::ess_tail(x))))
}
mx <- do.call(rbind, lapply(fits, function(f) {
  dr <- f$draws
  parts <- list(mu_b0 = array(dr$mu_b0, c(1, dim(dr$mu_b0))), mu_theta = array(dr$mu_theta, c(1, dim(dr$mu_theta))),
                B0 = dr$B0, beta_theta0 = dr$beta_theta0, theta0 = dr$theta0)
  do.call(rbind, lapply(names(parts), function(b) {
    st <- diag3(parts[[b]])
    data.frame(arm = f$arm, community = TRUE, replicate = f$replicate, block = b,
               element = seq_len(nrow(st)), draws = prod(dim(parts[[b]])[2:3]), st)
  }))
}))
earlier <- read.csv("dev/simstudy/b0-community-mean-20261016/mixing.csv")
earlier <- earlier[earlier$arm %in% c("base", "occupancy", "base_rare", "occupancy_rare",
                                      "base_cm", "occupancy_cm", "base_rare_cm", "occupancy_rare_cm"), ]
earlier$source <- "community-mean run"
mx$source <- "this run"
mx <- rbind(earlier[, names(mx)], mx)
write.csv(mx[mx$source == "this run", setdiff(names(mx), "source")], file.path(out, "mixing.csv"), row.names = FALSE)

arms <- c("base", "base_cm", "base_cm_both", "occupancy", "occupancy_cm", "occupancy_cm_both",
          "base_rare", "base_rare_cm", "base_rare_cm_both", "base_rare_cm_both_wide",
          "occupancy_rare", "occupancy_rare_cm", "occupancy_rare_cm_both", "occupancy_rare_cm_both_wide")
sm <- do.call(rbind, lapply(split(mx, list(mx$arm, mx$block), drop = TRUE), function(x) {
  data.frame(arm = x$arm[1], block = x$block[1], source = x$source[1], replicates = length(unique(x$replicate)),
             elements = nrow(x), rhat_median = median(x$rhat), rhat_max = max(x$rhat),
             share_rhat_over_1.01 = mean(x$rhat > 1.01), share_rhat_over_1.05 = mean(x$rhat > 1.05),
             replicates_rhat_under_1.05 = if (nrow(x) == length(unique(x$replicate))) sum(x$rhat < 1.05) else NA,
             ess_bulk_p10 = unname(quantile(x$ess_bulk, .1)), ess_bulk_median = median(x$ess_bulk),
             ess_tail_median = median(x$ess_tail))
}))
sm <- sm[order(match(sm$block, c("mu_b0", "mu_theta", "B0", "beta_theta0", "theta0")), match(sm$arm, arms)), ]
rownames(sm) <- NULL
write.csv(sm, file.path(out, "mixing-summary.csv"), row.names = FALSE)

# C16: how far the data pin the community means. Posterior sd of mu_theta and mu_b0 per replicate,
# against the hyperprior sd and against the floor sd / sqrt(S + sd^2 / hyper_sd^2) that each would
# have if every species' intercept were known exactly; and the within-fit correlation of the two.
idf <- do.call(rbind, lapply(fits, function(f) {
  a <- as.vector(f$draws$mu_b0); b <- as.vector(f$draws$mu_theta)
  S <- dim(f$draws$B0)[1]; hb <- f$intercept_prior$mu_b0_sd; ht <- f$collection_intercept_prior$mu_theta_sd
  data.frame(arm = f$arm, replicate = f$replicate, hyper_sd_b0 = hb, hyper_sd_theta = ht,
             mu_b0_sd = sd(a), mu_b0_floor = 1 / sqrt(S + 1 / hb^2),
             mu_theta_sd = sd(b), mu_theta_floor = 1 / sqrt(S + 1 / ht^2),
             cor_mu_b0_mu_theta = cor(a, b),
             cor_mean_B0_mean_theta0 = cor(as.vector(apply(f$draws$B0, 2:3, mean)),
                                           as.vector(apply(f$draws$beta_theta0, 2:3, mean))))
}))
write.csv(idf, file.path(out, "identification.csv"), row.names = FALSE)

# The collection rate on the probability scale: the mean over species of the posterior mean of
# plogis(intercept), against the true rate, per replicate. The intercept is on the fit's centred scale,
# so this is the rate at the collection covariate's mean.
rate_rows <- function(f, source) {
  tr <- f$rows[f$rows$block == "beta_theta" & f$rows$element %% 2 == 1, ]
  tr <- tr[order(tr$element), ]
  est <- apply(plogis(f$draws$beta_theta0), 1, mean)
  med <- apply(f$draws$beta_theta0, 1, median)
  data.frame(arm = f$arm, replicate = f$replicate, source = source,
             rate_true = mean(plogis(tr$truth)), rate_est = mean(est), rate_err = mean(est - plogis(tr$truth)),
             logit_mean_err = mean(tr$post_mean - tr$truth), logit_median_err = mean(med - tr$truth))
}
cr <- do.call(rbind, lapply(fits, rate_rows, source = "this run"))
if (nzchar(old_draws)) {
  old <- list.files(old_draws, "^(base|occupancy|base_rare|occupancy_rare)(_cm)?-[0-9]+[.]rds$", full.names = TRUE)
  cr <- rbind(do.call(rbind, lapply(old, function(p) rate_rows(readRDS(p), "community-mean run draws"))), cr)
}
crs <- do.call(rbind, lapply(split(cr, cr$arm), function(x) data.frame(
  arm = x$arm[1], replicates = nrow(x), rate_true = mean(x$rate_true), rate_est = mean(x$rate_est),
  rate_err = mean(x$rate_err), rate_err_se = sd(x$rate_err) / sqrt(nrow(x)),
  logit_mean_err = mean(x$logit_mean_err), logit_median_err = mean(x$logit_median_err))))
crs <- crs[order(match(crs$arm, arms)), ]; rownames(crs) <- NULL
write.csv(crs, file.path(out, "collection-rate.csv"), row.names = FALSE)

default_of <- c(base_cm_both = "base", occupancy_cm_both = "occupancy",
                base_rare_cm_both = "base_rare", occupancy_rare_cm_both = "occupancy_rare")
crit <- do.call(rbind, lapply(names(default_of), function(a) {
  one <- function(id, quantity, default, value, pass, rule)
    data.frame(criterion = id, arm = a, quantity = quantity, default = default, both = value, pass = pass, rule = rule)
  rows <- list()
  for (b in c("mu_b0", "mu_theta")) {
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
  i <- idf[idf$arm == a, ]
  rows[[length(rows) + 1]] <- one("C16", "median posterior sd of mu_theta (hyperprior sd 2)", NA,
                                  median(i$mu_theta_sd), NA, "descriptive; floor about 0.31")
  rows[[length(rows) + 1]] <- one("C16", "median within-fit correlation of mu_b0 and mu_theta", NA,
                                  median(i$cor_mu_b0_mu_theta), NA, "descriptive")
  do.call(rbind, rows)
}))
write.csv(crit, file.path(out, "criteria-mixing.csv"), row.names = FALSE)

options(width = 200)
cat("refits:", nrow(rp), " reproducing their stored rows exactly:", sum(rp$identical, na.rm = TRUE),
    " without stored rows yet:", sum(is.na(rp$identical)), "\n")
print(rp[is.na(rp$identical) | !rp$identical, ], row.names = FALSE)
cat("\n== mixing summary\n"); print(sm, digits = 3, row.names = FALSE)
cat("\n== community means by replicate\n"); print(mx[mx$block %in% c("mu_b0", "mu_theta") & mx$source == "this run", ], digits = 3, row.names = FALSE)
cat("\n== identification\n"); print(idf, digits = 3, row.names = FALSE)
print(aggregate(cbind(mu_b0_sd, mu_theta_sd, cor_mu_b0_mu_theta, cor_mean_B0_mean_theta0) ~ arm, idf, median), digits = 3)
cat("\n== collection rate\n"); print(crs, digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
