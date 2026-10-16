# Do the study's 1,000 + 1,000 chains give the same answers as 5,000 + 5,000, and do mu_b0 and
# mu_theta converge with longer chains?
# Run from the repository root: Rscript dev/simstudy/theta-community-mean-20261016/long.R
# Compares, replicate by replicate, the long refits (launch.sh long: replicates 1 to 5 of the six
# both-means arms) with the 1,000 + 1,000 refits of the same datasets (launch.sh draws), which
# reproduce the study's rows exactly.
# Writes long-chains.csv (per replicate) and long-chains-summary.csv.
here <- "dev/simstudy/theta-community-mean-20261016"
out <- Sys.getenv("OUT_DIR", here)
long_files <- list.files(file.path(here, "draws", "long"), "[.]rds$", full.names = TRUE)
one <- function(f) {
  L <- readRDS(f)
  S <- readRDS(file.path(here, "draws", "short", sprintf("%s-%02d.rds", L$arm, L$replicate)))
  stopifnot(identical(L$rows$truth, S$rows$truth))
  stat <- function(x) {
    b0 <- x$rows[x$rows$block == "B0", ]
    bt <- x$rows[x$rows$block == "beta_theta" & x$rows$element %% 2 == 1, ]
    mb <- x$draws$mu_b0; mt <- x$draws$mu_theta
    c(b0_error = mean(b0$post_mean - b0$truth), b0_coverage = mean(b0$covered),
      bt0_error = mean(bt$post_mean - bt$truth), bt0_coverage = mean(bt$covered),
      occ_error = mean(x$occupancy$occ_est) - mean(x$occupancy$occ_true),
      mu_b0 = mean(mb), mu_b0_sd = sd(as.vector(mb)),
      mu_b0_rhat = posterior::rhat(mb), mu_b0_ess = posterior::ess_bulk(mb),
      mu_theta = mean(mt), mu_theta_sd = sd(as.vector(mt)),
      mu_theta_rhat = posterior::rhat(mt), mu_theta_ess = posterior::ess_bulk(mt),
      mu_cor = cor(as.vector(mb), as.vector(mt)),
      b0_rhat_max = max(apply(x$draws$B0, 1, posterior::rhat)),
      b0_ess_median = median(apply(x$draws$B0, 1, posterior::ess_bulk)),
      bt0_rhat_max = max(apply(x$draws$beta_theta0, 1, posterior::rhat)),
      bt0_ess_median = median(apply(x$draws$beta_theta0, 1, posterior::ess_bulk)))
  }
  s <- stat(S); l <- stat(L)
  data.frame(arm = L$arm, replicate = L$replicate, quantity = names(s),
             short = unname(s), long = unname(l), difference = unname(l - s))
}
lc <- do.call(rbind, lapply(long_files, one))
write.csv(lc, file.path(out, "long-chains.csv"), row.names = FALSE)
arms <- c("base_cm_both", "occupancy_cm_both", "base_rare_cm_both", "base_rare_cm_both_wide",
          "occupancy_rare_cm_both", "occupancy_rare_cm_both_wide")
ls <- do.call(rbind, lapply(split(lc, list(lc$arm, lc$quantity), drop = TRUE), function(x)
  data.frame(arm = x$arm[1], quantity = x$quantity[1], replicates = nrow(x),
             short_mean = mean(x$short), long_mean = mean(x$long), long_max = max(x$long),
             long_min = min(x$long), mean_difference = mean(x$difference),
             max_abs_difference = max(abs(x$difference)))))
ls <- ls[order(match(ls$arm, arms), ls$quantity), ]; rownames(ls) <- NULL
write.csv(ls, file.path(out, "long-chains-summary.csv"), row.names = FALSE)
options(width = 160)
print(ls, digits = 3, row.names = FALSE)
