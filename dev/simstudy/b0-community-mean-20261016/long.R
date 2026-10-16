# Do the study's 1,000 + 1,000 chains give the same answers as 5,000 + 5,000?
# Run from the repository root: Rscript dev/simstudy/b0-community-mean-20261016/long.R
# Compares, replicate by replicate, the long refits of refit_long.R (replicates 1 to 5 of the
# four eDNA community arms and the two default rare arms) with the 1,000 + 1,000 refits of the
# same datasets from refit_draws.R, which reproduce the study's rows exactly.
# Writes long-chains.csv (per replicate) and long-chains-summary.csv.
here <- "dev/simstudy/b0-community-mean-20261016"
out <- Sys.getenv("OUT_DIR", here)
long_files <- list.files(file.path(here, "draws", "long"), "[.]rds$", full.names = TRUE)
one <- function(f) {
  L <- readRDS(f)
  S <- readRDS(file.path(here, "draws", sprintf("%s-%02d.rds", L$arm, L$replicate)))
  stopifnot(identical(L$rows$truth, S$rows$truth))
  stat <- function(x) {
    b0 <- x$rows[x$rows$block == "B0", ]
    bt <- x$rows[x$rows$block == "beta_theta" & x$rows$element %% 2 == 1, ]
    mu <- x$draws$mu_b0
    c(b0_error = mean(b0$post_mean - b0$truth), b0_coverage = mean(b0$covered),
      bt0_error = mean(bt$post_mean - bt$truth),
      occ_error = mean(x$occupancy$occ_est) - mean(x$occupancy$occ_true),
      mu_b0 = if (is.null(mu)) NA else mean(mu),
      mu_b0_rhat = if (is.null(mu)) NA else posterior::rhat(mu),
      mu_b0_ess = if (is.null(mu)) NA else posterior::ess_bulk(mu),
      b0_rhat_max = max(apply(x$draws$B0, 1, posterior::rhat)),
      b0_ess_median = median(apply(x$draws$B0, 1, posterior::ess_bulk)))
  }
  s <- stat(S); l <- stat(L)
  data.frame(arm = L$arm, replicate = L$replicate, quantity = names(s),
             short = unname(s), long = unname(l), difference = unname(l - s))
}
lc <- do.call(rbind, lapply(long_files, one))
write.csv(lc, file.path(out, "long-chains.csv"), row.names = FALSE)
arms <- c("base_rare", "base_rare_cm", "occupancy_rare", "occupancy_rare_cm", "base_cm", "occupancy_cm")
ls <- do.call(rbind, lapply(split(lc, list(lc$arm, lc$quantity), drop = TRUE), function(x)
  data.frame(arm = x$arm[1], quantity = x$quantity[1], replicates = nrow(x),
             short_mean = mean(x$short), long_mean = mean(x$long),
             mean_difference = mean(x$difference), max_abs_difference = max(abs(x$difference)))))
ls <- ls[order(match(ls$arm, arms), ls$quantity), ]; rownames(ls) <- NULL
write.csv(ls, file.path(out, "long-chains-summary.csv"), row.names = FALSE)
options(width = 160)
print(ls, digits = 3, row.names = FALSE)
