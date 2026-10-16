# Do the study's 1,000 + 1,000 chains give the same answers as 5,000 + 5,000 (D3 of README.md)?
# Run from the repository root: Rscript dev/simstudy/ordinary-m3-20261016/long.R
# Compares, replicate by replicate, the long refits (launch.sh long: replicates 1 to 5 of all four M = 3 arms)
# with the 1,000 + 1,000 refits of the same datasets (launch.sh draws), which reproduce the study's rows
# exactly, and repeats C8 and C9's community-minus-default comparison at both chain lengths.
# Writes long-chains.csv (per replicate), long-chains-summary.csv and long-cost.csv.
here <- "dev/simstudy/ordinary-m3-20261016"
out <- Sys.getenv("OUT_DIR", here)
long_files <- list.files(file.path(here, "draws", "long"), "[.]rds$", full.names = TRUE)
stat <- function(x) {
  b0 <- x$rows[x$rows$block == "B0", ]
  bt <- x$rows[x$rows$block == "beta_theta" & x$rows$element %% 2 == 1, ]
  B <- x$rows[x$rows$block == "B", ]
  mb <- x$draws$mu_b0; none <- is.null(mb)
  c(b0_error = mean(b0$post_mean - b0$truth), b0_coverage = mean(b0$covered),
    bt0_error = mean(bt$post_mean - bt$truth), bt0_coverage = mean(bt$covered), B_coverage = mean(B$covered),
    occ_error = mean(x$occupancy$occ_est) - mean(x$occupancy$occ_true),
    mu_b0 = if (none) NA else mean(mb), mu_b0_sd = if (none) NA else sd(as.vector(mb)),
    mu_b0_rhat = if (none) NA else posterior::rhat(mb), mu_b0_ess = if (none) NA else posterior::ess_bulk(mb),
    b0_rhat_max = max(apply(x$draws$B0, 1, posterior::rhat)),
    b0_ess_median = median(apply(x$draws$B0, 1, posterior::ess_bulk)),
    bt0_rhat_max = max(apply(x$draws$beta_theta0, 1, posterior::rhat)),
    bt0_ess_median = median(apply(x$draws$beta_theta0, 1, posterior::ess_bulk)))
}
one <- function(f) {
  L <- readRDS(f)
  S <- readRDS(file.path(here, "draws", "short", sprintf("%s-%02d.rds", L$arm, L$replicate)))
  stopifnot(identical(L$rows$truth, S$rows$truth))
  s <- stat(S); l <- stat(L)
  data.frame(arm = L$arm, replicate = L$replicate, quantity = names(s),
             short = unname(s), long = unname(l), difference = unname(l - s))
}
lc <- do.call(rbind, lapply(long_files, one))
write.csv(lc, file.path(out, "long-chains.csv"), row.names = FALSE)
arms <- c("base_m3", "base_cm_m3", "occupancy_m3", "occupancy_cm_m3")
lc <- lc[!is.na(lc$short) | !is.na(lc$long), ]
ls <- do.call(rbind, lapply(split(lc, list(lc$arm, lc$quantity), drop = TRUE), function(x)
  data.frame(arm = x$arm[1], quantity = x$quantity[1], replicates = nrow(x),
             short_mean = mean(x$short), long_mean = mean(x$long), long_max = max(x$long),
             long_min = min(x$long), mean_difference = mean(x$difference),
             max_abs_difference = max(abs(x$difference)))))
ls <- ls[order(match(ls$arm, arms), ls$quantity), ]; rownames(ls) <- NULL
write.csv(ls, file.path(out, "long-chains-summary.csv"), row.names = FALSE)

# The option's cost (community minus default) on the same replicates at each chain length.
lcost <- do.call(rbind, lapply(c("base", "occupancy"), function(st) do.call(rbind, lapply(
  c("b0_coverage", "b0_error", "occ_error", "bt0_error", "bt0_coverage", "B_coverage"), function(q) {
  a <- lc[lc$arm == paste0(st, "_cm_m3") & lc$quantity == q, ]; r <- lc[lc$arm == paste0(st, "_m3") & lc$quantity == q, ]
  m <- merge(r, a, by = "replicate")
  if (!nrow(m)) return(NULL)
  data.frame(stem = st, quantity = q, replicates = nrow(m),
             default_short = mean(m$short.x), community_short = mean(m$short.y), cost_short = mean(m$short.y - m$short.x),
             default_long = mean(m$long.x), community_long = mean(m$long.y), cost_long = mean(m$long.y - m$long.x))
}))))
write.csv(lcost, file.path(out, "long-cost.csv"), row.names = FALSE)
options(width = 180)
print(ls, digits = 3, row.names = FALSE)
cat("\n== community minus default at each chain length\n"); print(lcost, digits = 3, row.names = FALSE)
