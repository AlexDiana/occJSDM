# Mixing of mu_b0, B0 and the collection intercepts in the refits of replicates 1 to 10
# (refit_draws.R), community arms against their default counterparts on the same datasets.
# Run from the repository root: Rscript dev/simstudy/b0-community-mean-20261016/mixing.R
# First checks that every refit reproduces its arm's rows exactly: the community arms this
# run's, the default arms the R = 300 paper study's and the rare-dominated run's. A default
# refit that matches is this branch, option off, giving the earlier code's numbers.
# Writes reproduction.csv, mixing.csv (one row per element), mixing-summary.csv and
# criteria-mixing.csv (C12 and C13 of README.md).
here <- "dev/simstudy/b0-community-mean-20261016"
out <- Sys.getenv("OUT_DIR", here)
files <- list.files(file.path(here, "draws"), "[.]rds$", full.names = TRUE)
fits <- lapply(files, readRDS)
res <- list.files(file.path(here, "results"), pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
stored <- rbind(
  readRDS("dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds")$rows,
  readRDS("dev/simstudy/rare-dominated-20261015/results/simstudy-20261016-003851.rds")$rows,
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

diag3 <- function(arr) {
  # arr is [element, niter, nchain]; returns one row per element.
  t(apply(arr, 1, function(x) c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
                                ess_tail = posterior::ess_tail(x))))
}
mx <- do.call(rbind, lapply(fits, function(f) {
  dr <- f$draws
  parts <- list(
    mu_b0 = if (!is.null(dr$mu_b0)) array(dr$mu_b0, c(1, dim(dr$mu_b0))),
    B0 = dr$B0, beta_theta0 = dr$beta_theta0, theta0 = dr$theta0)
  do.call(rbind, lapply(names(parts), function(b) {
    a <- parts[[b]]; if (is.null(a)) return(NULL)
    st <- diag3(a)
    data.frame(arm = f$arm, community = !is.null(dr$mu_b0), replicate = f$replicate, block = b,
               element = seq_len(nrow(st)), draws = prod(dim(a)[2:3]), st)
  }))
}))
write.csv(mx, file.path(out, "mixing.csv"), row.names = FALSE)

sm <- do.call(rbind, lapply(split(mx, list(mx$arm, mx$block), drop = TRUE), function(x) {
  data.frame(arm = x$arm[1], block = x$block[1], replicates = length(unique(x$replicate)), elements = nrow(x),
             rhat_median = median(x$rhat), rhat_max = max(x$rhat),
             share_rhat_over_1.01 = mean(x$rhat > 1.01), share_rhat_over_1.05 = mean(x$rhat > 1.05),
             ess_bulk_p10 = unname(quantile(x$ess_bulk, .1)), ess_bulk_median = median(x$ess_bulk),
             ess_tail_median = median(x$ess_tail))
}))
arms <- c("base", "base_cm", "occupancy", "occupancy_cm", "binary", "binary_cm",
          "base_rare", "base_rare_cm", "occupancy_rare", "occupancy_rare_cm", "binary_rare", "binary_rare_cm")
sm <- sm[order(match(sm$block, c("mu_b0", "B0", "beta_theta0", "theta0")), match(sm$arm, arms)), ]
rownames(sm) <- NULL
write.csv(sm, file.path(out, "mixing-summary.csv"), row.names = FALSE)

pairs <- c(base_cm = "base", occupancy_cm = "occupancy", binary_cm = "binary",
           base_rare_cm = "base_rare", occupancy_rare_cm = "occupancy_rare", binary_rare_cm = "binary_rare")
crit <- do.call(rbind, lapply(names(pairs), function(a) {
  u <- mx[mx$arm == a & mx$block == "mu_b0", ]
  cb <- sm[sm$arm == a & sm$block == "B0", ]; db <- sm[sm$arm == pairs[[a]] & sm$block == "B0", ]
  if (!nrow(u) || !nrow(cb) || !nrow(db)) return(NULL)
  rbind(
    data.frame(criterion = "C12", arm = a, quantity = "replicates with mu_b0 R-hat under 1.05",
               default = NA, community = sum(u$rhat < 1.05), pass = sum(u$rhat < 1.05) >= 9, rule = "at least 9 of 10"),
    data.frame(criterion = "C12", arm = a, quantity = "mu_b0 median bulk ESS",
               default = NA, community = median(u$ess_bulk), pass = median(u$ess_bulk) >= 100, rule = "at least 100"),
    data.frame(criterion = "C13", arm = a, quantity = "B0 median bulk ESS",
               default = db$ess_bulk_median, community = cb$ess_bulk_median,
               pass = cb$ess_bulk_median >= .75 * db$ess_bulk_median, rule = "at least 75% of the default's"),
    data.frame(criterion = "C13", arm = a, quantity = "share of B0 with R-hat over 1.05",
               default = db$share_rhat_over_1.05, community = cb$share_rhat_over_1.05,
               pass = cb$share_rhat_over_1.05 <= db$share_rhat_over_1.05 + .05, rule = "no more than 0.05 higher"))
}))
write.csv(crit, file.path(out, "criteria-mixing.csv"), row.names = FALSE)

options(width = 180)
cat("refits:", nrow(rp), " reproducing their stored rows exactly:", sum(rp$identical, na.rm = TRUE),
    " without stored rows yet:", sum(is.na(rp$identical)), "\n")
print(rp[is.na(rp$identical) | !rp$identical, ], row.names = FALSE)
cat("\n== mixing summary\n"); print(sm, digits = 3, row.names = FALSE)
cat("\n== mu_b0 by replicate\n"); print(mx[mx$block == "mu_b0", ], digits = 3, row.names = FALSE)
cat("\n== criteria\n"); print(crit, digits = 3, row.names = FALSE)
