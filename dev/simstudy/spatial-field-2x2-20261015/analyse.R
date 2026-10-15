# Summarise the spatial 2x2 against base replicates 1 to 100 of the R = 300 paper study.
# Run from the repository root: Rscript dev/simstudy/spatial-field-2x2-20261015/analyse.R
# Reads the run's results file, or its checkpoint while the run is in flight.
# Writes coverage-by-arm.csv, paired-contrasts.csv, b-by-effect-size.csv and pairing-check.csv.
here <- "dev/simstudy/spatial-field-2x2-20261015"
out <- Sys.getenv("OUT_DIR", here)  # where the CSVs go; the run directory by default
base_file <- "dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds"
res <- list.files(file.path(here, "results"), pattern = "^simstudy-.*[.]rds$", full.names = TRUE)
new <- if (length(res)) readRDS(res[1])$rows else
  readRDS(list.files(file.path(here, "results"), pattern = "^checkpoint", full.names = TRUE)[1])$params
base <- readRDS(base_file)$rows
base <- base[base$scenario == "base" & base$replicate <= 100, ]
cols <- c("scenario", "replicate", "block", "element", "truth", "lower", "upper", "post_mean", "covered")
d <- rbind(base[, cols], new[new$block != "ERROR", cols])
d$width <- d$upper - d$lower; d$err <- d$post_mean - d$truth
arms <- c("base", "field_unfitted", "nofield_fitted", "no_spatial")
d <- d[d$scenario %in% arms, ]
reps <- Reduce(intersect, lapply(arms, function(a) unique(d$replicate[d$scenario == a])))
d <- d[d$replicate %in% reps, ]
cat(sprintf("replicates present in all four arms: %d\n", length(reps)))
blocks <- c("B", "B0", "G", "beta_theta", "theta0", "p", "q")
zq <- 2 * qnorm(0.975)

per_rep <- function(x, v) aggregate(x[[v]], list(replicate = x$replicate), mean)

# 1. Coverage and shape of the intervals by arm and block
cov_rows <- do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(blocks, function(b) {
  x <- d[d$scenario == a & d$block == b, ]; if (!nrow(x)) return(NULL)
  c_r <- per_rep(x, "covered")$x
  data.frame(arm = a, block = b, R = length(c_r), coverage = mean(c_r), se = sd(c_r) / sqrt(length(c_r)),
             mean_width = mean(x$width), rmse = sqrt(mean(x$err^2)),
             slope = unname(coef(lm(post_mean ~ truth, x))[2]),
             err_to_width = sqrt(mean((x$err / (x$width / zq))^2)))
}))))
write.csv(cov_rows, file.path(out, "coverage-by-arm.csv"), row.names = FALSE)

# 2. Paired contrasts on per-replicate coverage, width and error
contrasts <- list(
  c("fitting a spatial term, field present", "base", "field_unfitted"),
  c("fitting a spatial term, no field", "nofield_fitted", "no_spatial"),
  c("the field, spatial term fitted", "base", "nofield_fitted"),
  c("the field, no spatial term", "field_unfitted", "no_spatial"))
pc <- do.call(rbind, lapply(c("B", "B0"), function(b) do.call(rbind, lapply(contrasts, function(k) {
  a1 <- d[d$scenario == k[2] & d$block == b, ]; a0 <- d[d$scenario == k[3] & d$block == b, ]
  m <- merge(per_rep(a1, "covered"), per_rep(a0, "covered"), by = "replicate")
  w <- merge(per_rep(a1, "width"), per_rep(a0, "width"), by = "replicate")
  dd <- m$x.x - m$x.y
  data.frame(block = b, contrast = k[1], minuend = k[2], subtrahend = k[3], R = nrow(m),
             coverage_diff = mean(dd), coverage_diff_se = sd(dd) / sqrt(nrow(m)),
             width_ratio = mean(w$x.x) / mean(w$x.y),
             slope_minuend = unname(coef(lm(post_mean ~ truth, a1))[2]),
             slope_subtrahend = unname(coef(lm(post_mean ~ truth, a0))[2]))
}))))
# interaction: does fitting the spatial term matter more when a field is present?
inter <- do.call(rbind, lapply(c("B", "B0"), function(b) {
  g <- function(a) { r <- per_rep(d[d$scenario == a & d$block == b, ], "covered"); setNames(r$x, r$replicate) }
  v <- (g("base") - g("field_unfitted")) - (g("nofield_fitted") - g("no_spatial"))
  data.frame(block = b, contrast = "interaction", minuend = "(base - field_unfitted)", subtrahend = "(nofield_fitted - no_spatial)",
             R = length(v), coverage_diff = mean(v), coverage_diff_se = sd(v) / sqrt(length(v)),
             width_ratio = NA, slope_minuend = NA, slope_subtrahend = NA)
}))
pc <- rbind(pc, inter)
write.csv(pc, file.path(out, "paired-contrasts.csv"), row.names = FALSE)

# 3. B coverage by size of the true effect
bins <- c(0, 1, 2, 3, Inf); lab <- c("under 1", "1 to 2", "2 to 3", "over 3")
bb <- d[d$block == "B", ]; bb$bin <- cut(abs(bb$truth), bins, labels = lab, right = FALSE)
eb <- aggregate(cbind(covered, pull = abs(truth) - abs(post_mean)) ~ scenario + bin, bb, mean)
eb$n <- aggregate(covered ~ scenario + bin, bb, length)$covered
names(eb)[names(eb) == "covered"] <- "coverage"
write.csv(eb[order(eb$scenario, eb$bin), ], file.path(out, "b-by-effect-size.csv"), row.names = FALSE)

# 4. Pairing: are the truths element-for-element the same as base's?
pk <- do.call(rbind, lapply(setdiff(arms, "base"), function(a) do.call(rbind, lapply(blocks, function(b) {
  x <- d[d$scenario == "base" & d$block == b, ]; y <- d[d$scenario == a & d$block == b, ]
  m <- merge(x, y, by = c("replicate", "element"))
  if (!nrow(m)) return(NULL)
  data.frame(arm = a, block = b, n = nrow(m), max_abs_truth_diff = max(abs(m$truth.x - m$truth.y)),
             identical_post_means = mean(m$post_mean.x == m$post_mean.y))
}))))
write.csv(pk, file.path(out, "pairing-check.csv"), row.names = FALSE)

options(width = 150)
cat("\n== coverage by arm (B and B0)\n"); print(cov_rows[cov_rows$block %in% c("B", "B0"), ], digits = 3, row.names = FALSE)
cat("\n== paired contrasts\n"); print(pc, digits = 3, row.names = FALSE)
cat("\n== B by effect size\n"); print(eb[order(eb$bin, eb$scenario), ], digits = 3, row.names = FALSE)
cat("\n== pairing check\n"); print(pk, digits = 3, row.names = FALSE)
