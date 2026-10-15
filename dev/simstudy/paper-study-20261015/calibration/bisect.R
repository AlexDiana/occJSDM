# Summarise the bisection refits behind REPORT.md.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/paper-study-20261015/calibration/bisect.R \
#     <bisect directory> <July rds> <August 10 rds>
#
# The bisection directory holds one out-<sha>/ results folder per refitted
# commit, written by run-bisect.sh (copied beside this script). Every refit
# uses the production seeds, so replicate r has the same simulated data and
# truth as replicate r of the July, August and October runs. Only the code
# and the MCMC stream differ.

args <- commandArgs(trailingOnly = TRUE)
bis_dir <- args[1]
out_dir <- "dev/simstudy/paper-study-20261015/calibration"

read_rows <- function(f) {
  d <- readRDS(f)$rows
  d$width <- d$upper - d$lower
  d
}
latest <- function(dir) {
  f <- list.files(dir, pattern = "^simstudy-[0-9-]+[.]rds$", full.names = TRUE)
  if (!length(f)) return(NULL)
  sort(f, decreasing = TRUE)[1]
}

runs <- list(
  jul29 = read_rows(args[2]),
  aug10 = read_rows(args[3]))
# Main: 5382460 has the collection-alignment fix and serial RNG steps,
# 9a784a9 adds trait standardisation, cc1b4e6 merges the spatial PR (PR 8).
# Spatial branch: 54b5cca is its base (alignment fix and serial RNG steps, no
# spatial change), f1d2e12 changes the spatial basis, prior mean and range
# target, e9fd2d3 blocks the range and coefficient draws.
shas <- c("5382460", "9a784a9", "cc1b4e6", "54b5cca", "f1d2e12", "e9fd2d3")
for (s in shas) {
  f <- latest(file.path(bis_dir, paste0("out-", s)))
  if (!is.null(f)) runs[[paste0("c", s)]] <- read_rows(f)
}
runs$oct15 <- read_rows("dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds")
# The study's code refitted with the experimental half-Cauchy spatial
# amplitude prior (halfcauchy.R), on the same replicates.
hc <- file.path(bis_dir, "halfcauchy.rds")
if (file.exists(hc)) runs$oct15_halfcauchy <- read_rows(hc)

# Restrict every run to the replicates and cells the refits cover.
reps <- sort(unique(runs[[paste0("c", shas[1])]]$replicate))
cells <- c("base", "binary")
runs <- lapply(runs, function(d) d[d$replicate %in% reps & d$scenario %in% cells, ])

key <- function(d) paste(d$scenario, d$replicate, d$block, d$element)
ref <- runs$oct15

out <- do.call(rbind, lapply(names(runs), function(nm) {
  d <- runs[[nm]]
  # Score every run against the October truth: identical for B, B0 and
  # theta0, and the corrected scoring units for beta_theta.
  d$truth_oct <- ref$truth[match(key(d), key(ref))]
  d$cov_oct <- as.integer(d$lower <= d$truth_oct & d$truth_oct <= d$upper)
  do.call(rbind, lapply(cells, function(sc) {
    do.call(rbind, lapply(c("B0", "B", "beta_theta", "theta0"), function(bl) {
      x <- d[d$scenario == sc & d$block == bl, ]
      if (!nrow(x)) return(NULL)
      c_r <- tapply(x$cov_oct, x$replicate, mean)
      data.frame(run = nm, scenario = sc, block = bl, R = length(c_r),
                 coverage = mean(c_r), se = sd(c_r) / sqrt(length(c_r)),
                 mean_width = mean(x$width),
                 bias = mean(x$post_mean - x$truth_oct),
                 rmse = sqrt(mean((x$post_mean - x$truth_oct)^2)),
                 slope = unname(coef(lm(post_mean ~ truth_oct, x))[2]),
                 stringsAsFactors = FALSE)
    }))
  }))
}))
rownames(out) <- NULL
write.csv(out, file.path(out_dir, "bisect-base-binary.csv"), row.names = FALSE)
print(out, digits = 3)
