# Pre-run checks for the range-inside-grid diagnostic cells. No fitting.
#
# Usage, from the repository root:
#   Rscript dev/simstudy/range-inside-grid-20261015/check-pairing.R
#
# 1. draw_truth() returns identical output to the pre-change harness
#    (015d1dd) for every production cell, replicates 1 to 50.
# 2. Under the runner's L'Ecuyer-CMRG generator, the simulated B0 and B truths
#    of `base` and `binary` at replicates 1 to 3 match the R = 300 production
#    rows, so the seeds reproduce the production datasets.
# 3. Each diagnostic cell draws the same coordinates and JSDM coefficients as
#    its parent at replicates 1 to 50, and a different spatial field.
# 4. The simulated range on the fitter's standardised coordinates (raw range
#    over each axis's site SD), replicates 1 to 50, against the grid's 0.30.
#
# Writes standardised-range.csv beside this script.

run_dir <- "dev/simstudy/range-inside-grid-20261015"
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-simstudy.R")
new_env <- environment()

old_env <- new.env()
old_file <- tempfile(fileext = ".R")
writeLines(system2("git", c("show", "015d1dd:tests/testthat/helper-simstudy.R"),
                   stdout = TRUE), old_file)
sys.source(old_file, envir = old_env)

RNGkind("L'Ecuyer-CMRG")
sc <- simstudy_scenarios()
names(sc) <- vapply(sc, `[[`, "", "label")
production <- c("base", "binary", "d_overfit", "d_underfit", "low_information",
                "occupancy", "primers_3", "spatial_isolated", "species_20",
                "traits_isolated")
# Every cell that existed before is unchanged; the new list only appends.
old_sc <- old_env$simstudy_scenarios()
stopifnot(identical(old_sc, unname(sc[seq_along(old_sc)])),
          all(production %in% names(sc)[seq_along(old_sc)]))

# 1. Production truths unchanged.
for (lab in production) for (r in 1:50) {
  s <- sc[[lab]]
  seed <- simstudy_seed_for(s, r)
  if (!identical(old_env$draw_truth(s, seed), draw_truth(s, seed)))
    stop("draw_truth() changed for ", lab, " replicate ", r)
}
cat("1. draw_truth() identical to 015d1dd for 10 production cells x 50 replicates\n")

# 2. Seeds reproduce production datasets.
prod <- readRDS("dev/simstudy/paper-study-20261015/results/simstudy-20261015-172834.rds")$rows
for (lab in c("base", "binary")) for (r in 1:3) {
  s <- sc[[lab]]
  tp <- simstudy_simulate(draw_truth(s, simstudy_seed_for(s, r)))$true_params
  for (bl in c("B0", "B")) {
    p <- prod[prod$scenario == lab & prod$replicate == r & prod$block == bl, ]
    p <- p[order(p$element), ]
    d <- max(abs(p$truth - as.vector(tp$jsdmParams_true[[bl]])))
    if (!(d < 1e-12)) stop("production truth mismatch: ", lab, " ", r, " ", bl, " ", d)
  }
}
cat("2. base and binary B0, B truths match production rows, replicates 1 to 3\n")

# 3 and 4. Pairing of the diagnostic cells, and standardised ranges.
coords <- function(sim) {
  # Two-stage and occupancy data have one info row per sample; binary has
  # one per site.
  info <- sim$data_list$info
  if ("Site" %in% names(info)) info <- info[!duplicated(info$Site), ]
  as.matrix(info[c("Xs.1", "Xs.2")])
}
rows <- list()
for (lab in c("base_ls05", "binary_ls05", "base_ls075", "binary_ls075")) {
  s <- sc[[lab]]
  parent <- sc[[s$seed_label]]
  for (r in 1:50) {
    a <- simstudy_simulate(draw_truth(s, simstudy_seed_for(s, r)))
    b <- simstudy_simulate(draw_truth(parent, simstudy_seed_for(parent, r)))
    ja <- a$true_params$jsdmParams_true; jb <- b$true_params$jsdmParams_true
    if (!identical(ja$l_s, s$l_s)) stop(lab, " ", r, ": simulated l_s is not the cell's")
    for (nm in c("B0", "B", "G", "A", "C", "L", "U"))
      if (!identical(ja[[nm]], jb[[nm]])) stop(lab, " ", r, ": ", nm, " not paired")
    xa <- coords(a); xb <- coords(b)
    if (!identical(xa, xb)) stop(lab, " ", r, ": coordinates not paired")
    if (isTRUE(all.equal(ja$spatField, jb$spatField)))
      stop(lab, " ", r, ": spatial field unchanged")
    sds <- apply(xa, 2, sd)
    rows[[length(rows) + 1]] <- data.frame(
      scenario = lab, replicate = r, l_s_raw = s$l_s,
      std_axis1 = s$l_s / sds[1], std_axis2 = s$l_s / sds[2],
      parent_std_axis1 = 0.15 / sds[1], parent_std_axis2 = 0.15 / sds[2],
      field_sd = sd(as.vector(ja$spatField)),
      parent_field_sd = sd(as.vector(jb$spatField)))
  }
}
cat("3. diagnostic cells pair with their parents (B0, B, G, A, C, L, U, coordinates), field differs, replicates 1 to 50\n")
rng <- do.call(rbind, rows)
rownames(rng) <- NULL
write.csv(rng, file.path(run_dir, "standardised-range.csv"), row.names = FALSE)
cat("4. standardised range (min to max over both axes and 50 replicates):\n")
for (lab in unique(rng$scenario)) {
  x <- rng[rng$scenario == lab, ]
  v <- c(x$std_axis1, x$std_axis2); pv <- c(x$parent_std_axis1, x$parent_std_axis2)
  cat(sprintf("   %-12s raw %.3f -> %.3f to %.3f (parent 0.15 -> %.3f to %.3f); field sd %.3f (parent %.3f)\n",
              lab, x$l_s_raw[1], min(v), max(v), min(pv), max(pv),
              mean(x$field_sd), mean(x$parent_field_sd)))
}
