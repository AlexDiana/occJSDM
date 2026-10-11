# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-plot.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("config.R", "score.R", "plot.R")) source(file.path(D, f))
cfg <- stage2_config()
set.seed(5); rows <- list()
for (cell in cfg$cells) for (rep in 1:10) for (k in seq_len(nrow(cfg$arms))) {
  key <- arm_key(cfg$arms$package[k], cfg$arms$arm[k])
  rows[[length(rows) + 1L]] <- data.frame(dataset_id = sprintf("%s-rep%d", cell, rep), cell = cell, package = cfg$arms$package[k], arm = cfg$arms$arm[k], key = key,
    headroom_share = c(0, 0.3, 0.5, 0.45, 0.4)[k] + rnorm(1, 0, 0.1) - (if (key == "occJSDM/gp" && cell == "t2-n900" && rep == 1) 3 else 0),
    status = if (key == "gllvm/spatial2" && rep > 6) "unstable" else "pass", stringsAsFactors = FALSE)
}
agg <- do.call(rbind, rows)
agg$pattern <- NA_real_; agg$field_slope <- NA_real_; agg$field_correlation <- NA_real_; agg$rmse <- NA_real_; agg$rmse_env <- NA_real_; agg$rmse_oracle <- NA_real_; agg$skipped_species <- 0L
pairs <- paired_comparisons(agg, cfg)

test_that("both figures are drawn to non-empty PNG files", {
  f1 <- tempfile(fileext = ".png"); f2 <- tempfile(fileext = ".png")
  plot_headroom_share(agg, cfg, f1); plot_pairs(pairs, cfg, f2)
  expect_true(file.exists(f1) && file.size(f1) > 10000); expect_true(file.exists(f2) && file.size(f2) > 10000)
  keep <- Sys.getenv("STAGE2_PLOT_PREVIEW")
  if (nzchar(keep)) { file.copy(f1, file.path(keep, "headroom-share-preview.png"), overwrite = TRUE); file.copy(f2, file.path(keep, "pairs-preview.png"), overwrite = TRUE) }
})

test_that("results are copied with a runtime record", {
  study <- file.path(tempdir(), "plot-study"); dir.create(file.path(study, "scores"), recursive = TRUE)
  for (n in c("species", "aggregate", "pairs", "jobs")) write.csv(data.frame(x = 1), file.path(study, "scores", paste0(n, ".csv")), row.names = FALSE)
  writeLines("{}", file.path(study, "scores", "reading.json")); writeLines("{}", file.path(study, "timing-step1.json")); writeLines("{}", file.path(study, "microtiming.json"))
  cfg2 <- cfg; cfg2$packages <- c(occJSDM = "0.1.0", Hmsc = "3.3.7", gllvm = "2.0.15", posterior = "1.7.0"); cfg2$r_version <- "R 4.5.0"; cfg2$git <- "abc"
  cfg2$frozen <- list(support_points_900 = 900L, hmsc_method = "Full", hmsc_knots = NA_integer_, prediction_thin = 10L); saveRDS(cfg2, file.path(study, "config.rds"))
  out <- file.path(tempdir(), "plot-results"); copy_results(study, out)
  expect_true(all(file.exists(file.path(out, c("species.csv", "aggregate.csv", "pairs.csv", "jobs.csv", "reading.json", "timing-step1.json", "microtiming.json", "runtime.json")))))
  rt <- jsonlite::fromJSON(file.path(out, "runtime.json")); expect_equal(rt$packages$Hmsc, "3.3.7"); expect_match(rt$provisioning, "Route B")
})
cat("Plot tests passed.\n")
