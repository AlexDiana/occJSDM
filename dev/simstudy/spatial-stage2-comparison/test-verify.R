# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-verify.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
source(file.path(D, "verify.R"))
# score.R is sourced only to produce the tables that the independent code then has to reproduce.
for (f in c("config.R", "score.R")) source(file.path(D, f))
cfg <- stage2_config()

synthetic_species <- function(seed = 3) {
  set.seed(seed); rows <- list()
  for (cell in c("t1-r0.08-n900", "t2-n900")) for (rep in 1:10) for (k in seq_len(nrow(cfg$arms))) {
    key <- arm_key(cfg$arms$package[k], cfg$arms$arm[k]); level <- c(env = 0, gp = 0.40, spatial2 = 0.55, spatial6 = 0.50)[[ if (key == "occJSDM/env") "env" else if (key == "occJSDM/gp") "gp" else if (key == "Hmsc/spatial6") "spatial6" else "spatial2" ]]
    status <- if (key == "gllvm/spatial2" && rep > 7) "unstable" else "pass"
    rows[[length(rows) + 1L]] <- data.frame(dataset_id = sprintf("%s-rep%d", cell, rep), package = cfg$arms$package[k], arm = cfg$arms$arm[k], key = key,
      species = paste0("s", 1:4), occupancy = c(0.05, 0.25, 0.5, 0.75), status = status, rmse_env = 0.1, rmse_oracle = 0.05, rmse = 0.1,
      headroom_share = level + rnorm(4, 0, 0.05) + rnorm(1, 0, 0.03), pattern = 0.5, field_slope = NA_real_, field_correlation = NA_real_, stringsAsFactors = FALSE)
  }
  species <- do.call(rbind, rows); species$headroom_share[species$key == "occJSDM/env"] <- 0; species$headroom_share[3] <- NA; species
}

test_that("the independent recomputation reproduces score.R's pairs and flags a tampered row", {
  species <- synthetic_species(); agg <- aggregate_arm(species, cfg$gate_occupancy); pairs <- paired_comparisons(agg, cfg)
  ok <- verify_pairs(species, pairs, cfg)
  expect_equal(nrow(ok), 10L); expect_true(all(ok$passed))   # two synthetic cells x five comparisons; expect_true(any(grepl("gllvm", ok$check)))
  bad <- pairs; bad$estimate[3] <- bad$estimate[3] + 0.01
  expect_false(all(verify_pairs(species, bad, cfg)$passed))
  bad <- pairs; bad$verdict[5] <- "inconclusive"
  expect_false(all(verify_pairs(species, bad, cfg)$passed))
  bad <- pairs; bad$excluded[2] <- bad$excluded[2] + 1L
  expect_false(all(verify_pairs(species, bad, cfg)$passed))
})

test_that("the reading is recomputed from the pairs and compared with the stored reading", {
  species <- synthetic_species(); agg <- aggregate_arm(species, cfg$gate_occupancy); pairs <- paired_comparisons(agg, cfg)
  finite <- agg[is.finite(agg$headroom_share), ]   # the same definition of the largest median share that stage2_score uses
  reading <- development_reading(pairs, max(aggregate(headroom_share ~ cell + key, finite, median)$headroom_share))
  r <- verify_reading(species, pairs, jsonlite::fromJSON(jsonlite::toJSON(reading, auto_unbox = TRUE, na = "null", digits = NA)), cfg)
  expect_true(all(r$passed))
  wrong <- reading; wrong$reading <- if (identical(reading$reading, "trade-off")) "shared factors" else "trade-off"
  expect_false(all(verify_reading(species, pairs, jsonlite::fromJSON(jsonlite::toJSON(wrong, auto_unbox = TRUE, na = "null", digits = NA)), cfg)$passed))
})
cat("Verify tests passed.\n")
