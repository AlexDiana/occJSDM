# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-data.R
library(testthat)
D <- "dev/simstudy/spatial-stage2-comparison"
for (f in c("upstream/adapter-contract.R", "upstream/basis.R", "config.R", "data.R")) source(file.path(D, f))
PILOT <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"

test_that("config fixes arms, cells, budgets and thresholds", {
  cfg <- stage2_config()
  expect_equal(cfg$arms, data.frame(package = c("occJSDM", "occJSDM", "Hmsc", "Hmsc", "gllvm"),
                                    arm = c("env", "gp", "spatial2", "spatial6", "spatial2"), stringsAsFactors = FALSE))
  expect_setequal(cfg$cells, c("t1-r0.08-n300", "t1-r0.08-n900", "t2-n300", "t2-n900"))
  expect_equal(cfg$mcmc$initial[["nburn"]], 3000L); expect_equal(cfg$mcmc$extended[["niter"]], 16000L)
  expect_equal(cfg$delta, 0.10); expect_equal(cfg$limits$attempts, 400L); expect_equal(cfg$wall_hours, 48)
  expect_true(is.na(cfg$frozen$support_points_900)); expect_true(is.na(cfg$frozen$hmsc_method))
})

test_that("the dataset table lists 40 datasets with matching input hashes and oracle files", {
  d <- stage2_datasets(PILOT)
  expect_equal(nrow(d), 40L); expect_equal(as.vector(table(d$cell)), rep(10L, 4))
  expect_true(all(file.exists(d$input_file))); expect_true(all(file.exists(d$oracle_file)))
  man <- read.csv(file.path(PILOT, "input-manifest.csv"), stringsAsFactors = FALSE)
  expect_true(all(d$input_md5 == man$md5[match(basename(d$input_file), man$file)]))
  expect_true(all(grepl("jobs-resolution", d$oracle_file[d$dataset_id %in% c("t1-r0.08-n900-rep2", "t2-n900-rep3")])))
  expect_false(any(grepl("jobs-resolution", d$oracle_file[d$cell == "t1-r0.08-n300" & d$replicate == 1L])))
})

test_that("training, prediction and truth are built from a pilot community and pass the contract", {
  d <- stage2_datasets(PILOT); row <- d[d$dataset_id == "t2-n300-rep1", ]
  cm <- readRDS(row$input_file)
  tr <- load_training(cm, row$dataset_id, row$input_md5)
  expect_true(validate_training(tr))
  expect_equal(dim(tr$y), c(300L, 12L)); expect_equal(tr$species_id, cm$species$species)
  expect_equal(tr$xy, unname(cm$points[cm$index$uniform, ]), ignore_attr = TRUE)
  pr <- load_prediction(cm)
  expect_true(validate_prediction_data(pr, tr)); expect_equal(nrow(pr$xy), 1600L)
  expect_false(any(c("y", "psi", "field") %in% names(pr)))
  truth <- load_truth(cm)
  expect_equal(dim(truth$psi_lattice), c(1600L, 12L)); expect_equal(truth$occupancy, cm$species$occupancy)
  expect_equal(colMeans(truth$psi_lattice), cm$species$occupancy, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("the oracle reference loads with its status and hash", {
  d <- stage2_datasets(PILOT); row <- d[d$dataset_id == "t1-r0.08-n900-rep2", ]
  o <- load_oracle_reference(row$oracle_file)
  expect_equal(dim(o$probability_mean), c(1600L, 12L)); expect_true(all(o$probability_mean > 0 & o$probability_mean < 1))
  expect_true(all(o$screen$status %in% c("pass", "unresolved"))); expect_equal(nchar(o$md5), 32L)
  expect_identical(o$md5, row$oracle_md5)
})
cat("Data tests passed.\n")
