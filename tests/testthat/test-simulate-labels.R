# simulateOccJSDMData() labels its site and sample truths with the ids in
# info, so truth can be joined to observations by identifier.

test_that("two-stage truths are labelled by site, sample and species", {
  sim <- simulate_fixture(model = "two_stage", M = 2L, P = 1L, K = 2L)
  info <- sim$data_list$info
  species <- colnames(sim$data_list$OTU)
  truth <- sim$true_params

  site_ids <- as.character(sort(unique(info$Site)))
  sample_ids <- as.character(sort(unique(info$Sample)))

  expect_identical(dimnames(truth$z_true), list(site_ids, species))
  expect_identical(dimnames(truth$jsdmParams_true$eta), list(site_ids, species))
  expect_identical(dimnames(truth$w_true), list(sample_ids, species))
})

test_that("occupancy and binary truths are labelled too", {
  occ <- simulate_fixture(model = "occupancy", M = 2L)
  species <- colnames(occ$data_list$OTU)
  site_ids <- as.character(sort(unique(occ$data_list$info$Site)))
  expect_identical(dimnames(occ$true_params$z_true), list(site_ids, species))
  expect_null(occ$true_params$w_true)

  bin <- simulate_fixture(model = "binary")
  expect_identical(colnames(bin$true_params$jsdmParams_true$eta),
                   colnames(bin$data_list$OTU))
  expect_identical(rownames(bin$true_params$jsdmParams_true$eta),
                   as.character(seq_len(nrow(bin$data_list$OTU))))
})
