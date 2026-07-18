list_datasettings <- list(
  n = 5, S = 3, g = 2,
  M = rep(2, 5),
  P = 2,
  K = rep(3, 2 * sum(rep(2, 5))),
  ncov_psi = 2, ncov_theta = 2
)

list_params <- list(
  p = matrix(0.8, list_datasettings$P, list_datasettings$S),
  q = matrix(0.05, list_datasettings$P, list_datasettings$S),
  theta0 = rep(0.05, list_datasettings$S),
  theta_baseline = rep(0.7, list_datasettings$S)
)

list_jsdmParams <- list(
  gt = 2, d = 2, ds = 2,
  sigma_b = 1, sigma_bs = 1, sigma_ts = 1, sigma_h = 1, l_s = 1
)

test_that("toRunOccJSDMFormat returns info/OTU with expected dimensions and names", {
  set.seed(8214)

  sim <- simulateOccJSDMData(list_datasettings, list_params, list_jsdmParams)

  converted <- toRunOccJSDMFormat(
    sim,
    n = list_datasettings$n,
    M = list_datasettings$M,
    P = list_datasettings$P,
    K = list_datasettings$K
  )

  expect_named(converted, c("info", "OTU"))
  expect_s3_class(converted$info, "data.frame")

  N3 <- sum(list_datasettings$K)

  # One row per PCR replicate, matching in info and OTU
  expect_equal(nrow(converted$info), N3)
  expect_equal(nrow(converted$OTU), N3)

  # OTU has one column per species
  expect_equal(ncol(converted$OTU), list_datasettings$S)

  # info has id columns plus expanded covariates (2 X_psi + 1 X_theta,
  # since the intercept column is dropped by default)
  expect_true(all(c("Site", "Sample", "Primer") %in% names(converted$info)))
  expect_true(all(c("X_psi.1", "X_psi.2") %in% names(converted$info)))
  expect_true("X_theta" %in% names(converted$info))

  # Site/Sample/Primer take values in the expected ranges
  expect_true(all(converted$info$Site %in% seq_len(list_datasettings$n)))
  expect_true(all(converted$info$Primer %in% seq_len(list_datasettings$P)))
})

test_that("toRunOccJSDMFormat keeps the theta intercept column when requested", {
  set.seed(3067)

  sim <- simulateOccJSDMData(list_datasettings, list_params, list_jsdmParams)

  converted <- toRunOccJSDMFormat(
    sim,
    n = list_datasettings$n,
    M = list_datasettings$M,
    P = list_datasettings$P,
    K = list_datasettings$K,
    drop_theta_intercept = FALSE
  )

  expect_true(all(c("X_theta.1", "X_theta.2") %in% names(converted$info)))
  expect_true(all(converted$info$X_theta.1 == 1))
})

test_that("toRunOccJSDMFormat errors when n/M/P/K don't match sim", {
  set.seed(5921)

  sim <- simulateOccJSDMData(list_datasettings, list_params, list_jsdmParams)

  # A K vector inconsistent with n/M/P (and therefore with sim$data_list$y)
  # should error, whether inside createDataIdx() or the row-count check.
  expect_error(
    toRunOccJSDMFormat(
      sim,
      n = list_datasettings$n,
      M = list_datasettings$M,
      P = list_datasettings$P,
      K = rep(3, 2 * sum(rep(2, 5)) - 1)
    )
  )
})
