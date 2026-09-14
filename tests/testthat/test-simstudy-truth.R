# These tests check the quantities the simulation study scores. No MCMC is
# needed to distinguish raw generating coefficients from fitted-scale truth.
simstudy_truth_case <- function(beta_theta = NULL, collection_design = NULL) {
  list(fit = list(results_output = list(jsdm_output = list()),
                  infos = list(list_X_theta_mat = collection_design$list_matrix),
                  X_theta = collection_design$X),
       sim = list(true_params = list(jsdmParams_true = list(),
                                     beta_theta_true = beta_theta)),
       truth = list(params = list()))
}

test_that("simstudy scores recorded detections rather than latent read events", {
  x <- simstudy_truth_case()
  p <- matrix(c(.3, .5, .6, .8), 2, 2,
              dimnames = list(c("primer2", "primer1"), c("spB", "spA")))
  q <- p / 3
  x$sim$true_params$p_true <- p
  x$sim$true_params$q_true <- q
  blocks <- do.call(simstudy_param_blocks, x)

  # Independent tail probabilities at log(1.5), the read-rounding boundary
  # for threshold one; checking both detects use of one intensity for both.
  expect_equal(blocks$p$truth, p * 0.9999978314222251, tolerance = 1e-14)
  expect_equal(blocks$q$truth, q * 0.8631397734521575, tolerance = 1e-14)
  expect_identical(dimnames(blocks$q$truth), dimnames(q))
  expect_identical(x$sim$true_params$q_true, q)
})

test_that("simstudy detection truth honours overridden read-intensity settings", {
  x <- simstudy_truth_case()
  x$sim$true_params$p_true <- matrix(c(.4, .8), 1)
  x$sim$true_params$q_true <- matrix(c(.1, .3), 1)
  # A Normal centered at the read-rounding boundary retains exactly half
  # its draws. A mean one SD above it retains Phi(1).
  x$truth$params <- list(mu1 = log(1.5), sigma1 = 2,
                         mu0 = log(1.5) + 2, sigma0 = 2)
  blocks <- do.call(simstudy_param_blocks, x)
  expect_equal(blocks$p$truth, matrix(c(.2, .4), 1), tolerance = 1e-14)
  expect_equal(blocks$q$truth,
               matrix(c(.0841344746068543, .252403423820563), 1),
               tolerance = 1e-14)
})

test_that("simstudy collection truth preserves the raw-scale linear predictor", {
  raw <- data.frame(effort = c(2, 4, 6, 8), temperature = c(10, 14, 18, 22))
  design <- create_covariates_matrix(raw, spline_vars = FALSE,
                                      remove_intercept = FALSE)
  beta <- matrix(c(.3, 1, -.5, -.2, 0, 0), 3, 2,
                 dimnames = list(c("(Intercept)", "effort", "temperature"),
                                 c("spB", "spA")))
  # Lists may store named scales in a different order. Pair by covariate
  # identity, never by the incidental order of those named vectors.
  design$list_matrix$mean_df <- rev(design$list_matrix$mean_df)
  design$list_matrix$sd_df <- rev(design$list_matrix$sd_df)
  x <- simstudy_truth_case(beta, design)
  corrected <- do.call(simstudy_param_blocks, x)$beta_theta$truth
  raw_eta <- cbind(1, as.matrix(raw)) %*% beta
  expect_equal(unname(design$X %*% corrected), unname(raw_eta), tolerance = 1e-14)
  expect_equal(unname(corrected[1, ]), c(-2.7, -.2), tolerance = 1e-14)
  expect_identical(corrected[-1, 2], beta[-1, 2])
  expect_identical(dimnames(corrected), dimnames(beta))
  expect_identical(x$sim$true_params$beta_theta_true, beta)
})

test_that("simstudy rejects ambiguous collection coefficient and scale mappings", {
  raw <- data.frame(effort = c(2, 4, 6, 8))
  design <- create_covariates_matrix(raw, spline_vars = FALSE,
                                      remove_intercept = FALSE)
  beta <- matrix(c(.3, 1, -.2, 0), 2, 2)
  x <- simstudy_truth_case(beta, design)
  x$fit$infos$list_X_theta_mat$sd_df <- c(other = 2)
  expect_error(do.call(simstudy_param_blocks, x), "collection.*scal")

  x <- simstudy_truth_case(beta, design)
  rownames(x$sim$true_params$beta_theta_true) <- c("(Intercept)", "other")
  expect_error(do.call(simstudy_param_blocks, x), "collection.*mapping")
})

test_that("simstudy keeps absent detection blocks and intercept-only truth intact", {
  x <- simstudy_truth_case()
  blocks <- do.call(simstudy_param_blocks, x)
  expect_null(blocks$beta_theta$truth)
  expect_null(blocks$p$truth)
  expect_null(blocks$q$truth)

  beta <- matrix(c(-.8, .1), 1, dimnames = list("(Intercept)", c("spA", "spB")))
  x$sim$true_params$beta_theta_true <- beta
  blocks <- do.call(simstudy_param_blocks, x)
  expect_identical(blocks$beta_theta$truth, beta)
  expect_null(blocks$p$truth)
  expect_null(blocks$q$truth)
})
