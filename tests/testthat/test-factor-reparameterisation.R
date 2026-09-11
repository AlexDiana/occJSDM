# These checks catch factor-specific magnitude normalisation: it leaves U L
# unchanged while changing the loading Gram matrix and new-site variation.
factor_arrays <- function(U, L, niter = 1L, nchain = 1L) {
  list(U = array(U, c(nrow(U), ncol(U), niter, nchain)),
       L = array(L, c(nrow(L), ncol(L), niter, nchain)))
}

test_that("one-factor reparameterisation changes only the sign", {
  U <- matrix(c(-2, 0.5, 3), ncol = 1)
  for (loadings in list(c(-3, 2, -1), c(0, -3, 2), c(0, 0, 0))) {
    input <- factor_arrays(U, matrix(loadings, nrow = 1))
    out <- reparamFactorModel(input$U, input$L)
    expect_true(all(is.finite(out$L_output)))
    expect_true(all(is.finite(out$U_output)))
    expect_equal(abs(out$L_output), abs(input$L))
    expect_equal(abs(out$U_output), abs(input$U))
    expect_equal(out$U_output[, , 1, 1] %o% out$L_output[, , 1, 1],
                 U %*% matrix(loadings, nrow = 1))
    nonzero <- which(loadings != 0)
    if (length(nonzero)) expect_gt(out$L_output[1, nonzero[1], 1, 1], 0)
  }
})

test_that("orthogonal identification preserves each draw's predictor and covariance", {
  # Unequal scales and mixed signs; the old transform flips correlation 2,3
  # from +1/sqrt(10) to -1/sqrt(10).
  L <- rbind(c(2, 2, 1, -1), c(0, 1, -1, 3))
  U <- rbind(c(1, 2), c(-2, 0.5), c(3, -1))
  input <- factor_arrays(U, L, niter = 2, nchain = 2)
  input$L[, , 2, 1] <- 2 * L
  input$L[, , 1, 2] <- -L
  out <- reparamFactorModel(input$U, input$L)
  expect_equal(dim(out$U_output), dim(input$U))
  expect_equal(dim(out$L_output), dim(input$L))
  for (chain in 1:2) for (iter in 1:2) {
    before <- input$L[, , iter, chain]
    after <- out$L_output[, , iter, chain]
    expect_equal(out$U_output[, , iter, chain] %*% after, U %*% before)
    expect_equal(crossprod(after), crossprod(before))
    expect_equal(cov2cor(crossprod(after)), cov2cor(crossprod(before)))
    expect_true(all(diag(after) > 0))
    expect_equal(after[2, 1], 0, tolerance = 1e-12)
  }
  expect_equal(cov2cor(crossprod(out$L_output[, , 1, 1]))[2, 3], 1 / sqrt(10))
})

test_that("zero anchors and deficient or rectangular loadings keep their geometry", {
  loadings <- list(
    rbind(c(0, 2, -1), c(0, 1, 3)), # QR moves the zero anchor column
    rbind(c(1, 2, 3), c(2, 4, 6)), # deficient rank
    matrix(0, 2, 3),
    matrix(c(1, -2, 3), ncol = 1), # more factors than columns
    matrix(numeric(), 0, 3)       # no latent factors
  )
  for (L in loadings) {
    # One site also tests preservation of matrix shape after array slicing.
    U <- matrix(seq_len(nrow(L)), nrow = 1)
    input <- factor_arrays(U, L)
    out <- reparamFactorModel(input$U, input$L)
    after_L <- matrix(out$L_output, nrow(L), ncol(L))
    after_U <- matrix(out$U_output, nrow(U), ncol(U))
    expect_true(all(is.finite(after_L)))
    expect_true(all(is.finite(after_U)))
    expect_equal(dim(out$L_output), dim(input$L))
    expect_equal(after_U %*% after_L, U %*% L)
    expect_equal(crossprod(after_L), crossprod(L))
  }
})

test_that("stored factors and latent traits preserve their pre-rotation products", {
  captured <- list()
  original <- reparamFactorModel
  # Capture the real fit's inputs at its two post-processing boundaries.
  fitter <- runOccJSDM
  fit_env <- new.env(parent = environment(fitter))
  fit_env$reparamFactorModel <- function(U_output, L_output) {
    captured[[length(captured) + 1L]] <<- list(U = U_output, L = L_output)
    original(U_output, L_output)
  }
  environment(fitter) <- fit_env
  sim <- simulate_fixture(model = "continuous")
  fit <- suppressMessages(suppressWarnings(fitter(
    sim$data_list, listParams = list(n_factors = 2),
    occCovariates = fixture_occ_covariates(), MCMCparams = FIXTURE_MCMC)))
  expect_length(captured, 2)
  jsdm <- fit$results_output$jsdm_output
  for (i in 1:2) {
    U <- if (i == 1) jsdm$U_output else jsdm$A_output
    L <- if (i == 1) jsdm$L_output else jsdm$C_output
    for (chain in seq_len(dim(L)[4])) for (iter in seq_len(dim(L)[3])) {
      draw <- function(x) matrix(x[, , iter, chain], dim(x)[1], dim(x)[2])
      expect_equal(draw(U) %*% draw(L),
                   draw(captured[[i]]$U) %*% draw(captured[[i]]$L))
      expect_equal(crossprod(draw(L)), crossprod(draw(captured[[i]]$L)))
    }
  }
  raw_fit <- fit
  raw_fit$results_output$jsdm_output$L_output <- captured[[1]]$L
  expect_equal(returnResidualCorrelationMatrix(fit), returnResidualCorrelationMatrix(raw_fit))
  expect_equal(dim(returnOrdinationScores(fit)), c(3L, FIXTURE_N, 2L))
  expect_equal(dim(returnFactorLoadings(fit)), c(3L, 2L, FIXTURE_S))
  expect_s3_class(plotOrdinationScores(fit), "ggplot")
  expect_s3_class(plotFactorLoadings(fit), "ggplot")
  expect_null(plotBiplot(fit)$labels$caption)
  pred <- predictNewSites(fit, X_psi = sim$data_list$info[, fixture_occ_covariates()],
                          useSpatial = FALSE, verbose = FALSE)
  expect_true(all(is.finite(pred)))
})

test_that("new-site prediction retains the pre-rotation factor variance", {
  L <- rbind(c(2, 2, 1), c(0, 1, -1))
  input <- factor_arrays(matrix(c(1, 2), nrow = 1), L)
  rotated <- reparamFactorModel(input$U, input$L)
  set.seed(914)
  pred <- computeNewOutputs(
    X = matrix(0, 2, 1), B0_output = matrix(0, 3, 12000),
    B_output = array(0, c(1, 3, 12000)), Ks_all = array(0, c(0, 0, 0)),
    Bs_output = array(0, c(0, 0, 0)),
    L_output = array(rotated$L_output, c(2, 3, 12000)),
    sigmah_output = rep(0.75, 12000), idx_ls_output = rep(1, 12000),
    conflevels = c(0.025, 0.5, 0.975), useEnvCov = TRUE,
    useSpatial = FALSE, useBiotic = TRUE, model = "continuous", verbose = FALSE)
  # Analytic normal quantiles from sigma_h^2 L' L, independent of QR.
  expected <- qnorm(c(0.025, 0.5, 0.975)) %o% (0.75 * c(2, sqrt(5), sqrt(2)))
  expect_lt(max(abs(pred[, 1, ] - expected)), 0.2)
  expect_lt(max(abs(pred[, 2, ] - expected)), 0.2)
})
