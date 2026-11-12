# With an embedding, no b_prior is passed: naming embeddingVars must be enough.
fit_grouped <- function(sim, priors = if (embedding) list() else list(b_prior = "grouped"), embedding = TRUE,
                        params = list(), MCMCparams = FIXTURE_MCMC, spatCovariates = NULL) {
  dat <- sim$data_list
  occ <- grouped_occ_covariates(dat)
  if (embedding) params$embeddingVars <- dat$embeddingVars
  suppressMessages(suppressWarnings(runOccJSDM(dat, listParams = c(list(n_factors = 1L), params),
    listPriors = priors, occCovariates = occ, spatCovariates = spatCovariates, MCMCparams = MCMCparams)))
}

test_that("a grouped binary fit returns the new scales with the right shapes", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- fit_grouped(sim)
  jo <- fit$results_output$jsdm_output; mc <- fit$infos$MCMCparams
  expect_identical(unname(fit$infos$covariate_group), rep(c("measured", "embedding"), c(2L, 8L)))
  expect_identical(fit$infos$coef_prior$type, "grouped")
  expect_identical(dim(jo$lambda_output), c(8L, as.integer(mc$niter), as.integer(mc$nchain)))
  expect_identical(dimnames(jo$lambda_output)[[1]], names(fit$infos$covariate_group)[3:10])
  expect_identical(dim(jo$sigmaemb_output), dim(jo$sigmab_output))
  expect_identical(dim(jo$G_output)[1], ncol(sim$data_list$traits) + 1L)
  for (nm in c("lambda_output", "sigmaemb_output", "sigmab_output", "B_output"))
    expect_true(all(is.finite(jo[[nm]])), info = nm)
  expect_true(all(jo$lambda_output > 0) && all(jo$sigmaemb_output > 0))
  expect_gt(length(unique(as.vector(jo$sigmab_output))), 10L)   # no longer pinned
})

test_that("a grouped fit with no embedding columns learns only the measured scale", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- fit_grouped(sim, embedding = FALSE)
  jo <- fit$results_output$jsdm_output
  expect_null(jo$lambda_output); expect_null(jo$sigmaemb_output)
  expect_true(all(fit$infos$covariate_group == "measured"))
})

test_that("a fixed fit has none of the new elements", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- fit_grouped(sim, priors = list(), embedding = FALSE)
  expect_identical(fit$infos$coef_prior, list(type = "fixed"))
  expect_null(fit$results_output$jsdm_output$lambda_output)
  expect_identical(dim(fit$results_output$jsdm_output$G_output)[1], ncol(sim$data_list$traits))
})

test_that("bad embedding inputs stop before sampling", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  sim$data_list$info$X_emb.3 <- 1
  expect_error(fit_grouped(sim), "embeddingVars columns must vary across sites: X_emb.3")
  sim2 <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_error(fit_grouped(sim2, params = list(splineVars = "X_emb.1")),
               "embeddingVars must be numeric covariates without splines: X_emb.1")
})

test_that("a grouped fit runs on parallel chains", {
  # The chain workers read the prior, the trait matrix of the covariate side
  # and the starting scales from the setup list, and return the new draws for
  # runOccJSDM() to merge chain by chain. Fit identity across ncores is
  # test-parallel-chains.R's business and needs a reference BLAS.
  skip_on_cran()
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  set.seed(44)
  fit <- fit_grouped(sim, MCMCparams = list(nchain = 2L, nburn = 5L, niter = 6L,
                                            nthin = 1L, ncores = 2L))
  expect_identical(fit$infos$execution$workers, 2L)
  jo <- fit$results_output$jsdm_output
  expect_identical(dim(jo$lambda_output), c(8L, 6L, 2L))
  expect_identical(dim(jo$sigmaemb_output), c(6L, 2L))
  expect_true(all(is.finite(jo$lambda_output)) && all(is.finite(jo$sigmaemb_output)))
  expect_false(identical(jo$lambda_output[, , 1], jo$lambda_output[, , 2]))
})

test_that("the design columns record the covariate and kind each came from", {
  set.seed(5)
  df <- data.frame(a = rnorm(60), s = rnorm(60), f = factor(rep(c("x", "y", "z"), 20)),
                   few = rep(c(0, 1), 30), b = rnorm(60))
  built <- suppressMessages(occJSDM:::create_covariates_matrix(df, spline_vars = c("s", "few"),
                                                               remove_intercept = TRUE))
  got <- occJSDM:::occ_design_columns(built$X, built$list_matrix, built$X0)
  expect_identical(got$column, colnames(built$X))
  expect_identical(got$variable, rep(c("a", "s", "f", "few", "b"), c(1L, 5L, 2L, 1L, 1L)))
  # `few` was asked for as a spline and fell back to a linear term: still a spline.
  expect_identical(got$kind, rep(c("numeric", "spline", "factor", "spline", "numeric"),
                                 c(1L, 5L, 2L, 1L, 1L)))
  all_splines <- occJSDM:::create_covariates_matrix(df[c("a", "b")], spline_vars = TRUE,
                                                    remove_intercept = TRUE)
  expect_identical(occJSDM:::occ_design_columns(all_splines$X, all_splines$list_matrix,
                                                all_splines$X0)$kind, rep("spline", 10L))
  none <- occJSDM:::occ_design_columns(matrix(0, 60, 0), NULL, NULL)
  expect_identical(dim(none), c(0L, 3L))
  expect_error(occJSDM:::occ_design_columns(built$X[, -1], built$list_matrix, built$X0),
               "could not be matched")

  # A factor or a spline named as an embedding is rejected by the fit.
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  sim$data_list$info$X_emb.2 <- ifelse(sim$data_list$info$X_emb.2 > 0, "high", "low")
  expect_error(fit_grouped(sim), "embeddingVars must be numeric covariates without splines: X_emb.2")
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_error(fit_grouped(sim, params = list(splineVars = TRUE)),
               "embeddingVars must be numeric covariates without splines")
})

test_that("the starting scales leave measured columns at exactly 1", {
  cp <- occJSDM:::read_coef_prior(list(), list(embeddingVars = c("e1", "e2")), 10L, 2L,
    data.frame(column = c("m1", "e1", "e2"), variable = c("m1", "e1", "e2"), kind = "numeric"))
  expect_identical(occJSDM:::coef_scales_start(cp, 3L),
                   list(sigma_emb = 1, lambda = c(1, 1 / sqrt(2), 1 / sqrt(2))))
  expect_identical(occJSDM:::coef_scales_start(list(type = "fixed"), 3L),
                   list(sigma_emb = NA_real_, lambda = c(1, 1, 1)))
  expect_identical(occJSDM:::coef_prior_groups(list(type = "fixed"),
                     data.frame(column = character(0))), stats::setNames(character(0), character(0)))
})

test_that("the prior standard deviations are checked before they reach the draws", {
  msg <- "must be 2 finite positive numbers"
  expect_null(occJSDM:::check_coef_sds(c(1, 0.5), 2L))
  expect_null(occJSDM:::check_coef_sds(numeric(0), 0L))
  for (bad in list(NULL, 1, c(1, 1, 1), c(1, NA), c(1, Inf), c(1, 0), c(1, -1), c("1", "1")))
    expect_error(occJSDM:::check_coef_sds(bad, 2L), msg)

  # A sampler state without them stops in R, not inside the C++ draw.
  original_update <- occJSDM:::update_jSDMcoef
  local_mocked_bindings(update_jSDMcoef = function(list_data, list_params, ...) {
    list_params["sigma_b"] <- list(NULL)   # kept as an element: `$` would match sigma_bs
    original_update(list_data, list_params, ...)
  }, .package = "occJSDM")
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_error(fit_grouped(sim, priors = list(), embedding = FALSE), "must be 10 finite positive numbers")
})
