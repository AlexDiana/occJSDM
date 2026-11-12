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
  # The scales are redrawn every sweep and the draws are stored. If the grouped
  # update were skipped, or its results dropped, every scale would sit at its
  # starting value: 1 for sigma_b, sigma_emb_scale for sigma_emb and axis_scale
  # for each embedding column's covariate scale.
  cp <- fit$infos$coef_prior
  distinct <- function(x) length(unique(as.vector(x)))
  for (k in seq_len(nrow(jo$lambda_output)))
    expect_gt(distinct(jo$lambda_output[k, , 1]), 10L, label = rownames(jo$lambda_output)[k])
  expect_gt(distinct(jo$sigmaemb_output[, 1]), 10L)
  expect_gt(distinct(jo$sigmab_output[, 1]), 10L)
  expect_false(all(jo$lambda_output == cp$axis_scale))
  expect_false(all(jo$sigmaemb_output == cp$sigma_emb_scale))
  expect_false(all(jo$sigmab_output == 1))
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
                   list(sigma_emb = 1, lambda = c(1, 0.7 / 2, 0.7 / 2)))   # the default, 0.7 / K
  expect_identical(occJSDM:::coef_scales_start(cp, 3L)$lambda[2:3], rep(cp$axis_scale, 2))
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

grouped_fit <- function() cached_fixture("grouped_binary", function()
  fit_grouped(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)))

test_that("the scale and community accessors return one row per quantity", {
  fit <- grouped_fit()
  sc <- returnCovariateScales(fit)
  expect_identical(names(sc), c("scale", "covariate", "mean", "lower", "median", "upper"))
  expect_identical(sc$scale, c("measured group", "embedding leftover", rep("axis", 8)))
  expect_identical(sc$covariate[3:10], names(fit$infos$covariate_group)[3:10])
  expect_true(all(sc$lower <= sc$median & sc$median <= sc$upper))
  cr <- returnCommunityResponse(fit)
  expect_identical(cr$covariate, names(fit$infos$covariate_group))
  expect_identical(cr$group, unname(fit$infos$covariate_group))
})

test_that("blends are sign-fixed and reproduce the low-rank part", {
  fit <- grouped_fit(); bl <- returnEmbeddingBlends(fit, n_blends = 2L)
  expect_identical(dim(bl$weights), c(8L, 2L))
  expect_identical(nrow(bl$site_scores), nrow(fit$X_psi))
  expect_true(all(apply(bl$weights, 2, function(v) v[which.max(abs(v))]) > 0))
  expect_error(returnEmbeddingBlends(fit, n_blends = 99L), "n_blends cannot exceed the number of latent traits")
})

test_that("the embedding share lies between 0 and 1", {
  sh <- returnEmbeddingShare(grouped_fit())
  expect_identical(nrow(sh), 12L)
  expect_true(all(sh$mean >= 0 & sh$mean <= 1))
})

test_that("existing downstream functions work on a grouped fit", {
  fit <- grouped_fit(); sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_identical(nrow(returnTraitsCoeff(fit)) > 0, TRUE)
  thin <- thinOutput(fit, thin = 2)
  expect_identical(dim(thin$results_output$jsdm_output$lambda_output)[2],
                   length(seq(1, fit$infos$MCMCparams$niter, by = 2)))
  expect_no_error(suppressWarnings(utils::capture.output(computeDiagnostics(fit$results_output))))
  newX <- unique(sim$data_list$info[, grouped_occ_covariates(sim$data_list)])[1:5, ]
  pred <- predictNewSites(fit, X_psi = newX, useSpatial = FALSE, useBiotic = FALSE, verbose = FALSE)
  expect_false(anyNA(unlist(pred)))
})

test_that("the accessors refuse fits they do not apply to", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fixed <- fit_grouped(sim, priors = list(), embedding = FALSE)
  expect_error(returnCovariateScales(fixed), "needs a fit made with the grouped prior")
  expect_error(returnCommunityResponse(fixed), "needs a fit made with the grouped prior")
  expect_error(returnEmbeddingBlends(fixed), "needs a fit with embeddingVars and at least one latent trait")
})

test_that("a grouped fit with no traits behaves like any fit with no traits", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  sim$data_list$traits <- NULL
  fit <- fit_grouped(sim)
  expect_identical(fit$infos$coef_prior$type, "grouped")
  expect_error(returnTraitsCoeff(fit), "No Traits Present")
  expect_error(plotTraitsCoefficients(fit, covName = "X_psi.1"), "No Traits Present")
  cr <- returnCommunityResponse(fit)
  expect_identical(nrow(cr), 10L)
})

test_that("the default number of blends is bounded by the embedding columns and species", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 1L, S = 12L, gt = 2L)
  fit <- fit_grouped(sim, params = list(n_lattrait = 2L))
  expect_identical(dim(fit$results_output$jsdm_output$A_output)[2], 2L)
  bl <- returnEmbeddingBlends(fit)
  expect_identical(dim(bl$weights), c(1L, 1L))
  expect_identical(ncol(bl$site_scores), 1L)
  expect_error(returnEmbeddingBlends(fit, n_blends = 2L),
               "n_blends cannot exceed the number of embedding columns or species", fixed = TRUE)
})

test_that("n_blends must be a single positive whole number", {
  fit <- grouped_fit()
  for (bad in list(0L, 1.5, c(1L, 2L)))
    expect_error(returnEmbeddingBlends(fit, n_blends = bad),
                 "n_blends must be a single positive whole number", fixed = TRUE)
})

test_that("axes with no effect are shrunk more than axes with an effect", {
  skip_on_cran()
  sim <- simulate_grouped_fixture(n = 150L, S = 30L, p_meas = 2L, p_emb = 10L, gt = 2L, seed = 4L,
    embedding = list(n_blends = 2L, n_active = 3L, blend_sd = 1.2, sigma_emb = 0.1))
  fit <- fit_grouped(sim, MCMCparams = list(nchain = 2, nburn = 1500, niter = 1500, nthin = 1))
  lam <- apply(fit$results_output$jsdm_output$lambda_output, 1, median)
  active <- rowSums(abs(sim$true_params$jsdmParams_true$B[3:12, ])) > 0
  expect_gt(min(lam[active]), max(lam[!active]))
})

# A grouped fit used the prior and kept updating its scales. Checks the prior
# type, the group of every design column (measured first, then embedding, the
# embedding names matching the rows of lambda_output), and more than 10
# distinct values in the first chain for every embedding column's covariate
# scale and for the embedding scale. A fit that ignored the prior has no
# scales to look at; one that skipped their update holds every scale at its
# starting value.
expect_grouped_scales_move <- function(fit, p_meas, p_emb, arm) {
  jo <- fit$results_output$jsdm_output
  distinct <- function(x) length(unique(as.vector(x)))
  expect_identical(fit$infos$coef_prior$type, "grouped", info = arm)
  group <- fit$infos$covariate_group
  expect_identical(unname(group), rep(c("measured", "embedding"), c(p_meas, p_emb)), info = arm)
  expect_identical(names(group)[p_meas + seq_len(p_emb)], dimnames(jo$lambda_output)[[1]], info = arm)
  expect_gt(min(apply(jo$lambda_output[, , 1, drop = FALSE], 1, distinct)), 10L,
            label = paste0(arm, ": fewest distinct draws of any covariate scale"))
  expect_gt(distinct(jo$sigmaemb_output[, 1]), 10L,
            label = paste0(arm, ": distinct draws of the embedding scale"))
}

test_that("grouped fits run for every supported model type and layout", {
  skip_on_cran()
  for (model in c("occupancy", "two_stage")) {
    sim <- simulate_grouped_fixture(model = model, p_emb = 6L, S = 8L)
    # Replicate rows of one site must carry that site's embedding values.
    info <- sim$data_list$info
    spread <- vapply(info[sim$data_list$embeddingVars], function(v)
      max(tapply(v, info$Site, function(x) diff(range(x)))), numeric(1))
    expect_true(nrow(info) > length(unique(info$Site)) && all(spread == 0), info = model)
    fit <- fit_grouped(sim)
    expect_true(all(is.finite(fit$results_output$jsdm_output$lambda_output)), info = model)
    expect_grouped_scales_move(fit, 2L, 6L, model)
  }
  wide <- fit_grouped(simulate_grouped_fixture(n = 30L, p_emb = 40L, S = 8L))   # more columns than sites
  expect_true(all(is.finite(wide$results_output$jsdm_output$B_output)))
  expect_grouped_scales_move(wide, 2L, 40L, "more columns than sites")
  noblend <- simulate_grouped_fixture(p_emb = 6L, S = 8L, gt = 0L)
  fit0 <- fit_grouped(noblend, params = list(n_lattrait = 0L))
  expect_true(all(is.finite(fit0$results_output$jsdm_output$lambda_output)))
  expect_grouped_scales_move(fit0, 2L, 6L, "no latent traits")
  # Values are not compared between the two: fits are identical across worker
  # counts only under a reference BLAS, which test-parallel-chains.R covers.
  par <- fit_grouped(simulate_grouped_fixture(p_emb = 6L, S = 8L),
                     MCMCparams = c(FIXTURE_MCMC, list(ncores = 2)))
  ser <- fit_grouped(simulate_grouped_fixture(p_emb = 6L, S = 8L))
  expect_identical(par$infos$execution$workers, 2L)
  expect_identical(dim(par$results_output$jsdm_output$lambda_output),
                   dim(ser$results_output$jsdm_output$lambda_output))
  expect_grouped_scales_move(par, 2L, 6L, "two workers")
  expect_grouped_scales_move(ser, 2L, 6L, "serial")
})

test_that("a grouped fit runs with a spatial field", {
  skip_on_cran()
  sim <- simulate_fixture(model = "binary", useSpatField = TRUE, seed = 6L)
  fit <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 1L, n_supportpoints = FIXTURE_KNOTS, embeddingVars = fixture_occ_covariates()[2]),
    listPriors = list(b_prior = "grouped"), occCovariates = fixture_occ_covariates(),
    spatCovariates = fixture_spat_covariates(), MCMCparams = FIXTURE_MCMC)))
  expect_true(all(is.finite(fit$results_output$jsdm_output$lambda_output)))
  expect_grouped_scales_move(fit, 1L, 1L, "spatial field")
})
