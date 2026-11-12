test_that("b_prior = 'fixed' reproduces fits made before the option", {
  ref <- readRDS(test_path("fixtures", "fixed-prior-reference.rds"))
  skip_if_not(identical(La_library(), ref$la) && identical(extSoftVersion()[["BLAS"]], ref$blas),
    "BLAS or LAPACK differs from the reference")
  for (model in names(ref$fits)) {
    expect_identical(fixed_prior_reference_fit(model), ref$fits[[model]], info = model)
  }
})

test_that("the simulator adds an embedding block without touching existing calls", {
  a <- simulate_fixture(model = "binary", useSpatField = FALSE, seed = 3L)
  b <- simulate_fixture(model = "binary", useSpatField = FALSE, seed = 3L)
  expect_identical(a, b)
  expect_null(a$data_list$embeddingVars)

  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_identical(sim$data_list$embeddingVars, paste0("X_emb.", 1:8))
  expect_true(all(sim$data_list$embeddingVars %in% names(sim$data_list$info)))
  expect_identical(sim$true_params$jsdmParams_true$covariate_group,
                   rep(c("measured", "embedding"), c(2L, 8L)))
  expect_no_error(simulate_grouped_fixture(gt = 0L))   # no latent traits must still simulate
  expect_identical(dim(sim$true_params$jsdmParams_true$B), c(10L, 12L))
})

test_that("the embedding truth is low rank on the active axes and zero elsewhere", {
  sim <- simulate_grouped_fixture(p_emb = 8L, S = 12L,
    embedding = list(n_blends = 2L, n_active = 3L, sigma_emb = 0))
  Bemb <- sim$true_params$jsdmParams_true$B[3:10, ]
  active <- rowSums(abs(Bemb)) > 0
  expect_identical(sum(active), 3L)
  expect_lte(qr(Bemb)$rank, 2L)
})

test_that("a supplied embedding matrix and coefficient matrix are used as given", {
  X <- matrix(rnorm(60 * 8), 60, 8); B <- matrix(0, 8, 12); B[1, ] <- 1
  sim <- simulate_grouped_fixture(n = 60L, p_emb = 8L, S = 12L, embedding = list(X = X, B = B))
  expect_equal(unname(as.matrix(unique(sim$data_list$info[, paste0("X_emb.", 1:8)]))), X)
  expect_identical(unname(sim$true_params$jsdmParams_true$B[3:10, ]), B)
  expect_error(simulate_grouped_fixture(n = 60L, p_emb = 8L, embedding = list(X = X[1:10, ])),
               "embedding\\$X must have one row per site and n_emb columns")
})

test_that("a supplied embedding matrix enters the truth standardised and is returned as given", {
  set.seed(9L)
  X <- matrix(rnorm(60 * 8), 60, 8); B <- matrix(0, 8, 12); B[1:2, ] <- c(1, -1)
  a <- simulate_grouped_fixture(n = 60L, p_emb = 8L, S = 12L, embedding = list(X = X, B = B))
  b <- simulate_grouped_fixture(n = 60L, p_emb = 8L, S = 12L, embedding = list(X = 5 * X + 2, B = B))
  expect_identical(a$data_list$OTU, b$data_list$OTU)
  expect_equal(a$true_params$jsdmParams_true$B, b$true_params$jsdmParams_true$B)
  expect_equal(a$true_params$jsdmParams_true$eta, b$true_params$jsdmParams_true$eta)
  emb <- paste0("X_emb.", 1:8)
  expect_false(isTRUE(all.equal(a$data_list$info[, emb], b$data_list$info[, emb])))

  X0 <- X; X0[, 3] <- 1
  expect_error(simulate_grouped_fixture(n = 60L, p_emb = 8L, S = 12L, embedding = list(X = X0, B = B)),
               "embedding$X columns must vary across sites", fixed = TRUE)
})

cols <- function(p_meas = 2L, p_emb = 3L, kind = "numeric") data.frame(
  column = c(paste0("m", seq_len(p_meas)), paste0("e", seq_len(p_emb))),
  variable = c(paste0("m", seq_len(p_meas)), paste0("e", seq_len(p_emb))),
  kind = kind)
rcp <- function(priors = list(), params = list(), S = 10L, gt = 2L, columns = cols())
  occJSDM:::read_coef_prior(priors, params, S, gt, columns)

test_that("read_coef_prior resolves defaults", {
  expect_identical(rcp(), list(type = "fixed"))
  expect_identical(rcp(NULL, NULL), list(type = "fixed"))
  got <- rcp(list(), list(embeddingVars = paste0("e", 1:3)))   # embeddingVars alone switches it on
  expect_identical(got$type, "grouped")
  expect_identical(rcp(list(b_prior = "grouped"), list(embeddingVars = paste0("e", 1:3))), got)
  expect_identical(rcp(list(sigma_b_scale = 2), list(embeddingVars = "e1"))$sigma_b_scale, 2)
  expect_identical(got$is_emb, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_identical(got$sigma_b_scale, 1)
  expect_identical(got$sigma_emb_scale, 1)
  expect_equal(got$axis_scale, 1 / sqrt(3))
  expect_equal(got$c_base, c(1, 1, rep(1 / sqrt(2), 3)))
  expect_equal(rcp(list(b_prior = "grouped"), list(embeddingVars = "e1"), gt = 0L)$c_base,
               rep(1, 5))
  none <- rcp(list(b_prior = "grouped"))
  expect_false(any(none$is_emb)); expect_true(is.na(none$axis_scale))
})

test_that("read_coef_prior rejects settings that do not apply", {
  g <- list(b_prior = "grouped"); e <- list(embeddingVars = "e1")
  for (bad in list("Grouped", NA_character_, c("fixed", "grouped"), 1, TRUE))
    expect_error(rcp(list(b_prior = bad)), "b_prior must be 'fixed' or 'grouped'")
  expect_error(rcp(list(b_prior = "fixed"), e), "embeddingVars cannot be used with b_prior = 'fixed'")
  expect_error(rcp(g, e, S = 1L), "the grouped prior needs more than one species")
  expect_error(rcp(list(), e, S = 1L), "the grouped prior needs more than one species")
  expect_error(rcp(g, list(embeddingVars = c("e1", "zz"))),
               "embeddingVars must name occupancy covariates: zz")
  expect_error(rcp(g, e, columns = cols(kind = c("numeric", "numeric", "spline", "numeric", "numeric"))),
               "embeddingVars must be numeric covariates without splines: e1")
  expect_error(rcp(g, e, columns = cols(kind = c("numeric", "numeric", "factor", "numeric", "numeric"))),
               "embeddingVars must be numeric covariates without splines: e1")
  expect_error(rcp(list(sigma_b_scale = 2)), "sigma_b_scale requires the grouped prior")
  expect_error(rcp(c(g, axis_scale = 1)), "axis_scale requires embeddingVars")
  expect_error(rcp(c(g, sigma_emb_scale = 1)), "sigma_emb_scale requires embeddingVars")
  for (nm in c("sigma_b_scale", "sigma_emb_scale", "axis_scale")) for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "a"))
    expect_error(rcp(c(g, stats::setNames(list(bad), nm)), e),
                 paste(nm, "must be a finite positive number"))
})

test_that("coef_prior_sds builds the per-covariate prior scales", {
  f <- occJSDM:::coef_prior_sds(list(type = "fixed"), 0.3, NA, NULL, g_b = 2L, p = 4L)
  expect_identical(f, list(b = rep(0.3, 4), g = matrix(1, 2, 4), c = rep(1, 4)))
  cp <- rcp(list(b_prior = "grouped"), list(embeddingVars = paste0("e", 1:3)))
  lam <- c(1, 1, 0.5, 2, 0.1)
  s <- occJSDM:::coef_prior_sds(cp, 0.3, 0.2, lam, g_b = 2L, p = 5L)
  expect_equal(s$b, c(0.3, 0.3, 0.2 * 0.5, 0.2 * 2, 0.2 * 0.1))
  expect_equal(s$g, matrix(lam, 2, 5, byrow = TRUE))
  expect_equal(s$c, cp$c_base * lam)
})
