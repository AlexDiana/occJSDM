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
