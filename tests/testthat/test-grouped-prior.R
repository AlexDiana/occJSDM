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
  expect_equal(got$axis_scale, 0.7 / 3)                         # 0.7 / number of embedding columns
  expect_equal(rcp(list(), list(embeddingVars = "e1"))$axis_scale, 0.7)
  expect_identical(rcp(list(axis_scale = 0.2), list(embeddingVars = paste0("e", 1:3)))$axis_scale, 0.2)
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

bbsl_inputs <- function(n = 40L, S = 5L, p = 4L, d = 1L, g = 1L, gt = 1L) {
  set.seed(11)
  list(k = matrix(rbinom(n * S, 1, .4) - .5, n, S), X = matrix(rnorm(n * p), n, p),
       Tr = matrix(rnorm(S * g), S, g), U = matrix(rnorm(n * d), n, d),
       G = matrix(rnorm(g * p), g, p), A = matrix(rnorm(S * gt), S, gt),
       C = matrix(rnorm(gt * p), gt, p), Omega = matrix(rexp(n * S) / 4, n, S),
       empty = list(Gs = matrix(0, g, 0), As = matrix(0, S, gt), Cs = matrix(0, gt, 0),
                    Ks = matrix(0, n, 0), Xs_centers = matrix(0, n, 0)))
}
draw_bbsl <- function(f, x, sd_b, Tr_b = NULL, seed = 5L) {
  seed_both(seed)
  f(x$k, x$X, x$Tr, x$U, x$G, x$A, x$C, sd_b, x$empty$Gs, x$empty$As, x$empty$Cs, 1,
    x$empty$Ks, x$empty$Xs_centers, x$Omega, "binary", 1, 0, Tr_b)
}

test_that("a recycled scalar and a constant vector give the same joint draw", {
  x <- bbsl_inputs()
  for (f in list(occJSDM:::sample_BBsL_cpp, occJSDM:::sample_BBsL_parallel))
    expect_identical(draw_bbsl(f, x, 0.3), draw_bbsl(f, x, rep(0.3, 4)))
})

test_that("each covariate's prior scale is applied to that covariate", {
  x <- bbsl_inputs()
  M <- t(x$Tr %*% x$G + x$A %*% x$C)
  for (f in list(occJSDM:::sample_BBsL_cpp, occJSDM:::sample_BBsL_parallel)) {
    for (tight in 1:4) {
      sd_b <- rep(1, 4); sd_b[tight] <- 1e-6
      B <- draw_bbsl(f, x, sd_b)$B
      expect_equal(B[tight, ], M[tight, ], tolerance = 1e-4)
      expect_gt(max(abs(B[-tight, ] - M[-tight, ])), 1e-2)
    }
  }
  expect_identical(draw_bbsl(occJSDM:::sample_BBsL_cpp, x, c(.2, .5, 1, 2)),
                   draw_bbsl(occJSDM:::sample_BBsL_parallel, x, c(.2, .5, 1, 2)))
})

test_that("the covariate side can use its own trait matrix", {
  x <- bbsl_inputs()
  Tr_b <- cbind(1, x$Tr); G_b <- rbind(rnorm(4), x$G); x2 <- x; x2$G <- G_b
  out <- draw_bbsl(occJSDM:::sample_BBsL_cpp, x2, rep(1e-6, 4), Tr_b = Tr_b)
  expect_equal(out$B, t(Tr_b %*% G_b + x$A %*% x$C), tolerance = 1e-4)
  expect_equal(out$Bt, t(out$B) - Tr_b %*% G_b - x$A %*% x$C)
})

test_that("the joint draw rejects a prior scale vector of the wrong length", {
  x <- bbsl_inputs()
  for (f in list(occJSDM:::sample_BBsL_cpp, occJSDM:::sample_BBsL_parallel))
    for (bad in list(c(.3, .3), rep(.3, 5), numeric(0)))
      expect_error(draw_bbsl(f, x, bad), "sd_b must have length 1 or one entry per covariate")
})

test_that("the R reference draw takes the same prior scales and trait matrix", {
  # sample_BBsL() draws through the same C++ kernel, so with both generators
  # seeded alike it must agree with sample_BBsL_cpp() up to rounding.
  x <- bbsl_inputs()
  Tr_b <- cbind(1, x$Tr); x2 <- x; x2$G <- rbind(rnorm(4), x$G)
  cases <- list(list(x, 0.3, NULL), list(x, c(.2, .5, 1, 2), NULL),
                list(x2, c(.2, .5, 1, 2), Tr_b))
  for (case in cases) {
    cpp <- draw_bbsl(occJSDM:::sample_BBsL_cpp, case[[1]], case[[2]], Tr_b = case[[3]])
    ref <- if (is.null(case[[3]])) {
      # Tr_b left at its default, which is Tr.
      seed_both(5L)
      with(case[[1]], occJSDM:::sample_BBsL(k, X, Tr, U, G, A, C, case[[2]], empty$Gs, empty$As,
        empty$Cs, 1, empty$Ks, empty$Xs_centers, Omega, "binary", 1, 0))
    } else {
      draw_bbsl(occJSDM:::sample_BBsL, case[[1]], case[[2]], Tr_b = case[[3]])
    }
    cpp$B0 <- drop(cpp$B0)
    expect_equal(ref, cpp, tolerance = 1e-12)
  }
})

test_that("the range log-weights use one prior scale per covariate", {
  # Dense reference: for each species, 0.5 * h' Q^-1 h - 0.5 * log det Q, with
  # Q = Z' diag(omega) Z + P and h = Z' kappa + P m, summed over species, plus
  # the Gamma(a_l_s, b_l_s) log density of the range at that grid point.
  sim <- simulate_fixture(model = "binary", useSpatField = TRUE, seed = 9L)
  setup <- range_reference_setup(sim)          # in helper-grouped-prior.R
  sd_b <- c(0.2, 1.5)
  got <- with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, sd_b, sigma_bs,
    kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1))
  expect_equal(got, range_reference(setup, sd_b), tolerance = 1e-8)
  expect_identical(
    with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, 0.3, sigma_bs,
      kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1)),
    with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, c(0.3, 0.3), sigma_bs,
      kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1)))
})

test_that("the grouped range log-weights use one prior scale per covariate", {
  # Two visits per site share that site's spatial basis row, the case the
  # grouped crossproducts exist for. The dense reference does not group.
  sim <- simulate_fixture(model = "occupancy", useSpatField = TRUE, seed = 9L)
  setup <- range_reference_setup(sim)
  expect_gt(anyDuplicated(setup$location), 0L)
  sd_b <- c(0.2, 1.5)
  grouped <- function(sd_b) with(setup, occJSDM:::spatial_range_logweights_grouped(X, U, M_B,
    M_Bs, sd_b, sigma_bs, kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1, location))
  got <- grouped(sd_b)
  expect_equal(got, range_reference(setup, sd_b), tolerance = 1e-8)
  # Repeated labels route the dense entry point to the grouped function.
  expect_identical(with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, sd_b,
    sigma_bs, kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1, location = location)), got)
  expect_identical(grouped(0.3), grouped(c(0.3, 0.3)))
})

scale_inputs <- function() {
  set.seed(1); S <- 7L; gb <- 2L; gt <- 2L
  list(S = S, gb = gb, gt = gt, R = matrix(rnorm(S * 5), S, 5),
       G = matrix(rnorm(gb * 5), gb, 5), C = matrix(rnorm(gt * 5), gt, 5),
       cp = list(type = "grouped", is_emb = c(FALSE, FALSE, TRUE, TRUE, TRUE),
                 sigma_b_scale = 1, sigma_emb_scale = 1, axis_scale = 1 / sqrt(3),
                 c_base = c(1, 1, rep(1 / sqrt(gt), 3))),
       sigma_b = 0.7, sigma_emb = 0.4, lambda = c(1, 1, 0.5, 0.2, 0.9))
}
hc <- function(ss, n, sigma, scale) {
  aux <- occJSDM:::rinvgamma_cpp(1, 1 / sigma^2 + 1 / scale^2)
  sqrt(occJSDM:::rinvgamma_cpp((n + 1) / 2, ss / 2 + 1 / aux))
}

test_that("coefficient scales are drawn from their half-Cauchy conditionals", {
  x <- scale_inputs()
  seed_both(42)
  want_b <- hc(sum(x$R[, 1:2]^2), x$S * 2, x$sigma_b, 1)
  want_e <- hc(sum(sweep(x$R[, 3:5], 2, x$lambda[3:5], "/")^2), x$S * 3, x$sigma_emb, 1)
  want_l <- vapply(3:5, function(k) hc(
    sum(x$R[, k]^2) / want_e^2 + sum(x$G[, k]^2) + sum(x$C[, k]^2) / x$cp$c_base[k]^2,
    x$S + x$gb + x$gt, x$lambda[k], x$cp$axis_scale), numeric(1))
  seed_both(42)
  got <- with(x, occJSDM:::sample_coef_scales(R, G, C, cp, sigma_b, sigma_emb, lambda))
  expect_equal(got, list(sigma_b = want_b, sigma_emb = want_e, lambda = c(1, 1, want_l)))
})

test_that("with no embedding columns only the measured scale moves", {
  x <- scale_inputs(); x$cp$is_emb[] <- FALSE; x$cp$c_base[] <- 1
  seed_both(42); want <- hc(sum(x$R^2), x$S * 5, x$sigma_b, 1)
  seed_both(42)
  got <- with(x, occJSDM:::sample_coef_scales(R, G, C, cp, sigma_b, NA_real_, rep(1, 5)))
  expect_equal(got, list(sigma_b = want, sigma_emb = NA_real_, lambda = rep(1, 5)))
})

test_that("an axis scale tracks the size of its column", {
  x <- scale_inputs(); x$R[, 4] <- 0; x$G[, 4] <- 0; x$C[, 4] <- 0
  x$R[, 5] <- 30 * x$R[, 5]
  draws <- replicate(400, with(x, occJSDM:::sample_coef_scales(R, G, C, cp, sigma_b, sigma_emb, lambda))$lambda)
  expect_lt(median(draws[4, ]), 0.2)
  expect_gt(median(draws[5, ]), 5 * median(draws[4, ]))
})

test_that("sample_GC uses a prior and a noise scale per covariate", {
  set.seed(2); S <- 9L; p <- 3L; gb <- 2L; gt <- 2L
  B <- matrix(rnorm(p * S), p, S); Tr_b <- cbind(1, rnorm(S)); A <- matrix(rnorm(S * gt), S, gt)
  sd_b <- c(.3, 1, 2); sd_g <- matrix(c(.5, 1, 2), gb, p, byrow = TRUE); sd_c <- c(.2, .7, 1.5)
  seed_both(7); got <- occJSDM:::sample_GC(B, Tr_b, A, sd_b, sd_g, sd_c)
  seed_both(7)
  for (k in 1:p) {
    want <- occJSDM:::sampleBuniv(cbind(Tr_b, A), diag(c(sd_g[, k]^2, rep(sd_c[k]^2, gt))),
                                  rep(0, gb + gt), B[k, ], sd_b[k])
    expect_equal(c(got$G[, k], got$C[, k]), as.vector(want))
  }
  expect_equal(got$Bt, t(B) - Tr_b %*% got$G - A %*% got$C)
})

test_that("sample_A weights covariates by their noise scale", {
  set.seed(3); S <- 4L; p <- 6L; gt <- 2L
  B <- matrix(rnorm(p * S), p, S); C <- matrix(rnorm(gt * p), gt, p)
  Tr_b <- matrix(1, S, 1); G <- matrix(rnorm(p), 1, p); sd_b <- c(.2, .2, 1, 1, 3, 3)
  Btilde <- B - t(Tr_b %*% G)
  seed_both(8); got <- occJSDM:::sample_A(B, C, Tr_b, G, sd_b)
  seed_both(8)
  for (s in 1:S) expect_equal(got$A[s, ], as.vector(occJSDM:::sampleBuniv(
    t(C) / sd_b, diag(gt), rep(0, gt), Btilde[, s] / sd_b, 1)))
  seed_both(8); same <- occJSDM:::sample_A(B, C, Tr_b, G, rep(.4, p))
  seed_both(8)
  for (s in 1:S) expect_identical(same$A[s, ], as.vector(occJSDM:::sampleBuniv(
    t(C), diag(gt), rep(0, gt), Btilde[, s], .4)))
})
