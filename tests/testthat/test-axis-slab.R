# The cap on the embedding scales of the grouped prior (slab_scale).
# Design: dev/superpowers/specs/2026-11-25-axis-slab-design.md.
# axis_move_state() and axis_move_log_target() are in helper-scale-moves.R.

test_that("the capped standard deviations follow the design's formulas", {
  x <- axis_move_state(); cp <- x$cp; cp$slab_scale <- 0.8; emb <- 3:5
  sds <- occJSDM:::coef_prior_sds(cp, x$sigma_b, x$sigma_emb, x$lambda, 2L, 5L)
  lt <- 1 / sqrt(1 / x$lambda[emb]^2 + 1 / 0.8^2)
  om <- 1 / sqrt(1 / (x$sigma_emb * x$lambda[emb])^2 + 1 / 0.8^2)
  expect_equal(sds$b, c(rep(x$sigma_b, 2), om))
  expect_equal(sds$g, matrix(c(1, 1, lt), 2L, 5L, byrow = TRUE))
  expect_equal(sds$c, cp$c_base * c(1, 1, lt))
  eff <- occJSDM:::coef_effective_scales(cp, x$sigma_emb, x$lambda)
  expect_equal(eff$axis, c(1, 1, lt)); expect_equal(eff$leftover, c(NA, NA, om))
})

test_that("an infinite cap, or a record without one, gives the present standard deviations exactly", {
  x <- axis_move_state()
  want <- list(b = ifelse(x$cp$is_emb, x$sigma_emb, x$sigma_b) * x$lambda,
               g = matrix(x$lambda, 2L, 5L, byrow = TRUE), c = x$cp$c_base * x$lambda)
  expect_identical(occJSDM:::coef_prior_sds(x$cp, x$sigma_b, x$sigma_emb, x$lambda, 2L, 5L), want)
  cp <- x$cp; cp$slab_scale <- Inf
  expect_identical(occJSDM:::coef_prior_sds(cp, x$sigma_b, x$sigma_emb, x$lambda, 2L, 5L), want)
  expect_identical(occJSDM:::coef_effective_scales(cp, x$sigma_emb, x$lambda)$axis, x$lambda)
  expect_identical(occJSDM:::coef_C_start(cp, x$C, x$lambda), occJSDM:::coef_C_start(x$cp, x$C, x$lambda))
})

test_that("the capped scales stay finite, positive and under the cap at extreme axis scales", {
  cp <- axis_move_state()$cp; cp$slab_scale <- 2
  for (l in c(1e-150, 1e-12, 1e12, 1e150)) {
    sds <- occJSDM:::coef_prior_sds(cp, 0.7, 0.4, c(1, 1, l, l, l), 2L, 5L)
    v <- c(sds$b[3:5], sds$g[, 3:5], sds$c[3:5])
    expect_true(all(is.finite(v) & v > 0), info = format(l))
    expect_true(all(sds$b[3:5] <= 2) && all(sds$g[, 3:5] <= 2), info = format(l))
    if (l < 1) expect_equal(sds$b[3], 0.4 * l, tolerance = 1e-10)
    if (l > 1) expect_equal(sds$b[3], 2, tolerance = 1e-10)
  }
})

test_that("slab_scale is read and validated", {
  cols <- data.frame(column = c("m1", "e1", "e2"), variable = c("m1", "e1", "e2"), kind = "numeric")
  rd <- function(priors, params = list(embeddingVars = c("e1", "e2")), ...)
    occJSDM:::read_coef_prior(priors, params, 3L, 2L, cols, ...)
  expect_identical(rd(list())$slab_scale, Inf)          # the default until Task 5
  expect_identical(rd(list(slab_scale = 3))$slab_scale, 3)
  expect_identical(rd(list(slab_scale = Inf))$slab_scale, Inf)
  for (bad in list(0, -1, NA_real_, "2", c(1, 2)))
    expect_error(rd(list(slab_scale = bad)), "slab_scale must be a positive number or Inf")
  expect_error(rd(list(slab_scale = 2, b_prior = "grouped"), params = list()),
               "slab_scale requires embeddingVars")
  expect_error(rd(list(slab_scale = 2), params = list()), "slab_scale requires the grouped prior")
  expect_true(is.na(rd(list(b_prior = "grouped"), params = list())$slab_scale))
})
