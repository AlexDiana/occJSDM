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

# The scale draws under a cap (Task 2). Under a finite cap sigma_emb and the
# axis scales are drawn by univariate slice sampling on the log scale, from
# the densities f_sigma and f_k of the design's "Mathematics" section.

# The leftovers of axis_move_state(): B minus its trait and latent-trait parts.
slab_leftovers <- function(x) t(x$B) - x$Tr_b %*% x$G - x$A %*% x$C

# The quartiles of the density whose log is f, by numerical integration of
# exp(f) on a fine grid of the log scale (trapezoid rule), interpolated
# linearly between grid points.
slab_log_quartiles <- function(f, lo = -15, hi = 25, n = 40001L) {
  u <- seq(lo, hi, length.out = n)
  lf <- vapply(u, f, numeric(1))
  w <- exp(lf - max(lf))
  cdf <- c(0, cumsum((w[-1] + w[-n]) / 2))
  target <- c(.25, .5, .75) * cdf[n]
  i <- findInterval(target, cdf)
  u[i] + (target - cdf[i]) / (cdf[i + 1L] - cdf[i]) * (u[2] - u[1])
}

test_that("the slice draw samples a known density", {
  seed_both(3L); u <- 0; draws <- numeric(20000L)
  for (i in seq_along(draws))
    draws[i] <- u <- occJSDM:::slice_sample_log(function(u) dnorm(u, 1, 2, log = TRUE), u)
  share <- vapply(qnorm(c(.05, .25, .5, .75, .95), 1, 2), function(q) mean(draws < q), numeric(1))
  expect_lt(max(abs(share - c(.05, .25, .5, .75, .95))), 0.03)
})

test_that("the log densities of the capped scale draws are the target's", {
  x <- axis_move_state(); x$cp$slab_scale <- 0.8
  R <- t(x$B) - x$Tr_b %*% x$G - x$A %*% x$C; S <- nrow(R); emb <- 3:5; ls <- c(0.01, 0.3, 2, 40)
  for (k in emb) {
    target <- vapply(ls, function(l) { y <- x; y$lambda[k] <- l; axis_move_log_target(y) }, numeric(1))
    f <- vapply(ls, function(l) occJSDM:::axis_scale_logdens(
      log(l), sum(x$G[, k]^2) + sum(x$C[, k]^2) / x$cp$c_base[k]^2, nrow(x$G) + nrow(x$C),
      sum(R[, k]^2), S, x$sigma_emb, x$cp$axis_scale, 0.8), numeric(1))
    expect_equal(diff(f), diff(target), tolerance = 1e-8, info = format(k))
  }
  target <- vapply(ls, function(s) { y <- x; y$sigma_emb <- s; axis_move_log_target(y) }, numeric(1))
  f <- vapply(ls, function(s) occJSDM:::sigma_emb_logdens(
    log(s), colSums(R[, emb]^2), S, x$lambda[emb], x$cp$sigma_emb_scale, 0.8), numeric(1))
  expect_equal(diff(f), diff(target), tolerance = 1e-8)
})

test_that("the log densities stay finite far out on the log scale", {
  for (u in c(-350, -100, 0, 100, 350)) {
    expect_true(is.finite(occJSDM:::axis_scale_logdens(u, 0.3, 4L, 2, 7L, 0.4, 0.3, 2)), info = format(u))
    expect_true(is.finite(occJSDM:::sigma_emb_logdens(u, c(0.1, 2), 7L, c(0.5, 0.05), 1, 2)), info = format(u))
  }
})

test_that("a chain of capped axis scale draws has the quartiles of its conditional", {
  x <- axis_move_state(); R <- slab_leftovers(x); k <- 5L
  # Column 5's leftovers times 10: their standard deviation is then about 3.6,
  # far above what a cap of 0.8 allows, so the cap binds.
  ss_r <- 100 * sum(R[, k]^2)
  f <- function(u) occJSDM:::axis_scale_logdens(
    u, sum(x$G[, k]^2) + sum(x$C[, k]^2) / x$cp$c_base[k]^2, nrow(x$G) + nrow(x$C),
    ss_r, nrow(R), x$sigma_emb, x$cp$axis_scale, 0.8)
  q <- slab_log_quartiles(f)
  expect_gt(x$sigma_emb * exp(q[1]), 0.8)    # the leftover's scale is above the cap
  seed_both(4L); u <- log(x$lambda[k]); draws <- numeric(20000L)
  for (i in seq_along(draws)) draws[i] <- u <- occJSDM:::slice_sample_log(f, u)
  share <- vapply(q, function(v) mean(draws < v), numeric(1))
  # Tolerance. Successive draws are autocorrelated, so the Monte Carlo error
  # was measured: over 24 seeds of this chain the standard deviation of a
  # share was at most 0.0042 (0.0033, 0.0036, 0.0042 for the three), the
  # largest error in any seed was 0.011, and over 121 seeds the mean errors
  # were under 0.001. The tolerance is five of those standard deviations. With
  # the cap binding, G and C barely inform this scale, so the test is mostly of
  # the slice draw on a heavy tail; the form of f_k is pinned by the exact
  # test above. Dropping the Jacobian term u from f_k moves the shares by 0.15
  # to 0.30.
  expect_lt(max(abs(share - c(.25, .5, .75))), 0.021)
})

test_that("under a cap the scale updates leave the half-Cauchy priors of the scales unchanged", {
  # The successive-conditional check of test-grouped-prior.R ("the scale
  # updates leave the half-Cauchy priors of the scales unchanged"), with a cap
  # of 1 in the priors and the axis scale at 0.5, so that the cap binds often
  # (an axis scale is above 1 about 30% of the time, and the leftover's
  # sigma_emb * lambda more often still). The cap changes the prior of the
  # coefficients given the scales, not the prior of the scales themselves, so
  # each uncapped scale divided by its prior scale must still be half-Cauchy(1).
  # 20,000 sweeps, not the 50,000 of the uncapped test, to keep the run near
  # five seconds: the slice draws are slower than the inverse-gamma draws.
  S <- 3L; p <- 3L; gt <- 2L; g_b <- 2L; N <- 20000L; burn <- 1000L
  seed_both(1L)
  Tr_b <- cbind(1, rnorm(S))
  cp <- occJSDM:::read_coef_prior(
    list(sigma_b_scale = 0.8, sigma_emb_scale = 1.3, axis_scale = 0.5, slab_scale = 1),
    list(embeddingVars = c("e1", "e2")), S, gt,
    data.frame(column = c("m1", "e1", "e2"), variable = c("m1", "e1", "e2"), kind = "numeric"))
  expect_identical(cp$slab_scale, 1)
  start <- occJSDM:::coef_scales_start(cp, p)
  sigma_b <- cp$sigma_b_scale; sigma_emb <- start$sigma_emb; lambda <- start$lambda
  G <- matrix(rnorm(g_b * p), g_b, p); C <- matrix(rnorm(gt * p), gt, p)
  A <- matrix(rnorm(S * gt), S, gt)
  keep <- matrix(NA_real_, N + burn, 5L,
                 dimnames = list(NULL, c("sigma_b", "sigma_emb", "lambda_e1", "lambda_e2", "lambda_m1")))
  done <- tryCatch({
    for (i in seq_len(N + burn)) {
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      centre <- Tr_b %*% G + A %*% C
      B <- t(centre + matrix(rnorm(S * p), S, p) * rep(sds$b, each = S))
      sc <- occJSDM:::sample_coef_scales(t(B) - centre, G, C, cp, sigma_b, sigma_emb, lambda)
      sigma_b <- sc$sigma_b; sigma_emb <- sc$sigma_emb; lambda <- sc$lambda
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      gc <- occJSDM:::sample_GC(B, Tr_b, A, sds$b, sds$g, sds$c); G <- gc$G; C <- gc$C
      A <- occJSDM:::sample_A(B, C, Tr_b, G, sds$b)$A
      keep[i, ] <- c(sigma_b / cp$sigma_b_scale, sigma_emb / cp$sigma_emb_scale,
                     lambda[2:3] / cp$axis_scale, lambda[1])
    }
    TRUE
  }, error = function(e) conditionMessage(e))
  expect_identical(done, TRUE, label = "the chain ran to the end")
  if (!isTRUE(done)) return(invisible(NULL))
  keep <- keep[-seq_len(burn), ]
  expect_true(all(keep[, "lambda_m1"] == 1))
  expect_true(all(is.finite(keep) & keep > 0))
  quartiles <- tan(pi / 2 * c(.25, .5, .75))
  share <- vapply(quartiles, function(x) colMeans(keep[, 1:4] < x), numeric(4))
  # Tolerances. The Monte Carlo error was measured: over 24 seeds of this
  # chain the standard deviation of a share was at most 0.0164 for sigma_b,
  # 0.0146 for sigma_emb and 0.0181 for an axis scale, with no bias (no mean
  # error above 0.005), and the largest error of any of the twelve shares in
  # any seed was 0.056. Each tolerance is five of those standard deviations.
  # The uncapped draws run under the capped coefficient prior fail by far
  # more: they move the sigma_emb shares by 0.14 to 0.23 and the axis shares
  # by 0.24 to 0.37 (and this test then fails, as it did before the capped
  # draws were written).
  tolerance <- c(sigma_b = 0.085, sigma_emb = 0.075, lambda_e1 = 0.091, lambda_e2 = 0.091)
  for (nm in rownames(share)) for (j in 1:3)
    expect_lt(abs(share[nm, j] - c(.25, .5, .75)[j]), tolerance[[nm]],
              label = sprintf("error of the share of %s below its prior quartile %d", nm, j))
})

test_that("under a cap sigma_b is drawn as now, then sigma_emb, then each axis scale", {
  x <- axis_move_state(); cp <- x$cp; cp$slab_scale <- 0.8; R <- slab_leftovers(x); S <- nrow(R)
  emb <- 3:5; m_gc <- nrow(x$G) + nrow(x$C)
  axis_draws <- function(sigma_emb, lambda) {
    for (k in emb) lambda[k] <- exp(occJSDM:::slice_sample_log(function(u) occJSDM:::axis_scale_logdens(
      u, sum(x$G[, k]^2) + sum(x$C[, k]^2) / cp$c_base[k]^2, m_gc, sum(R[, k]^2), S,
      sigma_emb, cp$axis_scale, 0.8), log(lambda[k])))
    lambda
  }
  seed_both(42L)
  sigma_b <- occJSDM:::sample_scale_half_cauchy(sum(R[, 1:2]^2), S * 2, x$sigma_b, 1)
  sigma_emb <- exp(occJSDM:::slice_sample_log(function(v) occJSDM:::sigma_emb_logdens(
    v, colSums(R[, emb]^2), S, x$lambda[emb], 1, 0.8), log(x$sigma_emb)))
  want <- list(sigma_b = sigma_b, sigma_emb = sigma_emb, lambda = axis_draws(sigma_emb, x$lambda))
  seed_both(42L)
  expect_identical(occJSDM:::sample_coef_scales(R, x$G, x$C, cp, x$sigma_b, x$sigma_emb, x$lambda), want)
})

test_that("a held sigma_emb is not drawn under a cap", {
  x <- axis_move_state(); cp <- x$cp; cp$slab_scale <- 0.8; R <- slab_leftovers(x); S <- nrow(R)
  seed_both(42L)
  got <- occJSDM:::sample_coef_scales(R, x$G, x$C, cp, x$sigma_b, x$sigma_emb, x$lambda,
                                      hold_sigma_emb = TRUE)
  expect_identical(got$sigma_emb, x$sigma_emb)
  # No random number is drawn for it: the axis scales take the draws that
  # follow sigma_b's directly.
  seed_both(42L)
  sigma_b <- occJSDM:::sample_scale_half_cauchy(sum(R[, 1:2]^2), S * 2, x$sigma_b, 1)
  k <- 3L
  lambda3 <- exp(occJSDM:::slice_sample_log(function(u) occJSDM:::axis_scale_logdens(
    u, sum(x$G[, k]^2) + sum(x$C[, k]^2) / cp$c_base[k]^2, nrow(x$G) + nrow(x$C),
    sum(R[, k]^2), S, x$sigma_emb, cp$axis_scale, 0.8), log(x$lambda[k])))
  expect_identical(got$sigma_b, sigma_b)
  expect_identical(got$lambda[3], lambda3)
})

test_that("the capped draws return finite positive scales at extreme columns", {
  x <- axis_move_state(); cp <- x$cp; cp$slab_scale <- 0.8; R <- slab_leftovers(x)
  G <- x$G; C <- x$C
  R[, 4] <- 0; G[, 4] <- 0; C[, 4] <- 0     # a column with nothing in it
  R[, 5] <- 1e6 * R[, 5]                    # a column far above the cap
  seed_both(5L)
  draws <- replicate(200L, unlist(occJSDM:::sample_coef_scales(
    R, G, C, cp, x$sigma_b, x$sigma_emb, x$lambda)[c("sigma_emb", "lambda")]))
  expect_true(all(is.finite(draws) & draws > 0))
})

test_that("the capped draws run with no measured column and with no latent trait", {
  x <- axis_move_state(); cp <- x$cp; cp$slab_scale <- 0.8; R <- slab_leftovers(x)
  ok <- function(sc) all(is.finite(unlist(sc)) & unlist(sc) > 0)
  seed_both(6L)
  all_emb <- cp; all_emb$is_emb[] <- TRUE
  expect_true(ok(occJSDM:::sample_coef_scales(R, x$G, x$C, all_emb, x$sigma_b, x$sigma_emb, x$lambda)))
  no_trait <- cp; no_trait$c_base[] <- 1
  expect_true(ok(occJSDM:::sample_coef_scales(R, x$G, x$C[0, , drop = FALSE], no_trait,
                                              x$sigma_b, x$sigma_emb, x$lambda)))
})

test_that("with an infinite cap the scale draws are the present ones", {
  x <- axis_move_state(); R <- slab_leftovers(x)
  cp <- x$cp; cp$slab_scale <- Inf
  seed_both(42L)
  want <- occJSDM:::sample_coef_scales(R, x$G, x$C, x$cp, x$sigma_b, x$sigma_emb, x$lambda)
  seed_both(42L)
  expect_identical(occJSDM:::sample_coef_scales(R, x$G, x$C, cp, x$sigma_b, x$sigma_emb, x$lambda), want)
  seed_both(42L)
  held <- occJSDM:::sample_coef_scales(R, x$G, x$C, x$cp, x$sigma_b, x$sigma_emb, x$lambda, TRUE)
  seed_both(42L)
  expect_identical(occJSDM:::sample_coef_scales(R, x$G, x$C, cp, x$sigma_b, x$sigma_emb, x$lambda, TRUE), held)
})
