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
  expect_identical(rd(list())$slab_scale, 2)
  expect_identical(rd(list(), model = "continuous")$slab_scale, Inf)
  expect_identical(rd(list(slab_scale = 5), model = "continuous")$slab_scale, 5)
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

test_that("the slice draw stays exact when its stepping-out limit binds", {
  # With max_steps = 2 the interval can step out by one width in all, on one
  # side chosen at random, while N(0, 2^2) needs several widths, so the limit
  # binds on most draws. The random split of the limit between the two ends
  # (Neal 2003, Figure 3) is what keeps the draw exact then. Tolerances, one
  # per quantile: over 72 seeds of 200,000 draws the standard deviation of a
  # share was 0.0020, 0.0045, 0.0050, 0.0040 and 0.0017, with no bias (no
  # mean error above 1.2 standard errors), and each tolerance is five of
  # them, the larger of the two at symmetric quantiles. A fixed cap on each
  # side instead of the split (one step or two) moves the shares by about
  # 0.02 at the 5% and 95% points and 0.04 at the quartiles.
  seed_both(8L); u <- 0; draws <- numeric(200000L)
  for (i in seq_along(draws))
    draws[i] <- u <- occJSDM:::slice_sample_log(function(u) dnorm(u, 0, 2, log = TRUE), u,
                                                max_steps = 2L)
  p <- c(.05, .25, .5, .75, .95)
  share <- vapply(qnorm(p, 0, 2), function(q) mean(draws < q), numeric(1))
  tolerance <- c(0.010, 0.023, 0.025, 0.023, 0.010)
  for (j in seq_along(p))
    expect_lt(abs(share[j] - p[j]), tolerance[j],
              label = sprintf("error of the share below the %g quantile", p[j]))
})

test_that("the vector slice draw samples each element's density, also where stepping out binds", {
  # Three elements with different Normal densities, drawn together as
  # slice_sample_log_vec() draws the axis scales. Each element must have its
  # own distribution, as the scalar draw has (the two tests above, same
  # tolerances): with the default stepping out, and with max_steps = 2, where
  # the limit binds on most draws and only the random split keeps the draw
  # exact.
  mu <- c(1, -3, 0); sd <- c(2, 0.3, 2)
  f <- function(u, idx) dnorm(u, mu[idx], sd[idx], log = TRUE)
  p <- c(.05, .25, .5, .75, .95)
  for (case in list(list(n = 20000L, max_steps = 50L, tol = rep(0.03, 5)),
                    list(n = 200000L, max_steps = 2L, tol = c(0.010, 0.023, 0.025, 0.023, 0.010)))) {
    seed_both(9L); u <- c(0, 0, 0); draws <- matrix(NA_real_, case$n, 3L)
    for (i in seq_len(case$n))
      draws[i, ] <- u <- occJSDM:::slice_sample_log_vec(f, u, max_steps = case$max_steps)
    for (j in 1:3) {
      share <- vapply(qnorm(p, mu[j], sd[j]), function(q) mean(draws[, j] < q), numeric(1))
      expect_true(all(abs(share - p) < case$tol),
                  label = sprintf("element %d at max_steps = %d: shares %s", j, case$max_steps,
                                  paste(format(share, digits = 3), collapse = " ")))
    }
  }
  expect_error(occJSDM:::slice_sample_log_vec(function(u, idx) ifelse(idx == 2L, -Inf, 0), c(0, 0)),
               "a slice draw must start at a point where the log density is finite")
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
  # The axis scales are drawn together, one slice draw each
  # (slice_sample_log_vec()), given the sigma_emb passed in.
  axis_draws <- function(sigma_emb, lambda) {
    ss_gc <- colSums(x$G[, emb]^2) + colSums(x$C[, emb]^2) / cp$c_base[emb]^2
    ss_r <- colSums(R[, emb]^2)
    lambda[emb] <- exp(occJSDM:::slice_sample_log_vec(function(u, idx) occJSDM:::axis_scale_logdens(
      u, ss_gc[idx], m_gc, ss_r[idx], S, sigma_emb, cp$axis_scale, 0.8), log(lambda[emb])))
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
  emb <- 3:5
  ss_gc <- colSums(x$G[, emb]^2) + colSums(x$C[, emb]^2) / cp$c_base[emb]^2
  ss_r <- colSums(R[, emb]^2)
  lambda_emb <- exp(occJSDM:::slice_sample_log_vec(function(u, idx) occJSDM:::axis_scale_logdens(
    u, ss_gc[idx], nrow(x$G) + nrow(x$C), ss_r[idx], S, x$sigma_emb, cp$axis_scale, 0.8),
    log(x$lambda[emb])))
  expect_identical(got$sigma_b, sigma_b)
  expect_identical(got$lambda[emb], lambda_emb)
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

# The two embedding moves under a cap (Task 3): the design's "The axis scale
# move" and "The leftover scale move". Their exact-ratio tests are in
# test-axis-scale-move.R and test-sigma-emb-move.R.

test_that("with an infinite cap both moves are the present ones", {
  # A record with slab_scale = Inf, with the new arguments given, must return
  # what a record without the element returns without them, and leave the
  # random number stream where it leaves it.
  for (seed in 1:5) {
    x <- axis_move_state(seed); R <- slab_leftovers(x)
    cp <- x$cp; cp$slab_scale <- Inf
    for (constant in c(FALSE, TRUE)) {
      if (constant) { x$kappa[] <- 0; x$Omega[] <- 0 }   # most moves accepted
      label <- paste("seed", seed, if (constant) "constant likelihood")
      seed_both(seed)
      want <- occJSDM:::axis_scale_move(x$B, x$G, x$C, x$lambda, x$X, x$kappa, x$Omega, x$eta,
                                        x$cp, 0.5)
      want_next <- runif(1)
      seed_both(seed)
      got <- occJSDM:::axis_scale_move(x$B, x$G, x$C, x$lambda, x$X, x$kappa, x$Omega, x$eta,
                                       cp, 0.5, R = R, sigma_emb = x$sigma_emb)
      expect_identical(got, want, info = label)
      expect_identical(runif(1), want_next, info = label)
      seed_both(seed)
      want <- occJSDM:::sigma_emb_move(x$B, R, x$sigma_emb, x$X, x$kappa, x$Omega, x$eta, x$cp, 0.5)
      want_next <- runif(1)
      seed_both(seed)
      got <- occJSDM:::sigma_emb_move(x$B, R, x$sigma_emb, x$X, x$kappa, x$Omega, x$eta, cp, 0.5,
                                      lambda = x$lambda)
      expect_identical(got, want, info = label)
      expect_identical(runif(1), want_next, info = label)
    }
  }
})

test_that("under a cap a move at a column far above the cap returns a finite ratio", {
  # Column 5 far above the cap of 0.8, and column 4 far below it: the ratios
  # of the capped scales must be computed without overflow or cancellation.
  x <- axis_move_state(); x$cp$slab_scale <- 0.8
  x$lambda[5] <- 1e6; x$lambda[4] <- 1e-12
  R <- slab_leftovers(x)
  for (s in 1:20) {
    seed_both(s)
    mv <- with(x, occJSDM:::axis_scale_move(B, G, C, lambda, X, kappa, Omega, eta, cp, 2,
                                            R = R, sigma_emb = sigma_emb))
    expect_true(all(is.finite(mv$log_ratio[3:5])), info = format(s))
    expect_true(all(is.finite(mv$B)) && all(is.finite(mv$G)) && all(is.finite(mv$C)), info = format(s))
    sg <- with(x, occJSDM:::sigma_emb_move(B, R, sigma_emb, X, kappa, Omega, eta, cp, 2, lambda = lambda))
    expect_true(is.finite(sg$log_ratio), info = format(s))
    expect_true(all(is.finite(sg$B)) && all(is.finite(sg$eta)), info = format(s))
  }
  # An accepted move scales each part by the ratio of its capped scales, which
  # far above the cap is within about 1e-10 of 1 (row 5 of B barely moves)
  # and far below it within rounding of exp(delta) (row 4 scales as without
  # a cap).
  x$kappa[] <- 0; x$Omega[] <- 0                     # a constant likelihood
  eff0 <- occJSDM:::coef_effective_scales(x$cp, x$sigma_emb, x$lambda)
  seen <- c(FALSE, FALSE)
  for (s in 1:20) {
    seed_both(s)
    mv <- with(x, occJSDM:::axis_scale_move(B, G, C, lambda, X, kappa, Omega, eta, cp, 2,
                                            R = R, sigma_emb = sigma_emb))
    eff1 <- occJSDM:::coef_effective_scales(x$cp, x$sigma_emb, mv$lambda)
    rho_g <- eff1$axis / eff0$axis; rho_e <- eff1$leftover / eff0$leftover
    for (k in which(mv$accepted)) {
      expect_equal(mv$B[k, ], rho_g[k] * (x$B[k, ] - R[, k]) + rho_e[k] * R[, k],
                   tolerance = 1e-12, info = paste(s, k))
      expect_equal(mv$G[, k], rho_g[k] * x$G[, k], tolerance = 1e-12, info = paste(s, k))
    }
    if (mv$accepted[5]) expect_lt(max(abs(mv$B[5, ] / x$B[5, ] - 1)), 1e-8)
    if (mv$accepted[4]) expect_equal(mv$B[4, ], mv$lambda[4] / x$lambda[4] * x$B[4, ],
                                     tolerance = 1e-12, info = format(s))
    seen <- seen | mv$accepted[4:5]
  }
  expect_true(all(seen))                             # both checks were made
})

# Exact replays of the two moves under a cap, for one seed: the expected log
# ratio, decision and state, from the written-out target density as in the
# exact-ratio tests of the two moves.
slab_axis_replay <- function(x, seed, step = 0.5) {
  emb <- which(x$cp$is_emb)
  seed_both(seed)
  delta <- rnorm(length(emb), 0, step); logu <- log(runif(length(emb)))
  want <- x; log_ratio <- rep(NA_real_, ncol(x$X)); accepted <- rep(FALSE, ncol(x$X))
  for (j in seq_along(emb)) {
    k <- emb[j]; prop <- want
    eff0 <- occJSDM:::coef_effective_scales(want$cp, want$sigma_emb, want$lambda)
    prop$lambda[k] <- want$lambda[k] * exp(delta[j])
    eff1 <- occJSDM:::coef_effective_scales(prop$cp, prop$sigma_emb, prop$lambda)
    rho_g <- eff1$axis[k] / eff0$axis[k]; rho_e <- eff1$leftover[k] / eff0$leftover[k]
    prop$G[, k] <- want$G[, k] * rho_g; prop$C[, k] <- want$C[, k] * rho_g
    R <- t(want$B) - want$Tr_b %*% want$G - want$A %*% want$C; R[, k] <- R[, k] * rho_e
    prop$B <- t(want$Tr_b %*% prop$G + want$A %*% prop$C + R)
    log_ratio[k] <- axis_move_log_target(prop) - axis_move_log_target(want) +
      (nrow(want$G) + nrow(want$C)) * log(rho_g) + ncol(want$B) * log(rho_e)
    if (logu[j] < log_ratio[k]) { want <- prop; accepted[k] <- TRUE }
  }
  list(log_ratio = log_ratio, accepted = accepted, state = want)
}

slab_sigma_replay <- function(x, seed, step = 0.5) {
  emb <- which(x$cp$is_emb)
  seed_both(seed)
  delta <- rnorm(1, 0, step); logu <- log(runif(1))
  prop <- x; prop$sigma_emb <- x$sigma_emb * exp(delta)
  rho <- occJSDM:::coef_effective_scales(x$cp, prop$sigma_emb, x$lambda)$leftover[emb] /
    occJSDM:::coef_effective_scales(x$cp, x$sigma_emb, x$lambda)$leftover[emb]
  R <- slab_leftovers(x); R[, emb] <- R[, emb] * rep(rho, each = nrow(R))
  prop$B <- t(x$Tr_b %*% x$G + x$A %*% x$C + R)
  log_ratio <- axis_move_log_target(prop) - axis_move_log_target(x) + ncol(x$B) * sum(log(rho))
  accepted <- logu < log_ratio
  list(log_ratio = log_ratio, accepted = accepted, state = if (accepted) prop else x)
}

test_that("under a cap the moves run with no measured column and no latent trait", {
  # axis_move_state() with every column an embedding column, and C and A of
  # no rows and no columns: B is the community part plus the leftovers.
  x <- axis_move_state(); x$cp$slab_scale <- 0.8
  R <- slab_leftovers(x)
  x$cp$is_emb[] <- TRUE; x$cp$c_base[] <- 1; x$lambda <- c(0.3, 1.2, 0.5, 0.05, 0.9)
  x$C <- x$C[0, , drop = FALSE]; x$A <- x$A[, 0, drop = FALSE]
  x$B <- t(x$Tr_b %*% x$G + R); x$eta <- x$rest + x$X %*% x$B
  seed_both(7L)
  mv <- with(x, occJSDM:::axis_scale_move(B, G, C, lambda, X, kappa, Omega, eta, cp, 0.5,
                                          R = slab_leftovers(x), sigma_emb = sigma_emb))
  want <- slab_axis_replay(x, 7L)
  expect_true(all(is.finite(mv$log_ratio)))
  expect_equal(mv$log_ratio, want$log_ratio, tolerance = 1e-8)
  expect_identical(mv$accepted, want$accepted)
  expect_equal(mv[c("B", "G", "C", "lambda")], want$state[c("B", "G", "C", "lambda")])
  expect_identical(dim(mv$C), c(0L, 5L))
  seed_both(8L)
  sg <- with(x, occJSDM:::sigma_emb_move(B, slab_leftovers(x), sigma_emb, X, kappa, Omega, eta, cp,
                                         0.5, lambda = lambda))
  want <- slab_sigma_replay(x, 8L)
  expect_true(is.finite(sg$log_ratio))
  expect_equal(sg$log_ratio, want$log_ratio, tolerance = 1e-8)
  expect_identical(sg$accepted, want$accepted)
  expect_equal(sg[c("B", "sigma_emb")], want$state[c("B", "sigma_emb")])
})

test_that("under a cap the scale updates with all three moves leave the half-Cauchy priors unchanged", {
  # The test of test-sigma-b-move.R of the same kind ("the scale updates with
  # all three scale moves leave the half-Cauchy priors unchanged"), with a cap
  # of 1 and the moves given what they need under it: the sigma_emb move the
  # axis scales, the axis scale move the leftovers of the B the sigma_emb move
  # returns and the current sigma_emb. With a constant likelihood every step
  # leaves the prior invariant, and the cap changes the prior of the
  # coefficients given the scales, not that of the scales, so each uncapped
  # scale divided by its prior scale must still be half-Cauchy(1). 20,000
  # sweeps, as in the test of the capped draws above, for run time.
  S <- 3L; p <- 3L; gt <- 2L; g_b <- 2L; n <- 4L; N <- 20000L; burn <- 1000L
  seed_both(1L)
  Tr_b <- cbind(1, rnorm(S)); X <- matrix(rnorm(n * p), n, p)
  zero <- matrix(0, n, S)
  cp <- occJSDM:::read_coef_prior(
    list(sigma_b_scale = 0.8, sigma_emb_scale = 1.3, axis_scale = 0.5, slab_scale = 1),
    list(embeddingVars = c("e1", "e2")), S, gt,
    data.frame(column = c("m1", "e1", "e2"), variable = c("m1", "e1", "e2"), kind = "numeric"))
  start <- occJSDM:::coef_scales_start(cp, p)
  sigma_b <- cp$sigma_b_scale; sigma_emb <- start$sigma_emb; lambda <- start$lambda
  G <- matrix(rnorm(g_b * p), g_b, p); C <- matrix(rnorm(gt * p), gt, p)
  A <- matrix(rnorm(S * gt), S, gt)
  keep <- matrix(NA_real_, N + burn, 4L,
                 dimnames = list(NULL, c("sigma_b", "sigma_emb", "lambda_e1", "lambda_e2")))
  accepted <- 0L; accepted_b <- c(0L, 0L); accepted_axis <- 0L
  done <- tryCatch({
    for (i in seq_len(N + burn)) {
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      centre <- Tr_b %*% G + A %*% C
      B <- t(centre + matrix(rnorm(S * p), S, p) * rep(sds$b, each = S))
      sc <- occJSDM:::sample_coef_scales(t(B) - centre, G, C, cp, sigma_b, sigma_emb, lambda)
      sigma_b <- sc$sigma_b; sigma_emb <- sc$sigma_emb; lambda <- sc$lambda
      sb <- occJSDM:::sigma_b_move(B, t(B) - centre, sigma_b, X, zero, zero, zero, cp,
                                   c(0.5, occJSDM:::SIGMA_B_MOVE_LARGE_STEP))
      B <- sb$B; sigma_b <- sb$sigma_b
      accepted_b <- accepted_b + sb$accepted
      sg <- occJSDM:::sigma_emb_move(B, t(B) - centre, sigma_emb, X, zero, zero, zero, cp, 0.5,
                                     lambda = lambda)
      B <- sg$B; sigma_emb <- sg$sigma_emb
      accepted <- accepted + sg$accepted
      mv <- occJSDM:::axis_scale_move(B, G, C, lambda, X, zero, zero, sg$eta, cp, 0.5,
                                      R = t(B) - centre, sigma_emb = sigma_emb)
      B <- mv$B; G <- mv$G; C <- mv$C; lambda <- mv$lambda
      accepted_axis <- accepted_axis + sum(mv$accepted)
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      gc <- occJSDM:::sample_GC(B, Tr_b, A, sds$b, sds$g, sds$c); G <- gc$G; C <- gc$C
      A <- occJSDM:::sample_A(B, C, Tr_b, G, sds$b)$A
      keep[i, ] <- c(sigma_b / cp$sigma_b_scale, sigma_emb / cp$sigma_emb_scale,
                     lambda[2:3] / cp$axis_scale)
    }
    TRUE
  }, error = function(e) conditionMessage(e))
  expect_identical(done, TRUE, label = "the chain ran to the end (it does not when the scales run away)")
  if (!isTRUE(done)) return(invisible(NULL))
  keep <- keep[-seq_len(burn), ]
  expect_true(all(is.finite(keep) & keep > 0))
  expect_gt(accepted / (N + burn), 0.5)                # with no data most moves are accepted
  expect_true(all(accepted_b / (N + burn) > 0.5))      # both sigma_b proposals too
  expect_gt(accepted_axis / (2 * (N + burn)), 0.5)     # and the axis scale moves
  quartiles <- tan(pi / 2 * c(.25, .5, .75))
  share <- vapply(quartiles, function(x) colMeans(keep < x), numeric(4))
  # Tolerances. The Monte Carlo error was measured: over 24 seeds of this
  # chain (25 November 2026, on mains power) the standard deviation of a
  # share was at most 0.0078 for sigma_b, 0.0098 for sigma_emb and 0.0128 for
  # an axis scale, with no bias (no mean error above 1.5 standard errors),
  # and the largest error of any share in any seed was 0.026. Each tolerance
  # is five of those standard deviations. Running the two embedding moves in
  # their uncapped form under the capped prior moves the axis shares by up to
  # 0.04 to 0.08 (3 seeds), so this test catches that only in part; the form
  # of the two moves is pinned exactly by their exact-ratio tests.
  tolerance <- c(sigma_b = 0.039, sigma_emb = 0.049, lambda_e1 = 0.064, lambda_e2 = 0.064)
  for (nm in rownames(share)) for (j in 1:3)
    expect_lt(abs(share[nm, j] - c(.25, .5, .75)[j]), tolerance[[nm]],
              label = sprintf("error of the share of %s below its prior quartile %d", nm, j))
})

test_that("a grouped binary fit with a cap keeps every stored axis scale under it", {
  fit <- sigma_move_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L),
                        listPriors = list(slab_scale = 0.5))
  expect_identical(fit$infos$coef_prior$slab_scale, 0.5)
  jo <- fit$results_output$jsdm_output
  expect_identical(dim(jo$lambda_output), c(8L, 20L, 2L))
  expect_true(all(is.finite(jo$lambda_output) & jo$lambda_output > 0))
  expect_lte(max(jo$lambda_output), 0.5)
  for (nm in c("sigmaemb_output", "sigmab_output", "B_output", "G_output", "C_output"))
    expect_true(all(is.finite(jo[[nm]])), info = nm)
  # The three moves ran and recorded their acceptance, as without a cap.
  am <- fit$infos$axis_move
  expect_identical(dim(am$acceptance), c(8L, 2L))
  expect_true(all(am$acceptance >= 0 & am$acceptance <= 1))
  expect_gt(mean(am$acceptance), 0)
  expect_length(am$sigma_emb_acceptance, 2L)
  expect_true(all(am$sigma_emb_acceptance >= 0 & am$sigma_emb_acceptance <= 1))
  expect_identical(dim(am$sigma_b_acceptance), c(2L, 2L))
  expect_true(all(am$sigma_b_acceptance >= 0 & am$sigma_b_acceptance <= 1))
  expect_false(all(am$step == 0.5)); expect_false(all(am$sigma_emb_step == 0.5))
})

test_that("grouped occupancy and two-stage fits run with a cap", {
  for (model in c("occupancy", "two_stage")) {
    dat <- simulate_grouped_fixture(model = model, p_emb = 6L, S = 8L)$data_list
    fit <- suppressMessages(suppressWarnings(runOccJSDM(dat,
      listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
      listPriors = list(slab_scale = 0.5),
      occCovariates = grouped_occ_covariates(dat), MCMCparams = FIXTURE_MCMC)))
    expect_identical(fit$infos$coef_prior$slab_scale, 0.5, info = model)
    expect_false("axis_move" %in% names(fit$infos), info = model)   # the moves are binary only
    jo <- fit$results_output$jsdm_output
    expect_true(all(is.finite(jo$lambda_output) & jo$lambda_output > 0), info = model)
    expect_lte(max(jo$lambda_output), 0.5)
    expect_true(all(is.finite(jo$B_output)), info = model)
  }
})
