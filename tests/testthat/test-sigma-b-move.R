# The scale-expansion move on sigma_b, the measured group's scale under the
# grouped prior (sigma_b_move() in R/coef-prior.R). Design:
# dev/superpowers/specs/2026-11-22-sigma-b-move-design.md.
# axis_move_state() and axis_move_log_target() are in helper-scale-moves.R.

sigma_b_leftovers <- function(x, B = x$B) t(B) - x$Tr_b %*% x$G - x$A %*% x$C

run_sigma_b_move <- function(x, step = 0.5) with(x, occJSDM:::sigma_b_move(
  B, sigma_b_leftovers(x), sigma_b, X, kappa, Omega, eta, cp, step))

# The state after scaling the measured leftovers and sigma_b by r.
sigma_b_scaled <- function(x, r) {
  meas <- which(!x$cp$is_emb)
  R <- sigma_b_leftovers(x); R[, meas] <- R[, meas] * r
  x$B <- t(x$Tr_b %*% x$G + x$A %*% x$C + R); x$sigma_b <- x$sigma_b * r
  x$eta <- x$rest + x$X %*% x$B
  x
}

test_that("the sigma_b move accepts on the exact ratio of the target density", {
  # Replays the move's two random numbers and decides from the written-out
  # target density. The ratio is that of the target densities times r to the
  # power of the number of leftovers scaled (the Jacobian).
  n_accept <- 0L; n_reject <- 0L
  for (seed in 1:40) {
    x <- axis_move_state(seed)
    m <- ncol(x$B) * sum(!x$cp$is_emb)
    step <- c(0.5, 0.05)[1L + seed %% 2L]
    seed_both(300L + seed); got <- run_sigma_b_move(x, step)
    seed_both(300L + seed); delta <- rnorm(1, 0, step); logu <- log(runif(1))
    prop <- sigma_b_scaled(x, exp(delta))
    log_ratio <- axis_move_log_target(prop) - axis_move_log_target(x) + m * delta
    expect_equal(got$log_ratio, log_ratio, tolerance = 1e-8, info = paste("seed", seed))
    expect_identical(got$accepted, logu < log_ratio, info = paste("seed", seed))
    if (logu < log_ratio) {
      expect_equal(got[c("B", "sigma_b", "eta")], prop[c("B", "sigma_b", "eta")])
      expect_equal(got$R, sigma_b_leftovers(prop))
      n_accept <- n_accept + 1L
    } else {
      # A rejected move returns the state bit for bit.
      expect_identical(got[c("B", "sigma_b", "eta")], x[c("B", "sigma_b", "eta")])
      n_reject <- n_reject + 1L
    }
  }
  expect_gt(n_accept, 5L); expect_gt(n_reject, 5L)   # both outcomes were exercised
})

test_that("two proposals are made in order, the second from the state the first leaves", {
  n_first <- 0L
  for (seed in 1:30) {
    x <- axis_move_state(seed); m <- ncol(x$B) * sum(!x$cp$is_emb)
    step <- c(0.05, 1)
    seed_both(500L + seed); got <- run_sigma_b_move(x, step)
    seed_both(500L + seed); delta <- rnorm(2, 0, step); logu <- log(runif(2))
    cur <- x; want_ratio <- numeric(2); want_acc <- logical(2)
    for (j in 1:2) {
      prop <- sigma_b_scaled(cur, exp(delta[j]))
      want_ratio[j] <- axis_move_log_target(prop) - axis_move_log_target(cur) + m * delta[j]
      want_acc[j] <- logu[j] < want_ratio[j]
      if (want_acc[j]) cur <- prop
    }
    expect_equal(got$log_ratio, want_ratio, tolerance = 1e-8, info = paste("seed", seed))
    expect_identical(got$accepted, want_acc, info = paste("seed", seed))
    expect_equal(got[c("B", "sigma_b", "eta")], cur[c("B", "sigma_b", "eta")], tolerance = 1e-10)
    expect_equal(got$R, sigma_b_leftovers(cur), tolerance = 1e-10)
    n_first <- n_first + want_acc[1]
  }
  expect_gt(n_first, 5L)                               # the second proposal did start from a moved state
})

test_that("an accepted sigma_b move scales every measured leftover and nothing else", {
  x <- axis_move_state(3L); x$kappa[] <- 0; x$Omega[] <- 0   # a constant likelihood
  got <- NULL
  for (s in 1:20) { seed_both(s); got <- run_sigma_b_move(x); if (got$accepted) break }
  expect_true(got$accepted)
  r <- got$sigma_b / x$sigma_b
  expect_false(r == 1)
  R0 <- sigma_b_leftovers(x); R1 <- sigma_b_leftovers(x, got$B)
  expect_equal(R1[, 1:2], r * R0[, 1:2], tolerance = 1e-10)
  expect_identical(got$B[3:5, ], x$B[3:5, ])           # embedding rows bit for bit
  expect_equal(got$eta, x$rest + x$X %*% got$B, tolerance = 1e-12)
})

test_that("the sigma_b move does nothing and draws nothing without measured columns", {
  x <- axis_move_state(4L); x$cp$is_emb[] <- TRUE
  seed_both(5L); before <- .Random.seed
  got <- run_sigma_b_move(x, c(0.5, 1))
  expect_identical(.Random.seed, before)
  expect_identical(got[c("B", "sigma_b", "eta")], x[c("B", "sigma_b", "eta")])
  expect_identical(got$accepted, c(FALSE, FALSE))
  expect_identical(got$log_ratio, c(NA_real_, NA_real_))
})

test_that("with sigma_b and its leftovers near zero the ratio is finite and upward moves are accepted", {
  x <- sigma_b_scaled(axis_move_state(6L), 1e-12)     # far below any value a chain was seen at
  up <- 0L; total <- 0L
  for (s in 1:60) {
    seed_both(700L + s); got <- run_sigma_b_move(x, 1)
    seed_both(700L + s); delta <- rnorm(1, 0, 1)
    expect_true(is.finite(got$log_ratio))
    expect_equal(got$log_ratio, delta, tolerance = 1e-6)   # the likelihood and prior terms vanish
    if (delta > 0) { total <- total + 1L; up <- up + got$accepted }
  }
  expect_identical(up, total)                          # every upward proposal accepted
})
