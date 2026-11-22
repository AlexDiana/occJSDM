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

test_that("the scale updates with all three scale moves leave the half-Cauchy priors unchanged", {
  # The loop of the test of the same kind in test-sigma-emb-move.R, with the
  # sigma_b move added where update_jSDMcoef() runs it: after the scale draws
  # and before the sigma_emb move, with both of its proposals. With a constant
  # likelihood every step leaves the prior invariant, so each scale divided by
  # its prior scale must be half-Cauchy(1). Same seed, sizes and tolerances.
  S <- 3L; p <- 3L; gt <- 2L; g_b <- 2L; n <- 4L; N <- 50000L; burn <- 1000L
  seed_both(1L)
  Tr_b <- cbind(1, rnorm(S)); X <- matrix(rnorm(n * p), n, p)
  zero <- matrix(0, n, S)
  cp <- occJSDM:::read_coef_prior(
    list(sigma_b_scale = 0.8, sigma_emb_scale = 1.3, axis_scale = 0.5),
    list(embeddingVars = c("e1", "e2")), S, gt,
    data.frame(column = c("m1", "e1", "e2"), variable = c("m1", "e1", "e2"), kind = "numeric"))
  start <- occJSDM:::coef_scales_start(cp, p)
  sigma_b <- cp$sigma_b_scale; sigma_emb <- start$sigma_emb; lambda <- start$lambda
  G <- matrix(rnorm(g_b * p), g_b, p); C <- matrix(rnorm(gt * p), gt, p)
  A <- matrix(rnorm(S * gt), S, gt)
  keep <- matrix(NA_real_, N + burn, 4L,
                 dimnames = list(NULL, c("sigma_b", "sigma_emb", "lambda_e1", "lambda_e2")))
  accepted <- 0L; accepted_b <- c(0L, 0L)
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
      sg <- occJSDM:::sigma_emb_move(B, t(B) - centre, sigma_emb, X, zero, zero, zero, cp, 0.5)
      B <- sg$B; sigma_emb <- sg$sigma_emb
      accepted <- accepted + sg$accepted
      mv <- occJSDM:::axis_scale_move(B, G, C, lambda, X, zero, zero, sg$eta, cp, 0.5)
      B <- mv$B; G <- mv$G; C <- mv$C; lambda <- mv$lambda
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
  quartiles <- tan(pi / 2 * c(.25, .5, .75))
  share <- vapply(quartiles, function(x) colMeans(keep < x), numeric(4))
  # Tolerances: those of test-sigma-emb-move.R; see the comment at its end.
  tolerance <- c(sigma_b = 0.045, sigma_emb = 0.04, lambda_e1 = 0.05, lambda_e2 = 0.05)
  for (nm in rownames(share)) for (j in 1:3)
    expect_lt(abs(share[nm, j] - c(.25, .5, .75)[j]), tolerance[[nm]],
              label = sprintf("error of the share of %s below its prior quartile %d", nm, j))
})

test_that("a grouped binary fit runs the sigma_b move and records its acceptance", {
  for (mcmc in list(FIXTURE_MCMC, c(FIXTURE_MCMC, list(ncores = 2L)))) {
    fit <- sigma_move_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L), mcmc)
    am <- fit$infos$axis_move
    expect_length(am$sigma_b_step, 2L)                 # one per chain
    expect_true(all(am$sigma_b_step >= 0.01 & am$sigma_b_step <= 5))
    expect_false(all(am$sigma_b_step == 0.5))          # adapted during burn-in
    expect_identical(am$sigma_b_large_step, 1)
    for (nm in c("sigma_b_acceptance", "sigma_b_acceptance_burnin")) {
      expect_identical(dim(am[[nm]]), c(2L, 2L))
      expect_identical(rownames(am[[nm]]), c("adapted", "large"))
      expect_true(all(am[[nm]] >= 0 & am[[nm]] <= 1))
    }
    expect_true(all(is.finite(fit$results_output$jsdm_output$sigmab_output)))
  }
})

test_that("the coefficient update runs the sigma_b move first and hands its state on", {
  # Spies on the sigma_b and sigma_emb moves inside update_jSDMcoef(), in a
  # real fit, as the last test of test-sigma-emb-move.R does. The sigma_b move
  # must be given the leftovers and linear predictor of the B just drawn; the
  # sigma_emb move must be given the B, leftovers and linear predictor the
  # sigma_b move returns; and the sweep's sigma_b must be the move's.
  original_b <- occJSDM:::sigma_b_move
  original_emb <- occJSDM:::sigma_emb_move
  original_update <- occJSDM:::update_jSDMcoef
  seen_b <- list(); seen_emb <- list(); sweeps <- list()
  update <- original_update
  update_env <- new.env(parent = environment(update))
  update_env$sigma_b_move <- function(B, R, sigma_b, X, kappa, Omega, eta, coef_prior, step) {
    out <- original_b(B, R, sigma_b, X, kappa, Omega, eta, coef_prior, step)
    seen_b[[length(seen_b) + 1L]] <<- list(B = B, R = R, X = X, eta = eta, step = step, out = out)
    out
  }
  update_env$sigma_emb_move <- function(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step) {
    seen_emb[[length(seen_emb) + 1L]] <<- list(B = B, R = R, eta = eta)
    original_emb(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step)
  }
  environment(update) <- update_env
  local_mocked_bindings(update_jSDMcoef = function(list_data, list_params, ...) {
    out <- update(list_data, list_params, ...)
    sweeps[[length(sweeps) + 1L]] <<- list(Tr_b = list_data$Tr_b, before = list_params, out = out)
    out
  }, .package = "occJSDM")
  sigma_move_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L),
                 MCMCparams = list(nchain = 1L, nburn = 3L, niter = 3L, nthin = 1L))
  expect_identical(length(seen_b), 6L)
  expect_identical(length(seen_emb), 6L)
  for (i in seq_along(seen_b)) {
    sb <- seen_b[[i]]; sw <- sweeps[[i]]; label <- paste("sweep", i)
    leftovers <- t(sb$B) - sw$Tr_b %*% sw$before$G - sw$before$A %*% sw$before$C
    expect_equal(unname(sb$R), unname(leftovers), tolerance = 1e-10, info = label)
    eta <- matrix(sw$out$B0, nrow(sb$X), ncol(sb$B), byrow = TRUE) + sb$X %*% sb$B +
      sw$before$U %*% sw$out$L
    expect_equal(unname(sb$eta), unname(eta), tolerance = 1e-10, info = label)
    expect_identical(sb$step[2], occJSDM:::SIGMA_B_MOVE_LARGE_STEP, info = label)
    expect_identical(seen_emb[[i]]$B, sb$out$B, info = label)
    expect_identical(seen_emb[[i]]$R, sb$out$R, info = label)
    expect_identical(seen_emb[[i]]$eta, sb$out$eta, info = label)
    expect_identical(sw$out$sigma_b, sb$out$sigma_b, info = label)
  }
  # The adapted step starts at 0.5, is adapted after each of the three burn-in
  # sweeps from that sweep's first log ratio, and is then frozen. The large
  # step is never adapted.
  expect_identical(seen_b[[1]]$step, c(0.5, 1))
  for (i in 1:3) {
    want <- occJSDM:::axis_step_adapt(seen_b[[i]]$step[1], seen_b[[i]]$out$log_ratio[1], i)
    expect_identical(seen_b[[i + 1L]]$step, c(want, 1), info = paste("sweep", i))
  }
  expect_false(identical(seen_b[[4]]$step[1], 0.5))
  expect_identical(seen_b[[5]]$step, seen_b[[4]]$step)
  expect_identical(seen_b[[6]]$step, seen_b[[4]]$step)
})
