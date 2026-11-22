# The scale-expansion move on sigma_emb, the embedding leftover scale of the
# grouped prior (sigma_emb_move() in R/coef-prior.R): fix attempt 2 for the
# slow mixing recorded in dev/simstudy/grouped-mixing-20261118 (tag
# archive/grouped-mixing).
# axis_move_state() and axis_move_log_target() are in helper-scale-moves.R.

sigma_move_leftovers <- function(x, B = x$B) t(B) - x$Tr_b %*% x$G - x$A %*% x$C

run_sigma_move <- function(x, step = 0.5) with(x, occJSDM:::sigma_emb_move(
  B, sigma_move_leftovers(x), sigma_emb, X, kappa, Omega, eta, cp, step))

test_that("the sigma_emb move accepts on the exact ratio of the target density", {
  # Replays the move's two random numbers and decides from the written-out
  # target density. The move multiplies every leftover of the embedding block
  # by r = exp(delta) and adds delta to log sigma_emb, so the ratio is that of
  # the target densities times r to the power of the number of leftovers
  # scaled (species times embedding columns, the Jacobian).
  n_accept <- 0L; n_reject <- 0L
  for (seed in 1:40) {
    x <- axis_move_state(seed)
    emb <- which(x$cp$is_emb); m <- ncol(x$B) * length(emb)
    step <- c(0.5, 0.05)[1L + seed %% 2L]        # large steps are mostly rejected here, small ones accepted
    seed_both(300L + seed)
    got <- run_sigma_move(x, step)
    seed_both(300L + seed)
    delta <- rnorm(1, 0, step); logu <- log(runif(1)); r <- exp(delta)
    prop <- x
    prop$sigma_emb <- x$sigma_emb * r
    R <- sigma_move_leftovers(x); R[, emb] <- R[, emb] * r
    prop$B <- t(x$Tr_b %*% x$G + x$A %*% x$C + R)
    log_ratio <- axis_move_log_target(prop) - axis_move_log_target(x) + m * delta
    accepted <- logu < log_ratio
    expect_equal(got$log_ratio, log_ratio, tolerance = 1e-8, info = paste("seed", seed))
    expect_identical(got$accepted, accepted, info = paste("seed", seed))
    if (accepted) {
      expect_equal(got$B, prop$B, info = paste("seed", seed))
      expect_equal(got$sigma_emb, prop$sigma_emb, info = paste("seed", seed))
      expect_equal(got$eta, x$rest + x$X %*% prop$B, info = paste("seed", seed))
      n_accept <- n_accept + 1L
    } else {
      # A rejected move returns the state bit for bit.
      expect_identical(got[c("B", "sigma_emb", "eta")], x[c("B", "sigma_emb", "eta")],
                       info = paste("seed", seed))
      n_reject <- n_reject + 1L
    }
  }
  expect_gt(n_accept, 5L); expect_gt(n_reject, 5L)   # both outcomes were exercised
})

test_that("an accepted sigma_emb move scales every embedding leftover by one factor", {
  x <- axis_move_state(3L)
  x$kappa[] <- 0; x$Omega[] <- 0                     # a constant likelihood
  got <- NULL
  for (s in 1:20) {                                  # the first seed whose move is accepted
    seed_both(s); got <- run_sigma_move(x)
    if (got$accepted) break
  }
  expect_true(got$accepted)
  r <- got$sigma_emb / x$sigma_emb
  expect_false(r == 1)
  R0 <- sigma_move_leftovers(x); R1 <- sigma_move_leftovers(x, got$B)
  expect_equal(R1[, 3:5], r * R0[, 3:5], tolerance = 1e-10)
  # Measured rows of B come back bit for bit, so their leftovers are unchanged.
  expect_identical(got$B[1:2, ], x$B[1:2, ])
  # The linear predictor follows B.
  expect_equal(got$eta, x$rest + x$X %*% got$B, tolerance = 1e-12)
})

test_that("the sigma_emb move does nothing and draws nothing without embedding columns", {
  x <- axis_move_state(4L); x$cp$is_emb[] <- FALSE; x$lambda[] <- 1; x$sigma_emb <- NA_real_
  seed_both(5L); before <- .Random.seed
  got <- run_sigma_move(x)
  expect_identical(.Random.seed, before)
  expect_identical(got[c("B", "sigma_emb", "eta")], x[c("B", "sigma_emb", "eta")])
  expect_identical(got$accepted, FALSE)
  expect_identical(got$log_ratio, NA_real_)
})

test_that("the scale updates with both scale moves leave the half-Cauchy priors unchanged", {
  # The successive-conditional check of test-axis-scale-move.R with the
  # sigma_emb move added where update_jSDMcoef() runs it: after the scale
  # draws and before the axis scale move. With a constant likelihood every
  # step leaves the prior invariant, so each scale divided by its prior scale
  # must be half-Cauchy(1).
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
  accepted <- 0L
  done <- tryCatch({
    for (i in seq_len(N + burn)) {
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      centre <- Tr_b %*% G + A %*% C
      B <- t(centre + matrix(rnorm(S * p), S, p) * rep(sds$b, each = S))
      sc <- occJSDM:::sample_coef_scales(t(B) - centre, G, C, cp, sigma_b, sigma_emb, lambda)
      sigma_b <- sc$sigma_b; sigma_emb <- sc$sigma_emb; lambda <- sc$lambda
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
  quartiles <- tan(pi / 2 * c(.25, .5, .75))
  share <- vapply(quartiles, function(x) colMeans(keep < x), numeric(4))
  # Tolerances: see the comment at the end of this file.
  tolerance <- c(sigma_b = 0.045, sigma_emb = 0.04, lambda_e1 = 0.05, lambda_e2 = 0.05)
  for (nm in rownames(share)) for (j in 1:3)
    expect_lt(abs(share[nm, j] - c(.25, .5, .75)[j]), tolerance[[nm]],
              label = sprintf("error of the share of %s below its prior quartile %d", nm, j))
})

test_that("a grouped binary fit runs the sigma_emb move and records its acceptance", {
  fit <- sigma_move_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L))
  am <- fit$infos$axis_move
  expect_length(am$sigma_emb_step, 2L)                 # one per chain
  expect_true(all(am$sigma_emb_step >= 0.01 & am$sigma_emb_step <= 5))
  expect_false(all(am$sigma_emb_step == 0.5))          # adapted during burn-in
  expect_length(am$sigma_emb_acceptance, 2L)
  expect_length(am$sigma_emb_acceptance_burnin, 2L)
  expect_true(all(am$sigma_emb_acceptance >= 0 & am$sigma_emb_acceptance <= 1))
  jo <- fit$results_output$jsdm_output
  expect_true(all(is.finite(jo$sigmaemb_output) & jo$sigmaemb_output > 0))
  expect_true(all(is.finite(jo$B_output)))
})

test_that("the coefficient update runs the sigma_emb move on the current state, before the axis move", {
  # Spies on both moves inside update_jSDMcoef(), in a real fit. The sigma_emb
  # the sweep returns must be the move's (it is not redrawn afterwards); the
  # leftovers and linear predictor the move is given must be those implied by
  # the coefficients just drawn; and the axis scale move must be handed the B
  # and linear predictor the sigma_emb move returns.
  original_sigma <- occJSDM:::sigma_emb_move
  original_axis <- occJSDM:::axis_scale_move
  original_update <- occJSDM:::update_jSDMcoef
  seen_sigma <- list(); seen_axis <- list(); sweeps <- list()
  update <- original_update
  update_env <- new.env(parent = environment(update))
  update_env$sigma_emb_move <- function(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step) {
    out <- original_sigma(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step)
    seen_sigma[[length(seen_sigma) + 1L]] <<- list(B = B, R = R, X = X, kappa = kappa, eta = eta,
                                                   step = step, out = out)
    out
  }
  update_env$axis_scale_move <- function(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step) {
    seen_axis[[length(seen_axis) + 1L]] <<- list(B = B, eta = eta)
    original_axis(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step)
  }
  environment(update) <- update_env
  local_mocked_bindings(update_jSDMcoef = function(list_data, list_params, ...) {
    out <- update(list_data, list_params, ...)
    sweeps[[length(sweeps) + 1L]] <<- list(z = list_data$z, Tr_b = list_data$Tr_b, before = list_params, out = out)
    out
  }, .package = "occJSDM")
  sigma_move_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L),
                 MCMCparams = list(nchain = 1L, nburn = 3L, niter = 3L, nthin = 1L))
  expect_identical(length(seen_sigma), 6L)
  expect_identical(length(seen_axis), 6L)
  for (i in seq_along(seen_sigma)) {
    sg <- seen_sigma[[i]]; sw <- sweeps[[i]]; label <- paste("sweep", i)
    expect_identical(sw$out$sigma_emb, sg$out$sigma_emb, info = label)
    expect_identical(sg$kappa, sw$z - 0.5, info = label)
    eta <- matrix(sw$out$B0, nrow(sg$X), ncol(sg$B), byrow = TRUE) + sg$X %*% sg$B +
      sw$before$U %*% sw$out$L
    expect_equal(unname(sg$eta), unname(eta), tolerance = 1e-10, info = label)
    leftovers <- t(sg$B) - sw$Tr_b %*% sw$before$G - sw$before$A %*% sw$before$C
    expect_equal(unname(sg$R), unname(leftovers), tolerance = 1e-10, info = label)
    expect_identical(seen_axis[[i]]$B, sg$out$B, info = label)
    expect_identical(seen_axis[[i]]$eta, sg$out$eta, info = label)
  }
  # The step starts at 0.5, is adapted after each of the three burn-in sweeps
  # from that sweep's log ratio, and is then frozen.
  expect_identical(seen_sigma[[1]]$step, 0.5)
  for (i in 1:3) {
    want <- occJSDM:::axis_step_adapt(seen_sigma[[i]]$step, seen_sigma[[i]]$out$log_ratio, i)
    expect_identical(seen_sigma[[i + 1L]]$step, want, info = paste("sweep", i))
  }
  expect_false(identical(seen_sigma[[4]]$step, 0.5))
  expect_identical(seen_sigma[[5]]$step, seen_sigma[[4]]$step)
  expect_identical(seen_sigma[[6]]$step, seen_sigma[[4]]$step)
})

# Tolerances of the prior-invariance test, measured on 18 November 2026 with
# the test's own loop over seeds 1 to 8 (50,000 kept sweeps each, about 1,200
# to 1,800 effective draws of each log scale). The standard deviation over
# seeds of a share's error was at most 0.009 for sigma_b, 0.008 for sigma_emb
# and 0.012 for an axis scale, with no bias (largest mean error 0.007), and
# the largest error of any share in any seed was 0.019, 0.013 and 0.030. The
# tolerances of sigma_emb and of an axis scale are four or more of those
# standard deviations; sigma_b's is that of test-grouped-prior.R. A wrong
# acceptance ratio in the sigma_emb move fails by far more: leaving out the
# factor sigma' / sigma (the Jacobian of the log scale) moves the sigma_emb
# shares by 0.39 or sends the scales to where the regression draws fail, and
# a half-Cauchy scale wrong by a factor of 2 moves them by 0.20. About 88
# percent of sigma_emb moves are accepted in this run at a step of 0.5, and
# the effective draws of log sigma_emb about double (from about 800 to about
# 1,700 of 50,000).
