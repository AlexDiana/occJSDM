# The scale-expansion move on the axis scales of the grouped prior
# (axis_scale_move() in R/coef-prior.R): fix attempt 1 for the slow mixing
# recorded in dev/simstudy/grouped-mixing-20261118 (tag
# archive/grouped-mixing).

# axis_move_state() and axis_move_log_target() are in helper-scale-moves.R.

run_axis_move <- function(x, step = 0.5) with(x, occJSDM:::axis_scale_move(
  B, G, C, lambda, X, kappa, Omega, eta, cp, step))

test_that("the axis scale move accepts on the exact ratio of the target density", {
  # Replays the move's random numbers (all proposals, then all uniforms) and
  # decides each column from the written-out target density. The move
  # multiplies the column's leftovers, G and C by r = exp(delta) and adds
  # delta to the log scale, so the ratio is that of the target densities times
  # r to the power of the number of coefficients scaled (the Jacobian).
  n_accept <- 0L; n_reject <- 0L
  for (seed in 1:20) {
    x <- axis_move_state(seed)
    emb <- which(x$cp$is_emb); m <- ncol(x$B) + nrow(x$G) + nrow(x$C)
    seed_both(100L + seed)
    got <- run_axis_move(x)
    seed_both(100L + seed)
    delta <- rnorm(length(emb), 0, 0.5); logu <- log(runif(length(emb)))
    want <- x; log_ratio <- rep(NA_real_, 5L); accepted <- rep(FALSE, 5L)
    for (j in seq_along(emb)) {
      k <- emb[j]; r <- exp(delta[j])
      prop <- want
      prop$lambda[k] <- want$lambda[k] * r
      prop$G[, k] <- want$G[, k] * r; prop$C[, k] <- want$C[, k] * r
      R <- t(want$B) - want$Tr_b %*% want$G - want$A %*% want$C
      R[, k] <- R[, k] * r
      prop$B <- t(want$Tr_b %*% prop$G + want$A %*% prop$C + R)
      log_ratio[k] <- axis_move_log_target(prop) - axis_move_log_target(want) + m * delta[j]
      if (logu[j] < log_ratio[k]) { want <- prop; accepted[k] <- TRUE }
    }
    expect_equal(got$log_ratio, log_ratio, tolerance = 1e-8, info = seed)
    expect_identical(got$accepted, accepted, info = seed)
    expect_equal(got[c("B", "G", "C", "lambda")], want[c("B", "G", "C", "lambda")], info = seed)
    n_accept <- n_accept + sum(accepted); n_reject <- n_reject + sum(!accepted[emb])
  }
  expect_gt(n_accept, 5L); expect_gt(n_reject, 5L)   # both outcomes were exercised
})

test_that("an accepted move scales the leftovers, G and C of its column by one factor", {
  x <- axis_move_state(3L)
  x$kappa[] <- 0; x$Omega[] <- 0                     # a constant likelihood: most moves accepted
  seed_both(11L)
  got <- run_axis_move(x)
  leftovers <- function(B, G, C) t(B) - x$Tr_b %*% G - x$A %*% C
  R0 <- leftovers(x$B, x$G, x$C); R1 <- leftovers(got$B, got$G, got$C)
  r <- got$lambda / x$lambda
  expect_true(any(got$accepted))
  expect_identical(got$accepted, r != 1)
  for (k in 1:5) {
    expect_equal(got$B[k, ], r[k] * x$B[k, ], tolerance = 1e-12, info = k)
    expect_equal(R1[, k], r[k] * R0[, k], tolerance = 1e-10, info = k)
    expect_equal(got$G[, k], r[k] * x$G[, k], tolerance = 1e-12, info = k)
    expect_equal(got$C[, k], r[k] * x$C[, k], tolerance = 1e-12, info = k)
  }
  # Measured columns, and columns whose move was rejected, are returned bit for bit.
  same <- which(!got$accepted)
  expect_true(all(1:2 %in% same))
  expect_identical(got$B[same, ], x$B[same, ])
  expect_identical(got$G[, same], x$G[, same]); expect_identical(got$C[, same], x$C[, same])
  expect_identical(got$lambda[same], x$lambda[same])
})

test_that("the axis scale move does nothing and draws nothing without embedding columns", {
  x <- axis_move_state(4L); x$cp$is_emb[] <- FALSE; x$lambda[] <- 1
  seed_both(5L); before <- .Random.seed
  got <- run_axis_move(x)
  expect_identical(.Random.seed, before)
  expect_identical(got[c("B", "G", "C", "lambda")], x[c("B", "G", "C", "lambda")])
  expect_identical(got$accepted, rep(FALSE, 5L))
})

test_that("the scale updates with the axis scale move leave the half-Cauchy priors unchanged", {
  # The successive-conditional check of test-grouped-prior.R with the move
  # added where update_jSDMcoef() runs it: after the scale draws, before the G
  # and C draw. With a constant likelihood (kappa and Omega zero) every step
  # leaves the prior invariant, so each scale divided by its prior scale must
  # be half-Cauchy(1), with quartiles tan(pi / 8), 1 and tan(3 * pi / 8).
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
  keep <- matrix(NA_real_, N + burn, 5L,
                 dimnames = list(NULL, c("sigma_b", "sigma_emb", "lambda_e1", "lambda_e2", "lambda_m1")))
  accepted <- 0L
  done <- tryCatch({
    for (i in seq_len(N + burn)) {
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      centre <- Tr_b %*% G + A %*% C
      B <- t(centre + matrix(rnorm(S * p), S, p) * rep(sds$b, each = S))
      sc <- occJSDM:::sample_coef_scales(t(B) - centre, G, C, cp, sigma_b, sigma_emb, lambda)
      sigma_b <- sc$sigma_b; sigma_emb <- sc$sigma_emb; lambda <- sc$lambda
      mv <- occJSDM:::axis_scale_move(B, G, C, lambda, X, zero, zero, zero, cp, 0.5)
      B <- mv$B; G <- mv$G; C <- mv$C; lambda <- mv$lambda
      accepted <- accepted + sum(mv$accepted)
      sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, g_b, p)
      gc <- occJSDM:::sample_GC(B, Tr_b, A, sds$b, sds$g, sds$c); G <- gc$G; C <- gc$C
      A <- occJSDM:::sample_A(B, C, Tr_b, G, sds$b)$A
      keep[i, ] <- c(sigma_b / cp$sigma_b_scale, sigma_emb / cp$sigma_emb_scale,
                     lambda[2:3] / cp$axis_scale, lambda[1])
    }
    TRUE
  }, error = function(e) conditionMessage(e))
  expect_identical(done, TRUE, label = "the chain ran to the end (it does not when the scales run away)")
  if (!isTRUE(done)) return(invisible(NULL))
  keep <- keep[-seq_len(burn), ]
  expect_true(all(keep[, "lambda_m1"] == 1))
  expect_true(all(is.finite(keep) & keep > 0))
  expect_gt(accepted / (2 * (N + burn)), 0.5)          # with no data most moves are accepted
  quartiles <- tan(pi / 2 * c(.25, .5, .75))
  share <- vapply(quartiles, function(x) colMeans(keep[, 1:4] < x), numeric(4))
  # Tolerances: see the comment at the end of this file for how they were
  # measured and what a wrong acceptance ratio does to the shares.
  tolerance <- c(sigma_b = 0.045, sigma_emb = 0.06, lambda_e1 = 0.04, lambda_e2 = 0.04)
  for (nm in rownames(share)) for (j in 1:3)
    expect_lt(abs(share[nm, j] - c(.25, .5, .75)[j]), tolerance[[nm]],
              label = sprintf("error of the share of %s below its prior quartile %d", nm, j))
})

axis_fit <- function(sim, embedding = TRUE, MCMCparams = FIXTURE_MCMC) {
  dat <- sim$data_list
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = c(list(n_factors = 1L), if (embedding) list(embeddingVars = dat$embeddingVars)),
    occCovariates = grouped_occ_covariates(dat), MCMCparams = MCMCparams)))
}

test_that("a grouped binary fit runs the axis scale move and records its acceptance", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- axis_fit(sim)
  jo <- fit$results_output$jsdm_output
  am <- fit$infos$axis_move
  expect_identical(dim(am$step), c(8L, 2L))            # embedding columns by chains
  expect_identical(rownames(am$step), dimnames(jo$lambda_output)[[1]])
  expect_true(all(am$step >= 0.01 & am$step <= 5))
  expect_false(all(am$step == 0.5))                    # adapted during burn-in
  expect_identical(dim(am$acceptance), c(8L, 2L))
  expect_identical(rownames(am$acceptance), dimnames(jo$lambda_output)[[1]])
  expect_identical(dim(am$acceptance_burnin), c(8L, 2L))
  expect_true(all(am$acceptance >= 0 & am$acceptance <= 1))
  expect_gt(mean(am$acceptance), 0.1)
  expect_identical(dim(jo$lambda_output), c(8L, 20L, 2L))
  for (nm in c("lambda_output", "sigmaemb_output", "sigmab_output", "B_output", "G_output", "C_output"))
    expect_true(all(is.finite(jo[[nm]])), info = nm)
  expect_true(all(jo$lambda_output > 0))
})

test_that("the axis scale move is off outside grouped binary fits", {
  # First version: the binary model (one row per site) only.
  occ <- axis_fit(simulate_grouped_fixture(model = "occupancy", p_emb = 6L, S = 8L))
  expect_identical(occ$infos$coef_prior$type, "grouped")
  expect_false("axis_move" %in% names(occ$infos))
  fixed <- axis_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L), embedding = FALSE)
  expect_identical(fixed$infos$coef_prior, list(type = "fixed"))
  expect_false("axis_move" %in% names(fixed$infos))
})

test_that("the coefficient update hands the move the current state and keeps its result", {
  # A spy on axis_scale_move() inside update_jSDMcoef(), in a real fit. What
  # the move returns for B and the axis scales must be what the sweep returns
  # (neither is redrawn after the move), and the linear predictor it is given
  # must be the one implied by the coefficients it is given: the intercepts
  # and loadings just drawn, and the factor scores of the previous sweep.
  original_move <- occJSDM:::axis_scale_move
  original_update <- occJSDM:::update_jSDMcoef
  seen <- list(); sweeps <- list()
  update <- original_update
  update_env <- new.env(parent = environment(update))
  update_env$axis_scale_move <- function(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step) {
    out <- original_move(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step)
    seen[[length(seen) + 1L]] <<- list(B = B, lambda = lambda, X = X, kappa = kappa, eta = eta,
                                       step = step, out = out)
    out
  }
  environment(update) <- update_env
  local_mocked_bindings(update_jSDMcoef = function(list_data, list_params, ...) {
    out <- update(list_data, list_params, ...)
    sweeps[[length(sweeps) + 1L]] <<- list(z = list_data$z, U_before = list_params$U, out = out)
    out
  }, .package = "occJSDM")
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  axis_fit(sim, MCMCparams = list(nchain = 1L, nburn = 3L, niter = 3L, nthin = 1L))
  expect_identical(length(seen), 6L)
  expect_identical(length(sweeps), 6L)
  for (i in seq_along(seen)) {
    mv <- seen[[i]]; sw <- sweeps[[i]]
    expect_identical(sw$out$B, mv$out$B, info = paste("sweep", i))
    expect_identical(sw$out$lambda, mv$out$lambda, info = paste("sweep", i))
    expect_identical(mv$kappa, sw$z - 0.5, info = paste("sweep", i))
    eta <- matrix(sw$out$B0, nrow(mv$X), ncol(mv$B), byrow = TRUE) + mv$X %*% mv$B +
      sw$U_before %*% sw$out$L
    expect_equal(unname(mv$eta), unname(eta), tolerance = 1e-10, info = paste("sweep", i))
  }
  expect_true(any(vapply(seen, function(mv) any(mv$out$accepted), logical(1))))
  # The step of each embedding column starts at 0.5, is adapted after each of
  # the three burn-in sweeps from that sweep's log ratios, and is then frozen.
  expect_identical(seen[[1]]$step, rep(0.5, 8L))
  for (i in 1:3) {
    want <- occJSDM:::axis_step_adapt(seen[[i]]$step, seen[[i]]$out$log_ratio[3:10], i)
    expect_identical(seen[[i + 1L]]$step, want, info = paste("sweep", i))
  }
  expect_false(identical(seen[[4]]$step, rep(0.5, 8L)))
  expect_identical(seen[[5]]$step, seen[[4]]$step)
  expect_identical(seen[[6]]$step, seen[[4]]$step)
})

test_that("the burn-in adaptation moves each step toward the target acceptance rate", {
  # Robbins-Monro on the log step, as for the intercept shift (R/b0_shift.R):
  # log(step) moves by (acceptance probability - 0.44) / sqrt(t).
  step <- c(0.5, 0.5, 0.5, 2)
  log_ratio <- c(1.3, log(0.44), -3, -Inf)             # accepted for sure, on target, unlikely, never
  got <- occJSDM:::axis_step_adapt(step, log_ratio, 4L)
  expect_equal(got, step * exp((c(1, 0.44, exp(-3), 0) - 0.44) / 2))
  expect_gt(got[1], 0.5); expect_equal(got[2], 0.5); expect_lt(got[3], 0.5); expect_lt(got[4], 2)
  # Kept within loose bounds, and a column whose ratio is not finite counts as rejected.
  expect_identical(occJSDM:::axis_step_adapt(c(5, 0.01), c(10, -10), 1L), c(5, 0.01))
  expect_equal(occJSDM:::axis_step_adapt(1, NaN, 1L), exp(-0.44))
})

test_that("the axis scale move takes one step per embedding column", {
  x <- axis_move_state(2L)
  x$kappa[] <- 0; x$Omega[] <- 0
  seed_both(21L); got <- run_axis_move(x, step = c(0.1, 1, 3))
  seed_both(21L); delta <- rnorm(3, 0, c(0.1, 1, 3)); logu <- log(runif(3))
  # With a constant likelihood only the half-Cauchy ratio is left.
  log_ratio <- delta - log1p((x$lambda[3:5] * exp(delta) / 0.3)^2) + log1p((x$lambda[3:5] / 0.3)^2)
  expect_equal(got$log_ratio[3:5], log_ratio)
  expect_equal(got$lambda[3:5], x$lambda[3:5] * ifelse(logu < log_ratio, exp(delta), 1))
})

# Tolerances of the prior-invariance test, measured on 18 November 2026 with
# the test's own loop over seeds 1 to 8 (50,000 kept sweeps each, about 1,100
# to 1,500 effective draws of a log axis scale). The standard deviation over
# seeds of a share's error was at most 0.004 for sigma_b, 0.010 for sigma_emb
# and 0.010 for an axis scale, with no bias (largest mean error 0.006), and
# the largest error of any share in any seed was 0.008, 0.018 and 0.016. The
# axis tolerance of 0.04 is four of those standard deviations; the other two
# are those of test-grouped-prior.R. A wrong acceptance ratio fails by far
# more: leaving out the factor lambda' / lambda (the Jacobian of the log
# scale) moves the axis shares by 0.47 to 0.49, and a half-Cauchy scale wrong
# by a factor of 2 moves them by 0.30 to 0.34. About 88 percent of moves are
# accepted in this run at a step of 0.5.

# The axis scale move under a cap (slab_scale). Design:
# dev/superpowers/specs/2026-11-25-axis-slab-design.md, "Mathematics", "The
# axis scale move". Under a cap the move needs the leftovers of the B it is
# given and the current sigma_emb.

run_capped_axis_move <- function(x, step = 0.5) with(x, occJSDM:::axis_scale_move(
  B, G, C, lambda, X, kappa, Omega, eta, cp, step,
  R = t(B) - Tr_b %*% G - A %*% C, sigma_emb = sigma_emb))

test_that("under a cap the axis scale move accepts on the exact ratio of the target density", {
  # As the uncapped test above, with a cap of 0.8. Under the cap the move
  # multiplies column k of G and C by rho_g, the ratio of the effective axis
  # scales at the proposed and current lambda_k, and the column's leftovers by
  # rho_e, the same ratio for the leftover's capped standard deviation. The
  # Jacobian of that scaling is (g_b + q) log(rho_g) + S log(rho_e).
  n_accept <- 0L; n_reject <- 0L
  for (seed in 1:20) {
    x <- axis_move_state(seed); x$cp$slab_scale <- 0.8
    emb <- which(x$cp$is_emb)
    seed_both(100L + seed)
    got <- run_capped_axis_move(x)
    seed_both(100L + seed)
    delta <- rnorm(length(emb), 0, 0.5); logu <- log(runif(length(emb)))
    want <- x; log_ratio <- rep(NA_real_, 5L); accepted <- rep(FALSE, 5L)
    for (j in seq_along(emb)) {
      k <- emb[j]
      prop <- want
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
    expect_equal(got$log_ratio, log_ratio, tolerance = 1e-8, info = paste("seed", seed))
    expect_identical(got$accepted, accepted, info = paste("seed", seed))
    expect_equal(got[c("B", "G", "C", "lambda")], want[c("B", "G", "C", "lambda")],
                 info = paste("seed", seed))
    n_accept <- n_accept + sum(accepted); n_reject <- n_reject + sum(!accepted[emb])
  }
  expect_gt(n_accept, 5L); expect_gt(n_reject, 5L)   # both outcomes were exercised
})

test_that("under a cap an accepted axis move changes its column's parts by the two factors and nothing else", {
  x <- axis_move_state(3L); x$cp$slab_scale <- 0.8
  x$kappa[] <- 0; x$Omega[] <- 0                     # a constant likelihood: most moves accepted
  seed_both(11L)
  got <- run_capped_axis_move(x)
  leftovers <- function(B, G, C) t(B) - x$Tr_b %*% G - x$A %*% C
  R0 <- leftovers(x$B, x$G, x$C); R1 <- leftovers(got$B, got$G, got$C)
  eff0 <- occJSDM:::coef_effective_scales(x$cp, x$sigma_emb, x$lambda)
  eff1 <- occJSDM:::coef_effective_scales(x$cp, x$sigma_emb, got$lambda)
  rho_g <- eff1$axis / eff0$axis; rho_e <- eff1$leftover / eff0$leftover
  expect_true(any(got$accepted))
  expect_identical(got$accepted, got$lambda != x$lambda)
  # Under the cap the two factors differ, so the row of B is not scaled as a whole.
  expect_true(any(abs(rho_g - rho_e)[got$accepted] > 1e-3))
  for (k in which(got$accepted)) {
    expect_equal(got$G[, k], rho_g[k] * x$G[, k], tolerance = 1e-12, info = format(k))
    expect_equal(got$C[, k], rho_g[k] * x$C[, k], tolerance = 1e-12, info = format(k))
    expect_equal(R1[, k], rho_e[k] * R0[, k], tolerance = 1e-10, info = format(k))
    expect_equal(got$B[k, ], rho_g[k] * (x$B[k, ] - R0[, k]) + rho_e[k] * R0[, k],
                 tolerance = 1e-12, info = format(k))
  }
  # Measured columns, and columns whose move was rejected, are returned bit for bit.
  same <- which(!got$accepted)
  expect_true(all(1:2 %in% same))
  expect_identical(got$B[same, ], x$B[same, ])
  expect_identical(got$G[, same], x$G[, same]); expect_identical(got$C[, same], x$C[, same])
  expect_identical(got$lambda[same], x$lambda[same])
})

test_that("under a cap the coefficient update hands each move the current state", {
  # The spy test above with a cap. The sigma_emb move must be given the
  # current axis scales; the axis scale move must be given the B and sigma_emb
  # the sigma_emb move returns, and the leftovers of that B, not those from
  # before the sigma_emb move. What the axis scale move returns must be what
  # the sweep returns.
  original_sigma <- occJSDM:::sigma_emb_move
  original_move <- occJSDM:::axis_scale_move
  original_update <- occJSDM:::update_jSDMcoef
  seen_sigma <- list(); seen <- list(); sweeps <- list()
  update <- original_update
  update_env <- new.env(parent = environment(update))
  update_env$sigma_emb_move <- function(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step,
                                        lambda = NULL) {
    out <- original_sigma(B, R, sigma_emb, X, kappa, Omega, eta, coef_prior, step, lambda = lambda)
    seen_sigma[[length(seen_sigma) + 1L]] <<- list(R = R, lambda = lambda, out = out)
    out
  }
  update_env$axis_scale_move <- function(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step,
                                         R = NULL, sigma_emb = NULL) {
    out <- original_move(B, G, C, lambda, X, kappa, Omega, eta, coef_prior, step,
                         R = R, sigma_emb = sigma_emb)
    seen[[length(seen) + 1L]] <<- list(B = B, G = G, C = C, lambda = lambda, eta = eta, R = R,
                                       sigma_emb = sigma_emb, coef_prior = coef_prior, out = out)
    out
  }
  environment(update) <- update_env
  local_mocked_bindings(update_jSDMcoef = function(list_data, list_params, ...) {
    out <- update(list_data, list_params, ...)
    sweeps[[length(sweeps) + 1L]] <<- list(Tr_b = list_data$Tr_b, before = list_params, out = out)
    out
  }, .package = "occJSDM")
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  dat <- sim$data_list
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
    listPriors = list(slab_scale = 0.5),
    occCovariates = grouped_occ_covariates(dat),
    MCMCparams = list(nchain = 1L, nburn = 3L, niter = 3L, nthin = 1L))))
  expect_identical(length(seen_sigma), 6L)
  expect_identical(length(seen), 6L)
  expect_identical(length(sweeps), 6L)
  for (i in seq_along(seen)) {
    mv <- seen[[i]]; sg <- seen_sigma[[i]]; sw <- sweeps[[i]]; label <- paste("sweep", i)
    expect_identical(mv$coef_prior$slab_scale, 0.5, info = label)
    expect_identical(sg$lambda, mv$lambda, info = label)
    expect_identical(mv$B, sg$out$B, info = label)
    expect_identical(mv$eta, sg$out$eta, info = label)
    expect_identical(mv$sigma_emb, sg$out$sigma_emb, info = label)
    leftovers <- t(mv$B) - sw$Tr_b %*% mv$G - sw$before$A %*% mv$C
    expect_equal(unname(mv$R), unname(leftovers), tolerance = 1e-10, info = label)
    expect_identical(sw$out$B, mv$out$B, info = label)
    expect_identical(sw$out$lambda, mv$out$lambda, info = label)
    expect_identical(sw$out$sigma_emb, sg$out$sigma_emb, info = label)
  }
  # The check of the leftovers means something only if a sigma_emb move was
  # accepted, so that they differ from the ones it was given.
  moved <- vapply(seen_sigma, function(s) isTRUE(s$out$accepted), logical(1))
  expect_true(any(moved))
  i <- which(moved)[1]
  emb <- seen[[i]]$coef_prior$is_emb
  expect_gt(max(abs(seen[[i]]$R[, emb] - seen_sigma[[i]]$R[, emb])), 1e-6)
})
