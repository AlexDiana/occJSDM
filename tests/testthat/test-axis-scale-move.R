# The scale-expansion move on the axis scales of the grouped prior
# (axis_scale_move() in R/coef-prior.R). A prototype: fix attempt 1 for the
# slow mixing recorded in dev/simstudy/grouped-mixing-20261118.

# A state of the coefficient block drawn from its prior at the given scales,
# with an arbitrary likelihood state: `rest` stands for every term of the
# linear predictor other than X B, kappa and Omega for the response minus one
# half and the Polya-Gamma weights.
axis_move_state <- function(seed = 1L, n = 25L, S = 7L) {
  set.seed(seed); gb <- 2L; gt <- 2L; p <- 5L
  cp <- list(type = "grouped", is_emb = c(FALSE, FALSE, TRUE, TRUE, TRUE),
             sigma_b_scale = 1, sigma_emb_scale = 1, axis_scale = 0.3,
             c_base = c(1, 1, rep(1 / sqrt(gt), 3)))
  lambda <- c(1, 1, 0.5, 0.05, 0.9); sigma_b <- 0.7; sigma_emb <- 0.4
  sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, gb, p)
  X <- matrix(rnorm(n * p), n, p)
  Tr_b <- cbind(1, rnorm(S)); A <- matrix(rnorm(S * gt), S, gt)
  G <- matrix(rnorm(gb * p), gb, p) * sds$g
  C <- matrix(rnorm(gt * p), gt, p) * rep(sds$c, each = gt)
  R <- matrix(rnorm(S * p), S, p) * rep(sds$b, each = S)
  B <- t(Tr_b %*% G + A %*% C + R)
  rest <- matrix(rnorm(n * S, sd = 0.5), n, S)
  list(cp = cp, lambda = lambda, sigma_b = sigma_b, sigma_emb = sigma_emb, X = X, Tr_b = Tr_b,
       A = A, G = G, C = C, B = B, rest = rest, eta = rest + X %*% B,
       kappa = matrix(rbinom(n * S, 1, 0.4) - 0.5, n, S), Omega = matrix(rexp(n * S) / 4, n, S))
}

# The log of the density the sampler targets, written out term by term: the
# Normal priors of the leftovers, G and C at the standard deviations
# coef_prior_sds() gives, the half-Cauchy prior of each axis scale expressed
# as a density of the log scale (hence the added log lambda), and the
# likelihood given the Polya-Gamma weights. Terms the move cannot change
# (the priors of A and of the scales sigma_b and sigma_emb) are left out.
axis_move_log_target <- function(x) {
  S <- ncol(x$B); gt <- nrow(x$C); emb <- which(x$cp$is_emb)
  sds <- occJSDM:::coef_prior_sds(x$cp, x$sigma_b, x$sigma_emb, x$lambda, nrow(x$G), ncol(x$X))
  R <- t(x$B) - x$Tr_b %*% x$G - x$A %*% x$C
  eta <- x$rest + x$X %*% x$B
  sum(dnorm(R, 0, rep(sds$b, each = S), log = TRUE)) +
    sum(dnorm(x$G, 0, sds$g, log = TRUE)) +
    sum(dnorm(x$C, 0, rep(sds$c, each = gt), log = TRUE)) +
    sum(dcauchy(x$lambda[emb], 0, x$cp$axis_scale, log = TRUE) + log(x$lambda[emb])) +
    sum(x$kappa * eta - x$Omega * eta^2 / 2)
}

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
  expect_identical(am$step, 0.5)
  expect_identical(dim(am$acceptance), c(8L, 2L))      # embedding columns by chains
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
    expect_identical(sw$out$B, mv$out$B, info = i)
    expect_identical(sw$out$lambda, mv$out$lambda, info = i)
    expect_identical(mv$step, 0.5, info = i)
    expect_identical(mv$kappa, sw$z - 0.5, info = i)
    eta <- matrix(sw$out$B0, nrow(mv$X), ncol(mv$B), byrow = TRUE) + mv$X %*% mv$B +
      sw$U_before %*% sw$out$L
    expect_equal(mv$eta, eta, tolerance = 1e-10, info = i)
  }
  expect_true(any(vapply(seen, function(mv) any(mv$out$accepted), logical(1))))
})

# Tolerances of the prior-invariance test. To be filled in from measurement
# once the move exists (see the commit that adds it).
