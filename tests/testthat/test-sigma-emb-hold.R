# The staged burn-in of the grouped prior: sigma_emb is held small for the
# first half of burn-in (R/coef-prior.R, R/chains.R). A prototype: fix attempt
# 3 for the slow mixing recorded in dev/simstudy/grouped-mixing-20261118.

test_that("the hold lasts the first half of burn-in at a tenth of the prior scale", {
  expect_identical(vapply(c(0, 1, 2, 3, 20, 1500, 5001), occJSDM:::axis_hold_sweeps, integer(1)),
                   c(0L, 0L, 1L, 1L, 10L, 750L, 2500L))
  expect_equal(occJSDM:::axis_hold_value(list(sigma_emb_scale = 1)), 0.1)
  expect_equal(occJSDM:::axis_hold_value(list(sigma_emb_scale = 0.3)), 0.03)
})

test_that("a held sigma_emb is not drawn and the axis scales use the held value", {
  set.seed(1); S <- 7L; gb <- 2L; gt <- 2L
  R <- matrix(rnorm(S * 5), S, 5); G <- matrix(rnorm(gb * 5), gb, 5); C <- matrix(rnorm(gt * 5), gt, 5)
  cp <- list(type = "grouped", is_emb = c(FALSE, FALSE, TRUE, TRUE, TRUE),
             sigma_b_scale = 1, sigma_emb_scale = 1, axis_scale = 1 / sqrt(3),
             c_base = c(1, 1, rep(1 / sqrt(gt), 3)))
  lambda <- c(1, 1, 0.5, 0.2, 0.9)
  hc <- function(ss, n, sigma, scale) {
    aux <- occJSDM:::rinvgamma_cpp(1, 1 / sigma^2 + 1 / scale^2)
    sqrt(occJSDM:::rinvgamma_cpp((n + 1) / 2, ss / 2 + 1 / aux))
  }
  # The same random numbers as the usual update with the sigma_emb draw left out.
  seed_both(42)
  want_b <- hc(sum(R[, 1:2]^2), S * 2, 0.7, 1)
  want_l <- vapply(3:5, function(k) hc(
    sum(R[, k]^2) / 0.1^2 + sum(G[, k]^2) + sum(C[, k]^2) / cp$c_base[k]^2,
    S + gb + gt, lambda[k], cp$axis_scale), numeric(1))
  seed_both(42)
  got <- occJSDM:::sample_coef_scales(R, G, C, cp, 0.7, 0.1, lambda, hold_sigma_emb = TRUE)
  expect_identical(got$sigma_emb, 0.1)
  expect_equal(got, list(sigma_b = want_b, sigma_emb = 0.1, lambda = c(1, 1, want_l)))
  # Not held: sigma_emb is drawn, as before.
  seed_both(42)
  free <- occJSDM:::sample_coef_scales(R, G, C, cp, 0.7, 0.1, lambda)
  expect_false(free$sigma_emb == 0.1)
})

hold_fit <- function(MCMCparams, spy = TRUE) {
  seen <- new.env(); seen$sigma <- list(); seen$sweeps <- list()
  original_sigma <- occJSDM:::sigma_emb_move
  update <- occJSDM:::update_jSDMcoef
  update_env <- new.env(parent = environment(update))
  update_env$sigma_emb_move <- function(...) {
    out <- original_sigma(...)
    seen$sigma[[length(seen$sigma) + 1L]] <- list(sweep = length(seen$sweeps) + 1L, step = list(...)[[9]], out = out)
    out
  }
  environment(update) <- update_env
  local_mocked_bindings(update_jSDMcoef = function(list_data, list_params, ...) {
    out <- update(list_data, list_params, ...)
    seen$sweeps[[length(seen$sweeps) + 1L]] <- list(before = list_params$sigma_emb, after = out$sigma_emb)
    out
  }, .package = "occJSDM", .env = parent.frame())
  dat <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)$data_list
  fit <- suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
    occCovariates = grouped_occ_covariates(dat), MCMCparams = MCMCparams)))
  list(fit = fit, sigma = seen$sigma, sweeps = seen$sweeps)
}

test_that("a chain holds sigma_emb for the first half of burn-in and then releases it", {
  got <- hold_fit(list(nchain = 1L, nburn = 6L, niter = 2L, nthin = 1L))
  sw <- got$sweeps
  expect_identical(length(sw), 8L)
  # The chain starts at the held value and keeps it for sweeps 1 to 3.
  expect_identical(sw[[1]]$before, 0.1)
  for (i in 1:3) expect_identical(sw[[i]]$after, 0.1, info = paste("sweep", i))
  # From sweep 4 it is drawn again.
  for (i in 4:8) expect_false(sw[[i]]$after == 0.1, info = paste("sweep", i))
  # The scale move on sigma_emb is skipped while it is held: sweeps 4 to 8 only.
  expect_identical(vapply(got$sigma, function(s) s$sweep, integer(1)), 4:8)
  # Its step starts at 0.5 on release and is adapted in the burn-in sweeps that
  # remain (4, 5 and 6), counting from the release, and is then frozen.
  step <- lapply(got$sigma, function(s) s$step)
  expect_identical(step[[1]], 0.5)
  for (j in 1:3)
    expect_identical(step[[j + 1L]], occJSDM:::axis_step_adapt(step[[j]], got$sigma[[j]]$out$log_ratio, j),
                     info = paste("proposal", j))
  expect_identical(step[[5]], step[[4]])
  am <- got$fit$infos$axis_move
  expect_identical(am$hold_sweeps, 3L)
  expect_identical(am$hold_value, 0.1)
})

test_that("with too little burn-in to halve there is no hold", {
  got <- hold_fit(list(nchain = 1L, nburn = 1L, niter = 2L, nthin = 1L))
  expect_identical(got$sweeps[[1]]$before, 1)          # the usual start, the prior scale
  expect_identical(vapply(got$sigma, function(s) s$sweep, integer(1)), 1:3)
  expect_identical(got$fit$infos$axis_move$hold_sweeps, 0L)
})
