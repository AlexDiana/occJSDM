# Manufactured posterior draws permit analytic likelihood checks, without
# asserting that a stochastic fitted value equals a particular constant.
waic_test_fit <- function(y = matrix(c(0, 1), 2, 1), factors = 1L) {
  n <- nrow(y)
  S <- ncol(y)
  iterations <- 3L
  list(
    infos = list(model = "binary", ps = 0, n = n, S = S,
      n_factors = factors, siteNames = paste0("site", seq_len(n)),
      speciesNames = paste0("sp", seq_len(S)), OTU = y),
    X_psi = matrix(numeric(), n, 0),
    results_output = list(WAIC = -999, jsdm_output = list(
      B0_output = array(0, c(S, iterations, 1)),
      B_output = array(numeric(), c(0, S, iterations, 1)),
      L_output = array(0.8, c(factors, S, iterations, 1)),
      sigmah_output = matrix(1, iterations, 1))))
}

test_that("the default WAIC scores observed states rather than the saved latent score", {
  fit <- waic_test_fit()
  # Symmetry gives P(y = 0) = P(y = 1) = 1/2 exactly.
  expect_equal(extractWAIC(fit), -4 * log(0.5), tolerance = 1e-8)
  expect_equal(extractWAIC(fit, type = "legacy"), -999)
})

test_that("site integration preserves dependence between species", {
  fit <- waic_test_fit(matrix(c(1, 1, 0, 1), 2, 2, byrow = TRUE))
  expected <- vapply(seq_len(2), function(i) {
    log(integrate(function(u) {
      p <- plogis(0.8 * u)
      dbinom(fit$infos$OTU[i, 1], 1, p) *
        dbinom(fit$infos$OTU[i, 2], 1, p) * dnorm(u)
    }, -Inf, Inf, rel.tol = 1e-10)$value)
  }, numeric(1))
  result <- computeSiteWAIC(fit)
  expect_equal(unname(result$log_lik[1, ]), expected, tolerance = 1e-6)
  expect_gt(result$log_lik[1, 1], log(0.25))
  expect_lt(result$log_lik[1, 2], log(0.25))
  # Hidden fitted states must have no influence on the score.
  fit$results_output$z_output <- matrix(1, 2, 2)
  fit$results_output$jsdm_output$U_output <- array(1000, c(2, 1, 3, 1))
  expect_equal(computeSiteWAIC(fit)$log_lik, result$log_lik)
})

test_that("log likelihood summaries stay stable far below floating point underflow", {
  fit <- waic_test_fit(matrix(0, 2, 2), factors = 0L)
  fit$results_output$jsdm_output$B0_output[] <- 1000
  result <- computeSiteWAIC(fit)
  expect_equal(unname(result$log_lik), matrix(-2000, 3, 2))
  expect_equal(result$WAIC, 8000)
})

waic_two_stage_fit <- function() {
  fit <- waic_test_fit(matrix(c(2, NA, 0, 4, 1), 5, 1), factors = 0L)
  fit$infos$model <- "two_stage"
  fit$infos$n <- 2L
  fit$infos$siteNames <- c("one", "two")
  fit$infos$primerNames <- c("A", "B")
  fit$infos$threshold <- 2
  fit$infos$list_idx <- list(idx_z_w = c(1L, 1L, 2L),
    idx_w_k = c(1L, 1L, 2L, 3L, 3L), idx_p_k = c(2L, 2L, 1L, 2L, 1L))
  fit$X_psi <- matrix(numeric(), 2, 0)
  fit$X_theta <- cbind(1, c(-1, 0, 1))
  fit$results_output$jsdm_output$B0_output[] <- qlogis(0.4)
  fit$results_output$beta_theta_output <- array(rep(c(0.3, -0.7), 3), c(2, 1, 3, 1))
  fit$results_output$theta0_output <- array(0.05, c(1, 3, 1))
  fit$results_output$p_output <- array(rep(c(0.7, 0.9), 3), c(2, 1, 3, 1))
  fit$results_output$q_output <- array(rep(c(0.1, 0.2), 3), c(2, 1, 3, 1))
  fit
}

test_that("two-stage likelihood sums all hidden states with actual sample and primer identities", {
  fit <- waic_two_stage_fit()
  theta <- plogis(0.3 - 0.7 * c(-1, 0, 1))
  # Independent explicit enumeration of all z and w possibilities.
  expected <- vapply(1:2, function(site) {
    samples <- which(fit$infos$list_idx$idx_z_w == site)
    states <- expand.grid(rep(list(0:1), length(samples) + 1L))
    sum(apply(states, 1, function(state) {
      z <- state[1]
      w <- state[-1]
      value <- dbinom(z, 1, 0.4)
      for (m in seq_along(samples)) {
        value <- value * dbinom(w[m], 1, if (z == 1) theta[samples[m]] else 0.05)
        rows <- which(fit$infos$list_idx$idx_w_k == samples[m])
        for (r in rows) {
          if (is.na(fit$infos$OTU[r, 1])) next
          primer <- fit$infos$list_idx$idx_p_k[r]
          p <- if (w[m] == 1) c(0.7, 0.9)[primer] else c(0.1, 0.2)[primer]
          value <- value * dbinom(as.integer(fit$infos$OTU[r, 1] >= 2), 1, p)
        }
      }
      value
    }))
  }, numeric(1))
  result <- computeSiteWAIC(fit)
  expect_equal(unname(result$log_lik[1, ]), log(expected), tolerance = 1e-12)
  expect_error(computeSiteWAIC(fit, threshold = 1), "threshold.*fit")
  fit$infos$threshold <- NULL
  expect_error(computeSiteWAIC(fit), "threshold")
  expect_equal(computeSiteWAIC(fit, threshold = 2)$log_lik, result$log_lik)
})

test_that("all-missing sites are omitted and an empty survey cannot be scored", {
  fit <- waic_two_stage_fit()
  fit$infos$OTU[4:5, ] <- NA
  result <- computeSiteWAIC(fit)
  expect_equal(result$pointwise$site, "one")
  expect_equal(result$omitted_sites, "two")
  expect_true(is.na(result$SE))
  fit$infos$OTU[] <- NA
  expect_error(computeSiteWAIC(fit), "observed")
})

test_that("matched iteration and chain draws supply environmental and detection parameters", {
  fit <- waic_test_fit(factors = 0L)
  fit$X_psi <- matrix(c(-1, 2), 2, 1)
  fit$results_output$jsdm_output$B0_output <- array(seq(-0.2, 0.3, by = 0.1), c(1, 3, 2))
  fit$results_output$jsdm_output$B_output <- array(seq(0.1, 0.6, by = 0.1), c(1, 1, 3, 2))
  fit$results_output$jsdm_output$L_output <- array(numeric(), c(0, 1, 3, 2))
  fit$results_output$jsdm_output$sigmah_output <- matrix(1, 3, 2)
  fit$results_output$jsdm_output$B0_output[1, 1, 2] <- -0.4
  result <- computeSiteWAIC(fit, draws = c(4L, 2L))
  expected <- rbind(
    dbinom(c(0, 1), 1, plogis(-0.4 + c(-1, 2) * 0.4), log = TRUE),
    dbinom(c(0, 1), 1, plogis(-0.1 + c(-1, 2) * 0.2), log = TRUE))
  expect_equal(unname(result$log_lik), expected)
  expect_equal(result$draw_ids$iteration, c(1L, 2L))
  expect_equal(result$draw_ids$chain, c(2L, 1L))
  expect_equal(result$p_waic, sum(apply(expected, 2, var)))
  expect_error(computeSiteWAIC(fit, draws = 1), "two")
  expect_error(computeSiteWAIC(fit, draws = c(1, 1)), "unique")
})

test_that("unsupported targets and inaccurate integration fail explicitly", {
  fit <- waic_test_fit()
  fit$infos$ps <- 1
  expect_error(computeSiteWAIC(fit), "non-spatial")
  expect_equal(extractWAIC(fit, type = "legacy"), -999)
  fit$infos$ps <- 0
  fit$infos$model <- "continuous"
  expect_error(computeSiteWAIC(fit), "binary.*occupancy.*two_stage")
  fit <- waic_test_fit(matrix(c(1, 1, 0, 1), 2, 2, byrow = TRUE))
  expect_error(computeSiteWAIC(fit, quadrature = c(2, 3), tolerance = 1e-12), "quadrature")
  expect_error(computeSiteWAIC(fit, max_nodes = 2), "max_nodes")
})

test_that("WAIC comparisons use paired sites and reject different observations", {
  first <- computeSiteWAIC(waic_test_fit())
  fit <- waic_test_fit()
  fit$results_output$jsdm_output$B0_output[] <- 0.5
  second <- computeSiteWAIC(fit)
  difference <- first$pointwise$waic - second$pointwise$waic
  comparison <- compareSiteWAIC(first, second)
  expect_equal(comparison$difference, sum(difference))
  expect_equal(comparison$SE, sqrt(2 * var(difference)))
  fit$infos$OTU[1, 1] <- 1
  expect_error(compareSiteWAIC(first, computeSiteWAIC(fit)), "same.*observations")
})

test_that("occupancy likelihood integrates site states across field samples", {
  fit <- waic_two_stage_fit()
  fit$infos$model <- "occupancy"
  fit$infos$OTU <- matrix(c(1, 0, 1), 3, 1)
  fit$infos$threshold <- 1
  theta <- plogis(0.3 - 0.7 * c(-1, 0, 1))
  expected <- c(0.4 * theta[1] * (1 - theta[2]) + 0.6 * 0.05 * 0.95,
                0.4 * theta[3] + 0.6 * 0.05)
  expect_equal(unname(computeSiteWAIC(fit)$log_lik[1, ]), log(expected))
})

test_that("joint likelihood is invariant to rotations and uses the factor SD", {
  fit <- waic_test_fit(matrix(c(1, 1, 0, 1, 0, 1), 2, 3, byrow = TRUE), factors = 2L)
  loadings <- matrix(c(0.2, 0.7, -0.5, 0.1, 0.8, -0.4), 2, 3)
  fit$results_output$jsdm_output$L_output[] <- rep(loadings, 3)
  fit$results_output$jsdm_output$sigmah_output[] <- 0.6
  set.seed(3)
  rng <- .Random.seed
  score <- computeSiteWAIC(fit, tolerance = 1e-8)
  expect_identical(.Random.seed, rng)
  rotation <- matrix(c(cos(0.7), -sin(0.7), sin(0.7), cos(0.7)), 2, 2)
  fit$results_output$jsdm_output$L_output[] <- rep(rotation %*% loadings * 0.6, 3)
  fit$results_output$jsdm_output$sigmah_output[] <- 1
  expect_equal(computeSiteWAIC(fit, tolerance = 1e-8)$log_lik, score$log_lik, tolerance = 1e-7)
})

test_that("new fits retain thresholds and work with only summarized latent states", {
  for (model in c("binary", "occupancy", "two_stage")) {
    fit <- fit_fixture(simulate_fixture(model = model, useSpatField = FALSE), spatial = FALSE, d = 1L)
    score <- suppressWarnings(computeSiteWAIC(fit, draws = c(1L, 21L)))
    expect_true(is.finite(score$WAIC))
    expect_equal(nrow(score$pointwise), fit$infos$n)
    expect_equal(score$draw_ids$chain, c(1L, 2L))
    if (model != "binary") expect_equal(fit$infos$threshold, 1)
  }
})

test_that("uncertain and impossible likelihoods are reported rather than silently dropped", {
  fit <- waic_test_fit(factors = 0L)
  fit$results_output$jsdm_output$B0_output[] <- c(-4, 0, 4)
  expect_warning(score <- computeSiteWAIC(fit), "p_waic > 0.4")
  expect_gt(score$p_waic, 0.8)
  fit <- waic_two_stage_fit()
  fit$results_output$p_output[] <- 0
  fit$results_output$q_output[] <- 0
  expect_error(computeSiteWAIC(fit), "Nonfinite site log likelihood")
  fit <- waic_two_stage_fit()
  fit$infos$list_idx$idx_p_k[1] <- 3L
  expect_error(computeSiteWAIC(fit), "idx_p_k")
  fit <- waic_test_fit()
  fit$results_output$jsdm_output$L_output[1] <- NA
  expect_error(computeSiteWAIC(fit), "L_output")
})

test_that("binary comparisons retain real site identities even with repeated response rows", {
  fit <- waic_test_fit(matrix(1, 2, 1), factors = 0L)
  # Binary fits historically store siteNames=1:n, not data_info$Site.
  fit$infos$siteNames <- 1:2
  fit$infos$data_info <- data.frame(Site = c("A", "B"), environment = c(-1, 1))
  fit$X_psi <- matrix(c(-1, 1), 2, 1)
  fit$results_output$jsdm_output$B_output <- array(1, c(1, 1, 3, 1))
  first <- computeSiteWAIC(fit)
  expect_equal(first$pointwise$site, c("A", "B"))
  fit$infos$data_info <- fit$infos$data_info[2:1, ]
  fit$X_psi <- fit$X_psi[2:1, , drop = FALSE]
  expect_error(compareSiteWAIC(first, computeSiteWAIC(fit)), "same.*observations")
})
