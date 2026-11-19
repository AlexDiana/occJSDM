# The convergence check for fits made with the grouped prior and an embedding
# (R/grouped-convergence.R; design in
# dev/superpowers/specs/2026-11-19-grouped-convergence-check-design.md).

# ---- the two checks of the embedding's combined effect ----------------------

effect_from <- function(scores, loadings) scores %*% t(loadings)
set.seed(1); n <- 40L; S <- 12L
U <- qr.Q(qr(matrix(rnorm(n * 3), n, 3))); V <- matrix(rnorm(S * 3), S, 3)

test_that("identical chains agree fully and do not disagree", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E1), c(n, S, 3L)))
  expect_true(got$assessed)
  expect_identical(got$n_gradients, 2L)
  expect_equal(got$agreement, 1); expect_equal(got$disagreement, 0)
  expect_true(all(is.na(diag(got$agreement_pairs))))
})

test_that("chains with the same gradients mixed differently agree fully", {
  turn <- matrix(c(cos(0.7), sin(0.7), -sin(0.7), cos(0.7)), 2, 2)
  scores <- U[, 1:2] %*% diag(c(9, 5))
  E1 <- effect_from(scores, V[, 1:2])
  E2 <- effect_from(scores %*% turn, V[, 1:2] %*% turn)      # the same effect, other blends
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E2), c(n, S, 2L)))
  expect_equal(got$agreement, 1); expect_equal(got$disagreement, 0)
})

test_that("a chain with other gradients is seen by both checks", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  E3 <- effect_from(U[, c(1, 3)] %*% diag(c(9, 5)), V[, 1:2])   # second gradient replaced
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E3), c(n, S, 3L)))
  expect_lt(got$agreement, 0.90); expect_gt(got$disagreement, 0.15)
  expect_equal(got$agreement_pairs[1, 2], 1)
  expect_lt(got$agreement_pairs[1, 3], 0.90)
})

test_that("a small or absent effect, or a single chain, is not assessed", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  small <- occJSDM:::grouped_effect_checks(array(c(E1, E1) * 1e-3, c(n, S, 2L)))
  expect_false(small$assessed); expect_match(small$why_not, "small")
  zero <- occJSDM:::grouped_effect_checks(array(0, c(n, S, 2L)))
  expect_false(zero$assessed); expect_identical(zero$size, 0)
  expect_true(is.na(zero$agreement) && is.na(zero$disagreement))
  expect_false(any(vapply(zero[c("size", "agreement", "disagreement")], is.nan, logical(1))))
  one <- occJSDM:::grouped_effect_checks(array(E1, c(n, S, 1L)))
  expect_false(one$assessed); expect_match(one$why_not, "chain")
})

# ---- the state and R-hat checks, the verdict, the exported function ---------

# The arguments of grouped_convergence() for a made-up fit: `p_meas` measured
# columns followed by 3 embedding columns, 6 species, 200 draws per chain. The
# embedding coefficients are a rank-2 pattern plus a little noise per draw,
# the same pattern in every chain or, with gradients = "other", with its
# second direction replaced in the last chain. The measured coefficients are
# noise about fixed values, shifted by `shift_measured` in the last chain.
# sigma_emb is drawn about the given median in each chain.
made_fit <- function(nchain, sigma_emb_median, shift_measured = 0, gradients = "same", p_meas = 2L) {
  set.seed(7); S <- 6L; niter <- 200L; n <- 30L; K <- 3L; p <- p_meas + K
  X <- scale(matrix(rnorm(n * p), n, p))
  L <- qr.Q(qr(matrix(rnorm(K * K), K, K))); W <- matrix(rnorm(S * 2), S, 2)
  pattern <- L[, 1:2] %*% diag(c(2, 1.4)) %*% t(W)
  other <- L[, c(1, 3)] %*% diag(c(2, 1.4)) %*% t(W)
  meas_mean <- matrix(rnorm(p_meas * S), p_meas, S)
  B <- array(0, c(p, S, niter, nchain))
  for (ch in seq_len(nchain)) for (it in seq_len(niter)) {
    emb <- if (gradients == "other" && ch == nchain) other else pattern
    B[p_meas + 1:K, , it, ch] <- emb + rnorm(K * S, sd = 0.05)
    if (p_meas > 0L)
      B[seq_len(p_meas), , it, ch] <- meas_mean + rnorm(p_meas * S, sd = 0.1) +
        if (ch == nchain) shift_measured else 0
  }
  list(B_output = B,
       sigmab_output = matrix(exp(-0.5 + 0.1 * rnorm(niter * nchain)), niter, nchain),
       sigmaemb_output = vapply(sigma_emb_median, function(m) exp(log(m) + 0.2 * rnorm(niter)), numeric(niter)),
       lambda_output = array(exp(-2 + 0.2 * rnorm(K * niter * nchain)), c(K, niter, nchain)),
       X = X, is_emb = rep(c(FALSE, TRUE), c(p_meas, K)))
}
gc_of <- function(...) do.call(occJSDM:::grouped_convergence, made_fit(...))
check_row <- function(x, name) x$checks[x$checks$check == name, ]

test_that("a fit whose chains agree and found shared structure is ok", {
  got <- gc_of(3L, c(0.3, 0.3, 0.3))
  expect_identical(got$verdict, "ok")
  expect_identical(got$checks$check, c("state", "rhat_measured", "rhat_scales", "gradients",
                                       "disagreement", "disagreement_beyond"))
  expect_identical(names(got$checks), c("check", "value", "limit", "flagged", "meaning"))
  expect_false(any(got$checks$flagged))
  expect_identical(got$chains$state, rep("shared structure", 3L))
  expect_equal(got$chains$sigma_emb_median, c(0.3, 0.3, 0.3), tolerance = 0.1)
  expect_true(got$effect$assessed)
  expect_identical(check_row(got, "gradients")$limit, 0.90)
  expect_identical(check_row(got, "disagreement")$limit, 0.15)
})

test_that("chains in different states are flagged", {
  got <- gc_of(2L, c(0.3, 6))
  row <- check_row(got, "state")
  expect_true(row$flagged); expect_identical(row$value, 1)
  expect_match(row$meaning, "different states", fixed = TRUE)
  expect_identical(got$chains$state, c("shared structure", "leftovers dominate"))
  expect_identical(got$verdict, "review")
})

test_that("every chain with leftovers dominating is flagged, with both readings given", {
  got <- gc_of(2L, c(6, 7))
  row <- check_row(got, "state")
  expect_true(row$flagged); expect_identical(row$value, 2)
  expect_match(row$meaning, "respond independently", fixed = TRUE)
  expect_match(row$meaning, "5,000 burn-in", fixed = TRUE)
  expect_identical(got$verdict, "review")
})

test_that("a measured coefficient that differs between chains is flagged by R-hat", {
  got <- gc_of(2L, c(0.3, 0.3), shift_measured = 5)
  row <- check_row(got, "rhat_measured")
  expect_true(row$flagged); expect_gt(row$value, 1.1)
  expect_false(check_row(got, "rhat_scales")$flagged)
  expect_identical(got$verdict, "review")
})

test_that("a chain with other gradients is flagged", {
  got <- gc_of(3L, c(0.3, 0.3, 0.3), gradients = "other")
  expect_true(check_row(got, "gradients")$flagged)
  expect_lt(check_row(got, "gradients")$value, 0.90)
  expect_identical(got$verdict, "review")
})

test_that("with one chain only the state is assessed", {
  got <- gc_of(1L, 0.3)
  expect_identical(check_row(got, "state")$value, 0)
  for (nm in c("rhat_measured", "rhat_scales", "gradients", "disagreement", "disagreement_beyond")) {
    expect_true(is.na(check_row(got, nm)$value), info = nm)
    expect_false(check_row(got, nm)$flagged, info = nm)
  }
  expect_match(check_row(got, "gradients")$meaning, "two chains", fixed = TRUE)
  expect_identical(got$verdict, "ok")
})

test_that("with no measured covariate the R-hat check covers the scales only", {
  got <- gc_of(2L, c(0.3, 0.3), p_meas = 0L)
  expect_true(is.na(check_row(got, "rhat_measured")$value))
  expect_false(check_row(got, "rhat_measured")$flagged)
  expect_true(is.finite(check_row(got, "rhat_scales")$value))
})

test_that("how far the chains' coefficients differ is reported and never flagged", {
  # Beyond the fitted sites the chains can differ where the fitted effects
  # agree: a difference between chains' coefficients that lies in a direction
  # the sites do not vary in changes nothing at those sites.
  x <- made_fit(2L, c(0.3, 0.3))
  same <- do.call(occJSDM:::grouped_convergence, x)
  row <- check_row(same, "disagreement_beyond")
  expect_false(row$flagged); expect_true(is.na(row$limit))
  expect_match(row$meaning, "unlike the fitted", fixed = TRUE)
  expect_lt(row$value, 0.05)
  expect_match(check_row(same, "gradients")$meaning, "fitted sites", fixed = TRUE)
  expect_match(check_row(same, "disagreement")$meaning, "fitted sites", fixed = TRUE)
  # Remove one embedding direction from the sites, and move chain 2's
  # coefficients along it: nothing changes at the sites, the coefficients do.
  emb <- which(x$is_emb); v <- c(1, -1, 0.5); v <- v / sqrt(sum(v^2))
  x$X[, emb] <- x$X[, emb] - (x$X[, emb] %*% v) %*% t(v)
  before <- do.call(occJSDM:::grouped_convergence, x)
  x$B_output[emb, , , 2] <- x$B_output[emb, , , 2] + 3 * v
  after <- do.call(occJSDM:::grouped_convergence, x)
  expect_equal(check_row(after, "disagreement")$value, check_row(before, "disagreement")$value)
  expect_equal(check_row(after, "gradients")$value, check_row(before, "gradients")$value)
  expect_gt(check_row(after, "disagreement_beyond")$value, 10 * check_row(before, "disagreement_beyond")$value)
  expect_false(check_row(after, "disagreement_beyond")$flagged)
  expect_identical(after$verdict, before$verdict)
})

gc_fit <- function(sim, embedding = TRUE) {
  dat <- sim$data_list
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = c(list(n_factors = 1L), if (embedding) list(embeddingVars = dat$embeddingVars)),
    occCovariates = grouped_occ_covariates(dat), MCMCparams = FIXTURE_MCMC)))
}

test_that("returnGroupedConvergence() works on fits with embeddingVars and refuses others", {
  fit <- gc_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L))
  got <- returnGroupedConvergence(fit)
  expect_identical(names(got), c("verdict", "checks", "chains", "effect"))
  expect_true(got$verdict %in% c("ok", "review"))
  expect_identical(nrow(got$checks), 6L)
  expect_identical(nrow(got$chains), 2L)
  plain <- gc_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L), embedding = FALSE)
  expect_error(returnGroupedConvergence(plain), "returnGroupedConvergence() needs a fit with embeddingVars", fixed = TRUE)
  # Nothing in it depends on the model type.
  occ <- gc_fit(simulate_grouped_fixture(model = "occupancy", p_emb = 6L, S = 8L))
  expect_identical(nrow(returnGroupedConvergence(occ)$checks), 6L)
})

# ---- the automatic report ---------------------------------------------------

# A results_output for computeDiagnostics(), from made_fit(): rows 1 and 2 of
# B are measured, rows 3 to 5 embedding. One coefficient can be moved in
# chain 2 so that its R-hat is far over 1.1.
made_results <- function(shift_embedding = 0, shift_measured = 0) {
  x <- made_fit(2L, c(0.3, 0.3))
  x$B_output[4, 2, , 2] <- x$B_output[4, 2, , 2] + shift_embedding
  x$B_output[1, 3, , 2] <- x$B_output[1, 3, , 2] + shift_measured
  set.seed(8)
  list(jsdm_output = list(B0_output = array(rnorm(6 * 200 * 2), c(6L, 200L, 2L)), B_output = x$B_output,
                          sigmab_output = x$sigmab_output, sigmaemb_output = x$sigmaemb_output,
                          lambda_output = x$lambda_output))
}
said <- function(expr) {
  messages <- character(0); warnings <- character(0)
  withCallingHandlers(expr,
    message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(messages = messages, warnings = warnings)
}
has <- function(x, what) any(grepl(what, x, fixed = TRUE))
emb_rows <- c(FALSE, FALSE, TRUE, TRUE, TRUE)

test_that("computeDiagnostics() does not warn on a single embedding coefficient when told which rows they are", {
  results <- made_results(shift_embedding = 5)
  before <- said(computeDiagnostics(results))
  expect_true(has(before$warnings, "Convergence issue in 'B_output'"))
  expect_true(has(before$messages, "Parameter Set: B (Flattened N = 30)"))
  expect_false(has(before$messages, "B (measured)"))
  after <- said(computeDiagnostics(results, embedding_rows = emb_rows))
  expect_identical(after$warnings, character(0))
  expect_true(has(after$messages, "Parameter Set: B (measured) (Flattened N = 12)"))
  expect_true(has(after$messages, "Parameter Set: B (embedding) (Flattened N = 18)"))
  expect_true(has(after$messages, "Not a test of convergence under the grouped prior: see returnGroupedConvergence()"))
  for (nm in c("log sigmab", "log sigmaemb", "log lambda"))
    expect_true(has(after$messages, paste0("Parameter Set: ", nm)), info = nm)
})

test_that("computeDiagnostics() still warns on a measured coefficient", {
  results <- made_results(shift_measured = 5)
  after <- said(computeDiagnostics(results, embedding_rows = emb_rows))
  expect_true(has(after$warnings, "Convergence issue in 'B (measured)'"))
  expect_false(has(after$warnings, "embedding"))
})

test_that("computeDiagnostics() reports as before when not told the embedding rows", {
  results <- made_results()
  a <- said(computeDiagnostics(results)); b <- said(computeDiagnostics(results, embedding_rows = NULL))
  expect_identical(a, b)
  expect_true(has(a$messages, "Parameter Set: lambda")); expect_false(has(a$messages, "log lambda"))
  expect_false(has(a$messages, "Not a test of convergence"))
})

test_that("the summary of the checks names the verdict and explains every flag", {
  flagged <- said(occJSDM:::print_grouped_convergence(gc_of(2L, c(0.3, 6))))$messages
  expect_true(has(flagged, "Grouped prior: convergence checks"))
  expect_true(any(grepl("Verdict +\\| review", flagged)))
  expect_true(any(grepl("^  ! .*different states", flagged)))
  ok <- said(occJSDM:::print_grouped_convergence(gc_of(3L, c(0.3, 0.3, 0.3))))$messages
  expect_true(any(grepl("Verdict +\\| ok", ok)))
  expect_false(any(grepl("^  ! ", ok)))
  expect_true(any(grepl("Gradients", ok)) && any(grepl("beyond", ok)))
})

test_that("a fit with embeddingVars ends with the checks, and warns when they flag it", {
  dat <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)$data_list
  fit <- NULL
  grouped <- said(fit <- runOccJSDM(dat, listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
                                    occCovariates = grouped_occ_covariates(dat), MCMCparams = FIXTURE_MCMC))
  expect_identical(sum(grepl("Grouped prior: convergence checks", grouped$messages, fixed = TRUE)), 1L)
  expect_true(has(grouped$messages, "B (embedding)"))
  expect_identical(has(grouped$warnings, "The grouped prior's convergence checks flag this fit for review: see returnGroupedConvergence()."),
                   returnGroupedConvergence(fit)$verdict == "review")
  expect_false(has(grouped$warnings, "Convergence issue in 'B_output'"))
  plain <- said(runOccJSDM(dat, listParams = list(n_factors = 1L),
                           occCovariates = grouped_occ_covariates(dat), MCMCparams = FIXTURE_MCMC))
  expect_false(has(plain$messages, "Grouped prior: convergence checks"))
  expect_false(has(plain$messages, "B (embedding)"))
})
