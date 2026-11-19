# The convergence check for fits made with the grouped prior and an embedding
# (R/grouped-convergence.R; design in
# dev/superpowers/specs/2026-11-19-grouped-convergence-check-design.md).

# ---- the two checks of the embedding's combined effect ----------------------

effect_from <- function(scores, loadings) scores %*% t(loadings)
set.seed(1); n <- 40L; S <- 12L
U <- qr.Q(qr(matrix(rnorm(n * 3), n, 3))); V <- matrix(rnorm(S * 3), S, 3)

rms <- function(x) sqrt(mean(x^2))
set.seed(2); rest <- matrix(rnorm(n * S, sd = 0.35), n, S)     # what two gradients do not hold

test_that("identical chains agree fully and do not disagree", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E1), c(n, S, 3L)))
  expect_true(got$assessed)
  expect_identical(got$n_gradients, 2L)
  expect_equal(got$agreement, 1); expect_equal(got$disagreement, 0)
  expect_true(all(is.na(diag(got$agreement_pairs))))
  expect_true(all(is.na(diag(got$disagreement_pairs))))
})

test_that("identical chains agree fully whatever share of the effect their gradients hold", {
  # The gradients hold 90 percent or more of the effect, not all of it. The
  # agreement is the share of a chain's effect within another's gradients as a
  # share of what its own gradients hold, so chains that are the same score 1,
  # and the limit is not the score of a perfect fit.
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2]) + rest
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1), c(n, S, 2L)))
  d <- svd(E1)$d
  expect_lt(sum(d[1:2]^2) / sum(d^2), 0.99)
  expect_identical(got$n_gradients, 2L)
  expect_equal(got$agreement, 1); expect_equal(got$disagreement, 0)
})

test_that("chains with the same gradients and other responses agree on the gradients and differ in effect", {
  mix <- matrix(c(1, 0.4, -0.6, 1.3), 2, 2)                  # not a rotation: the effect changes
  scores <- U[, 1:2] %*% diag(c(9, 5))
  E1 <- effect_from(scores, V[, 1:2])
  E2 <- effect_from(scores %*% mix, V[, 1:2])                # the same span across sites
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E2), c(n, S, 2L)))
  expect_equal(got$agreement, 1)
  expect_equal(got$disagreement, rms(E1 - E2) / rms((E1 + E2) / 2))
})

test_that("the disagreement is the largest difference between two chains, however many agree", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  E3 <- effect_from(U[, c(1, 3)] %*% diag(c(9, 5)), V[, 1:2])
  two <- occJSDM:::grouped_effect_checks(array(c(E1, E3), c(n, S, 2L)))
  five <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E1, E1, E3), c(n, S, 5L)))
  apart <- rms(E1 - E3)
  expect_equal(two$disagreement, apart / rms((E1 + E3) / 2))
  expect_equal(five$disagreement, apart / rms((4 * E1 + E3) / 5))
  expect_equal(five$disagreement_pairs[1, 5], five$disagreement)
  expect_equal(five$disagreement_pairs[5, 1], five$disagreement)
  expect_equal(five$disagreement_pairs[1, 2], 0)
  expect_gt(five$disagreement, 0.8 * two$disagreement)       # one chain apart is not diluted
})

test_that("a chain with other gradients is seen by both checks", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  E3 <- effect_from(U[, c(1, 3)] %*% diag(c(9, 5)), V[, 1:2])   # second gradient replaced
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E3), c(n, S, 3L)))
  expect_lt(got$agreement, occJSDM:::GROUPED_LIMITS$gradients)
  expect_gt(got$disagreement, occJSDM:::GROUPED_LIMITS$disagreement)
  expect_equal(got$agreement_pairs[1, 2], 1)
  expect_lt(got$agreement_pairs[1, 3], occJSDM:::GROUPED_LIMITS$gradients)
})

test_that("an effect with missing values is not assessed and says why", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  E2 <- E1; E2[3, 4] <- NaN
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E2), c(n, S, 2L)))
  expect_false(got$assessed); expect_true(got$missing)
  expect_match(got$why_not, "missing", fixed = TRUE)
  expect_false(occJSDM:::grouped_effect_checks(array(c(E1, E1), c(n, S, 2L)))$missing)
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
  expect_identical(got$checks$check, c("state", "rhat_measured", "rhat_scales", "rhat_axis", "gradients",
                                       "disagreement", "disagreement_beyond"))
  expect_identical(names(got$checks), c("check", "value", "limit", "flagged", "meaning"))
  expect_false(any(got$checks$flagged))
  expect_identical(got$chains$state, rep("shared structure", 3L))
  expect_equal(got$chains$sigma_emb_median, c(0.3, 0.3, 0.3), tolerance = 0.1)
  expect_true(got$effect$assessed)
  expect_identical(check_row(got, "gradients")$limit, 0.90)
  expect_identical(check_row(got, "disagreement")$limit, 0.30)
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

test_that("the two scales shared by all columns decide the verdict, the axis scales do not", {
  # One axis scale differs between the chains: reported, with how many are
  # over the limit, and the verdict is unchanged. Single axis scales, like
  # single embedding coefficients, are not pinned down, and in full-length
  # fits whose chains agree a few of them fail R-hat.
  x <- made_fit(3L, c(0.3, 0.3, 0.3))
  x$lambda_output[2, , 3] <- x$lambda_output[2, , 3] * 20
  got <- do.call(occJSDM:::grouped_convergence, x)
  row <- check_row(got, "rhat_axis")
  expect_gt(row$value, 1.1); expect_false(row$flagged); expect_true(is.na(row$limit))
  expect_match(row$meaning, "1 of 3 axis scales", fixed = TRUE)
  expect_match(row$meaning, "not part of the verdict", fixed = TRUE)
  expect_false(check_row(got, "rhat_scales")$flagged)
  expect_identical(got$verdict, "ok")
  quiet <- check_row(do.call(occJSDM:::grouped_convergence, made_fit(3L, c(0.3, 0.3, 0.3))), "rhat_axis")
  expect_lt(quiet$value, 1.1); expect_match(quiet$meaning, "0 of 3 axis scales", fixed = TRUE)
  # sigma_emb differs between the chains, all of them under 1: flagged.
  got <- gc_of(2L, c(0.15, 0.8))
  expect_identical(check_row(got, "state")$value, 0)
  expect_true(check_row(got, "rhat_scales")$flagged)
  expect_identical(got$verdict, "review")
  # sigma_b differs between the chains: flagged.
  x <- made_fit(2L, c(0.3, 0.3)); x$sigmab_output[, 2] <- x$sigmab_output[, 2] * 5
  expect_true(check_row(do.call(occJSDM:::grouped_convergence, x), "rhat_scales")$flagged)
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
  for (nm in c("rhat_measured", "rhat_scales", "rhat_axis", "gradients", "disagreement", "disagreement_beyond")) {
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

test_that("how far the chains' coefficients differ is the largest difference between two chains", {
  x <- made_fit(3L, c(0.3, 0.3, 0.3), gradients = "other")
  got <- do.call(occJSDM:::grouped_convergence, x)
  emb <- which(x$is_emb)
  Bm <- lapply(1:3, function(ch) apply(x$B_output[emb, , , ch], c(1, 2), mean))
  pooled <- Reduce(`+`, Bm) / 3
  pairs <- c(rms(Bm[[1]] - Bm[[2]]), rms(Bm[[1]] - Bm[[3]]), rms(Bm[[2]] - Bm[[3]]))
  expect_equal(check_row(got, "disagreement_beyond")$value, max(pairs) / rms(pooled))
})

test_that("the checks' values are those of the pieces computed apart", {
  # The whole of B goes through one pass that copies none of it; the values
  # must be the ones the rows taken out and summarised by themselves give.
  x <- made_fit(3L, c(0.3, 0.3, 0.3), shift_measured = 0.05)
  got <- do.call(occJSDM:::grouped_convergence, x)
  meas <- x$B_output[!x$is_emb, , , , drop = FALSE]
  expect_equal(check_row(got, "rhat_measured")$value,
               max(occJSDM:::blockRhatESS(meas, 12L, 200L, 3L)$rhat))
  E <- vapply(1:3, function(ch) x$X[, x$is_emb] %*% apply(x$B_output[x$is_emb, , , ch], c(1, 2), mean),
              matrix(0, 30L, 6L))
  ref <- occJSDM:::grouped_effect_checks(E)
  expect_equal(got$effect$agreement_pairs, ref$agreement_pairs)
  expect_equal(got$effect$disagreement_pairs, ref$disagreement_pairs)
})

test_that("the checks do not copy the coefficient draws", {
  skip_if_not(capabilities("profmem"))
  # 2 measured and 30 embedding columns, 40 species, 200 draws, 2 chains: the
  # embedding rows of one chain are 1.9 MB and the measured rows 0.26 MB. No
  # single allocation may come near either.
  set.seed(11); p <- 32L; S <- 40L; niter <- 200L; nchain <- 2L; is_emb <- rep(c(FALSE, TRUE), c(2L, 30L))
  B <- array(rnorm(p * S * niter * nchain), c(p, S, niter, nchain))
  sb <- matrix(exp(rnorm(niter * nchain, sd = 0.1)), niter, nchain); se <- sb * 0.3
  lam <- array(exp(rnorm(30L * niter * nchain, sd = 0.1)), c(30L, niter, nchain))
  X <- matrix(rnorm(50L * p), 50L, p)
  big <- function(expr) {
    f <- tempfile(); on.exit(unlink(f))
    utils::Rprofmem(f, threshold = 2e5); force(expr); utils::Rprofmem(NULL)
    sum(grepl("^[0-9]", readLines(f)))                       # allocations, not R's "new page" lines
  }
  expect_identical(big(occJSDM:::grouped_convergence(B, sb, se, lam, X, is_emb)), 0L)
  results <- list(jsdm_output = list(B0_output = array(rnorm(S * niter * nchain), c(S, niter, nchain)),
                                     B_output = B, sigmab_output = sb, sigmaemb_output = se, lambda_output = lam))
  expect_identical(big(suppressMessages(suppressWarnings(computeDiagnostics(results, embedding_rows = is_emb)))), 0L)
})

test_that("missing values in the draws are flagged, not an error", {
  x <- made_fit(2L, c(0.3, 0.3)); x$sigmaemb_output[5, 2] <- NaN
  got <- do.call(occJSDM:::grouped_convergence, x)
  row <- check_row(got, "state")
  expect_true(row$flagged); expect_true(is.na(row$value))
  expect_match(row$meaning, "missing", fixed = TRUE)
  expect_identical(got$chains$state, c("shared structure", "not known"))
  expect_identical(got$verdict, "review")
  expect_silent(suppressMessages(occJSDM:::print_grouped_convergence(got)))

  x <- made_fit(2L, c(0.3, 0.3)); x$B_output[4, 2, 7, 1] <- NaN        # an embedding coefficient
  got <- do.call(occJSDM:::grouped_convergence, x)
  expect_true(check_row(got, "gradients")$flagged)
  expect_match(check_row(got, "gradients")$meaning, "missing", fixed = TRUE)
  expect_false(grepl("small", check_row(got, "gradients")$meaning))
  expect_identical(got$verdict, "review")

  x <- made_fit(2L, c(0.3, 0.3)); x$B_output[1, 2, 7, 1] <- NaN        # a measured coefficient
  got <- do.call(occJSDM:::grouped_convergence, x)
  expect_true(check_row(got, "rhat_measured")$flagged)
  expect_match(check_row(got, "rhat_measured")$meaning, "missing", fixed = TRUE)

  x <- made_fit(2L, c(0.3, 0.3)); x$lambda_output[2, 9, 1] <- NaN      # an axis scale
  got <- do.call(occJSDM:::grouped_convergence, x)
  expect_true(check_row(got, "rhat_axis")$flagged)
  expect_match(check_row(got, "rhat_axis")$meaning, "missing", fixed = TRUE)
  expect_false(check_row(got, "rhat_scales")$flagged)
  expect_identical(got$verdict, "review")

  x <- made_fit(2L, c(0.3, 0.3)); x$sigmab_output[3, 1] <- NaN         # the measured group's scale
  got <- do.call(occJSDM:::grouped_convergence, x)
  expect_true(check_row(got, "rhat_scales")$flagged)
  expect_match(check_row(got, "rhat_scales")$meaning, "missing", fixed = TRUE)
})

test_that("the state's value and limit are on the same scale", {
  got <- gc_of(2L, c(0.3, 6))
  expect_identical(check_row(got, "state")$limit, 0)
  expect_identical(check_row(got, "state")$value > check_row(got, "state")$limit, check_row(got, "state")$flagged)
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
  expect_identical(nrow(got$checks), 7L)
  expect_identical(nrow(got$chains), 2L)
  plain <- gc_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L), embedding = FALSE)
  expect_error(returnGroupedConvergence(plain), "returnGroupedConvergence() needs a fit with embeddingVars", fixed = TRUE)
  # Nothing in it depends on the model type.
  occ <- gc_fit(simulate_grouped_fixture(model = "occupancy", p_emb = 6L, S = 8L))
  expect_identical(nrow(returnGroupedConvergence(occ)$checks), 7L)
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
  expect_true(has(ok, "axis scales (log): 0 of 3 over 1.10"))
  expect_true(has(ok, "not part of the verdict"))
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
  # What is printed is what the function returns for the same fit.
  check <- returnGroupedConvergence(fit)
  expect_true(has(grouped$messages, paste(sprintf("%.2f", check$chains$sigma_emb_median), collapse = ", ")))
  expect_true(has(grouped$messages, paste("Verdict      |", check$verdict)))
  rhat <- check$checks$value[check$checks$check == "rhat_scales"]
  expect_true(has(grouped$messages, sprintf("sigma_emb and sigma_b (log), max: %.2f", rhat)))
  plain <- said(runOccJSDM(dat, listParams = list(n_factors = 1L),
                           occCovariates = grouped_occ_covariates(dat), MCMCparams = FIXTURE_MCMC))
  expect_false(has(plain$messages, "Grouped prior: convergence checks"))
  expect_false(has(plain$messages, "B (embedding)"))
})

# ---- the step at the end of a fit ------------------------------------------

jsdm_of <- function(x) x[c("B_output", "sigmab_output", "sigmaemb_output", "lambda_output")]

test_that("the end-of-fit step prints the checks and warns only on review", {
  x <- made_fit(2L, c(0.3, 6))
  review <- said(got <- occJSDM:::report_grouped_convergence(jsdm_of(x), x$X, x$is_emb))
  expect_identical(got$verdict, "review")
  expect_true(has(review$messages, "Grouped prior: convergence checks"))
  expect_identical(review$warnings, "The grouped prior's convergence checks flag this fit for review: see returnGroupedConvergence().")
  x <- made_fit(3L, c(0.3, 0.3, 0.3))
  ok <- said(got <- occJSDM:::report_grouped_convergence(jsdm_of(x), x$X, x$is_emb))
  expect_identical(got$verdict, "ok")
  expect_identical(ok$warnings, character(0))
})

test_that("an error in the checks is a warning at the end of a fit, and the fit is returned", {
  x <- made_fit(2L, c(0.3, 0.3))
  testthat::local_mocked_bindings(grouped_convergence = function(...) stop("boom"), .package = "occJSDM")
  got <- "not run"
  out <- said(got <- occJSDM:::report_grouped_convergence(jsdm_of(x), x$X, x$is_emb))
  expect_null(got)
  expect_identical(length(out$warnings), 1L)
  expect_match(out$warnings, "could not be computed", fixed = TRUE)
  expect_match(out$warnings, "boom", fixed = TRUE)
  expect_match(out$warnings, "returnGroupedConvergence()", fixed = TRUE)
  dat <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)$data_list
  fit <- NULL
  whole <- said(fit <- runOccJSDM(dat, listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
                                  occCovariates = grouped_occ_covariates(dat), MCMCparams = FIXTURE_MCMC))
  expect_false(is.null(fit$results_output$jsdm_output$B_output))
  expect_true(has(whole$warnings, "could not be computed"))
})

test_that("computeDiagnostics() refuses embedding_rows of the wrong length", {
  results <- made_results()
  expect_error(computeDiagnostics(results, embedding_rows = c(FALSE, TRUE)),
               "embedding_rows must have one entry per row of B_output (5), not 2", fixed = TRUE)
})
