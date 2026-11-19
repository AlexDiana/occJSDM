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
