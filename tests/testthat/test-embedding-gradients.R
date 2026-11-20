# Gradients of the embedding response and their comparison with covariates
# (R/embedding-gradients.R; design in
# dev/superpowers/specs/2026-11-20-embedding-gradients-design.md).

test_that("gradient_agreement() gives the pairs grouped_effect_checks() reports", {
  set.seed(1); n <- 40L; S <- 12L
  U <- qr.Q(qr(matrix(rnorm(n * 3), n, 3))); V <- matrix(rnorm(S * 3), S, 3)
  E1 <- U[, 1:2] %*% diag(c(9, 5)) %*% t(V[, 1:2])
  E2 <- U[, 2:3] %*% diag(c(9, 5)) %*% t(V[, 1:2])
  E <- array(c(E1, E2), c(n, S, 2L))
  got <- occJSDM:::gradient_agreement(E, 2L)
  expect_identical(dim(got), c(2L, 2L)); expect_true(all(is.na(diag(got))))
  expect_equal(got, occJSDM:::grouped_effect_checks(E)$agreement_pairs)
  expect_equal(occJSDM:::gradient_agreement(array(c(E1, E1), c(n, S, 2L)), 2L)[1, 2], 1)
})

test_that("embedding_effect() rebuilds each chain's effect from the chain means", {
  set.seed(2); n <- 10L; p <- 5L; S <- 3L; is_emb <- c(FALSE, TRUE, TRUE, FALSE, TRUE)
  X <- matrix(rnorm(n * p), n, p); B1 <- matrix(rnorm(p * S), p, S); B2 <- matrix(rnorm(p * S), p, S)
  got <- occJSDM:::embedding_effect(cbind(as.vector(B1), as.vector(B2)), X, is_emb, S)
  expect_identical(dim(got$B_mean), c(3L, S, 2L)); expect_identical(dim(got$E), c(n, S, 2L))
  expect_equal(got$B_mean[, , 2], B2[is_emb, ]); expect_equal(got$E[, , 1], X[, is_emb] %*% B1[is_emb, ])
})
