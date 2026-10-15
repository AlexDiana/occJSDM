# The weighted Gram matrix Z' diag(w) Z is symmetric, so it is formed with the
# symmetric BLAS routine (dsyrk, through crossprod(sqrt(w) * Z) in R and
# Zw.t() * Zw in Armadillo), which does half the arithmetic of the general
# product crossprod(Z, w * Z) and returns an exactly symmetric matrix.

test_that("weighted_gram() forms Z' diag(w) Z exactly symmetrically", {
  set.seed(2026101401)
  Z <- matrix(rnorm(60 * 7), 60, 7)
  w <- c(runif(58), 0, 0)
  Q <- weighted_gram(Z, w)
  expect_equal(Q, crossprod(Z, w * Z), tolerance = 1e-12)
  expect_identical(Q, t(Q))
  expect_identical(dim(weighted_gram(matrix(0, 60, 0), w)), c(0L, 0L))
  expect_equal(weighted_gram(Z, rep(0, 60)), matrix(0, 7, 7))
})

test_that("the C++ spatial precision uses the same symmetric routine as weighted_gram()", {
  # With the full basis layout XtOmegaX_SoR() scales the design by sqrt(Omega)
  # and calls dsyrk, exactly as weighted_gram() does, so the two agree bit for
  # bit; the general product it used before differed from this at round-off.
  set.seed(2026101402)
  n <- 50L; ps <- 6L
  X <- cbind(1, rnorm(n))
  Ks <- matrix(rnorm(n * ps), n, ps)
  centers <- matrix(rep(seq_len(ps), each = n), n, ps)
  Omega <- exp(rnorm(n))
  actual <- XtOmegaX_SoR(X, ps, Omega, centers, Ks)
  expect_identical(actual, weighted_gram(cbind(X, Ks), Omega))
  expect_identical(actual, t(actual))
})
