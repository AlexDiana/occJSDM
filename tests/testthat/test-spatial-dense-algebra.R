# Spatial slots carry 1-based coefficient indices. Dense BLAS arithmetic
# must preserve that interpretation for full, permuted and sparse layouts.
spatial_algebra_cases <- function() {
  Ks <- matrix(seq(-1.2, 1.1, length.out = 24), 6, 4)
  canonical <- matrix(rep(1:4, each = 6), 6, 4)
  permutation <- c(3, 1, 4, 2)
  row_permuted <- canonical
  row_permuted[2, ] <- permutation
  sparse <- matrix(c(1, 3, 2, 4, 1, 2, 4, 2, 3, 1, 3, 4), 6, 2)
  list(
    full = list(Ks = Ks, centers = canonical, knots = 4L),
    permuted = list(Ks = Ks[, permutation], centers = canonical[, permutation], knots = 4L),
    row_permuted = list(Ks = Ks, centers = row_permuted, knots = 4L),
    sparse = list(Ks = Ks[, 1:2], centers = sparse, knots = 4L),
    one_knot = list(Ks = Ks[, 1, drop = FALSE], centers = matrix(1, 6, 1), knots = 1L),
    no_knots = list(Ks = matrix(0, 6, 0), centers = matrix(0, 6, 0), knots = 0L)
  )
}

test_that("spatial precision and contributions match an explicit design matrix", {
  X <- cbind(1, c(-1, .2, .7, -2, .5, 1))
  Omega <- c(.4, 0, .9, 1.7, .2, 1.1)
  for (x in spatial_algebra_cases()) {
    # Independently expand each site's indexed slots into coefficient columns.
    H <- matrix(0, nrow(X), x$knots)
    for (i in seq_len(nrow(H))) for (j in seq_len(ncol(x$centers))) {
      H[i, x$centers[i, j]] <- H[i, x$centers[i, j]] + x$Ks[i, j]
    }
    design <- cbind(X, H)
    expected_precision <- crossprod(design, Omega * design)
    actual_precision <- XtOmegaX_SoR(X, x$knots, Omega, x$centers, x$Ks)
    expect_equal(actual_precision, unname(expected_precision), tolerance = 1e-12)
    B <- matrix(seq_len(x$knots * 3) / 7, x$knots, 3)
    expect_equal(KsBproduct(x$Ks, B, x$centers), H %*% B, tolerance = 1e-12)
    # No environmental columns is another live spatial-model configuration.
    expect_equal(XtOmegaX_SoR(matrix(0, 6, 0), x$knots, Omega, x$centers, x$Ks),
                 crossprod(H, Omega * H), tolerance = 1e-12)
  }
})

test_that("spatial conditional mean shifts match the Gaussian linear term", {
  X <- cbind(1, c(-1, .2, .7, -2, .5, 1))
  Omega <- c(.4, 0, .9, 1.7, .2, 1.1)
  k <- c(.7, -1, .4, 1.2, -.3, .6)
  for (x in spatial_algebra_cases()) {
    H <- matrix(0, nrow(X), x$knots)
    for (i in seq_len(nrow(H))) for (j in seq_len(ncol(x$centers))) {
      H[i, x$centers[i, j]] <- H[i, x$centers[i, j]] + x$Ks[i, j]
    }
    design <- cbind(X, H)
    prior_precision <- diag(seq(.5, 1.5, length.out = ncol(design)))
    b <- seq(-.4, .2, length.out = ncol(design))
    precision <- crossprod(design, Omega * design) + prior_precision
    expected <- solve(precision, crossprod(design, k))
    # Identical innovations cancel, isolating XtK_SoR through its real caller.
    setOccJSDMSeed(318L)
    shifted <- sampleB_SoR(X, prior_precision, b, k, Omega,
                           x$centers, x$Ks, x$knots)
    setOccJSDMSeed(318L)
    zero <- sampleB_SoR(X, prior_precision, b, rep(0, length(k)), Omega,
                        x$centers, x$Ks, x$knots)
    expect_equal(as.vector(shifted - zero), as.vector(expected), tolerance = 1e-12)
  }
})

test_that("dense and permuted representations produce the same conditional draws", {
  cases <- spatial_algebra_cases()
  X <- cbind(1, c(-1, .2, .7, -2, .5, 1))
  Omega <- c(.4, 0, .9, 1.7, .2, 1.1)
  k <- c(.7, -1, .4, 1.2, -.3, .6)
  sample <- function(x) {
    setOccJSDMSeed(179L)
    sampleB_SoR(X, diag(6), rep(.3, 6), k, Omega, x$centers, x$Ks, x$knots)
  }
  expect_equal(sample(cases$full), sample(cases$permuted), tolerance = 1e-12)
})
