# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-generator.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/kernel.R")
source("dev/simstudy/spatial-recoverability-pilot/seeds.R")
source("dev/simstudy/spatial-recoverability-pilot/generator.R")

test_that("species table crosses three classes with four occupancies", {
  sp <- pilot_species()
  expect_equal(nrow(sp), 12L)
  expect_equal(unname(table(sp$class, sp$occupancy)), matrix(1L, 3, 4), ignore_attr = TRUE)
  expect_equal(sort(unique(exp(sp$log_range))), c(0.03, 0.08, 0.15))
  expect_equal(sp$slope, rep(c(0.4, -0.4), 6))
})

test_that("arrangements have the declared counts, separations and axis balance", {
  for (n in c(100L, 300L, 900L)) {
    a <- make_arrangements(n, seed = 11L)
    expect_equal(nrow(a$uniform), n); expect_equal(nrow(a$clustered), n)
    expect_true(all(a$uniform >= 0 & a$uniform <= 1)); expect_true(all(a$clustered >= 0 & a$clustered <= 1))
    expect_equal(nrow(a$centres), n / 10L)
    d <- as.matrix(dist(a$centres)); diag(d) <- Inf
    expect_true(min(d) >= a$separation - 1e-12)
    expect_equal(a$separation, 0.2 * sqrt(10 / (n / 10)))
    for (m in list(a$uniform, a$clustered)) expect_lt(abs(sd(m[, 1]) / sd(m[, 2]) - 1), 0.1 + 1e-12)
  }
})

test_that("truth 2 community hits its occupancy targets on the lattice", {
  cm <- make_community("t2", n = 100L, replicate = 1L)
  lat <- cm$index$lattice
  expect_equal(unname(colMeans(cm$psi[lat, ])), cm$species$occupancy, tolerance = 1e-6)
  expect_equal(nrow(cm$points), 100L + 100L + 1600L); expect_false(anyDuplicated(cm$points) > 0)
  expect_true(all(cm$z %in% 0:1))
  expect_equal(cm$spatial_sd, rep(1, 12))
  expect_equal(unname(exp(cm$species$log_range)), cm$ranges)
  expect_true(all(cm$conditioning$route %in% c("cholesky", "eigen")))
  expect_true(all(cm$conditioning$error <= 1e-5))
})

test_that("truth 1 community has unit spatial variance and the 90/10 split in expectation", {
  vars <- sapply(1:12, function(r) {
    cm <- make_community("t1", n = 100L, replicate = 1000L + r, range = 0.08)
    c(total = mean(apply(cm$spatial[cm$index$lattice, ], 2, var)),
      shared = mean(apply(cm$shared_part[cm$index$lattice, ], 2, var)))
  })
  expect_equal(mean(vars["total", ]), 1, tolerance = 0.25)
  expect_equal(mean(vars["shared", ]) / mean(vars["total", ]), 0.9, tolerance = 0.1)
  cm <- make_community("t1", 100L, 1L, range = 0.08)
  expect_equal(unname(rowSums(cm$loadings^2)), rep(1, 12))
  expect_equal(unname(colMeans(cm$psi[cm$index$lattice, ])), cm$species$occupancy, tolerance = 1e-6)
})

test_that("communities reproduce from their replicate and differ across replicates", {
  a <- make_community("t2", 100L, 3L); b <- make_community("t2", 100L, 3L); c3 <- make_community("t2", 100L, 4L)
  expect_identical(a$z, b$z); expect_false(identical(a$z, c3$z))
  expect_false(identical(a$arrangements$uniform, c3$arrangements$uniform))
})

test_that("the cell table has 9 truth 1 cells and 3 truth 2 cells", {
  cells <- pilot_cells(); expect_equal(nrow(cells), 12L)
  expect_equal(sum(cells$truth == "t1"), 9L); expect_true(all(is.na(cells$range[cells$truth == "t2"])))
  expect_false(anyDuplicated(cells$cell) > 0)
})
cat("Generator tests passed.\n")
