# Run from the repository root: Rscript dev/simstudy/spatial-design-sweep/test-generator.R
library(testthat)
source("dev/simstudy/spatial-design-sweep/generator.R")

test_that("seeds are deterministic integers distinct across labels", {
  expect_identical(sweep_seed("spatial-design-data", 1L), sweep_seed("spatial-design-data", 1L))
  expect_true(sweep_seed("spatial-design-data", 1L) != sweep_seed("spatial-design-sites", 1L))
  expect_true(sweep_seed("spatial-design-data", 2L) == sweep_seed("spatial-design-data", 1L) + 1L)
})

test_that("arrangements have 100 sites each inside the unit square", {
  a <- make_arrangements(sweep_seed("spatial-design-sites", 1L))
  expect_identical(names(a), c("spread", "pairs", "clustered", "grid"))
  for (nm in names(a)) {
    expect_identical(dim(a[[nm]]), c(100L, 2L))
    expect_identical(colnames(a[[nm]]), c("x", "y"))
    expect_true(all(a[[nm]] >= 0 & a[[nm]] <= 1))
    expect_false(anyDuplicated(a[[nm]]) > 0)
    ratio <- sd(a[[nm]][, 1]) / sd(a[[nm]][, 2])
    expect_true(abs(ratio - 1) <= 0.1)
  }
})

test_that("pairs extend the first 80 spread sites with partners at distance 0.01", {
  a <- make_arrangements(sweep_seed("spatial-design-sites", 1L))
  expect_identical(a$pairs[1:80, ], a$spread[1:80, ])
  d <- as.matrix(dist(a$pairs))[81:100, 1:80]
  expect_equal(unname(apply(d, 1, min)), rep(0.01, 20), tolerance = 1e-12)
})

test_that("cluster centres are separated and points lie within 0.02 of a centre", {
  a <- make_arrangements(sweep_seed("spatial-design-sites", 1L))
  centres <- attr(a, "cluster_centres")
  expect_identical(dim(centres), c(10L, 2L))
  dc <- as.matrix(dist(centres)); diag(dc) <- Inf
  expect_true(min(dc) >= 0.2)
  for (k in 1:10) {
    pts <- a$clustered[(k - 1) * 10 + 1:10, ]
    d <- sqrt(colSums((t(pts) - centres[k, ])^2))
    expect_true(all(d <= 0.02 + 1e-12))
  }
})

test_that("the grid is the exact 10 by 10 lattice", {
  a <- make_arrangements(1L)
  g <- (1:10 - 0.5) / 10
  expect_equal(unname(a$grid), unname(as.matrix(expand.grid(x = g, y = g))))
})

test_that("the landscape has the right pieces and target prevalence on the lattice", {
  land <- make_landscape(1L)
  expect_identical(dim(land$points), c(1920L, 2L))
  expect_identical(land$index$lattice, 321:1920)
  expect_identical(land$index$pairs, c(1:80, 101:120))
  expect_identical(dim(land$field), c(1920L, 8L))
  expect_equal(unname(colMeans(land$psi[land$index$lattice, ])),
               c(.05, .05, .25, .25, .25, .75, .75, .75), tolerance = 1e-6)
  expect_equal(land$range, 0.03)
  expect_true(all(land$z %in% 0:1))
})

test_that("the generator is reproducible", {
  expect_identical(make_sweep_input(2L), make_sweep_input(2L))
})

test_that("surveys map rows to sites, samples and primers consistently", {
  input <- make_sweep_input(1L)
  for (arr in c("spread", "pairs", "clustered", "grid")) {
    s <- input$surveys[[arr]]
    expect_identical(dim(s$two_stage$OTU), c(2400L, 8L))
    expect_identical(nrow(s$two_stage$info), 2400L)
    expect_identical(unname(s$two_stage$info$Site), rep(1:100, each = 24))
    expect_identical(unname(s$two_stage$info$Sample), rep(1:200, each = 12))
    expect_identical(unname(s$two_stage$info$Primer), rep(rep(1:2, each = 6), times = 200))
    expect_identical(unname(s$binary$OTU), unname(input$landscape$z[input$landscape$index[[arr]], ]))
    expect_identical(colnames(s$binary$OTU), input$landscape$species)
    expect_equal(s$binary$info$longitude, input$landscape$points[input$landscape$index[[arr]], 1])
  }
})

test_that("design statistics are computed for every arrangement", {
  input <- make_sweep_input(1L)
  st <- input$statistics
  expect_identical(st$arrangement, c("spread", "pairs", "clustered", "grid"))
  expect_true(all(c("mean_nearest_neighbour", "fraction_with_half_neighbour",
                    "effective_rank", "standardised_range_x", "standardised_range_y",
                    "axis_sd_ratio", "environment_span", "in_grid") %in% names(st)))
  expect_true(st$fraction_with_half_neighbour[st$arrangement == "clustered"] >
              st$fraction_with_half_neighbour[st$arrangement == "spread"])
  expect_equal(st$mean_nearest_neighbour[st$arrangement == "grid"], 0.1)
})
cat("Generator tests passed.\n")
