# Stage 1 of the count model: likelihood pieces tested alone. Nothing here is
# called by the sampler loop yet. Design:
# dev/superpowers/specs/2026-11-12-count-model-design.md.

test_that("the count model has one name, counts", {
  eta <- matrix(c(0, log(2), log(3), 0, 0, 0), 3, 2)
  expect_equal(occJSDM:::SD(eta, "counts"), c(stats::sd(c(1, 2, 3)), 0))
  expect_error(occJSDM:::SD(eta, "count"), "Using unavailable model")
})
