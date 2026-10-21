# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-seeds.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/seeds.R")
source("dev/simstudy/spatial-design-sweep/generator.R")   # for sweep_seed()

test_that("seeds are positive integers and distinct across labels", {
  labels <- c("sites-t1-r0.03-n100-rep1", "sites-t1-r0.03-n100-rep2", "data-t2-n900-rep10",
              "oracle-t2-n300-rep4-clustered-species07")
  s <- vapply(labels, pilot_seed, integer(1))
  expect_true(all(s > 0L)); expect_false(anyDuplicated(s) > 0)
})

test_that("pilot seeds do not collide with the sweep's seed namespace", {
  sweep <- unlist(lapply(c("spatial-design-sites", "spatial-design-data", "spatial-design-detection",
                           "spatial-design-survey-spread", "spatial-design-survey-clustered"),
                         function(l) vapply(1:3, function(r) sweep_seed(l, r), integer(1))))
  pilot <- vapply(sprintf("sites-t1-r0.03-n100-rep%d", 1:10), pilot_seed, integer(1))
  expect_length(intersect(sweep, pilot), 0L)
})

test_that("pilot_label joins parts", {
  expect_identical(pilot_label("sites", "t2", "n100", "rep3"), "sites-t2-n100-rep3")
})
cat("Seed tests passed.\n")
