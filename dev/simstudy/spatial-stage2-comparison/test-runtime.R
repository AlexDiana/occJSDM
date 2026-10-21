# Run: dev/simstudy/spatial-stage2-comparison/run-r -f dev/simstudy/spatial-stage2-comparison/test-runtime.R
library(testthat)
test_that("the pinned runtime loads every package the study needs", {
  expect_identical(R.version$major, "4"); expect_identical(R.version$minor, "5.0")
  for (p in c("occJSDM", "Hmsc", "gllvm", "posterior", "jsonlite", "testthat"))
    expect_true(requireNamespace(p, quietly = TRUE), info = p)
  expect_identical(as.character(packageVersion("occJSDM")), "0.1.0")
  expect_true(grepl("spatial-package-comparison-20261008/library", find.package("occJSDM")))
})
test_that("copied upstream sources match their recorded hashes", {
  s <- read.csv("dev/simstudy/spatial-stage2-comparison/sources.csv", stringsAsFactors = FALSE)
  expect_true(all(unname(tools::md5sum(file.path("dev/simstudy/spatial-stage2-comparison", s$file))) == s$md5))
  expect_equal(nrow(s), 11L)
})
cat("Runtime tests passed.\n")
