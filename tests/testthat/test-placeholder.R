test_that("package loads", {
  expect_true(isNamespaceLoaded("occJSDM") || requireNamespace("occJSDM", quietly = TRUE))
})
