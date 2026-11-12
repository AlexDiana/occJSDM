test_that("b_prior = 'fixed' reproduces fits made before the option", {
  ref <- readRDS(test_path("fixtures", "fixed-prior-reference.rds"))
  skip_if_not(identical(La_library(), ref$la) && identical(extSoftVersion()[["BLAS"]], ref$blas),
    "BLAS or LAPACK differs from the reference")
  for (model in names(ref$fits)) {
    expect_identical(fixed_prior_reference_fit(model), ref$fits[[model]], info = model)
  }
})
