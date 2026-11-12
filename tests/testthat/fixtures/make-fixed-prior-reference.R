# Regenerate with: Rscript tests/testthat/fixtures/make-fixed-prior-reference.R
# Run from the package root, under the same R and BLAS the tests will use.
# Built from the code at commit 4f2ee6c, before the grouped prior existed.
devtools::load_all(quiet = TRUE)
invisible(testthat::source_test_helpers("tests/testthat"))
models <- c("binary", "continuous", "occupancy", "two_stage")
fits <- stats::setNames(lapply(models, fixed_prior_reference_fit), models)
saveRDS(list(fits = fits, la = La_library(), r = R.version.string, commit = "4f2ee6c"),
        "tests/testthat/fixtures/fixed-prior-reference.rds", compress = "xz")
