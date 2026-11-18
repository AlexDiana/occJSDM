# The reference fits of the first test in test-grouped-prior.R: what the code
# gave at commit 4f2ee6c, before the grouped prior existed.
#
# Do NOT run this script at the branch head. It would save what the present
# code gives under the stamp "4f2ee6c", and the test would then compare the
# code with itself. To regenerate honestly (for another R or BLAS, say), load
# the package from a checkout of 4f2ee6c. From the package root:
#   git worktree add --detach <dir> 4f2ee6c
#   mkdir -p <dir>/tests/testthat/fixtures
#   cp tests/testthat/fixtures/make-fixed-prior-reference.R <dir>/tests/testthat/fixtures/
#   cp tests/testthat/helper-grouped-prior.R <dir>/tests/testthat/
#   (cd <dir> && Rscript tests/testthat/fixtures/make-fixed-prior-reference.R)
#   cp <dir>/tests/testthat/fixtures/fixed-prior-reference.rds tests/testthat/fixtures/
#   git worktree remove --force <dir>
# Neither copied file exists at 4f2ee6c; helper-fixtures.R is the same there
# as here. Run under the same R and BLAS the tests will use. Checked on
# 12 November 2026: fits made this way are identical to the saved reference.
devtools::load_all(quiet = TRUE)
invisible(testthat::source_test_helpers("tests/testthat"))
models <- c("binary", "continuous", "occupancy", "two_stage")
fits <- stats::setNames(lapply(models, fixed_prior_reference_fit), models)
saveRDS(list(fits = fits, la = La_library(), blas = extSoftVersion()[["BLAS"]], r = R.version.string, commit = "4f2ee6c"),
        "tests/testthat/fixtures/fixed-prior-reference.rds", compress = "xz")
