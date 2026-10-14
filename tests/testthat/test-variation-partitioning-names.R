test_that("returnVariationPartitioning returns one row of fractions per species", {
  vp <- returnVariationPartitioning(variance_partitioning_plot_fixture())
  expect_identical(nrow(vp), 2L)
})

test_that("the old variance-partitioning names warn and match the new ones", {
  fit <- variance_partitioning_plot_fixture()
  expect_warning(old <- returnVariancePartitioning(fit),
                 "returnVariationPartitioning")
  expect_identical(old, returnVariationPartitioning(fit))
  local_mocked_bindings(plotVarPart = function(varPart_output, speciesNames)
    speciesNames)
  expect_warning(old_plot <- plotVariancePartitioning(fit),
                 "plotVariationPartitioning")
  expect_identical(old_plot, plotVariationPartitioning(fit))
})
