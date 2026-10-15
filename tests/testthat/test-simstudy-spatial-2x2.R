# The spatial 2x2 arms of the simulation study: base (field simulated and
# fitted), no_spatial (neither), field_unfitted and nofield_fitted. See the
# comment above those cells in helper-simstudy.R.

sc_2x2 <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

test_that("fit_spatial decides whether a cell's fit gets a spatial term", {
  sc <- sc_2x2()
  # Every cell that leaves fit_spatial unset keeps the historical rule.
  for (s in sc) if (is.null(s$fit_spatial))
    expect_identical(simstudy_fits_spatial(s), isTRUE(s$useSpatField),
                     info = s$label)
  expect_true(simstudy_fits_spatial(sc$base))
  expect_false(simstudy_fits_spatial(sc$no_spatial))
  # The two new arms decouple simulation from fitting.
  expect_true(isTRUE(sc$field_unfitted$useSpatField))
  expect_false(simstudy_fits_spatial(sc$field_unfitted))
  expect_false(isTRUE(sc$nofield_fitted$useSpatField))
  expect_true(simstudy_fits_spatial(sc$nofield_fitted))
})

test_that("the 2x2 arms hand the fitter a spatial term exactly when asked", {
  skip_on_cran()
  sc <- sc_2x2()
  args <- function(lab) {
    s <- sc[[lab]]
    simstudy_fit_args(simstudy_simulate(draw_truth(s, simstudy_seed_for(s, 1L))), s)
  }
  a_base <- args("base"); a_fu <- args("field_unfitted")
  a_nf <- args("nofield_fitted"); a_ns <- args("no_spatial")

  expect_identical(a_base$spatCovariates, c("Xs.1", "Xs.2"))
  expect_identical(a_base$listParams$n_supportpoints, 20L)
  expect_null(a_fu$spatCovariates)
  expect_null(a_fu$listParams$n_supportpoints)
  expect_identical(a_nf$spatCovariates, c("Xs.1", "Xs.2"))
  expect_identical(a_nf$listParams$n_supportpoints, 20L)
  expect_null(a_ns$spatCovariates)
  # The coordinates the no-field fit is pointed at exist in its data.
  expect_true(all(c("Xs.1", "Xs.2") %in% names(a_nf$data$info)))
  # Apart from the spatial term, paired arms pass the fitter the same thing.
  for (nm in c("data", "occCovariates", "collCovariates", "listPriors")) {
    expect_identical(a_fu[[nm]], a_base[[nm]], info = nm)
    expect_identical(a_nf[[nm]], a_ns[[nm]], info = nm)
  }
  expect_identical(a_fu$listParams[c("n_factors", "n_lattrait")],
                   a_base$listParams[c("n_factors", "n_lattrait")])
})

test_that("the 2x2 arms pair as documented", {
  skip_on_cran()
  # If this fails, the comparison between arms is no longer paired the way
  # the cells' comment says, and their report would over-claim.
  sc <- sc_2x2()
  sim <- function(lab) {
    s <- sc[[lab]]
    simstudy_simulate(draw_truth(s, simstudy_seed_for(s, 1L)))
  }
  b <- sim("base"); fu <- sim("field_unfitted")
  nf <- sim("nofield_fitted"); ns <- sim("no_spatial")

  # Same simulation settings and seed: identical datasets.
  expect_identical(fu$data_list, b$data_list)
  expect_identical(nf$data_list, ns$data_list)
  # Across the field and no-field pairs the whole JSDM layer still matches...
  for (nm in c("B0", "B", "G", "A", "C", "L"))
    expect_identical(b$true_params$jsdmParams_true[[nm]],
                     ns$true_params$jsdmParams_true[[nm]], info = nm)
  # ...but the field changes the latent occupancy, so the data differ.
  expect_false(identical(b$true_params$z_true, ns$true_params$z_true))
})
