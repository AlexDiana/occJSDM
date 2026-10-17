# The sigma_b prior arms: no_spatial_sbhc, no_spatial_sbhc_long,
# no_spatial_sbig2, no_spatial_sbig2_long and no_spatial_sbhc_gt4_long, each
# paired with an arm already run on no_spatial's datasets. See the comment
# above those cells in helper-simstudy.R and
# dev/simstudy/sigma-b-prior-20261017/README.md.

sc_sb <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

sb_arms <- c("no_spatial_sbhc", "no_spatial_sbhc_long", "no_spatial_sbig2",
             "no_spatial_sbig2_long", "no_spatial_sbhc_gt4_long")
hc <- list(sigma_b_prior = "half_cauchy", sigma_b_scale = 1)
ig2 <- list(a_sigmab = 2, b_sigmab = 1 / 9)

test_that("only the sigma_b prior arms set a sigma_b prior, and each sets what its name says", {
  sc <- sc_sb()
  default <- list(type = "inverse_gamma", shape = 10, rate = 1)
  sb_names <- c("sigma_b_prior", "sigma_b_scale", "a_sigmab", "b_sigmab")
  # Other cells may carry other priors (the M ladder's tight collection
  # priors do); none may touch sigma_b's.
  for (s in sc) if (!s$label %in% sb_arms) {
    expect_false(any(sb_names %in% names(s$listPriors)), info = s$label)
    expect_identical(occJSDM:::read_sigma_b_prior(s$listPriors), default, info = s$label)
  }
  expect_identical(sc$no_spatial_sbhc$listPriors, hc)
  expect_identical(sc$no_spatial_sbhc_long$listPriors, hc)
  expect_identical(sc$no_spatial_sbhc_gt4_long$listPriors, hc)
  expect_identical(sc$no_spatial_sbig2$listPriors, ig2)
  expect_identical(sc$no_spatial_sbig2_long$listPriors, ig2)
  # The weak inverse-gamma keeps the default's prior mean of sigma_b^2, 1/9.
  expect_equal(ig2$b_sigmab / (ig2$a_sigmab - 1), 1 / 9)
  long <- sc$no_spatial_long$MCMCparams
  expect_identical(sc$no_spatial_sbhc_long$MCMCparams, long)
  expect_identical(sc$no_spatial_sbig2_long$MCMCparams, long)
  expect_identical(sc$no_spatial_sbhc_gt4_long$MCMCparams, long)
  expect_null(sc$no_spatial_sbhc$MCMCparams)
  expect_null(sc$no_spatial_sbig2$MCMCparams)
  expect_identical(sc$no_spatial_sbhc_gt4_long$fit_gt, sc$no_spatial_gt4_long$fit_gt)
  for (a in setdiff(sb_arms, "no_spatial_sbhc_gt4_long")) expect_null(sc[[a]]$fit_gt, info = a)
  for (a in sb_arms) {
    expect_false(isTRUE(sc[[a]]$useSpatField), info = a)
    expect_false(simstudy_fits_spatial(sc[[a]]), info = a)
    expect_identical(sc[[a]]$seed_label, "base", info = a)
  }
})

test_that("the sigma_b prior arms fit exactly their comparison arm's data and settings", {
  skip_on_cran()
  sc <- sc_sb()
  args <- function(lab) {
    s <- sc[[lab]]
    simstudy_fit_args(simstudy_simulate(draw_truth(s, simstudy_seed_for(s, 1L))), s)
  }
  a_ns <- args("no_spatial"); a_gt4 <- args("no_spatial_gt4_long")
  a <- lapply(stats::setNames(sb_arms, sb_arms), args)
  for (lab in sb_arms) {
    expect_identical(a[[lab]]$data, a_ns$data, info = lab)
    for (nm in c("occCovariates", "collCovariates", "spatCovariates"))
      expect_identical(a[[lab]][[nm]], a_ns[[nm]], info = paste(lab, nm))
  }
  for (lab in setdiff(sb_arms, "no_spatial_sbhc_gt4_long"))
    expect_identical(a[[lab]]$listParams, a_ns$listParams, info = lab)
  expect_identical(a$no_spatial_sbhc_gt4_long$listParams, a_gt4$listParams)
  expect_identical(a_gt4$listParams, list(n_factors = 2L, n_lattrait = 4L))
  expect_identical(a_ns$listPriors, list())
  # The priors reach runOccJSDM()'s reader as the arms intend.
  expect_identical(occJSDM:::read_sigma_b_prior(a$no_spatial_sbhc$listPriors),
                   list(type = "half_cauchy", scale = 1))
  expect_identical(occJSDM:::read_sigma_b_prior(a$no_spatial_sbig2$listPriors),
                   list(type = "inverse_gamma", shape = 2, rate = 1 / 9))
})
