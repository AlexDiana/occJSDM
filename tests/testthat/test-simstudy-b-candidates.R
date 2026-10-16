# The non-spatial candidates for B's undercoverage: no_spatial_long,
# no_factors, no_factors_fit2, no_spatial_notraits, no_spatial_gt4 and
# no_spatial_gt4_long, each paired with no_spatial. See the comment above
# those cells in helper-simstudy.R.

sc_bc <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

bc_arms <- c("no_spatial_long", "no_factors", "no_factors_fit2",
             "no_spatial_notraits", "no_spatial_gt4", "no_spatial_gt4_long")

test_that("only the new arms set fit_gt or their own MCMC settings", {
  sc <- sc_bc()
  for (s in sc) if (!s$label %in% bc_arms) {
    expect_null(s$fit_gt, info = s$label)
    expect_null(s$MCMCparams, info = s$label)
  }
  expect_identical(sc$no_spatial_long$MCMCparams,
                   list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1))
  expect_identical(sc$no_spatial_gt4_long$MCMCparams, sc$no_spatial_long$MCMCparams)
  expect_identical(sc$no_spatial_gt4_long$fit_gt, sc$no_spatial_gt4$fit_gt)
  for (a in bc_arms) {
    expect_false(isTRUE(sc[[a]]$useSpatField), info = a)
    expect_false(simstudy_fits_spatial(sc[[a]]), info = a)
    expect_identical(sc[[a]]$seed_label, "base", info = a)
  }
})

test_that("a cell's own MCMC settings replace the caller's", {
  sc <- sc_bc()
  seen <- NULL
  # The helper is sourced, not part of the package, so swap the fit by hand.
  env <- environment(simstudy_replicate)
  real_fit <- get("simstudy_fit", envir = env)
  assign("simstudy_fit", function(sim, scenario, MCMCparams = NULL) {
    seen <<- MCMCparams
    stop("stub fit")
  }, envir = env)
  withr::defer(assign("simstudy_fit", real_fit, envir = env))
  short <- list(nchain = 2, nburn = 10, niter = 10, nthin = 1)
  expect_error(simstudy_replicate(sc$no_spatial_long, 1L, short), "stub fit")
  expect_identical(seen, sc$no_spatial_long$MCMCparams)
  expect_error(simstudy_replicate(sc$no_spatial, 1L, short), "stub fit")
  expect_identical(seen, short)
})

test_that("the new arms hand the fitter what their comment says", {
  skip_on_cran()
  sc <- sc_bc()
  args <- function(lab) {
    s <- sc[[lab]]
    simstudy_fit_args(simstudy_simulate(draw_truth(s, simstudy_seed_for(s, 1L))), s)
  }
  a_ns <- args("no_spatial")
  a <- lapply(stats::setNames(bc_arms, bc_arms), args)

  expect_identical(a_ns$listParams, list(n_factors = 2L, n_lattrait = 2L))
  expect_identical(a$no_spatial_long$listParams, a_ns$listParams)
  expect_identical(a$no_factors$listParams, list(n_factors = 0L, n_lattrait = 2L))
  expect_identical(a$no_factors_fit2$listParams, list(n_factors = 2L, n_lattrait = 2L))
  expect_identical(a$no_spatial_notraits$listParams, a_ns$listParams)
  expect_identical(a$no_spatial_gt4$listParams, list(n_factors = 2L, n_lattrait = 4L))
  expect_identical(a$no_spatial_gt4_long$listParams, list(n_factors = 2L, n_lattrait = 4L))

  # Traits are dropped only in the arm that says so.
  expect_null(a$no_spatial_notraits$data$traits)
  expect_false(is.null(a_ns$data$traits))
  expect_identical(a$no_spatial_notraits$data,
                   a_ns$data[setdiff(names(a_ns$data), "traits")])
  # The fit-only arms fit no_spatial's datasets; the two d = 0 arms share theirs.
  expect_identical(a$no_spatial_long$data, a_ns$data)
  expect_identical(a$no_spatial_gt4$data, a_ns$data)
  expect_identical(a$no_spatial_gt4_long$data, a_ns$data)
  expect_identical(a$no_factors_fit2$data, a$no_factors$data)
  for (lab in bc_arms) for (nm in c("occCovariates", "collCovariates",
                                    "spatCovariates", "listPriors"))
    expect_identical(a[[lab]][[nm]], a_ns[[nm]], info = paste(lab, nm))
})

test_that("the new arms pair with no_spatial as documented", {
  skip_on_cran()
  # If this fails, the comparison between arms is no longer paired the way
  # the cells' comment says, and their report would over-claim.
  sc <- sc_bc()
  sim <- function(lab) {
    s <- sc[[lab]]
    simstudy_simulate(draw_truth(s, simstudy_seed_for(s, 1L)))
  }
  ns <- sim("no_spatial"); d0 <- sim("no_factors")

  # B, B0 and the trait layer are no_spatial's in every arm, d = 0 included.
  for (lab in bc_arms) {
    x <- sim(lab)
    for (nm in c("B0", "B", "G", "A", "C"))
      expect_identical(x$true_params$jsdmParams_true[[nm]],
                       ns$true_params$jsdmParams_true[[nm]], info = paste(lab, nm))
  }
  # The truth has no factors in the d = 0 arms, and their data differ.
  expect_identical(dim(d0$true_params$jsdmParams_true$L), c(0L, 10L))
  expect_false(identical(d0$data_list, ns$data_list))
  expect_identical(sim("no_factors_fit2")$data_list, d0$data_list)
})
