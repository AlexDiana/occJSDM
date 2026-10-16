# The fit_b0_prior setting of the simulation study harness and the six arms
# that use it (base_cm, occupancy_cm, binary_cm and their rare counterparts).
# See the comment above those cells in helper-simstudy.R. No fits: simulation
# and fit arguments only, so these run in seconds.

sc_cm <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

cm_counterparts <- c(base_cm = "base", occupancy_cm = "occupancy",
                     binary_cm = "binary", base_rare_cm = "base_rare",
                     occupancy_rare_cm = "occupancy_rare",
                     binary_rare_cm = "binary_rare")

test_that("cells without fit_b0_prior hand the fitter the priors they always did", {
  sc <- sc_cm()
  sim <- simstudy_simulate(draw_truth(sc$base, simstudy_seed_for(sc$base, 1L)))
  # The both-means arms are tested in test-simstudy-theta-community.R and the
  # M = 3 arms in test-simstudy-rare-m3.R and test-simstudy-ordinary-m3.R.
  for (s in sc) if (!s$label %in% names(cm_counterparts) &&
                    is.null(s$fit_theta_intercept_prior) &&
                    !grepl("_m3$", s$label)) {
    expect_null(s$fit_b0_prior, label = s$label)
    expect_identical(simstudy_fit_args(sim, s)$listPriors, s$listPriors %||% list(),
                     label = s$label)
  }
})

test_that("each community arm fits its counterpart's data with only b0_prior changed", {
  sc <- sc_cm()
  for (arm in names(cm_counterparts)) {
    a <- sc[[arm]]; p <- sc[[cm_counterparts[[arm]]]]
    expect_identical(a$fit_b0_prior, "community", label = arm)
    expect_identical(simstudy_seed_for(a, 7L), simstudy_seed_for(p, 7L), label = arm)
    # Everything in the specification except the label, the flag and the seed
    # label (checked through the seeds above) agrees.
    drop <- c("label", "fit_b0_prior", "seed_label")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = arm)
    sim_a <- simstudy_simulate(draw_truth(a, simstudy_seed_for(a, 1L)))
    sim_p <- simstudy_simulate(draw_truth(p, simstudy_seed_for(p, 1L)))
    expect_identical(sim_a, sim_p, label = arm)
    args_a <- simstudy_fit_args(sim_a, a); args_p <- simstudy_fit_args(sim_p, p)
    expect_identical(args_a$listPriors, c(args_p$listPriors, list(b0_prior = "community")),
                     label = arm)
    expect_identical(args_a[names(args_a) != "listPriors"],
                     args_p[names(args_p) != "listPriors"], label = arm)
  }
})

test_that("fit_b0_prior is added to a cell's own listPriors, not in place of it", {
  s <- utils::modifyList(SIMSTUDY_BASE, list(label = "x", fit_b0_prior = "community",
                                             listPriors = list(b_betatheta_slope_var = .5)))
  sim <- simstudy_simulate(draw_truth(sc_cm()$base, simstudy_seed_for(sc_cm()$base, 1L)))
  expect_identical(simstudy_fit_args(sim, s)$listPriors,
                   list(b_betatheta_slope_var = .5, b0_prior = "community"))
})
