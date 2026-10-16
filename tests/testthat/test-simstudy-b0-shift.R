# The b0_mean setting of the simulation study harness and the rare-dominated
# arms that use it (base_rare, occupancy_rare, binary_rare). See the comment
# above those cells and simstudy_shifted_simulator() in helper-simstudy.R.
# No fits: simulation only, so these run in seconds.

sc_b0 <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

sim_at <- function(s, replicate = 1L) {
  simstudy_simulate(draw_truth(s, simstudy_seed_for(s, replicate)))
}

rare_parents <- c(base_rare = "base", occupancy_rare = "occupancy",
                  binary_rare = "binary")

test_that("cells without b0_mean simulate exactly as before", {
  sc <- sc_b0()
  # The rare arms' community-mean refits share their b0_mean
  # (test-simstudy-b0-community.R).
  shifted <- c(names(rare_parents), paste0(names(rare_parents), "_cm"))
  for (s in sc) if (!s$label %in% shifted)
    expect_null(s$b0_mean, label = s$label)

  # draw_truth() sets the seed, so calling it again replays the stream.
  truth <- draw_truth(sc$base, simstudy_seed_for(sc$base, 1L))
  expect_false("b0_mean" %in% names(truth))
  direct <- simulateOccJSDMData(truth$datasettings, truth$params,
                                truth$jsdmParams, model = truth$model)
  truth <- draw_truth(sc$base, simstudy_seed_for(sc$base, 1L))
  expect_identical(simstudy_simulate(truth), direct)
})

test_that("a zero shift reproduces the package simulator bit for bit", {
  # The shifted simulator is a copy of the package's own; with nothing to
  # add it must give the same output from the same seed, or the copy is wrong.
  sc <- sc_b0()
  for (lab in c("base", "binary")) {
    truth <- draw_truth(sc[[lab]], simstudy_seed_for(sc[[lab]], 1L))
    direct <- simulateOccJSDMData(truth$datasettings, truth$params,
                                  truth$jsdmParams, model = truth$model)
    truth <- draw_truth(sc[[lab]], simstudy_seed_for(sc[[lab]], 1L))
    copied <- simstudy_shifted_simulator(0)(truth$datasettings, truth$params,
                                            truth$jsdmParams, model = truth$model)
    expect_identical(copied, direct, label = lab)
  }
})

test_that("the shifted simulator leaves the package namespace alone", {
  before <- get("simulateData", envir = asNamespace("occJSDM"))
  f <- simstudy_shifted_simulator(-2.75)
  expect_identical(get("simulateData", envir = asNamespace("occJSDM")), before)
  expect_false(identical(environment(f), asNamespace("occJSDM")))
  expect_error(simstudy_shifted_simulator("a"))
  expect_error(simstudy_shifted_simulator(c(-1, -2)))
})

test_that("rare arms pair with their parents, differing only in B0", {
  sc <- sc_b0()
  for (rare in names(rare_parents)) {
    r <- sc[[rare]]; p <- sc[[rare_parents[[rare]]]]
    expect_identical(r$b0_mean, -2.75)
    expect_identical(simstudy_seed_for(r, 7L), simstudy_seed_for(p, 7L))
    # draw_truth() consumes the same RNG and records the shift, nothing else.
    tr <- draw_truth(r, simstudy_seed_for(r, 1L))
    tp <- draw_truth(p, simstudy_seed_for(p, 1L))
    expect_identical(tr$b0_mean, -2.75)
    tr$b0_mean <- NULL
    expect_identical(tr, tp, label = rare)

    sr <- sim_at(r); sp <- sim_at(p)
    jr <- sr$true_params$jsdmParams_true; jp <- sp$true_params$jsdmParams_true
    expect_identical(jr$B0, jp$B0 - 2.75, label = rare)
    expect_equal(jr$eta, jp$eta - 2.75, tolerance = 1e-12, label = rare)
    same <- setdiff(names(jp), c("B0", "eta", "varPart"))
    for (nm in same) expect_identical(jr[[nm]], jp[[nm]], label = paste(rare, nm))
    # Covariates, coordinates, collection design and traits match...
    expect_identical(sr$data_list$info, sp$data_list$info, label = rare)
    expect_identical(sr$data_list$traits, sp$data_list$traits, label = rare)
    for (nm in c("beta_theta_true", "p_true", "q_true"))
      expect_identical(sr$true_params[[nm]], sp$true_params[[nm]],
                       label = paste(rare, nm))
    # ...and the shift does what it is for: fewer occupied sites.
    z <- function(s) if (is.null(s$true_params$z_true)) s$data_list$OTU else
      s$true_params$z_true
    expect_lt(mean(z(sr)), mean(z(sp)))

    # Fitted exactly like the parent, apart from the data.
    ar <- simstudy_fit_args(sr, r); ap <- simstudy_fit_args(sp, p)
    for (nm in c("listParams", "listPriors", "occCovariates",
                 "collCovariates", "spatCovariates"))
      expect_identical(ar[[nm]], ap[[nm]], label = paste(rare, nm))
  }
})
