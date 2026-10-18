# The 100-species arms with B0's community mean (base_cm_s100 and
# base_rare_cm_s100). See the comment above those cells in helper-simstudy.R.
# No fits: simulation and fit arguments only.

sc_s100 <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

sim_s100 <- function(s, replicate = 1L) {
  simstudy_simulate(draw_truth(s, simstudy_seed_for(s, replicate)))
}

# Each 100-species arm and its 10-species counterpart.
s100_arms <- c(base_cm_s100 = "base_cm_m3", base_rare_cm_s100 = "base_rare_cm_m3")

test_that("the 100-species arms are these two", {
  expect_setequal(grep("_s100$", names(sc_s100()), value = TRUE), names(s100_arms))
})

test_that("each 100-species arm is its 10-species M = 3 counterpart with S = 100", {
  sc <- sc_s100()
  for (arm in names(s100_arms)) {
    a <- sc[[arm]]; p <- sc[[s100_arms[[arm]]]]
    expect_identical(a$S, 100L, label = arm)
    expect_identical(p$S, 10L, label = arm)
    expect_identical(a$seed_label, "s100", label = arm)
    expect_identical(a$fit_b0_prior, "community", label = arm)
    drop <- c("label", "S", "seed_label")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = arm)
  }
})

test_that("the two 100-species arms build 100-species datasets that differ only through B0", {
  sc <- sc_s100()
  o <- sc$base_cm_s100; r <- sc$base_rare_cm_s100
  expect_identical(r[setdiff(names(r), c("label", "b0_mean"))],
                   o[setdiff(names(o), c("label", "b0_mean"))])
  expect_null(o$b0_mean)
  expect_identical(r$b0_mean, -2.75)
  sim_o <- sim_s100(o); sim_r <- sim_s100(r)
  for (sim in list(sim_o, sim_r)) {
    expect_identical(ncol(sim$data_list$OTU), 100L)
    expect_identical(nrow(sim$data_list$traits), 100L)
    expect_length(sim$true_params$jsdmParams_true$B0, 100L)
    expect_identical(dim(sim$true_params$beta_theta_true), c(2L, 100L))
    # 3 field samples x 2 PCR primers x 3 replicates at each of 100 sites.
    expect_identical(as.vector(table(sim$data_list$info$Site)), rep(18L, 100L))
  }
  jo <- sim_o$true_params$jsdmParams_true; jr <- sim_r$true_params$jsdmParams_true
  expect_equal(jr$B0, jo$B0 - 2.75, tolerance = 1e-12)
  for (k in c("B", "A", "C", "U", "L", "spatField"))
    expect_identical(jr[[k]], jo[[k]], label = k)
  expect_identical(sim_r$true_params$beta_theta_true, sim_o$true_params$beta_theta_true)
  # The shift does what it is for: far fewer occupied sites.
  expect_lt(mean(sim_r$true_params$z_true), mean(sim_o$true_params$z_true) - 0.2)
})

test_that("the 100-species arms hand the fitter only b0_prior = 'community'", {
  sc <- sc_s100()
  for (arm in names(s100_arms)) {
    args <- simstudy_fit_args(sim_s100(sc[[arm]]), sc[[arm]])
    expect_identical(args$listPriors, list(b0_prior = "community"), label = arm)
    expect_identical(args$listParams$n_factors, 2L, label = arm)
    expect_identical(args$collCovariates, "X_theta", label = arm)
  }
})
