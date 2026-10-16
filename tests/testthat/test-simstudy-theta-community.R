# The fit_theta_intercept_prior setting of the simulation study harness and
# the six arms that use it with fit_b0_prior (base_cm_both, occupancy_cm_both,
# their rare counterparts and the two _wide arms). See the comment above those
# cells in helper-simstudy.R. No fits: simulation and fit arguments only.

sc_both <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

# Each arm against its default-prior counterpart, and the extra listPriors it adds.
both_counterparts <- list(
  base_cm_both = list("base", list()),
  occupancy_cm_both = list("occupancy", list()),
  base_rare_cm_both = list("base_rare", list()),
  occupancy_rare_cm_both = list("occupancy_rare", list()),
  base_rare_cm_both_wide = list("base_rare", list(mu_b0_sd = 5, mu_theta_sd = 5)),
  occupancy_rare_cm_both_wide = list("occupancy_rare", list(mu_b0_sd = 5, mu_theta_sd = 5)))

test_that("only the both-means arms set fit_theta_intercept_prior", {
  sc <- sc_both()
  setting <- names(sc)[!vapply(sc, function(s) is.null(s$fit_theta_intercept_prior), logical(1))]
  expect_setequal(setting, names(both_counterparts))
})

test_that("each both-means arm fits its counterpart's data with only the priors changed", {
  sc <- sc_both()
  for (arm in names(both_counterparts)) {
    a <- sc[[arm]]; p <- sc[[both_counterparts[[arm]][[1]]]]
    extra <- both_counterparts[[arm]][[2]]
    expect_identical(a$fit_b0_prior, "community", label = arm)
    expect_identical(a$fit_theta_intercept_prior, "community", label = arm)
    expect_identical(simstudy_seed_for(a, 7L), simstudy_seed_for(p, 7L), label = arm)
    # The B0-only arm of the same counterpart has the same seeds too.
    b0_only <- sc[[paste0(both_counterparts[[arm]][[1]], "_cm")]]
    expect_identical(simstudy_seed_for(a, 7L), simstudy_seed_for(b0_only, 7L), label = arm)
    # Everything in the specification except the label, the flags, the seed
    # label (checked through the seeds above) and the hyperpriors agrees.
    drop <- c("label", "fit_b0_prior", "fit_theta_intercept_prior", "seed_label", "listPriors")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = arm)
    expect_null(p$listPriors, label = arm)
    sim_a <- simstudy_simulate(draw_truth(a, simstudy_seed_for(a, 1L)))
    sim_p <- simstudy_simulate(draw_truth(p, simstudy_seed_for(p, 1L)))
    expect_identical(sim_a, sim_p, label = arm)
    args_a <- simstudy_fit_args(sim_a, a); args_p <- simstudy_fit_args(sim_p, p)
    expect_identical(args_a$listPriors,
                     c(extra, list(b0_prior = "community", theta_intercept_prior = "community")),
                     label = arm)
    expect_identical(args_a[names(args_a) != "listPriors"],
                     args_p[names(args_p) != "listPriors"], label = arm)
  }
})

test_that("fit_theta_intercept_prior is added to a cell's own listPriors", {
  s <- utils::modifyList(SIMSTUDY_BASE, list(label = "x", fit_theta_intercept_prior = "community",
                                             listPriors = list(b_betatheta_slope_var = .5)))
  sim <- simstudy_simulate(draw_truth(sc_both()$base, simstudy_seed_for(sc_both()$base, 1L)))
  expect_identical(simstudy_fit_args(sim, s)$listPriors,
                   list(b_betatheta_slope_var = .5, theta_intercept_prior = "community"))
})

test_that("mu_theta rows compare the draws with the realised mean collection intercept", {
  sc <- sc_both()
  s <- sc$base_cm_both
  truth <- draw_truth(s, simstudy_seed_for(s, 1L))
  sim <- simstudy_simulate(truth)
  # A stand-in fit: the scaling the harness reads, and draws for two blocks.
  niter <- 6L; nchain <- 2L; S <- s$S
  X <- sim$data_list$info$X_theta
  fit <- list(
    infos = list(list_X_theta_mat = list(names_df = "X_theta",
                                         mean_df = c(X_theta = mean(X)),
                                         sd_df = c(X_theta = stats::sd(X)))),
    X_theta = cbind("(Intercept)" = 1, X_theta = 0),
    results_output = list(
      jsdm_output = list(),
      beta_theta_output = array(0, c(2L, S, niter, nchain)),
      mu_theta_output = array(seq_len(niter * nchain) / 10, c(niter, nchain))))
  blocks <- simstudy_param_blocks(fit, sim, truth)
  bt <- blocks$beta_theta$truth
  expect_identical(blocks$mu_theta$truth, mean(bt[1, ]))
  # Without the option there is no mu_theta row.
  fit$results_output$mu_theta_output <- NULL
  expect_null(simstudy_param_blocks(fit, sim, truth)$mu_theta$post)
})
