# The rare-dominated arms at 3 field samples per site (base_rare_m3,
# occupancy_rare_m3 and their _cm_m3 and _cm_both_m3 refits). See the comment
# above those cells in helper-simstudy.R. No fits: simulation and fit
# arguments only.

sc_m3 <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

# Each M = 3 arm, the M = 2 arm with the same priors, and the priors it adds.
m3_arms <- list(
  base_rare_m3 = list("base_rare", list()),
  occupancy_rare_m3 = list("occupancy_rare", list()),
  base_rare_cm_m3 = list("base_rare_cm", list(b0_prior = "community")),
  occupancy_rare_cm_m3 = list("occupancy_rare_cm", list(b0_prior = "community")),
  base_rare_cm_both_m3 = list("base_rare_cm_both",
                              list(b0_prior = "community", theta_intercept_prior = "community")),
  occupancy_rare_cm_both_m3 = list("occupancy_rare_cm_both",
                                   list(b0_prior = "community", theta_intercept_prior = "community")))

test_that("each M = 3 arm is its M = 2 counterpart with 3 field samples per site", {
  sc <- sc_m3()
  expect_setequal(grep("_m3$", names(sc), value = TRUE), names(m3_arms))
  for (arm in names(m3_arms)) {
    a <- sc[[arm]]; p <- sc[[m3_arms[[arm]][[1]]]]
    expect_identical(a$M, 3L, label = arm)
    expect_identical(p$M, 2L, label = arm)
    expect_identical(simstudy_seed_for(a, 7L), simstudy_seed_for(p, 7L), label = arm)
    drop <- c("label", "M")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = arm)
  }
})

test_that("the three priors at M = 3 fit identical datasets", {
  sc <- sc_m3()
  for (stem in c("base_rare", "occupancy_rare")) {
    arms <- paste0(stem, c("_m3", "_cm_m3", "_cm_both_m3"))
    sims <- lapply(arms, function(k) simstudy_simulate(draw_truth(sc[[k]], simstudy_seed_for(sc[[k]], 1L))))
    expect_identical(sims[[2]], sims[[1]], label = arms[2])
    expect_identical(sims[[3]], sims[[1]], label = arms[3])
    args <- Map(function(k, sim) simstudy_fit_args(sim, sc[[k]]), arms, sims)
    for (i in seq_along(arms)) {
      expect_identical(args[[i]]$listPriors, m3_arms[[arms[i]]][[2]], label = arms[i])
      expect_identical(args[[i]][names(args[[i]]) != "listPriors"],
                       args[[1]][names(args[[1]]) != "listPriors"], label = arms[i])
    }
  }
})

test_that("M = 3 datasets share the M = 2 arms' occupancy truth but not their field samples", {
  sc <- sc_m3()
  for (stem in c("base_rare", "occupancy_rare")) {
    a <- sc[[paste0(stem, "_m3")]]; p <- sc[[stem]]
    sim_a <- simstudy_simulate(draw_truth(a, simstudy_seed_for(a, 1L)))
    sim_p <- simstudy_simulate(draw_truth(p, simstudy_seed_for(p, 1L)))
    ja <- sim_a$true_params$jsdmParams_true; jp <- sim_p$true_params$jsdmParams_true
    expect_identical(ja$B0, jp$B0, label = stem)
    expect_identical(ja$eta, jp$eta, label = stem)
    expect_identical(sim_a$true_params$z_true, sim_p$true_params$z_true, label = stem)
    expect_identical(sim_a$true_params$beta_theta_true[1, ], sim_p$true_params$beta_theta_true[1, ],
                     label = stem)
    # Three field samples per site where the M = 2 arm has two.
    expect_identical(nrow(sim_a$data_list$OTU) / nrow(sim_p$data_list$OTU), 1.5, label = stem)
    expect_identical(as.vector(table(sim_a$data_list$info$Site)) / as.vector(table(sim_p$data_list$info$Site)),
                     rep(1.5, a$n), label = stem)
  }
})
