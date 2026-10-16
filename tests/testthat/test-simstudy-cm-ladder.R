# The field-sample ladder with B0's community mean (M5_cm, M10_cm, M20_cm)
# and the rare-dominated ladder under both priors (base_rare_m5 to _m20 and
# base_rare_cm_m5 to _m20). See the comment above those cells in
# helper-simstudy.R. No fits: simulation and fit arguments only.

sc_ladder <- function() {
  sc <- simstudy_scenarios()
  names(sc) <- vapply(sc, `[[`, character(1), "label")
  sc
}

sim_ladder <- function(s, replicate = 1L) {
  simstudy_simulate(draw_truth(s, simstudy_seed_for(s, replicate)))
}

ladder_M <- c(5L, 10L, 20L)
ord_cm_arms <- setNames(paste0("M", ladder_M), paste0("M", ladder_M, "_cm"))
rare_arms <- paste0("base_rare_m", ladder_M)
rare_cm_arms <- setNames(rare_arms, paste0("base_rare_cm_m", ladder_M))

test_that("the ladder arms are these nine", {
  sc <- sc_ladder()
  expect_setequal(grep("^M[0-9]+_cm$|^base_rare(_cm)?_m([0-9]+)$", names(sc), value = TRUE),
                  c(names(ord_cm_arms), rare_arms, names(rare_cm_arms), "base_rare_m3", "base_rare_cm_m3"))
})

test_that("each ordinary community-mean arm refits its ladder arm's datasets", {
  sc <- sc_ladder()
  for (arm in names(ord_cm_arms)) {
    a <- sc[[arm]]; p <- sc[[ord_cm_arms[[arm]]]]
    expect_identical(a$fit_b0_prior, "community", label = arm)
    expect_identical(a$seed_label, "mladder", label = arm)
    drop <- c("label", "fit_b0_prior")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = arm)
    sim_a <- sim_ladder(a); sim_p <- sim_ladder(p)
    expect_identical(sim_a, sim_p, label = arm)
    args_a <- simstudy_fit_args(sim_a, a); args_p <- simstudy_fit_args(sim_p, p)
    expect_identical(args_p$listPriors, list(), label = arm)
    expect_identical(args_a$listPriors, list(b0_prior = "community"), label = arm)
    expect_identical(args_a[names(args_a) != "listPriors"],
                     args_p[names(args_p) != "listPriors"], label = arm)
  }
})

test_that("each rare ladder arm is base_rare_m3 with more field samples per site", {
  sc <- sc_ladder()
  for (i in seq_along(ladder_M)) {
    a <- sc[[rare_arms[i]]]; p <- sc$base_rare_m3
    expect_identical(a$M, ladder_M[i], label = rare_arms[i])
    expect_identical(simstudy_seed_for(a, 7L), simstudy_seed_for(sc$base_rare, 7L),
                     label = rare_arms[i])
    drop <- c("label", "M")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = rare_arms[i])
  }
})

test_that("rare ladder datasets share base_rare's occupancy truth but not its field samples", {
  sc <- sc_ladder()
  sim_p <- sim_ladder(sc$base_rare)
  jp <- sim_p$true_params$jsdmParams_true
  for (i in seq_along(ladder_M)) {
    arm <- rare_arms[i]
    sim_a <- sim_ladder(sc[[arm]])
    ja <- sim_a$true_params$jsdmParams_true
    expect_identical(ja$B0, jp$B0, label = arm)
    expect_identical(ja$eta, jp$eta, label = arm)
    expect_identical(sim_a$true_params$z_true, sim_p$true_params$z_true, label = arm)
    expect_identical(sim_a$true_params$beta_theta_true[1, ], sim_p$true_params$beta_theta_true[1, ],
                     label = arm)
    expect_identical(as.vector(table(sim_a$data_list$info$Site)) /
                       as.vector(table(sim_p$data_list$info$Site)),
                     rep(ladder_M[i] / 2, sc[[arm]]$n), label = arm)
  }
})

test_that("each rare community-mean arm refits its default arm's datasets", {
  sc <- sc_ladder()
  for (arm in names(rare_cm_arms)) {
    a <- sc[[arm]]; p <- sc[[rare_cm_arms[[arm]]]]
    expect_identical(a$fit_b0_prior, "community", label = arm)
    drop <- c("label", "fit_b0_prior")
    expect_identical(a[setdiff(names(a), drop)], p[setdiff(names(p), drop)], label = arm)
    sim_a <- sim_ladder(a); sim_p <- sim_ladder(p)
    expect_identical(sim_a, sim_p, label = arm)
    args_a <- simstudy_fit_args(sim_a, a); args_p <- simstudy_fit_args(sim_p, p)
    expect_identical(args_p$listPriors, list(), label = arm)
    expect_identical(args_a$listPriors, list(b0_prior = "community"), label = arm)
    expect_identical(args_a[names(args_a) != "listPriors"],
                     args_p[names(args_p) != "listPriors"], label = arm)
  }
})
