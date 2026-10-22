# Distinct values identify the parameter, species, iteration and chain without
# running MCMC. Pooling chains, selecting a wrong slot or permuting axes changes
# the literal values and per-chain means checked below.
posterior_draws_fit <- function() {
  list(
    infos = list(speciesNames = c("OTU_1", "OTU_2", "OTU_3"),
                 primerNames = c(10L, 20L)),
    X_psi = matrix(0, 1, 2,
                   dimnames = list(NULL, c("X_psi.a", "X_psi.b"))),
    X_theta = matrix(0, 1, 2,
                     dimnames = list(NULL, c("(Intercept)", "X_theta"))),
    results_output = list(
      beta_theta_output = array(as.numeric(1:48), c(2, 3, 4, 2)),
      p_output = array(100 + 1:48, c(2, 3, 4, 2)),
      q_output = array(200 + 1:48, c(2, 3, 4, 2)),
      theta0_output = array(1000 + 1:24, c(3, 4, 2)),
      jsdm_output = list(
        B0_output = array(2000 + 1:24, c(3, 4, 2)),
        B_output = array(3000 + 1:48, c(2, 3, 4, 2))
      )
    )
  )
}

test_that("posterior draws retain each parameter's exact species and chain order", {
  fit <- posterior_draws_fit()
  cases <- list(
    list(parameter = "beta_theta", extent = 2L,
         first = c(1, 7, 13, 19), second = c(25, 31, 37, 43)),
    list(parameter = "beta_psi", extent = 2L,
         first = c(3001, 3007, 3013, 3019), second = c(3025, 3031, 3037, 3043)),
    list(parameter = "p", extent = 2L,
         first = c(101, 107, 113, 119), second = c(125, 131, 137, 143)),
    list(parameter = "q", extent = 2L,
         first = c(201, 207, 213, 219), second = c(225, 231, 237, 243)),
    list(parameter = "theta0", extent = 1L,
         first = c(1001, 1004, 1007, 1010), second = c(1013, 1016, 1019, 1022)),
    list(parameter = "beta0_psi", extent = 1L,
         first = c(2001, 2004, 2007, 2010), second = c(2013, 2016, 2019, 2022))
  )
  for (case in cases) {
    draws <- returnPosteriorDraws(fit, case$parameter)
    expect_equal(dim(draws), c(case$extent, 3L, 4L, 2L), info = case$parameter)
    expect_equal(as.numeric(draws[1, 1, , 1]), case$first, info = case$parameter)
    expect_equal(as.numeric(draws[1, 1, , 2]), case$second, info = case$parameter)
    expect_identical(dimnames(draws)[[2]], c("OTU_1", "OTU_2", "OTU_3"))
  }
})

test_that("named selections preserve dimensions and requested ordering", {
  fit <- posterior_draws_fit()
  draws <- returnPosteriorDraws(fit, "beta_theta", covariate = "X_theta",
                                species = c("OTU_3", "OTU_1"))
  expect_equal(dim(draws), c(1L, 2L, 4L, 2L))
  expect_identical(dimnames(draws)[[1]], "X_theta")
  expect_identical(dimnames(draws)[[2]], c("OTU_3", "OTU_1"))
  expect_equal(as.numeric(draws[1, 1, , 2]), c(30, 36, 42, 48))
  expect_equal(as.numeric(draws[1, 2, , 2]), c(26, 32, 38, 44))

  slopes <- returnPosteriorDraws(fit, "beta_psi",
                                 covariate = c("X_psi.b", "X_psi.a"),
                                 species = "OTU_2")
  expect_equal(as.numeric(slopes[, 1, 1, 1]), c(3004, 3003))
})

test_that("primer selections match identifiers rather than their positions", {
  fit <- posterior_draws_fit()
  numeric_id <- returnPosteriorDraws(fit, "p", primer = 20, species = "OTU_3")
  character_id <- returnPosteriorDraws(fit, "p", primer = "20", species = "OTU_3")
  expect_identical(numeric_id, character_id)
  expect_equal(dim(numeric_id), c(1L, 1L, 4L, 2L))
  expect_identical(dimnames(numeric_id)[[1]], "20")
  expect_equal(as.numeric(numeric_id[1, 1, , 2]), c(130, 136, 142, 148))
  expect_error(returnPosteriorDraws(fit, "p", primer = 2), "primer.*not found")
})

test_that("posterior draws allow per-chain summaries without pooling", {
  draws <- returnPosteriorDraws(posterior_draws_fit(), "beta_theta",
                                covariate = "X_theta", species = "OTU_3")
  means <- apply(draws, c(2, 4), mean)
  expect_equal(unname(means), matrix(c(15, 39), nrow = 1L))
})

test_that("single species and chain keep the iteration and chain axes", {
  fit <- posterior_draws_fit()
  fit$infos$speciesNames <- "OTU_3"
  fit$results_output$theta0_output <- array(c(1003, 1006, 1009, 1012), c(1, 4, 1))
  draws <- returnPosteriorDraws(fit, "theta0", species = "OTU_3")
  expect_equal(dim(draws), c(1L, 1L, 4L, 1L))
  expect_equal(as.numeric(draws), c(1003, 1006, 1009, 1012))
  expect_identical(dimnames(draws)[[1]], "theta0")
  expect_identical(dimnames(draws)[[2]], "OTU_3")
})

test_that("unavailable or unknown parameters produce actionable errors", {
  fit <- posterior_draws_fit()
  fit$results_output$p_output <- NULL
  expect_error(returnPosteriorDraws(fit, "p"), "p.*no saved draws")
  expect_error(returnPosteriorDraws(fit, "unknown"), "parameter.*beta0_psi")
  expect_error(returnPosteriorDraws(fit, c("p", "q")), "parameter.*single")
  expect_error(returnPosteriorDraws(NULL, "p"), "fitmodel")
})

test_that("unknown labels and inappropriate selectors are rejected", {
  fit <- posterior_draws_fit()
  expect_error(returnPosteriorDraws(fit, "q", species = "absent"), "species.*not found")
  expect_error(returnPosteriorDraws(fit, "beta_psi", covariate = "absent"),
               "covariate.*not found")
  expect_error(returnPosteriorDraws(fit, "theta0", covariate = "X_theta"), "covariate")
  expect_error(returnPosteriorDraws(fit, "beta_theta", primer = "10"), "primer")
  expect_error(returnPosteriorDraws(fit, "p", covariate = "X_theta"), "covariate")
  expect_error(returnPosteriorDraws(fit, "p", species = character()), "species")
  expect_error(returnPosteriorDraws(fit, "beta_psi", covariate = NA_character_), "covariate")
})

test_that("traceplots use accessor labels and preserve all chain values", {
  draws <- returnPosteriorDraws(posterior_draws_fit(), "beta_theta",
                                covariate = "X_theta", species = "OTU_3")
  plot <- plotTraceplot(draws)
  expect_identical(unique(plot$data$label1), "X_theta")
  expect_identical(unique(plot$data$label2), "OTU_3")
  expect_equal(as.numeric(plot$data$value[plot$data$chain == "2"]), c(30, 36, 42, 48))
  expect_identical(unique(plot$data$iter), 1:4)
})

test_that("traceplot labels default to array names and allow explicit overrides", {
  draws <- array(as.numeric(1:12), c(3, 4, 1),
                 dimnames = list(c("species_a", "species_b", "species_c"), NULL, NULL))
  plot <- plotTraceplot(draws)
  expect_identical(unique(plot$data$label1), c("species_a", "species_b", "species_c"))
  explicit <- plotTraceplot(draws, dimnames1 = c("A", "B", "C"))
  expect_identical(unique(explicit$data$label1), c("A", "B", "C"))
})
