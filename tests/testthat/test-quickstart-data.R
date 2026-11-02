# The Quickstart's shipped example: a simulated river survey and its fit.

test_that("quickstartdata has the agreed survey design", {
  d <- quickstartdata
  expect_named(d, c("info", "OTU", "traits"))
  expect_equal(nrow(d$info), 1800L)
  expect_equal(length(unique(d$info$Site)), 100L)
  expect_true(all(table(unique(d$info[c("Site", "Sample")])$Site) == 3L))
  expect_equal(sort(unique(d$info$Primer)), 1:2)
  expect_true(all(table(interaction(d$info$Sample, d$info$Primer, drop = TRUE)) == 3L))
  # "not_protected" first, so the coefficient reads "protected compared with not protected".
  expect_identical(levels(d$info$protected), c("not_protected", "protected"))
  expect_true(all(c("elevation", "hmi", "water_volume") %in% names(d$info)))
  expect_equal(nrow(d$OTU), nrow(d$info))
  expect_true(ncol(d$OTU) >= 18L && ncol(d$OTU) <= 22L)
  expect_true(all(grepl("^[A-Z][a-z]+_[A-Z][a-z]+_[0-9]+$", colnames(d$OTU))))
  expect_identical(rownames(d$traits), colnames(d$OTU))
  expect_identical(colnames(d$traits), c("log10_body_mass", "red_listed"))
  expect_true(all(d$traits[, "red_listed"] %in% c(0, 1)))
  expect_true(any(d$traits[, "red_listed"] == 1) && any(d$traits[, "red_listed"] == 0))
})

test_that("quickstartresults is the agreed non-spatial fit with a handful of flags", {
  f <- quickstartresults
  expect_identical(f$infos$model, "two_stage")
  expect_equal(f$infos$n, 100L)
  expect_equal(f$infos$S, ncol(quickstartdata$OTU))
  expect_equal(f$infos$n_factors, 2)
  rda <- test_path("../../data/quickstartresults.rda")  # absent under R CMD check
  if (file.exists(rda)) expect_lt(file.size(rda), 25e6)
  diag <- returnConvergenceDiagnostics(f)
  n_flag <- sum(diag$rhat > 1.01 | diag$ess < 400)
  expect_gte(n_flag, 2L)
  expect_lte(n_flag, 10L)
})

test_that("every function the Quickstart calls runs on the slimmed quickstartresults", {
  f <- quickstartresults
  expect_no_error({
    computePredictiveOccupancyProbs(f); computeConditionalOccupancyProbs(f)
    returnOccupancyRates(f); returnOccupancyCovariates(f); returnTraitsCoeff(f)
    returnCollectionCovariates(f); returnResidualCorrelationMatrix(f)
    returnOrdinationScores(f); returnVariationPartitioning(f)
    returnCovariateEffect(f, covName = "hmi")
    returnPosteriorDraws(f, parameter = "theta0")
  })
  plots <- list(
    plotOccupancyRates(f), plotDetectionRates(f), plotCollectionRates(f),
    plotStage1FPRates(f), plotStage2FPRates(f), plotFPTPStage2Rates(f),
    plotOccupancyCovariates(f, covName = "hmi"),
    plotOccupancyGradient(f, covName = "elevation", idx_species = 1:2),
    plotCovariateEffect(f, covNames = "hmi", idx_species = 1:2)[[1]],  # a list, one plot per covariate
    plotCollectionCovariates(f, covName = "water_volume"),
    plotTraitsCoefficients(f, covName = "hmi"),
    plotResidualCorrelationMatrix(f), plotBiplot(f), plotOrdinationScores(f),
    plotCumulativeSpeciesDetections(f, M = 6, K = 6, primer = 1)
  )
  # ggplot2 4 plots are S7 objects, so expect_s3_class() does not apply.
  for (p in plots) expect_true(inherits(p, "ggplot"))
})

test_that("new-site predictions keep the factor's levels and show the planted protected effect", {
  f <- quickstartresults
  new_sites <- data.frame(
    elevation = rep(c(1000, 3000), 2),
    hmi = 0.2,
    protected = factor(rep(c("not_protected", "protected"), each = 2),
                       levels = levels(quickstartdata$info$protected))
  )
  prediction <- predictNewSites(f, X_psi = new_sites, useBiotic = FALSE, verbose = FALSE)
  red <- which(quickstartdata$traits[f$infos$speciesNames, "red_listed"] == 1)
  # Median occupancy: protected (rows 3, 4) above not protected (rows 1, 2) for Red Listed species.
  med <- prediction[2, , red, drop = FALSE]
  expect_gt(mean(med[1, 3:4, ] - med[1, 1:2, ]), 0)
})
