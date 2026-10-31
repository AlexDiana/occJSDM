# inferDataModel() decides between the no-stage, one-stage (occupancy) and
# two-stage model types from the shape of data$info alone. These tests pin
# that decision down for the two ways a survey can number its field samples:
# unique across the whole survey, or restarting at 1 within each site.
#
# The site-local case was misclassified until the check moved from the Sample
# column alone to (Site, Sample) pairs: a one-stage survey that numbers its
# samples 1, 2, 3 at every site repeats each Sample ID across sites, which
# looked like repeated observations within a sample and so like PCR
# replicates.

# Build data$info for n sites with M field samples per site and K PCR
# replicates per sample (K = 1 gives one-stage data, one row per sample).
# Sample IDs are either unique across the survey or restart at each site.
make_info <- function(n = 4L, M = 3L, K = 1L, site_local = FALSE) {
  info <- expand.grid(k = seq_len(K), m = seq_len(M), site = seq_len(n),
                      KEEP.OUT.ATTRS = FALSE)
  info <- info[order(info$site, info$m, info$k), ]
  sample_id <- if (site_local) info$m else (info$site - 1L) * M + info$m
  data.frame(Site = info$site,
             Sample = sample_id,
             Primer = 1L,
             X_psi.1 = stats::rnorm(nrow(info)))
}

make_data <- function(info, seed = 1L) {
  set.seed(seed)
  OTU <- matrix(stats::rbinom(nrow(info) * 2L, 1, 0.5), ncol = 2L,
                dimnames = list(NULL, c("OTU_1", "OTU_2")))
  list(info = info, OTU = OTU)
}

infer_quietly <- function(data) suppressMessages(inferDataModel(data))

test_that("one-stage data with survey-wide sample IDs is inferred as occupancy", {
  data <- make_data(make_info(K = 1L, site_local = FALSE))
  expect_false(any(duplicated(data$info$Sample)))
  expect_identical(infer_quietly(data), "occupancy")
})

test_that("one-stage data with site-local sample IDs is inferred as occupancy", {
  # Sample runs 1, 2, 3 at every site, so Sample alone is duplicated across
  # sites but each (Site, Sample) pair appears exactly once.
  data <- make_data(make_info(K = 1L, site_local = TRUE))
  expect_true(any(duplicated(data$info$Sample)))
  expect_false(any(duplicated(data$info[, c("Site", "Sample")])))
  expect_identical(infer_quietly(data), "occupancy")
})

test_that("two-stage data with survey-wide sample IDs is inferred as two_stage", {
  data <- make_data(make_info(K = 2L, site_local = FALSE))
  expect_identical(infer_quietly(data), "two_stage")
})

test_that("two-stage data with site-local sample IDs is inferred as two_stage", {
  data <- make_data(make_info(K = 2L, site_local = TRUE))
  expect_true(any(duplicated(data$info[, c("Site", "Sample")])))
  expect_identical(infer_quietly(data), "two_stage")
})

test_that("one row per site is inferred as a no-stage model", {
  data <- make_data(make_info(M = 1L, K = 1L))
  expect_identical(infer_quietly(data), "binary")
})

test_that("a one-stage fit accepts site-local sample IDs", {
  # The fit groups by (Site, Sample) everywhere, so site-local IDs must run
  # end to end, and the fitted object must record the one-stage model.
  info <- make_info(n = FIXTURE_N, M = 2L, K = 1L, site_local = TRUE)
  info$X_psi.EnvCov.1 <- stats::rnorm(nrow(info))
  info$X_psi.EnvCov.2 <- stats::rnorm(nrow(info))
  data <- make_data(info)
  data$OTU <- cbind(data$OTU,
                    OTU_3 = stats::rbinom(nrow(info), 1, 0.5),
                    OTU_4 = stats::rbinom(nrow(info), 1, 0.5))
  fit <- suppressMessages(suppressWarnings(
    runOccJSDM(data,
               occCovariates = fixture_occ_covariates(),
               listParams = list(n_factors = 2L),
               MCMCparams = FIXTURE_MCMC)))
  expect_identical(fit$infos$model, "occupancy")
  expect_identical(fit$infos$n, FIXTURE_N)
})

test_that("NA outside the two-stage model stops with a clean message", {
  data <- make_data(make_info(K = 1L))
  data$OTU[1, 1] <- NA
  err <- expect_error(suppressMessages(runOccJSDM(data)),
                      "allowed only in the two-stage model")
  expect_false(grepl("FALSE", conditionMessage(err)))
})

test_that("a one-stage fit accepts unequal samples per site", {
  # A lost field sample is deleted, not recorded as NA, which leaves some
  # sites with fewer samples than others.
  info <- make_info(n = FIXTURE_N, M = 3L, K = 1L)
  info <- info[!(info$Site <= FIXTURE_N %/% 2L & info$Sample %% 3L == 0L), ]
  info$X_psi.EnvCov.1 <- stats::rnorm(nrow(info))
  info$X_psi.EnvCov.2 <- stats::rnorm(nrow(info))
  data <- make_data(info)
  fit <- suppressMessages(suppressWarnings(
    runOccJSDM(data,
               occCovariates = fixture_occ_covariates(),
               MCMCparams = FIXTURE_MCMC)))
  expect_identical(fit$infos$model, "occupancy")
  expect_identical(sort(unique(as.integer(table(info$Site)))), c(2L, 3L))
})
