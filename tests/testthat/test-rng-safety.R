# Random draws must not depend on TBB scheduling, even when deterministic
# probability calculations are allowed several threads.
with_rng_threads <- function(n, code) {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS") else
      Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old)
  }, add = TRUE)
  RcppParallel::setThreadOptions(numThreads = n)
  force(code)
}

test_that("PG draws follow one advancing stream regardless of requested threads", {
  # Equal predictors must receive distinct draws, rather than replaying the
  # same stream from each TBB worker (whose OpenMP thread ID is always zero).
  eta <- matrix(0, 256, 32)
  draw <- function(threads) with_rng_threads(threads, {
    setOccJSDMSeed(4242L)
    samplePGvariables_parallel(eta)
  })
  serial <- draw(1L)
  parallel_requested <- draw(4L)
  expect_identical(parallel_requested, serial)
  expect_length(unique(as.vector(parallel_requested)), length(eta))
  next_draw <- with_rng_threads(4L, samplePGvariables_parallel(eta))
  expect_false(identical(next_draw, parallel_requested))
})

test_that("collection draws agree with the serial sampler under several threads", {
  # This exercises both PG draws and the R/Armadillo normal draw reached by
  # sample_beta_nocov_cpp_TS(), including the collection slope.
  n <- 64L
  S <- 32L
  args <- list(w = matrix(rep(c(0, 1, 1, 0), n * S / 4), n, S),
               z = matrix(1, n, S), beta_theta = matrix(0, 2, S),
               idx_z = seq_len(n), X_theta = cbind(1, seq(-1, 1, length.out = n)),
               b_betatheta = c(0, 0), B_betatheta = diag(2))
  draw <- function(fun, threads) with_rng_threads(threads, {
    set.seed(42)
    setOccJSDMSeed(4242L)
    value <- do.call(fun, args)
    list(value = value, seed = .Random.seed)
  })
  reference <- draw(sample_betatheta_cpp, 1L)
  actual <- draw(sample_betatheta_cpp_parallel, 4L)
  expect_identical(actual$value, reference$value)
  expect_identical(actual$seed, reference$seed)
})

test_that("public fits reproduce across thread requests for every model", {
  for (model in c("binary", "continuous", "occupancy", "two_stage")) {
    sim <- simulate_fixture(model = model, P = if (model == "two_stage") 2L else 1L)
    fit <- function(threads) with_rng_threads(threads, {
      set.seed(932L)
      value <- suppressMessages(suppressWarnings(runOccJSDM(
        sim$data_list, occCovariates = fixture_occ_covariates(),
        collCovariates = if (model %in% c("occupancy", "two_stage")) "X_theta" else NULL,
        spatCovariates = fixture_spat_covariates(),
        listParams = list(n_factors = 2, n_supportpoints = FIXTURE_KNOTS),
        MCMCparams = list(nchain = 2, nburn = 3, niter = 5, nthin = 1))))
      expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS"), as.character(threads))
      list(output = value$results_output, seed = .Random.seed)
    })
    expect_identical(fit(4L), fit(1L), info = model)
    expect_identical(fit(4L), fit(4L), info = paste(model, "repeat"))
  }
})
