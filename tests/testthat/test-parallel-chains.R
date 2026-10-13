# Parallel MCMC chains. Design: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.
# Rules of this file: at most two workers (R CMD check's limit), structural
# assertions in the always-on tier, fit identity only where the spec allows it.

alive <- function(pids) vapply(pids, function(p) isTRUE(tools::pskill(p, signal = 0L)), logical(1))

test_that("MCMCparams$ncores is validated", {
  for (bad in list(0, 1.5, "2", c(1, 2), NA_real_, Inf)) {
    expect_error(occJSDM:::resolve_ncores(list(ncores = bad)), "single whole number")
  }
  expect_identical(occJSDM:::resolve_ncores(list()), 1L)
  expect_identical(occJSDM:::resolve_ncores(list(ncores = 2)), 2L)      # a double is accepted
  expect_identical(occJSDM:::resolve_ncores(list(ncores = 3L)), 3L)
  # The detected count is injected: detectCores() can return NA on some
  # systems, and the warning must not depend on what the operating system says.
  expect_warning(occJSDM:::resolve_ncores(list(ncores = 5), available = 4L), "exceeds the 4 cores")
  expect_silent(occJSDM:::resolve_ncores(list(ncores = 4), available = 4L))
  expect_silent(occJSDM:::resolve_ncores(list(ncores = 64), available = NA_integer_))
  expect_identical(occJSDM:::resolve_ncores(list(ncores = 64), available = NA_integer_), 64L)
})

test_that("the core budget splits into workers and threads per worker", {
  b <- occJSDM:::plan_chain_budget
  pick <- function(x) x[c("workers", "threads_per_worker", "thread_setting")]
  expect_identical(pick(b(4L, 1L)), list(workers = 1L, threads_per_worker = 1L, thread_setting = "ambient"))
  expect_identical(pick(b(4L, 4L)), list(workers = 4L, threads_per_worker = 1L, thread_setting = "budget"))
  expect_identical(pick(b(4L, 8L)), list(workers = 4L, threads_per_worker = 2L, thread_setting = "budget"))
  expect_identical(pick(b(2L, 3L)), list(workers = 2L, threads_per_worker = 1L, thread_setting = "budget"))
  expect_identical(pick(b(1L, 4L)), list(workers = 1L, threads_per_worker = 4L, thread_setting = "budget"))
  expect_identical(b(4L, 8L)$ncores, 8L)
})

test_that("the ambient RcppParallel count follows the environment variable", {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS") else Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old))
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "3")
  expect_identical(occJSDM:::ambient_rcpp_threads(), 3L)
  Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
  expect_identical(occJSDM:::ambient_rcpp_threads(), as.integer(RcppParallel::defaultNumThreads()))
})

test_that("per-chain WAIC accumulators merge to the single-pass values", {
  set.seed(11)
  loglik <- matrix(rnorm(7 * 12, mean = -1), nrow = 7, ncol = 12)   # 7 observations, 12 draws
  single <- occJSDM:::create_waic_quantities(7)
  for (t in 1:12) single <- occJSDM:::update_waic_summary(loglik[, t], single, t)
  accumulate <- function(cols) {
    part <- occJSDM:::create_waic_quantities(7)
    for (i in seq_along(cols)) part <- occJSDM:::update_waic_summary(loglik[, cols[i]], part, i)
    c(part, list(n = length(cols)))
  }
  merged <- occJSDM:::merge_waic_parts(list(accumulate(1:5), accumulate(6:12)))   # unequal parts on purpose
  expect_equal(merged$n, 12L)
  expect_equal(merged$mean_log, single$mean_log, tolerance = 1e-12)
  expect_equal(merged$mean_lik, single$mean_lik, tolerance = 1e-12)
  expect_equal(merged$M2, single$M2, tolerance = 1e-12)
  expect_equal(occJSDM:::compute_waic(merged, 12), occJSDM:::compute_waic(single, 12), tolerance = 1e-12)
  three <- occJSDM:::merge_waic_parts(list(accumulate(1:4), accumulate(5:8), accumulate(9:12)))
  expect_equal(three$M2, single$M2, tolerance = 1e-12)
  expect_null(occJSDM:::merge_waic_parts(list(NULL, NULL)))
  expect_equal(occJSDM:::merge_waic_parts(list(NULL, accumulate(1:12)))$M2, single$M2, tolerance = 1e-12)
})

test_that("latent-state sums add elementwise and tolerate absent parts", {
  a <- list(z = matrix(1, 2, 2), psi = matrix(2, 2, 2), theta = NULL, w = NULL)
  b <- list(z = matrix(3, 2, 2), psi = matrix(4, 2, 2), theta = NULL, w = NULL)
  s <- occJSDM:::add_latent_sums(NULL, a)
  expect_identical(s, a)
  s <- occJSDM:::add_latent_sums(s, b)
  expect_equal(s$z, matrix(4, 2, 2))
  expect_equal(s$psi, matrix(6, 2, 2))
  expect_null(s$theta)
  expect_null(occJSDM:::add_latent_sums(NULL, NULL))
})

test_that("latent-state sums return the other side for an absent part", {
  a <- list(z = matrix(1, 2, 2), psi = matrix(2, 2, 2), theta = NULL, w = NULL)
  expect_identical(occJSDM:::add_latent_sums(a, NULL), a)
  expect_identical(occJSDM:::add_latent_sums(NULL, a), a)
})

test_that("a latent-state name missing or NULL on either side sums to NULL", {
  total <- list(z = matrix(1, 2, 2), psi = matrix(2, 2, 2))
  part <- list(z = matrix(3, 2, 2), theta = matrix(5, 2, 2))
  s <- occJSDM:::add_latent_sums(total, part)
  expect_equal(s$z, matrix(4, 2, 2))
  expect_null(s$psi)       # in total, absent from part: not left as total's unsummed matrix
  expect_null(s$theta)     # in part, absent from total: not silently dropped into a partial sum
  expect_setequal(names(s), c("z", "psi", "theta"))
  # The same holds when the name is present on one side as an explicit NULL.
  s2 <- occJSDM:::add_latent_sums(list(z = matrix(1, 2, 2), psi = matrix(2, 2, 2)),
                                  list(z = matrix(3, 2, 2), psi = NULL))
  expect_equal(s2$z, matrix(4, 2, 2))
  expect_null(s2$psi)
  s3 <- occJSDM:::add_latent_sums(list(z = matrix(1, 2, 2), psi = NULL),
                                  list(z = matrix(3, 2, 2), psi = matrix(2, 2, 2)))
  expect_null(s3$psi)
})

test_that("latent-state sums fold over three parts and keep their names", {
  mk <- function(zv, pv) list(z = matrix(zv, 2, 2), psi = matrix(pv, 2, 2), theta = NULL, w = NULL)
  total <- NULL
  for (part in list(mk(1, 2), mk(3, 4), mk(5, 6))) total <- occJSDM:::add_latent_sums(total, part)
  expect_equal(total$z, matrix(9, 2, 2))
  expect_equal(total$psi, matrix(12, 2, 2))
  expect_null(total$theta)
  expect_null(total$w)
  expect_named(total, c("z", "psi", "theta", "w"))   # NULL entries survive a second and a third add
})

test_that("chain seeds are a 2 by nchain integer matrix drawn from R's RNG", {
  set.seed(3); a <- occJSDM:::draw_chain_seeds(3L)
  set.seed(3); b <- occJSDM:::draw_chain_seeds(3L)
  expect_identical(a, b)
  expect_identical(dim(a), c(2L, 3L))
  expect_identical(rownames(a), c("r", "cpp"))
  expect_type(a, "integer")
  expect_false(anyDuplicated(as.vector(a)) > 0)
  set.seed(4); expect_false(identical(occJSDM:::draw_chain_seeds(3L), a))
})

test_that("a chain's seeds do not depend on how many chains follow it", {
  set.seed(8); one <- occJSDM:::draw_chain_seeds(1L)
  set.seed(8); two <- occJSDM:::draw_chain_seeds(2L)
  expect_identical(one[, 1], two[, 1])
})

test_that("restoring the session RNG puts the value and the kinds back", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(9); invisible(runif(2))
  saved <- occJSDM:::save_session_rng()
  expected <- runif(3)
  occJSDM:::seed_chain_rng(123L, 456L)       # pins Mersenne-Twister etc., as a chain does
  invisible(rnorm(5))
  occJSDM:::restore_session_rng(saved)
  expect_identical(RNGkind()[1:3], c("L'Ecuyer-CMRG", "Inversion", "Rejection"))
  expect_identical(runif(3), expected)
})

test_that("chain 1 of a two-chain fit equals the one-chain fit under the same seed", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(515); one <- fit_fixture(sim, MCMCparams = list(nchain = 1, nburn = 20, niter = 20, nthin = 1))
  set.seed(515); two <- fit_fixture(sim, MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1))
  expect_identical(one$results_output$jsdm_output$B0_output[, , 1], two$results_output$jsdm_output$B0_output[, , 1])
  expect_identical(one$results_output$p_output[, , , 1], two$results_output$p_output[, , , 1])
  expect_identical(one$infos$chain_seeds[, 1], two$infos$chain_seeds[, 1])
  expect_false(identical(two$results_output$jsdm_output$B0_output[, , 1], two$results_output$jsdm_output$B0_output[, , 2]))
})

test_that("every model type fits through run_one_chain() in process", {
  # Direct coverage of the extraction for all four models, so that a missing
  # setup field (such as the z or w state the parent sets from y) is caught
  # here and not only by the wider suite.
  for (model in c("binary", "continuous", "occupancy", "two_stage")) {
    sim <- simulate_fixture(model = model,
                            M = if (model %in% c("occupancy", "two_stage")) 2L else 1L,
                            P = if (model == "two_stage") 2L else 1L,
                            K = if (model == "two_stage") 2L else 1L)
    set.seed(1)
    fit <- fit_fixture(sim, spatial = FALSE, MCMCparams = list(nchain = 2, nburn = 5, niter = 5, nthin = 1))
    expect_identical(fit$infos$model, model, info = model)
    expect_identical(dim(fit$results_output$jsdm_output$B0_output), c(FIXTURE_S, 5L, 2L), info = model)
    expect_true(is.finite(fit$results_output$WAIC), info = model)
    expect_identical(dim(fit$infos$chain_seeds), c(2L, 2L), info = model)
  }
})

test_that("a fit leaves the session RNG just after its seed draws and clears the setup", {
  # The state after a fit is the state after the chain seeds were drawn, so a
  # third chain's seeds, had one been requested, are exactly the next draws.
  # Comparing two fits rather than replaying the draws by hand keeps the test
  # true even if the data preparation before the loop uses R's RNG.
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(21); two <- fit_fixture(sim, MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1))
  next_pair <- sample.int(.Machine$integer.max, 2L)
  set.seed(21); three <- fit_fixture(sim, MCMCparams = list(nchain = 3, nburn = 20, niter = 20, nthin = 1))
  expect_identical(unname(three$infos$chain_seeds[, 3]), next_pair)
  expect_identical(two$infos$chain_seeds, three$infos$chain_seeds[, 1:2])
  expect_false(exists("setup", envir = occJSDM:::chain_env, inherits = FALSE))
})

test_that("the new infos fields record the sampling settings, seeds and execution", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(6); fit <- fit_fixture(sim, MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1))
  expect_identical(fit$infos$MCMCparams, list(nchain = 2L, nburn = 20L, niter = 20L, nthin = 1L))
  expect_identical(dim(fit$infos$chain_seeds), c(2L, 2L))
  ex <- fit$infos$execution
  expect_identical(ex$ncores, 1L); expect_identical(ex$workers, 1L)
  expect_identical(ex$thread_setting, "ambient")
  expect_identical(ex$threads_per_worker, occJSDM:::ambient_rcpp_threads())
  expect_null(ex$progress_file); expect_null(ex$blas_threads_requested)
  expect_identical(ex$blas_declared, nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")))
  expect_length(ex$elapsed, 2L)
  expect_null(occJSDM:::drop_execution(fit)$infos$execution)
})

test_that("summarised latent means equal the mean over all retained draws", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  mc <- list(nchain = 2, nburn = 20, niter = 20, nthin = 1)
  set.seed(31); summarised <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 2, n_supportpoints = FIXTURE_KNOTS),
    occCovariates = fixture_occ_covariates(), spatCovariates = fixture_spat_covariates(),
    MCMCparams = mc, summarisedLatentPresences = TRUE)))
  set.seed(31); full <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 2, n_supportpoints = FIXTURE_KNOTS),
    occCovariates = fixture_occ_covariates(), spatCovariates = fixture_spat_covariates(),
    MCMCparams = mc, summarisedLatentPresences = FALSE)))
  z_mean <- apply(full$results_output$z_output, c(1, 2), mean)
  expect_equal(summarised$results_output$z_output, z_mean, tolerance = 1e-12)
  expect_identical(summarised$results_output$WAIC, full$results_output$WAIC)
  expect_identical(summarised$results_output$jsdm_output$B0_output, full$results_output$jsdm_output$B0_output)
})
