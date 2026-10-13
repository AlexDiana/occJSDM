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
  psi_mean <- apply(full$results_output$psi_output, c(1, 2), mean)
  # The running sum keeps psi's dimnames, the per-draw array does not.
  expect_equal(unname(summarised$results_output$psi_output), psi_mean, tolerance = 1e-12)
  expect_identical(summarised$results_output$theta_output, full$results_output$theta_output)
  expect_identical(summarised$results_output$w_output, full$results_output$w_output)
  expect_identical(summarised$results_output$WAIC, full$results_output$WAIC)
  expect_identical(summarised$results_output$jsdm_output$B0_output, full$results_output$jsdm_output$B0_output)
})

test_that("every WAIC part of every chain records the draws it accumulated", {
  # The parent checks each non-NULL part's draw count, w and y included, so a
  # chain that accumulated a different number of draws cannot be merged with
  # the wrong weight.
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  mc <- list(nchain = 2, nburn = 2, niter = 3, nthin = 1)
  original <- occJSDM:::run_chain_batch
  seen <- list()
  local_mocked_bindings(run_chain_batch = function(handle, batch, seeds, task = occJSDM:::run_chain_task) {
    res <- original(handle, batch, seeds)
    seen <<- c(seen, res)
    res
  }, .package = "occJSDM")
  set.seed(2); fit_fixture(sim, MCMCparams = mc)
  expect_length(seen, 2L)
  for (res in seen) for (nm in c("jsdm", "w", "y")) expect_equal(res$waic[[nm]]$n, 3, info = nm)

  local_mocked_bindings(run_chain_batch = function(handle, batch, seeds, task = occJSDM:::run_chain_task) {
    res <- original(handle, batch, seeds)
    for (i in seq_along(res)) if (res[[i]]$chain == 2L) res[[i]]$waic$y$n <- res[[i]]$waic$y$n - 1
    res
  }, .package = "occJSDM")
  set.seed(2); expect_error(fit_fixture(sim, MCMCparams = mc), "Current wAIC iter wrong")
  expect_false(exists("setup", envir = occJSDM:::chain_env, inherits = FALSE))
})

test_that("a chain that returns no result stops the fit with a clear message", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  mc <- list(nchain = 2, nburn = 2, niter = 3, nthin = 1)
  original <- occJSDM:::run_chain_batch
  # A batch that comes back without one of its chains.
  local_mocked_bindings(run_chain_batch = function(handle, batch, seeds, task = occJSDM:::run_chain_task) {
    res <- original(handle, batch, seeds)
    Filter(function(r) r$chain != 2L, res)
  }, .package = "occJSDM")
  set.seed(2); expect_error(fit_fixture(sim, MCMCparams = mc), "No result was returned for chain(s) 2", fixed = TRUE)
  # A chain whose result arrives without its WAIC parts.
  local_mocked_bindings(run_chain_batch = function(handle, batch, seeds, task = occJSDM:::run_chain_task) {
    res <- original(handle, batch, seeds)
    for (i in seq_along(res)) if (res[[i]]$chain == 2L) res[[i]]["waic"] <- list(NULL)
    res
  }, .package = "occJSDM")
  set.seed(2); expect_error(fit_fixture(sim, MCMCparams = mc), "No result was returned for chain(s) 2", fixed = TRUE)
  expect_false(exists("setup", envir = occJSDM:::chain_env, inherits = FALSE))
})

test_that("the R-code fingerprint matches across fresh workers and changes when a body changes", {
  skip_on_cran()
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  parent <- occJSDM:::package_fingerprint()
  workers <- parallel::clusterCall(handle$cl, occJSDM:::package_fingerprint)
  for (w in workers) expect_identical(w, parent)
  expect_true(occJSDM:::check_worker_builds(handle))
  expect_identical(handle$worker_threads, c("1", "1"))
  expect_identical(length(handle$pids), 2L)
  expect_true(all(alive(handle$pids)))
  # The name order does not follow the session's collation: a parent that
  # changed LC_COLLATE still matches the fresh workers. testthat sets
  # LC_COLLATE = C, in the environment too, so the workers start in C and
  # the parent switches to a collation that orders the names differently.
  old_collate <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old_collate), add = TRUE)
  worker_collate <- parallel::clusterCall(handle$cl, Sys.getlocale, "LC_COLLATE")[[1]]
  new_collate <- if (identical(worker_collate, "C")) "en_US.UTF-8" else "C"
  skip_if(identical(suppressWarnings(Sys.setlocale("LC_COLLATE", new_collate)), ""),
          paste(new_collate, "is not available"))
  for (w in workers) expect_identical(w, occJSDM:::package_fingerprint())
  expect_true(occJSDM:::check_worker_builds(handle))
  # negative case on a toy namespace-like environment
  fp <- function(env) {
    nms <- sort(ls(env, all.names = TRUE))
    text <- unlist(lapply(nms, function(n) c(n, deparse(get(n, envir = env)))))
    tmp <- tempfile(); on.exit(unlink(tmp)); writeLines(text, tmp); unname(tools::md5sum(tmp))
  }
  e1 <- new.env(); e1$f <- function(x) x + 1
  e2 <- new.env(); e2$f <- function(x) x + 2
  expect_false(identical(fp(e1), fp(e2)))
})

test_that("an error or an interrupt during worker start-up leaves no worker alive", {
  skip_on_cran()
  # The injected bootstraps run in the workers with globalenv() as their
  # environment, so they learn where to record their PID, and which process
  # to interrupt, from environment variables, which the workers inherit at
  # process start.
  pid_dir <- tempfile("worker-pids-"); dir.create(pid_dir)
  Sys.setenv(OCCJSDM_TEST_PIDDIR = pid_dir, OCCJSDM_TEST_PARENT_PID = Sys.getpid())
  on.exit(Sys.unsetenv(c("OCCJSDM_TEST_PIDDIR", "OCCJSDM_TEST_PARENT_PID")), add = TRUE)
  failing_boot <- function(libpaths, threads, dev_path) {
    d <- Sys.getenv("OCCJSDM_TEST_PIDDIR"); writeLines(as.character(Sys.getpid()), file.path(d, as.character(Sys.getpid())))
    stop("boot failed")
  }
  expect_error(occJSDM:::start_chain_cluster(2L, 1L, tempfile(fileext = ".log"), bootstrap = failing_boot), "boot failed")
  Sys.sleep(0.5)
  pids <- as.integer(list.files(pid_dir))
  expect_length(pids, 2L)
  expect_false(any(alive(pids)))
  unlink(file.path(pid_dir, list.files(pid_dir)))

  skip_on_os("windows")
  # Interrupts the parent itself once its PID is recorded, so the interrupt
  # always lands while the parent waits on the bootstrap.
  slow_boot <- function(libpaths, threads, dev_path) {
    d <- Sys.getenv("OCCJSDM_TEST_PIDDIR"); writeLines(as.character(Sys.getpid()), file.path(d, as.character(Sys.getpid())))
    system(sprintf("kill -INT %s", Sys.getenv("OCCJSDM_TEST_PARENT_PID")))
    Sys.sleep(10); list(pid = Sys.getpid(), threads = "1")
  }
  outcome <- tryCatch({
    occJSDM:::start_chain_cluster(2L, 1L, tempfile(fileext = ".log"), bootstrap = slow_boot); "completed"
  }, interrupt = function(c) "interrupted")
  expect_identical(outcome, "interrupted")
  Sys.sleep(0.5)
  pids <- as.integer(list.files(pid_dir))
  expect_length(pids, 2L)
  expect_false(any(alive(pids)))
})

test_that("a build mismatch is reported with both paths and the differing field", {
  parent <- list(path = "/a/occJSDM", version = "0.1.1", dll_md5 = "x", code_md5 = "y")
  expect_true(occJSDM:::compare_fingerprints(parent, parent, 1L))
  worker <- parent; worker$code_md5 <- "z"; worker$path <- "/b/occJSDM"
  expect_error(occJSDM:::compare_fingerprints(parent, worker, 2L), "Worker 2 .*path.*/a/occJSDM.*/b/occJSDM")
  worker <- parent; worker$dll_md5 <- "q"
  expect_error(occJSDM:::compare_fingerprints(parent, worker, 1L), "dll_md5")
})

test_that("a chain that fails in a worker is reported by number and no worker survives", {
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  pids <- handle$pids
  occJSDM:::clear_chain_setup()                      # the real task will fail: no setup in the workers
  res <- occJSDM:::run_chain_batch(handle, 1:2, matrix(1:4, 2, dimnames = list(c("r", "cpp"), NULL)))
  expect_s3_class(res[[2]], "occJSDM_chain_failure")
  expect_error(occJSDM:::stop_if_chain_failed(res[[2]]), "^Chain 2 failed: ")
  occJSDM:::stop_chain_cluster(handle)
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("cleanup kills every worker even when stopCluster() itself errors", {
  skip_on_cran()
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  pids <- handle$pids
  close(handle$cl[[1]]$con)                          # a broken connection makes stopCluster() throw
  expect_true(alive(pids[2]))
  expect_silent(occJSDM:::stop_chain_cluster(handle))
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("abandoning running workers leaves no worker process alive", {
  skip_on_cran(); skip_on_os("windows")
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  pids <- handle$pids
  system(sprintf("(sleep 1; kill -INT %d) > /dev/null 2>&1 &", Sys.getpid()))   # timer starts after the cluster is up
  t0 <- Sys.time()
  outcome <- tryCatch({
    occJSDM:::run_chain_batch(handle, 1:2, matrix(1:4, 2), task = function(k, seeds) { Sys.sleep(10); k })
    "completed"
  }, interrupt = function(c) "interrupted")
  waited <- as.numeric(Sys.time() - t0, units = "secs")
  expect_identical(outcome, "interrupted")
  expect_lt(waited, 8)                               # the workers were still busy when cleanup begins
  expect_true(all(alive(pids)))
  occJSDM:::stop_chain_cluster(handle)
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("an interrupt during shutdown still kills every worker", {
  skip_on_cran(); skip_on_os("windows")
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  pids <- handle$pids
  # Abandon both workers in a 30-second task, as an interrupted fit does.
  system(sprintf("(sleep 1; kill -INT %d) > /dev/null 2>&1 &", Sys.getpid()))
  abandoned <- tryCatch({
    occJSDM:::run_chain_batch(handle, 1:2, matrix(1:4, 2), task = function(k, seeds) { Sys.sleep(30); k })
    FALSE
  }, interrupt = function(c) TRUE)
  expect_true(abandoned)
  expect_true(all(alive(pids)))
  # stopCluster() now waits up to 5 s for the busy workers to close their
  # sockets; a second interrupt arrives 1 s into that wait.
  system(sprintf("(sleep 1; kill -INT %d) > /dev/null 2>&1 &", Sys.getpid()))
  returned <- FALSE
  outcome <- tryCatch({
    occJSDM:::stop_chain_cluster(handle)
    returned <- TRUE
    Sys.sleep(5)                                     # absorbs the signal should stop_chain_cluster() return first
    "completed"
  }, interrupt = function(c) "interrupted")
  expect_identical(outcome, "interrupted")
  expect_false(returned)                             # the interrupt landed inside stop_chain_cluster()
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("the BLAS thread variables are set for the workers only", {
  skip_on_cran()
  vars <- names(occJSDM:::blas_thread_variables())
  before <- Sys.getenv(vars, unset = NA_character_, names = TRUE)
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  expect_identical(Sys.getenv(vars, unset = NA_character_, names = TRUE), before)
  seen <- parallel::clusterCall(handle$cl, function(v) Sys.getenv(v, names = TRUE), vars)
  for (s in seen) expect_identical(unname(s), unname(occJSDM:::blas_thread_variables()))
  expect_identical(handle$blas_threads_requested, occJSDM:::blas_thread_variables())
})

test_that("the RcppParallel stack size survives the in-process thread pin and reaches the workers", {
  skip_on_cran()
  vars <- c("RCPP_PARALLEL_NUM_THREADS", "RCPP_PARALLEL_STACK_SIZE")
  old <- Sys.getenv(vars, unset = NA_character_, names = TRUE)
  on.exit(for (nm in names(old)) if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(setNames(old[[nm]], nm))), add = TRUE)
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "3", RCPP_PARALLEL_STACK_SIZE = "8388608")
  restore <- occJSDM:::pin_session_threads(2L)
  expect_identical(unname(Sys.getenv(vars)), c("2", "8388608"))
  restore()
  expect_identical(unname(Sys.getenv(vars)), c("3", "8388608"))
  # cleanup after an error, registered the way runOccJSDM() registers it
  f <- function() { restore <- occJSDM:::pin_session_threads(2L); on.exit(restore(), add = TRUE); stop("boom") }
  expect_error(f(), "boom")
  expect_identical(unname(Sys.getenv(vars)), c("3", "8388608"))
  # an unset thread count comes back unset, the stack size untouched
  Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
  restore <- occJSDM:::pin_session_threads(2L); restore()
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_), NA_character_)
  expect_identical(Sys.getenv("RCPP_PARALLEL_STACK_SIZE"), "8388608")
  # workers get the budget's thread count and inherit the stack size
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile(fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  seen <- parallel::clusterCall(handle$cl, function(v) Sys.getenv(v), vars)
  for (s in seen) expect_identical(unname(s), c("1", "8388608"))
})

test_that("the reset-then-restore helper leaves the same state after an error in either driver", {
  # Helper-level: it checks save_session_rng() and restore_session_rng()
  # around a failing batch in both drivers. The production orchestration,
  # which registers the same helper with on.exit() in runOccJSDM(), is
  # exercised with a real interrupt in Task 5's tests.
  skip_on_cran()
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  run_guarded <- function(driver) {
    saved <- occJSDM:::save_session_rng()
    on.exit(occJSDM:::restore_session_rng(saved))
    driver()
  }
  seeds <- matrix(1:4, 2, dimnames = list(c("r", "cpp"), NULL))
  failing <- function(k, seeds) occJSDM:::chain_failure(k, "boom")
  serial_error <- function() for (r in occJSDM:::run_chain_batch(NULL, 1:2, seeds, task = failing)) occJSDM:::stop_if_chain_failed(r)
  parallel_error <- function() {
    handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile(fileext = ".log"))
    on.exit(occJSDM:::stop_chain_cluster(handle))
    for (r in occJSDM:::run_chain_batch(handle, 1:2, seeds, task = failing)) occJSDM:::stop_if_chain_failed(r)
  }
  next_draws <- function() list(rnorm(3), runif(3), sample(100, 3))
  for (normal_kind in c("Inversion", "Box-Muller")) {
    RNGkind("Mersenne-Twister", normal_kind, "Rejection")
    set.seed(77); invisible(rnorm(1))                 # under Box-Muller this leaves a cached deviate pending
    expect_error(run_guarded(serial_error), "Chain 1 failed: boom"); after_serial <- next_draws()
    set.seed(77); invisible(rnorm(1))
    expect_error(run_guarded(parallel_error), "Chain 1 failed: boom"); after_parallel <- next_draws()
    expect_identical(after_serial, after_parallel)
  }
})

test_that("cleanup closes the sockets that stopCluster() leaves open when it gives up", {
  skip_on_cran()
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  close(handle$cl[[1]]$con)                          # stopCluster() throws here
  occJSDM:::stop_chain_cluster(handle)
  # Node 2's socket is closed now, not by the garbage collector later, with a warning.
  expect_false(any(vapply(handle$cl, function(n) occJSDM:::is_live_connection(n$con), logical(1))))
  # A second call is a no-op: a node's stale connection object must not write
  # into a connection that has since taken its number.
  numbers <- vapply(handle$cl, function(n) as.integer(n$con), integer(1))
  opened <- list()
  repeat {                                           # R hands out the lowest free number
    path <- tempfile()
    opened <- c(opened, list(file(path, "wb")))
    if (as.integer(opened[[length(opened)]]) %in% numbers) break
  }
  on.exit(for (con in opened) close(con), add = TRUE)
  occJSDM:::stop_chain_cluster(handle)
  flush(opened[[length(opened)]])
  expect_identical(file.size(path), 0)
  # Cleanup closes only a connection that is still the node's own: R reuses a
  # closed connection's number, and close() on the stale object would close
  # whichever connection holds that number now.
  a <- file(tempfile(), "wb")
  expect_true(occJSDM:::is_live_connection(a))
  close(a)
  b <- file(tempfile(), "wb")
  on.exit(close(b), add = TRUE)
  expect_identical(as.integer(b), as.integer(a))     # b has taken a's number
  expect_false(occJSDM:::is_live_connection(a))
  expect_true(occJSDM:::is_live_connection(b))
})
