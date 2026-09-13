# Focused RNG diagnostics, independent of posterior coverage/recovery.
# Run from the package root. Outputs belong outside the source checkout.
# Rscript dev/simstudy/validate_rng_safety.R --out=/absolute/output/directory
# Add --source=/path/to/instrumented/package for the C++ thread audit.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
source_dir <- normalizePath(option("source", "."))
out_dir <- option("out", "")
if (!nzchar(out_dir)) stop("Supply --out with a dedicated output directory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(source_dir, "src", "functions.cpp")) &&
    any(grepl("rngAuditCounts", readLines(file.path(source_dir, "src", "functions.cpp"))))) {
  Rcpp::compileAttributes(source_dir)
}
devtools::load_all(source_dir, quiet = TRUE)
source("tests/testthat/helper-fixtures.R")
old_threads <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
RcppParallel::setThreadOptions(numThreads = 4)

moment_summary <- function(x, name, expected_mean, expected_var) {
  pairs <- matrix(as.vector(x), ncol = 2, byrow = TRUE)
  data.frame(distribution = name, draws = length(x), mean = mean(x),
             expected_mean, mean_error = mean(x) - expected_mean,
             mean_mcse = sqrt(expected_var / length(x)), variance = var(as.vector(x)),
             expected_var, adjacent_correlation = cor(pairs[, 1], pairs[, 2]),
             duplicate_draws = length(x) - length(unique(as.vector(x))))
}

# Independently known moments of PG(1, 0): mean 1/4, variance 1/24.
set.seed(6123)
setOccJSDMSeed(6182L)
pg <- samplePGvariables_parallel(matrix(0, 256, 128))
moments <- moment_summary(pg, "PG(1,0)", 1 / 4, 1 / 24)

# A zero design contributes no likelihood precision: the conditional draw
# is exactly its Normal(0.7, variance 1.6) prior. This tests the R-backed
# collection normal leaf separately from the C++ normal leaf below.
S <- 32768L
collection <- sample_betatheta_cpp_parallel(
  matrix(0, 1, S), matrix(1, 1, S), matrix(0, 1, S), 1L,
  matrix(0, 1, 1), 0.7, matrix(1.6, 1, 1))
moments <- rbind(moments, moment_summary(collection, "collection Normal", 0.7, 1.6))

normal <- replicate(S, sampleB_SoR(
  matrix(0, 1, 1), matrix(1 / 1.6, 1, 1), 0.7, 0, 1,
  matrix(0, 1, 0), matrix(0, 1, 0), 0L))
moments <- rbind(moments, moment_summary(normal, "C++ Normal", 0.7, 1.6))
stopifnot(all(abs(moments$mean_error) < 5 * moments$mean_mcse),
          all(abs(moments$variance / moments$expected_var - 1) < 0.06),
          all(abs(moments$adjacent_correlation) < 0.04),
          all(moments$duplicate_draws == 0))
write.csv(moments, file.path(out_dir, "moments.csv"), row.names = FALSE)
print(moments, digits = 6)

# Keep every advertised sampler configuration available. The basic matrix
# covers all four data models with spatial fields, traits and factors; the
# extras cover splines and zero-factor/intercept-only branches.
configs <- list(binary = list(model = "binary"), continuous = list(model = "continuous"),
                occupancy = list(model = "occupancy"), two_stage = list(model = "two_stage"),
                splines = list(model = "two_stage", splines = TRUE),
                intercept_only = list(model = "two_stage", intercept_only = TRUE))
fit_checks <- lapply(names(configs), function(label) {
  cfg <- configs[[label]]
  sim <- simulate_fixture(model = cfg$model, P = if (cfg$model == "two_stage") 2L else 1L)
  empty <- isTRUE(cfg$intercept_only)
  if (empty) sim$data_list$traits <- NULL
  fit_once <- function(threads) {
    RcppParallel::setThreadOptions(numThreads = threads)
    set.seed(12831)
    fit <- suppressMessages(runOccJSDM(
      sim$data_list,
      occCovariates = if (empty) NULL else fixture_occ_covariates(),
      collCovariates = if (empty || !cfg$model %in% c("occupancy", "two_stage")) NULL else "X_theta",
      spatCovariates = if (empty) NULL else fixture_spat_covariates(),
      listParams = list(n_factors = if (empty) 0 else 2, n_supportpoints = FIXTURE_KNOTS,
                        splineVars = isTRUE(cfg$splines)),
      MCMCparams = list(nchain = 2, nburn = 5, niter = 8, nthin = 1)))
    stopifnot(identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS"), as.character(threads)))
    list(results = fit$results_output, seed = .Random.seed)
  }
  one <- fit_once(1L)
  four <- fit_once(4L)
  repeat_four <- fit_once(4L)
  stopifnot(identical(one, four), identical(four, repeat_four))
  data.frame(configuration = label, same_1_4 = identical(one, four),
             repeat_4 = identical(four, repeat_four))
})
fit_checks <- do.call(rbind, fit_checks)
capture.output(warnings(), file = file.path(out_dir, "fit-warnings.txt"))
write.csv(fit_checks, file.path(out_dir, "fit-checks.csv"), row.names = FALSE)
print(fit_checks)

if (exists("rngAuditCounts", mode = "function")) {
  counts <- unlist(rngAuditCounts())
  stopifnot(all(c("cpp runif", "cpp rnorm", "R::rbeta", "R::rbinom",
                  "R::rgamma", "arma::randn") %in% names(counts)))
  write.csv(data.frame(leaf = names(counts), calls = as.numeric(counts)),
            file.path(out_dir, "main-thread-counts.csv"), row.names = FALSE)
  print(counts)
}
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
if (is.na(old_threads)) Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS") else
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old_threads)
