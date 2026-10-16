# Are fits that do not use theta_intercept_prior = "community" bit-identical
# before and after the change? Covers default fits and B0-only community fits.
# Usage: Rscript check_default_identity.R <old package root> <new package root> <output csv>
# Fits the same seeded models under each build, each in its own R process, and
# compares results_output and infos (less infos$execution, which records
# timings and process ids). The roots are exports of the branch base and of
# the change, for example: git archive cf9c133 | tar -x -C <old root>
# The old code ignores listPriors names it does not know, so the explicit
# theta_intercept_prior = "fixed" configuration is also its default fit, and
# control_theta_community, which uses the new option, must differ: it shows
# that the comparison can fail.
args <- commandArgs(trailingOnly = TRUE)
old_root <- normalizePath(args[1], mustWork = TRUE)
new_root <- normalizePath(args[2], mustWork = TRUE)
out_csv <- args[3]

fit_all <- function(root, file) {
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('suppressMessages(devtools::load_all("%s", quiet = TRUE))', root),
    sprintf('source(file.path("%s", "tests", "testthat", "helper-fixtures.R"))', root),
    'Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1"); RcppParallel::setThreadOptions(numThreads = 1)',
    'mc <- list(nchain = 2L, nburn = 30L, niter = 30L, nthin = 1L)',
    'fit <- function(model, spatial, priors = list(), P = 1L) {',
    '  sim <- simulate_fixture(model = model, P = P)',
    '  set.seed(99)',
    '  suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,',
    '    occCovariates = fixture_occ_covariates(),',
    '    collCovariates = if (model %in% c("occupancy", "two_stage")) "X_theta" else NULL,',
    '    spatCovariates = if (spatial) fixture_spat_covariates() else NULL,',
    '    listParams = c(list(n_factors = 2L), if (spatial) list(n_supportpoints = FIXTURE_KNOTS)),',
    '    listPriors = priors, MCMCparams = mc)))',
    '}',
    'fits <- list(',
    '  two_stage_spatial = fit("two_stage", TRUE, P = 2L),',
    '  two_stage_plain = fit("two_stage", FALSE),',
    '  occupancy_spatial = fit("occupancy", TRUE),',
    '  occupancy_plain = fit("occupancy", FALSE),',
    '  binary_spatial = fit("binary", TRUE),',
    '  binary_plain = fit("binary", FALSE),',
    '  continuous_plain = fit("continuous", FALSE),',
    '  two_stage_sigma_b0_2 = fit("two_stage", TRUE, list(sigma_b0 = 2)),',
    '  two_stage_theta_fixed = fit("two_stage", TRUE, list(theta_intercept_prior = "fixed")),',
    '  two_stage_b0_community = fit("two_stage", TRUE, list(b0_prior = "community"), P = 2L),',
    '  occupancy_b0_community = fit("occupancy", FALSE, list(b0_prior = "community", mu_b0_sd = 3)),',
    '  binary_b0_community = fit("binary", TRUE, list(b0_prior = "community")),',
    '  control_theta_community = fit("two_stage", TRUE, list(theta_intercept_prior = "community")))',
    'fits <- lapply(fits, function(f) { f$infos$execution <- NULL; f })',
    sprintf('saveRDS(fits, "%s")', file)), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (status != 0) stop("fits failed under ", root)
  readRDS(file)
}
old <- fit_all(old_root, tempfile(fileext = ".rds"))
new <- fit_all(new_root, tempfile(fileext = ".rds"))
res <- data.frame(
  config = names(old),
  results_identical = vapply(names(old), function(k)
    identical(old[[k]]$results_output, new[[k]]$results_output), logical(1)),
  infos_identical = vapply(names(old), function(k)
    identical(old[[k]]$infos, new[[k]]$infos), logical(1)),
  whole_fit_identical = vapply(names(old), function(k) identical(old[[k]], new[[k]]), logical(1)),
  row.names = NULL)
print(res)
write.csv(res, out_csv, row.names = FALSE)
