# Refit a shipped example fit at the default intercept prior, which since
# 1b4a2cd is b0_prior = "community" for these two-stage fits. Seeds, data,
# arguments and MCMC settings are those of the previous build scripts:
#   sampleresults:     dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R
#   quickstartresults: dev/release/2026-11-01-quickstart-data/fit-quickstartresults.R
# b0_prior is deliberately not passed. The reference-BLAS identity gate of the
# 13 October script is dropped: the laptop now runs Accelerate (14 October).
#
# Run from the repository root:
#   Rscript dev/release/2026-11-04-community-prior-refit/refit-shipped.R sampleresults|quickstartresults LIB
# Replaces data/<name>.rda and writes <name>-diagnostics.csv and
# <name>-provenance.rds beside this script.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[[1]] %in% c("sampleresults", "quickstartresults"))
name <- args[[1]]
lib <- normalizePath(args[[2]], mustWork = TRUE)
here <- "dev/release/2026-11-04-community-prior-refit"
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

source_files <- system2("git", c("ls-files", "R", "src", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)
source_hashes <- stats::setNames(tools::md5sum(source_files), source_files)
rda <- file.path("data", paste0(name, ".rda"))
old_fit_md5 <- unname(tools::md5sum(rda))

if (name == "sampleresults") {
  input_file <- "data/sampledata.rda"
  seed <- 20261007L
  input <- new.env(parent = emptyenv()); load(input_file, envir = input)
  call_args <- list(
    data = input$sampledata,
    listParams = list(n_factors = 2), threshold = 1,
    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
    collCovariates = c("X_theta.1", "X_theta.2"),
    spatCovariates = c("Xs.1", "Xs.2"),
    summarisedLatentPresences = TRUE, listPriors = list(),
    MCMCparams = list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1, ncores = 2))
} else {
  input_file <- "data/quickstartdata.rda"
  seed <- 20261101L
  input <- new.env(parent = emptyenv()); load(input_file, envir = input)
  call_args <- list(
    data = input$quickstartdata,
    listParams = list(n_factors = 2), threshold = 1,
    occCovariates = c("elevation", "hmi", "protected"),
    collCovariates = c("water_volume"),
    MCMCparams = list(nchain = 2, nburn = 2000, niter = 2000, nthin = 2, ncores = 2),
    summarisedLatentPresences = TRUE)
}

warnings_seen <- character()
started <- Sys.time()
set.seed(seed)
fit <- withCallingHandlers(
  do.call(runOccJSDM, call_args),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
finished <- Sys.time()
stopifnot(identical(fit$infos$model, "two_stage"),
          !is.null(fit$results_output$jsdm_output$mu_b0_output),
          identical(source_hashes, stats::setNames(tools::md5sum(source_files), source_files)))

assign(name, fit)
candidate <- tempfile(paste0(".", name, "-"), tmpdir = "data", fileext = ".rda")
save(list = name, file = candidate, compress = "xz")
stopifnot(file.rename(candidate, rda))

diagnostics <- returnConvergenceDiagnostics(fit)
write.csv(diagnostics, file.path(here, paste0(name, "-diagnostics.csv")), row.names = FALSE)
provenance <- list(
  revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  dirty = system2("git", c("status", "--porcelain", "R", "src"), stdout = TRUE),
  source_hashes = source_hashes,
  input_md5 = unname(tools::md5sum(input_file)),
  previous_fit_md5 = old_fit_md5,
  fit_md5 = unname(tools::md5sum(rda)),
  seed = seed, rng_kind = RNGkind(),
  settings = call_args[setdiff(names(call_args), "data")],
  intercept_prior = fit$infos$intercept_prior, b0_shift = fit$infos$b0_shift,
  chain_seeds = fit$infos$chain_seeds,
  execution = fit$infos$execution,
  started = started, finished = finished,
  minutes = as.numeric(difftime(finished, started, units = "mins")),
  warnings = unique(warnings_seen),
  power = system2("pmset", c("-g", "batt"), stdout = TRUE),
  low_power_mode = grep("lowpowermode", system2("pmset", "-g", stdout = TRUE), value = TRUE),
  blas = extSoftVersion()[["BLAS"]], la_library = La_library(),
  r_version = R.version.string, installed_library = lib)
saveRDS(provenance, file.path(here, paste0(name, "-provenance.rds")))

flagged <- diagnostics[!is.finite(diagnostics$rhat) | diagnostics$rhat > 1.01 |
                         !is.finite(diagnostics$ess) | diagnostics$ess < 400,
                       c("param", "label1", "label2", "rhat", "ess")]
cat("minutes:", round(provenance$minutes, 2), "\n")
cat("rows:", nrow(diagnostics), " flagged:", nrow(flagged), "\n")
print(flagged[order(-flagged$rhat), ], row.names = FALSE)
cat("mu_b0:\n"); print(diagnostics[diagnostics$param == "mu_b0", ], row.names = FALSE)
cat("warnings:\n"); cat(unique(warnings_seen), sep = "\n")
cat("REFIT COMPLETE\n")
