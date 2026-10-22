# Reproduce the beta quickstart's shipped fit from an isolated installation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
repo <- normalizePath(args[[1]], mustWork = TRUE)
lib <- normalizePath(args[[2]], mustWork = TRUE)
evidence <- normalizePath(args[[3]], mustWork = TRUE)
stopifnot(!file.exists(file.path(evidence, "provenance.rds")))
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
RcppParallel::setThreadOptions(numThreads = 1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seed <- 20261007L
settings <- list(
  listParams = list(n_factors = 2),
  threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = c("X_theta.1", "X_theta.2"),
  spatCovariates = c("Xs.1", "Xs.2"),
  MCMCparams = list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1),
  summarisedLatentPresences = TRUE,
  listPriors = list()
)
source_files <- system2("git", c("-C", shQuote(repo), "ls-files", "R", "src", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)
source_hashes <- tools::md5sum(file.path(repo, source_files))
names(source_hashes) <- source_files
old_fit_hash <- unname(tools::md5sum(file.path(repo, "data/sampleresults.rda")))
input <- new.env(parent = emptyenv())
load(file.path(repo, "data/sampledata.rda"), envir = input)
stopifnot(is.list(input$sampledata), ncol(input$sampledata$OTU) == 10L)
warnings_seen <- character()
started <- Sys.time()
set.seed(seed)
sampleresults <- withCallingHandlers(
  do.call(runOccJSDM, c(list(data = input$sampledata), settings)),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
  }
)
finished <- Sys.time()
stopifnot(
  identical(sampleresults$infos$model, "two_stage"),
  sampleresults$infos$n == 100L,
  sampleresults$infos$S == 10L,
  sampleresults$infos$n_factors == 2L,
  sampleresults$infos$threshold == 1,
  sampleresults$infos$intercept_prior$sd == 1,
  identical(dim(sampleresults$results_output$jsdm_output$B_output), c(2L, 10L, 5000L, 2L)),
  all(is.finite(sampleresults$results_output$psi_output)),
  all(sampleresults$results_output$psi_output >= 0 & sampleresults$results_output$psi_output <= 1),
  identical(source_hashes, setNames(tools::md5sum(file.path(repo, source_files)), source_files))
)
diagnostics <- returnConvergenceDiagnostics(sampleresults)
write.csv(diagnostics, file.path(evidence, "diagnostics.csv"), row.names = FALSE)
candidate <- tempfile(".sampleresults-", tmpdir = file.path(repo, "data"), fileext = ".rda")
save(sampleresults, file = candidate, compress = "xz")
stopifnot(file.rename(candidate, file.path(repo, "data/sampleresults.rda")))
provenance <- list(
  revision = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE),
  source_hashes = source_hashes,
  input_md5 = unname(tools::md5sum(file.path(repo, "data/sampledata.rda"))),
  previous_fit_md5 = old_fit_hash,
  fit_md5 = unname(tools::md5sum(file.path(repo, "data/sampleresults.rda"))),
  seed = seed, rng_kind = RNGkind(), settings = settings,
  requested_threads = list(RcppParallel = 1L, environment = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"))),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
  warnings = warnings_seen, installed_library = lib, session = sessionInfo()
)
saveRDS(provenance, file.path(evidence, "provenance.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(evidence, "session-info.txt"))
flagged <- diagnostics[!is.finite(diagnostics$rhat) | diagnostics$rhat > 1.01 | !is.finite(diagnostics$ess) | diagnostics$ess < 400, ]
write.csv(flagged, file.path(evidence, "diagnostics-flagged.csv"), row.names = FALSE)
cat("REFIT COMPLETE: ", nrow(diagnostics), " coefficient diagnostics; ", nrow(flagged), " flagged at Rhat > 1.01 or ESS < 400.\n", sep = "")
cat("Elapsed seconds: ", provenance$elapsed_seconds, "\n", sep = "")
