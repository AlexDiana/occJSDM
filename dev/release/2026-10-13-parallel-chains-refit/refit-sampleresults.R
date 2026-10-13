# Refit the quickstart's shipped example under per-chain seeding, twice (one
# core and two cores), require identity apart from the execution record, and
# ship the two-core object. Rscript refit-sampleresults.R REPO LIB EVIDENCE BLAS_EVIDENCE_RDS
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
repo <- normalizePath(args[[1]], mustWork = TRUE)
lib <- normalizePath(args[[2]], mustWork = TRUE)
evidence <- normalizePath(args[[3]], mustWork = TRUE)
blas_json <- normalizePath(args[[4]], mustWork = TRUE)
stopifnot(!file.exists(file.path(evidence, "provenance.rds")), nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")))
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seed <- 20261007L
settings <- list(
  listParams = list(n_factors = 2), threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = c("X_theta.1", "X_theta.2"),
  spatCovariates = c("Xs.1", "Xs.2"),
  summarisedLatentPresences = TRUE, listPriors = list())
mcmc <- list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1)
source_files <- system2("git", c("-C", shQuote(repo), "ls-files", "R", "src", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)
source_hashes <- tools::md5sum(file.path(repo, source_files)); names(source_hashes) <- source_files
old_fit_hash <- unname(tools::md5sum(file.path(repo, "data/sampleresults.rda")))
input <- new.env(parent = emptyenv()); load(file.path(repo, "data/sampledata.rda"), envir = input)
stopifnot(is.list(input$sampledata), ncol(input$sampledata$OTU) == 10L)
fit_at <- function(ncores) {
  warnings_seen <- character(); started <- Sys.time()
  set.seed(seed)
  fit <- withCallingHandlers(
    do.call(runOccJSDM, c(list(data = input$sampledata), settings, list(MCMCparams = c(mcmc, list(ncores = ncores))))),
    warning = function(w) warnings_seen <<- c(warnings_seen, conditionMessage(w)))
  list(fit = fit, warnings = warnings_seen, started = started, finished = Sys.time(),
       seed_state = get(".Random.seed", envir = globalenv()))
}
one <- fit_at(1L)
two <- fit_at(2L)
stopifnot(identical(occJSDM:::drop_execution(one$fit), occJSDM:::drop_execution(two$fit)),
          identical(one$seed_state, two$seed_state),
          two$fit$infos$execution$workers == 2L, two$fit$infos$execution$threads_per_worker == 1L,
          one$fit$infos$execution$workers == 1L)
sampleresults <- two$fit
stopifnot(
  identical(sampleresults$infos$model, "two_stage"), sampleresults$infos$n == 100L,
  sampleresults$infos$S == 10L, sampleresults$infos$n_factors == 2L, sampleresults$infos$threshold == 1,
  identical(dim(sampleresults$results_output$jsdm_output$B_output), c(2L, 10L, 5000L, 2L)),
  all(is.finite(sampleresults$results_output$psi_output)),
  identical(source_hashes, stats::setNames(tools::md5sum(file.path(repo, source_files)), source_files)))
diagnostics <- returnConvergenceDiagnostics(sampleresults)
write.csv(diagnostics, file.path(evidence, "diagnostics.csv"), row.names = FALSE)
candidate <- tempfile(".sampleresults-", tmpdir = file.path(repo, "data"), fileext = ".rda")
save(sampleresults, file = candidate, compress = "xz")
stopifnot(file.rename(candidate, file.path(repo, "data/sampleresults.rda")))
file.copy(blas_json, file.path(evidence, "blas-evidence.rds"), overwrite = TRUE)
provenance <- list(
  revision = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE),
  source_hashes = source_hashes,
  input_md5 = unname(tools::md5sum(file.path(repo, "data/sampledata.rda"))),
  previous_fit_md5 = old_fit_hash,
  fit_md5 = unname(tools::md5sum(file.path(repo, "data/sampleresults.rda"))),
  seed = seed, rng_kind = RNGkind(), settings = settings, MCMCparams = c(mcmc, list(ncores = 2L)),
  chain_seeds = sampleresults$infos$chain_seeds,
  execution_one_core = one$fit$infos$execution, execution_two_cores = two$fit$infos$execution,
  identical_apart_from_execution = TRUE,
  blas_evidence_file = "blas-evidence.rds", blas_declaration_basis = "blas-declaration-basis.txt",
  elapsed_seconds = c(one_core = as.numeric(difftime(one$finished, one$started, units = "secs")),
                      two_cores = as.numeric(difftime(two$finished, two$started, units = "secs"))),
  warnings = list(one_core = one$warnings, two_cores = two$warnings),
  installed_library = lib, session = sessionInfo())
saveRDS(provenance, file.path(evidence, "provenance.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(evidence, "session-info.txt"))
flagged <- diagnostics[!is.finite(diagnostics$rhat) | diagnostics$rhat > 1.01 | !is.finite(diagnostics$ess) | diagnostics$ess < 400, ]
write.csv(flagged, file.path(evidence, "diagnostics-flagged.csv"), row.names = FALSE)
cat("REFIT COMPLETE: ", nrow(diagnostics), " coefficient diagnostics; ", nrow(flagged), " flagged.\n", sep = "")
