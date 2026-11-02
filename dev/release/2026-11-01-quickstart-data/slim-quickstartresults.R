# Turn a fit saved by fit-quickstartresults.R into data/quickstartresults.rda,
# keeping every THIN-th stored draw so the shipped file stays small.
#
# Run from the repository root:
#   Rscript dev/release/2026-11-01-quickstart-data/slim-quickstartresults.R FIT.rds THIN
# Also writes diagnostics.csv, warnings.txt and provenance.rds beside this script.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
saved <- readRDS(args[[1]])
thin <- as.integer(args[[2]])
here <- "dev/release/2026-11-01-quickstart-data"
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

fit <- saved$fit
niter <- fit$infos$MCMCparams$niter
nchain <- fit$infos$MCMCparams$nchain
keep <- seq(thin, niter, by = thin)

# Posterior draw arrays end in (iteration, chain); thin the iteration axis.
thin_draws <- function(x) {
  d <- dim(x)
  if (is.null(d) || length(d) < 2L || d[length(d)] != nchain || d[length(d) - 1L] != niter) return(x)
  idx <- rep(list(TRUE), length(d))
  idx[[length(d) - 1L]] <- keep
  do.call(`[`, c(list(x), idx, list(drop = FALSE)))
}
if (thin > 1L) {
  ro <- fit$results_output
  for (nm in names(ro)) if (nm != "jsdm_output") ro[[nm]] <- thin_draws(ro[[nm]])
  ro$jsdm_output <- lapply(ro$jsdm_output, thin_draws)
  fit$results_output <- ro
  fit$infos$MCMCparams$niter <- length(keep)
  fit$infos$MCMCparams$shipped_thinning <- thin
}

quickstartresults <- fit
candidate <- tempfile(".quickstartresults-", tmpdir = "data", fileext = ".rda")
save(quickstartresults, file = candidate, compress = "xz")
stopifnot(file.rename(candidate, "data/quickstartresults.rda"))

write.csv(returnConvergenceDiagnostics(quickstartresults), file.path(here, "diagnostics.csv"), row.names = FALSE)
writeLines(unique(saved$warnings), file.path(here, "warnings.txt"))
provenance <- saved$provenance
provenance$session <- NULL  # sessionInfo() objects do not survive other R versions well
provenance$shipped_thinning <- thin
saveRDS(provenance, file.path(here, "provenance.rds"))
cat("size MB:", round(file.size("data/quickstartresults.rda") / 1e6, 1), "\n")
