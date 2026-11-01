# Fit the Quickstart's example survey (quickstartdata) without a spatial field
# and with the default priors, the Quickstart's live-fit call, and save the
# fit with its provenance.
#
# Run from the repository root:
#   Rscript dev/release/2026-11-01-quickstart-data/fit-quickstartresults.R NBURN NITER NTHIN OUT.rds
# OUT.rds receives list(fit, warnings, provenance). slim-quickstartresults.R
# turns a chosen fit into data/quickstartresults.rda.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
nburn <- as.integer(args[[1]]); niter <- as.integer(args[[2]]); nthin <- as.integer(args[[3]])
out <- args[[4]]
stopifnot(!file.exists(out))

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
load("data/quickstartdata.rda")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

warnings_seen <- character()
started <- Sys.time()
set.seed(20261101)
fit <- withCallingHandlers(
  runOccJSDM(
    data = quickstartdata,
    listParams = list(n_factors = 2),
    threshold = 1,
    occCovariates = c("elevation", "hmi", "protected"),
    collCovariates = c("water_volume"),
    MCMCparams = list(nchain = 2, nburn = nburn, niter = niter, nthin = nthin, ncores = 2),
    summarisedLatentPresences = TRUE
  ),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
finished <- Sys.time()

provenance <- list(
  revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  dirty = system2("git", c("status", "--porcelain", "R", "src", "data"), stdout = TRUE),
  settings = list(nburn = nburn, niter = niter, nthin = nthin, nchain = 2, ncores = 2, seed = 20261101L),
  started = started, finished = finished,
  minutes = as.numeric(difftime(finished, started, units = "mins")),
  power = system2("pmset", c("-g", "batt"), stdout = TRUE),
  low_power_mode = grep("lowpowermode", system2("pmset", "-g", stdout = TRUE), value = TRUE),
  blas = extSoftVersion()[["BLAS"]], la_library = La_library(),
  session = utils::sessionInfo())
saveRDS(list(fit = fit, warnings = warnings_seen, provenance = provenance), out)

diag <- returnConvergenceDiagnostics(fit)
flagged <- diag[diag$rhat > 1.01 | diag$ess < 400, c("param", "label1", "label2", "rhat", "ess")]
cat("minutes:", round(provenance$minutes, 1), "\n")
cat("rows:", nrow(diag), " flagged:", nrow(flagged), "\n")
print(flagged[order(-flagged$rhat), ], row.names = FALSE)
cat("warnings:\n"); cat(unique(warnings_seen), sep = "\n")
