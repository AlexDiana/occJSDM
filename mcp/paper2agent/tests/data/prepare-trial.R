# Prepare synthetic inputs and grader references through the existing R API.
# This is trial preparation, not an MCP server or a model evaluation.
# Usage: Rscript prepare-trial.R NEW_OUTPUT_DIRECTORY PINNED_R_LIBRARY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply a new output directory and the pinned R library.")
lib <- normalizePath(args[[2]])
expected_package <- file.path(lib, "occJSDM")
if (!dir.exists(expected_package)) stop("Pinned occJSDM installation is missing from the requested library.")
.libPaths(c(lib, .libPaths()))
if (!identical(normalizePath(find.package("occJSDM")), normalizePath(expected_package))) {
  stop("occJSDM did not resolve to the pinned library.")
}
stopifnot(requireNamespace("jsonlite", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))
script <- normalizePath(sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())]))
root <- normalizePath(file.path(dirname(script), "../../../.."))
source_files <- c(list.files(file.path(root, "R"), full.names = TRUE),
  list.files(file.path(root, "src"), pattern = "\\.(cpp|h)$|^Makevars", full.names = TRUE),
  file.path(root, "DESCRIPTION"), file.path(root, "NAMESPACE"),
  file.path(root, "tests/testthat/helper-fixtures.R"))
hashes <- setNames(lapply(source_files, function(path) digest::digest(file = path, algo = "sha256")),
                  substring(source_files, nchar(root) + 2L))
build_path <- file.path(lib, "occJSDM-source.json")
if (!file.exists(build_path)) stop("Pinned build provenance is missing; run install-reference.R first.")
build <- jsonlite::read_json(build_path)
if (!identical(build$source_sha256, hashes)) stop("Pinned source hashes do not match this checkout.")
for (name in names(build$installed_sha256)) {
  artifact <- file.path(expected_package, name)
  if (!file.exists(artifact) || !identical(digest::digest(file = artifact, algo = "sha256"),
                                         build$installed_sha256[[name]])) {
    stop("Pinned installed package artifact changed: ", name)
  }
}
library(occJSDM)
out <- normalizePath(args[[1]], mustWork = FALSE)
if (file.exists(out)) stop("Output directory already exists; choose a new directory.")
student <- file.path(out, "student")
grader <- file.path(out, "grader")
dir.create(student, recursive = TRUE)
dir.create(grader, recursive = TRUE)

fixtures <- new.env()
source(file.path(root, "tests/testthat/helper-fixtures.R"), local = fixtures)
sim <- fixtures$simulate_fixture(useSpatField = FALSE, seed = 1701L)
data <- sim$data_list
saveRDS(data, file.path(student, "data.rds"))
no_traits <- data
no_traits$traits <- NULL
saveRDS(no_traits, file.path(student, "no-traits.rds"))
malformed <- data
malformed$OTU <- malformed$OTU[-nrow(malformed$OTU), , drop = FALSE]
saveRDS(malformed, file.path(student, "malformed.rds"))
counts <- fixtures$simulate_fixture(model = "binary", useSpatField = FALSE,
                                    seed = 1701L)$data_list
counts$OTU <- counts$OTU * 5 + 1
saveRDS(counts, file.path(student, "unreplicated-counts.rds"))

settings <- list(listParams = list(n_factors = 2L), threshold = 1L,
                 occCovariates = fixtures$fixture_occ_covariates(),
                 collCovariates = "X_theta", spatCovariates = NULL,
                 MCMCparams = fixtures$FIXTURE_MCMC,
                 summarisedLatentPresences = TRUE)
set.seed(1702L)
start <- proc.time()
fit <- do.call(runOccJSDM, c(list(data = data), settings))
elapsed <- unname((proc.time() - start)[["elapsed"]])
saveRDS(fit, file.path(student, "fit-short.rds"))
diagnostics <- returnConvergenceDiagnostics(fit)
write.csv(diagnostics, file.path(grader, "diagnostics.csv"), row.names = FALSE, na = "")
saveRDS(diagnostics, file.path(grader, "diagnostics.rds"))
draws <- returnOccupancyRates(fit)
summary <- data.frame(species = colnames(draws), mean = colMeans(draws),
                      q2.5 = apply(draws, 2, quantile, 0.025),
                      q97.5 = apply(draws, 2, quantile, 0.975), row.names = NULL)
write.csv(summary, file.path(grader, "baseline-occupancy.csv"), row.names = FALSE)

# A genuine one-chain fit supplies an unavailable-Rhat interpretation exercise.
# This is not a fitting request accepted by the planned two-chain MCP interface.
one_chain_settings <- settings
one_chain_settings$MCMCparams$nchain <- 1L
set.seed(1704L)
one_chain <- do.call(runOccJSDM, c(list(data = data), one_chain_settings))
unavailable <- returnConvergenceDiagnostics(one_chain)
write.csv(unavailable, file.path(student, "diagnostics-unavailable.csv"),
          row.names = FALSE, na = "")
stopifnot(all(is.na(unavailable$rhat)))

reference <- list(schema_version = 1L, status = "local_references_prepared",
  source_commit = build$source_commit,
  source_sha256 = hashes, package_version = as.character(packageVersion("occJSDM")),
  package_path = find.package("occJSDM"), R_version = R.version.string,
  build_provenance_sha256 = digest::digest(file = build_path, algo = "sha256"),
  build_configuration = build$build_commands,
  thread_environment = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "RCPP_PARALLEL_NUM_THREADS"), unset = NA)),
  RcppParallel_default_threads = RcppParallel::defaultNumThreads(),
  BLAS = extSoftVersion()[["BLAS"]],
  generator_sha256 = digest::digest(file = script, algo = "sha256"),
  input_sha256 = digest::digest(file = file.path(student, "data.rds"), algo = "sha256"),
  simulation_seed = 1701L, fit_seed = 1702L, unavailable_fit_seed = 1704L,
  fit_settings = settings, unavailable_fit_settings = one_chain_settings,
  fit_elapsed_seconds = elapsed, scientific_interpretation = FALSE,
  warning = "20 burn-in and 20 kept iterations are only for execution tests, not scientific interpretation.",
  dimensions = list(rows = nrow(data$info), species = ncol(data$OTU)),
  inferred_model = getFromNamespace("inferDataModel", "occJSDM")(data),
  diagnostics = list(rows = nrow(diagnostics), rhat_flagged = sum(diagnostics$rhat > 1.01, na.rm = TRUE),
                     ess_flagged = sum(diagnostics$ess < 400, na.rm = TRUE),
                     rhat_unavailable = sum(is.na(diagnostics$rhat)),
                     ess_unavailable = sum(is.na(diagnostics$ess))),
  cloud_trials = "not_run", model_trials = "not_run")
jsonlite::write_json(reference, file.path(grader, "reference.json"), auto_unbox = TRUE,
                     pretty = TRUE, null = "null", na = "null", digits = NA)
writeLines(capture.output(sessionInfo()), file.path(grader, "session-info.txt"))

# Export an explicit source/lesson inventory, not the whole development repo.
# The trial protocol, expected answers, tests and grader outputs stay outside it.
project <- file.path(out, "student-project")
dir.create(project)
inventory <- c("DESCRIPTION", "NAMESPACE", "README.md", "LICENSE", "CITATION.cff",
  "occJSDM.Rproj", ".Rbuildignore",
  file.path("R", list.files(file.path(root, "R"), pattern = "\\.R$", full.names = FALSE)),
  file.path("src", list.files(file.path(root, "src"), pattern = "\\.(cpp|h)$|^Makevars", full.names = FALSE)),
  file.path("man", list.files(file.path(root, "man"), pattern = "\\.Rd$", full.names = FALSE)),
  "data/sampledata.rda", "data/sampleresults.rda", "vignettes/occJSDM.Rmd",
  "vignettes/lesson-links.R", "vignettes/teaching.css")
inventory <- inventory[file.exists(file.path(root, inventory))]
for (rel in inventory) {
  dest <- file.path(project, rel)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(root, rel), dest)) stop("Could not copy source file: ", rel)
}
dir.create(file.path(project, "trial"))
if (!all(file.copy(list.files(student, full.names = TRUE), file.path(project, "trial")))) {
  stop("Could not copy trial inputs.")
}
dir.create(file.path(project, "teaching"))
for (name in c("assistant-instructions.md", "glossary.md")) {
  if (!file.copy(file.path(root, "mcp/paper2agent/class", name), file.path(project, "teaching", name))) {
    stop("Could not copy teaching file: ", name)
  }
}
project_files <- list.files(project, recursive = TRUE, full.names = TRUE, all.files = TRUE)
project_files <- project_files[!dir.exists(project_files)]
project_hashes <- setNames(lapply(project_files, function(path) digest::digest(file = path, algo = "sha256")),
                           substring(project_files, nchar(project) + 2L))
jsonlite::write_json(list(status = "draft_preflight_project_pending_teaching_review",
  files_sha256 = project_hashes, excluded = c("tests", "dev", ".git", ".superpowers", "mcp", "grader")),
  file.path(grader, "student-project-manifest.json"), auto_unbox = TRUE, pretty = TRUE)

cat("Prepared synthetic trial inputs and native R references in", out, "\n")
