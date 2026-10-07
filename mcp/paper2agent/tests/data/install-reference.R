# Install from this checkout and bind the installed artifacts to source hashes.
# Usage: Rscript install-reference.R NEW_R_LIBRARY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new, empty R library directory.")
stopifnot(requireNamespace("jsonlite", quietly = TRUE), requireNamespace("digest", quietly = TRUE))
script <- normalizePath(sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())]))
root <- normalizePath(file.path(dirname(script), "../../../.."))
lib <- args[[1]]
if (file.exists(lib)) stop("Choose a new library directory; existing environments are preserved.")
dir.create(lib, recursive = TRUE)
lib <- normalizePath(lib)
source_files <- c(list.files(file.path(root, "R"), full.names = TRUE),
  list.files(file.path(root, "src"), pattern = "\\.(cpp|h)$|^Makevars", full.names = TRUE),
  file.path(root, "DESCRIPTION"), file.path(root, "NAMESPACE"),
  file.path(root, "tests/testthat/helper-fixtures.R"))
source_hashes <- function() setNames(lapply(source_files, function(path)
  digest::digest(file = path, algo = "sha256")), substring(source_files, nchar(root) + 2L))
before <- source_hashes()
log <- file.path(lib, "install.log")
status <- system2(file.path(R.home("bin"), "R"),
  shQuote(c("CMD", "INSTALL", "--preclean", "--no-multiarch", paste0("--library=", lib), root)),
  stdout = log, stderr = log)
if (status != 0L) stop("R CMD INSTALL failed; see ", log)
if (!identical(before, source_hashes())) stop("Source changed during installation; no provenance was published.")
package <- normalizePath(file.path(lib, "occJSDM"))
artifacts <- list.files(package, recursive = TRUE, full.names = TRUE, all.files = TRUE)
artifacts <- artifacts[!dir.exists(artifacts)]
artifact_hashes <- setNames(lapply(artifacts, function(path) digest::digest(file = path, algo = "sha256")),
  substring(artifacts, nchar(package) + 2L))
record <- list(source_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  source_sha256 = before, installed_sha256 = artifact_hashes, R_version = R.version.string,
  build_commands = list(install = "R CMD INSTALL --preclean --no-multiarch",
    compiler = system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE),
    openmp_flags = system2(file.path(R.home("bin"), "R"), c("CMD", "config", "SHLIB_OPENMP_CXXFLAGS"), stdout = TRUE)))
jsonlite::write_json(record, file.path(lib, "occJSDM-source.json"), auto_unbox = TRUE,
                     pretty = TRUE, null = "null")
cat("Installed pinned reference package and recorded build provenance in", lib, "\n")
