# Explicit activation for Rscript --vanilla. Callers must bind the project path.
p2a_project <- Sys.getenv("P2A_R_PROJECT", unset = "")
if (!nzchar(p2a_project)) {
  p2a_root <- Sys.getenv("P2A_PROJECT_ROOT", unset = "")
  if (nzchar(p2a_root)) p2a_project <- file.path(p2a_root,"r-runtime")
  else p2a_project <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE))
}
p2a_project <- normalizePath(p2a_project, mustWork = TRUE)
Sys.setenv(RENV_PROJECT = p2a_project, RENV_PATHS_ROOT = file.path(p2a_project,"state"), RENV_PATHS_LIBRARY = file.path(p2a_project,"library"), RENV_PATHS_CACHE = file.path(p2a_project,"cache"), RENV_PATHS_SANDBOX = file.path(p2a_project,"sandbox"), RENV_CONFIG_CACHE_ENABLED = "FALSE")
p2a_boot <- file.path(p2a_project,"bootstrap-library")
if (!dir.exists(file.path(p2a_boot,"renv"))) stop("Project-owned renv bootstrap is missing")
.libPaths(c(p2a_boot, .Library))
source(file.path(p2a_project,"renv","activate.R"))
p2a_lib <- renv::paths$library(project = p2a_project)
if (!dir.exists(file.path(p2a_lib,"occJSDM"))) warning("Pinned occJSDM source has not yet been installed")
if (!startsWith(find.package("renv"),p2a_project)) stop("renv resolved outside the selected project")
if (any(grepl("/Users/[^/]+/Library/R", .libPaths()))) stop("User R library fallback is prohibited")
rm(p2a_boot, p2a_lib)
