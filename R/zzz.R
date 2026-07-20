.onLoad <- function(libname, pkgname) {
  # Limit default parallelism for CRAN compliance
  options(mc.cores = min(2L, parallel::detectCores()))
}


