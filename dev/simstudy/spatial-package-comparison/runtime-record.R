record_toolchain <- function(run_root,code) {
  private_R <- file.path(R.home(),"bin/R")
  config_names <- c("CC","CXX","CXXFLAGS","SHLIB_OPENMP_CXXFLAGS","LDFLAGS")
  flags <- setNames(lapply(config_names,function(nm) system2(private_R,c("CMD","config",nm),stdout=TRUE,stderr=TRUE)),config_names)
  compiler <- system2("/usr/bin/clang","--version",stdout=TRUE,stderr=TRUE)
  python <- Sys.getenv("RETICULATE_PYTHON"); stopifnot(nzchar(python))
  py <- system2(python,shQuote(file.path(code,"runtime-probe.py")),stdout=TRUE,stderr=FALSE)
  if(!is.null(attr(py,"status"))) stop("Python runtime probe failed")
  py <- jsonlite::fromJSON(paste(py,collapse="\n")); stopifnot(!py$cuda_available)
  binary <- c(file.path(R.home(),"bin/exec/R"),file.path(R.home(),"lib/libR.dylib"))
  files <- c(setNames(md5(binary),normalizePath(binary)),unlist(py$files))
  records <- list(schema_version=1L,R=R.version.string,R_home=R.home(),platform=Sys.info(),compiler=compiler,flags=flags,python=py,files=files,
    note="Rebuilt-runtime identity collected after acceptance; not a claim about the lost historical compiler or binaries")
  write_once(records,file.path(run_root,"runtime/toolchain-record.rds"))
  priors <- list(occJSDM=list(revision=study_config()$occ_revision,
    intercept=occJSDM:::read_intercept_prior(list()),spatial=occJSDM:::read_spatial_sd_prior(list(),TRUE),
    coefficient_variance=list(shape=10,rate=1),traits=0L,biotic_factors=0L,
    defaults_source_hash=md5(file.path(run_root,"sources/occJSDM-9af8597/R/runOccJSDM.R"))),
    Hmsc="Full native prior object retained in every FitRecord; nfMin=nfMax=2 for spatial2; no random level in fixed-basis arms",
    gllvm="VA likelihood; no row effect; fixed coefficients/loadings/ranges held fixed in reported conditional intervals; Hessian not requested",
    sjSDM=study_config()$sjsdm)
  write_once(priors,file.path(run_root,"runtime/prior-record.rds")); invisible(records)
}
verify_toolchain <- function(run_root) {
  path <- file.path(run_root,"runtime/toolchain-record.rds")
  if(!file.exists(path)) {
    if(isTRUE(readRDS(file.path(run_root,"config.rds"))$smoke_only)) return(invisible(TRUE))
    stop("Missing external runtime identity")
  }
  r <- readRDS(path)
  if(!identical(r$R,R.version.string)) stop("R runtime changed")
  if(any(!file.exists(names(r$files)))||any(md5(names(r$files))!=unname(r$files))) stop("External runtime binary/source hash mismatch")
  invisible(TRUE)
}
