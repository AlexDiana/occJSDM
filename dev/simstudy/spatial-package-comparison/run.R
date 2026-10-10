# Invoke through ./run-r so recursive R processes use the frozen compatible runtime.
args <- commandArgs(trailingOnly=TRUE)
options <- list()
for(arg in args) {
  if(!grepl("^--[^=]+=",arg)) stop("Use --name=value arguments")
  key <- sub("^--([^=]+)=.*$","\\1",arg); value <- sub("^--[^=]+=","",arg)
  if(!is.null(options[[key]])) stop("Duplicate argument: ",key)
  options[[key]] <- value
}
# Keep base::options available; this list is only the parsed CLI.
cli <- options; rm(options)
all_args <- commandArgs(); filearg <- grep("^--file=",all_args,value=TRUE)
script <- if(length(filearg)) sub("^--file=","",filearg[1]) else all_args[match("-f",all_args)+1L]
if(is.na(script)) stop("Cannot locate script; invoke run-r -f /absolute/path/run.R")
code <- dirname(normalizePath(script,mustWork=TRUE))
for(f in c("config.R","archive.R","runtime-record.R","basis.R","adapter-contract.R","adapter-occjsdm.R","adapter-sjsdm.R","adapter-hmsc.R","adapter-gllvm.R","oracle-reference.R","diagnostics.R","select.R","runner.R","score.R","attribution.R","verify.R","native-metadata.R","verify-teaching.R","export-teaching.R","plot.R")) source(file.path(code,f))
mode <- cli$mode
if(is.null(mode)) stop("--mode is required")
if(mode=="worker") { worker_task(cli$task); quit(status=0L) }
root <- normalizePath(cli$run,mustWork=TRUE)
if(mode %in% c("preflight","fit","oracle","select","score","export","smoke")) check_input_mode(cli[["input-mode"]])
if(mode=="preflight") {
  repo <- normalizePath(cli$repo,mustWork=TRUE)
  input_path <- file.path(root,"input-record.rds")
  if(!file.exists(input_path)) {
    if(cli[["input-mode"]]=="reconstructed") ir <- reconstruct_inputs(repo,root) else {
      if(is.null(cli$archive)) stop("Restored mode requires --archive")
      ir <- verify_archive(cli$archive,repo); ir$original_archive <- normalizePath(cli$archive)
      for(i in seq_len(nrow(ir$files))) {
        target <- file.path(root,ir$files$relative_path[i]); dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
        if(file.exists(target)) { if(md5(target)!=ir$files$md5[i]) stop("Existing different original input preserved") } else stopifnot(file.copy(ir$files$path[i],target))
        ir$files$path[i] <- target
      }
      write_once(ir,input_path)
    }
  }
  ir <- readRDS(input_path); stopifnot(identical(ir$mode,cli[["input-mode"]]))
  import_inputs(ir,root)
  write_once(study_config(),file.path(root,"config.rds"))
  write_once(comparison_manifest(study_config(),ir$mode),file.path(root,"comparison-manifest.rds"))
  write_once(oracle_manifest(study_config(),ir$mode),file.path(root,"oracle-manifest.rds"))
  files <- resolve_manifest(ir$files,root); stopifnot(all(ir$evidence$checks$passed))
  evidence <- audit_inputs(lapply(files$path,readRDS),repo)
  resolve_manifest(readRDS(file.path(root,"import-manifest.rds")),root)
  stopifnot(identical(readRDS(file.path(root,"config.rds")),study_config()))
  stopifnot(identical(readRDS(file.path(root,"comparison-manifest.rds")),comparison_manifest(study_config(),ir$mode)),
    identical(readRDS(file.path(root,"oracle-manifest.rds")),oracle_manifest(study_config(),ir$mode)))
  prepare_bases(root); prepare_oracles(root)
  verify_toolchain(root)
  verify_package_identity(readRDS(file.path(root,"runtime/package-records.rds")))
  verify_execution_ready(root,code); verify_storage(root)
  record <- freeze_sources(code,repo,root)
  write_once(list(schema_version=1L,snapshot_hash=md5(file.path(root,"snapshot-record.rds")),input_mode=ir$mode,
    checked_inputs=files$md5,audit_checks=nrow(evidence$checks)),file.path(root,"preflight.rds"))
  cat("PREFLIGHT PASSED",ir$mode,"3 communities; 192 comparison cells; 96 oracle jobs\n")
} else if(mode=="verify") {
  print(verify_results(root,phase=if(is.null(cli$phase)) "full" else cli$phase))
} else if(mode=="score") {
  if(!isTRUE(readRDS(file.path(root,"config.rds"))$smoke_only)) {
    verification <- readRDS(file.path(root,paste0("verification-",if(is.null(cli$phase)) "full" else cli$phase,".rds")))
    stopifnot(identical(verification$snapshot_hash,md5(file.path(root,"snapshot-record.rds"))),
      all(verification$checks$status %in% c("passed","unavailable")))
  }
  score_run(root,phase=if(is.null(cli$phase)) "full" else cli$phase)
} else if(mode=="export") {
  repo <- normalizePath(cli$repo,mustWork=TRUE)
  bundle <- build_teaching_bundle(root,repo)
  path <- if(is.null(cli$bundle)) file.path(repo,"vignettes/teaching-data/spatial-package-comparison.rds") else cli$bundle
  save_teaching_bundle(bundle,path)
  if(!is.null(cli$figures)) plot_spatial_comparison(bundle,cli$figures)
  cat("VERIFIED TEACHING EXPORT",normalizePath(path),"\n")
} else if(mode=="checkpoint") {
  result <- timing_checkpoint(root)
  cat("Remaining conservative wall hours:",result$remaining_hours,"; limit:",result$limit_hours,"; passed:",result$passed,"\n")
} else if(mode=="queue") {
  phase <- if(is.null(cli$phase)) "pilot" else cli$phase; kind <- if(is.null(cli$kind)) "all" else cli$kind
  queue <- plan_queue(root,phase,kind,if(is.null(cli$jobs)) NULL else strsplit(cli$jobs,",",fixed=TRUE)[[1]])
  atomic_json(queue,cli$output)
} else if(mode %in% c("fit","oracle","smoke")) {
  verify_snapshot(root)
  if(!file.exists(file.path(root,"preflight.rds"))) stop("Preflight is required")
  config <- readRDS(file.path(root,"config.rds")); preflight <- readRDS(file.path(root,"preflight.rds"))
  if(!isTRUE(config$smoke_only) && !identical(preflight$input_mode,cli[["input-mode"]])) stop("Input mode cannot switch after preflight")
  if(mode=="smoke" && !isTRUE(config$smoke_only)) stop("Smoke mode requires a synthetic archive")
  phase <- if(is.null(cli$phase)) "pilot" else cli$phase
  validate_phase(phase,readRDS(file.path(root,"config.rds")))
  workers <- if(is.null(cli$workers)) 4L else as.integer(cli$workers)
  stopifnot(workers %in% 1:4)
  snap <- file.path(root,"snapshot/study")
  status <- system2("/usr/bin/python3",c(shQuote(file.path(snap,"coordinate.py")),shQuote(file.path(snap,"run-r")),
    shQuote(file.path(snap,"run.R")),shQuote(root),phase,workers,if(mode=="oracle") "oracle" else "all",shQuote(if(is.null(cli$jobs)) "" else cli$jobs)))
  if(status!=0L) stop("Dispatcher stopped; inspect owned run logs and statuses")
} else if(mode=="select") {
  q <- plan_queue(root,if(is.null(cli$phase)) "pilot" else cli$phase)
  cat(q$status,"\n")
} else stop("Mode not yet implemented: ",mode)
