# Administrative recovery for the 2026-10-08 loaded-version spelling defect.
# Statistical functions remain those in the original immutable snapshot.
install_runtime_repair <- function(root, control) {
  manifest <- jsonlite::read_json(file.path(control,"manifest.json"),simplifyVector=TRUE)
  files <- file.path(control,names(manifest$files))
  if(any(!file.exists(files)) || any(md5(files)!=unname(manifest$files))) stop("Runtime repair source mismatch")
  verify_snapshot(root); verify_execution_ready(root); verify_toolchain(root)
  old <- readLines(file.path(root,"snapshot/study/runner.R"),warn=FALSE)
  corrected <- sub("as.character(getNamespaceVersion(r$package))","as.character(package_version(getNamespaceVersion(r$package)))",old,fixed=TRUE)
  if(sum(old!=corrected)!=1L || !identical(corrected,readLines(file.path(control,"runner.R"),warn=FALSE))) stop("Repair exceeds the one-line version correction")
  sys.source(file.path(control,"runner.R"),envir=.GlobalEnv)
  original_worker <- worker_task
  worker_task <<- function(path) {
    request <- readRDS(path)
    write_once(list(amendment="E12",repair_manifest=md5(file.path(control,"manifest.json")),original_snapshot=md5(file.path(root,"snapshot-record.rds"))),
      file.path(request$directory,"runtime-validation-repair.rds"))
    original_worker(path)
  }
  invisible(TRUE)
}

recover_hmsc_attempt <- function(path, repair_identity) {
  request <- readRDS(file.path(path,"request.rds"))
  if(request$mode!="fit" || request$job$package!="Hmsc") stop("Recovery is restricted to Hmsc fits")
  exit <- jsonlite::read_json(file.path(path,"exit.json"),simplifyVector=TRUE)
  error <- paste(readLines(file.path(path,"stderr.log"),warn=FALSE),collapse="\n")
  expected <- paste("Error in verify_package_identity(packages) : ","  Loaded runtime package mismatch: Hmsc","Calls: worker_task -> verify_package_identity","Execution halted",sep="\n")
  if(!identical(exit$exit_status,1L) || !identical(error,expected)) stop("Not the isolated Hmsc version validation failure")
  required <- file.path(path,c("request.rds","fit.rds","diagnostics.rds","training-prediction.rds"))
  if(any(!file.exists(required))) stop("Incomplete native artifacts cannot be recovered")
  receipt_path <- file.path(path,"version-recovery/receipt.rds")
  has_selection <- file.exists(file.path(dirname(dirname(path)),"selection.rds"))
  if(has_selection && !file.exists(receipt_path)) stop("Recovery cannot revise a completed selection")
  if(file.exists(receipt_path)) {
    receipt <- readRDS(receipt_path)
    if(!identical(receipt$repair_identity,repair_identity) || any(md5(required)!=receipt$artifacts)) stop("Recovery identity/artifacts changed")
    for(n in names(receipt$original_records)) {
      archive <- file.path(path,"version-recovery",paste0("original-",n))
      if(!file.exists(archive)||md5(archive)!=receipt$original_records[[n]]) stop("Original administrative record changed")
    }
    completed <- file.path(path,"complete.rds")
    original_hash <- if("complete.rds" %in% names(receipt$original_records)) receipt$original_records[["complete.rds"]] else NULL
    original_manifest <- file.exists(completed) && !is.null(original_hash) && md5(completed)==original_hash
    if(file.exists(completed) && !original_manifest && resume_state(path,request$identity)$status=="complete") return(invisible(FALSE))
  } else {
    existing <- file.path(path,c("summary.rds","complete.rds")); existing <- existing[file.exists(existing)]
    if(file.exists(file.path(path,"complete.rds"))) {
      old <- resume_state(path,request$identity)
      if(old$status!="interrupted" || !identical(old$summary$reason,"Worker exit 1")) stop("Existing record is not the administrative failure")
    } else if(length(existing)) stop("Unexpected partial completion record")
    receipt <- list(amendment="E12",repair_identity=repair_identity,identity=request$identity,artifacts=md5(required),
      original_records=setNames(md5(existing),basename(existing)),recovered_at=Sys.time(),
      reason="Hmsc namespace version 3.3-7 equals package version 3.3.7; native fitting, diagnostics and training prediction finished before the text-comparison error")
    dir.create(dirname(receipt_path),recursive=TRUE,showWarnings=FALSE)
    for(p in existing) {
      target <- file.path(dirname(receipt_path),paste0("original-",basename(p)))
      if(file.exists(target)) { if(md5(target)!=md5(p)) stop("Different original archive exists") } else stopifnot(file.copy(p,target))
    }
    write_once(receipt,receipt_path)
  }
  if(file.exists(file.path(dirname(dirname(path)),"selection.rds"))) stop("Recovery cannot revise a completed selection")
  fit <- readRDS(required[2])$fit; diagnostics <- readRDS(required[3]); probabilities <- readRDS(required[4])
  stopifnot(identical(fit$package,"Hmsc"),nrow(diagnostics)>0L,is.logical(diagnostics$failed),!anyNA(diagnostics$failed),
    is.matrix(probabilities$p_full),all(is.finite(probabilities$p_full)),all(probabilities$p_full>=0 & probabilities$p_full<=1))
  elapsed <- exit$elapsed_seconds
  observed <- file.path(path,"held-worker-finished.json")
  if(file.exists(observed)) elapsed <- jsonlite::read_json(observed,simplifyVector=TRUE)$elapsed_seconds_upper_bound
  summary <- list(attempt_id=request$settings$attempt_id,status="complete",acceptable=!any(diagnostics$failed),objective=NA_real_,stationarity=NA_real_,
    training_probability=probabilities$p_full,diagnostics=diagnostics,fit_seconds=fit$elapsed_seconds,elapsed_seconds=elapsed)
  # Preserve the original administrative records byte for byte before superseding
  # them. No scientific selection existed when this recovery was authorised.
  for(n in names(receipt$original_records)) {
    current <- file.path(path,n)
    if(file.exists(current) && md5(current)==receipt$original_records[[n]]) unlink(current)
  }
  complete_attempt(path,request$identity,summary)
  stopifnot(resume_state(path,request$identity)$status=="complete")
  invisible(TRUE)
}

recover_hmsc_run <- function(root,control) {
  verify_snapshot(root); verify_toolchain(root); verify_storage(root)
  loadNamespace("Hmsc")
  verify_package_identity(readRDS(file.path(root,"runtime/package-records.rds")))
  identity <- c(original_snapshot=md5(file.path(root,"snapshot-record.rds")),repair_manifest=md5(file.path(control,"manifest.json")))
  completed <- file.path(root,"runtime-recovery.rds")
  if(file.exists(completed)) {
    receipt <- readRDS(completed)
    if(!identical(receipt$identity,identity) || any(md5(file.path(receipt$attempts,"complete.rds"))!=receipt$completion_hashes)) stop("Completed recovery changed")
    for(path in receipt$attempts) recover_hmsc_attempt(path,identity)
    return(invisible(TRUE))
  }
  jobs <- list.dirs(file.path(root,"jobs"),recursive=FALSE,full.names=TRUE)
  jobs <- jobs[grepl("^Hmsc-",basename(jobs))]
  paths <- unlist(lapply(jobs,function(job)list.dirs(file.path(job,"attempts"),recursive=FALSE,full.names=TRUE)),use.names=FALSE)
  stopifnot(length(paths)==16L)
  for(path in paths) recover_hmsc_attempt(path,identity)
  write_once(list(amendment="E12",identity=identity,attempts=paths,completion_hashes=setNames(md5(file.path(paths,"complete.rds")),paths)),file.path(root,"runtime-recovery.rds"))
  cat("Recovered",length(paths),"Hmsc attempts without refitting or revising any selection\n")
}
