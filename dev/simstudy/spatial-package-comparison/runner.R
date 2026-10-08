attempt_settings <- function(package,job_id,attempt_id,config=study_config(),source_hashes=character()) {
  oracle <- package=="oracle"
  list(seed=attempt_seed(job_id,attempt_id,if(oracle) "oracle" else "fit"),
    prediction_seed=attempt_seed(job_id,attempt_id,if(oracle) "oracle_prediction" else "prediction"),
    diagnostic_seed=attempt_seed(job_id,attempt_id,"diagnostic"),job_id=as.integer(job_id),attempt_id=as.integer(attempt_id),
    source_hashes=source_hashes,mcmc=config$mcmc[[if(attempt_id==1L) "initial" else "extended"]],sjsdm=config$sjsdm,gllvm=config$gllvm)
}
complete_attempt <- function(path,identity,summary) {
  write_once(summary,file.path(path,"summary.rds"))
  files <- list.files(path,pattern="\\.(rds|so|dylib|cpp|h)$",full.names=TRUE,recursive=TRUE)
  files <- files[basename(files)!="complete.rds"]
  write_once(list(schema_version=1L,identity=identity,status=summary$status,
    files=setNames(md5(files),substring(files,nchar(path)+2L))),file.path(path,"complete.rds"))
}
resume_state <- function(path,expected_hashes) {
  done <- file.path(path,"complete.rds")
  if(!file.exists(done)) return(list(status=if(dir.exists(path)) "interrupted" else "not_started"))
  record <- readRDS(done)
  if(!identical(record$identity,expected_hashes)) stop("Completion identity mismatch: ",path)
  files <- file.path(path,names(record$files))
  if(any(!file.exists(files))||any(md5(files)!=record$files)) stop("Completion file hash mismatch: ",path)
  list(status=record$status,summary=readRDS(file.path(path,"summary.rds")),completion=record)
}
source_manifest <- function(path) {
  files <- list.files(path,recursive=TRUE,full.names=TRUE,all.files=TRUE)
  files <- files[!dir.exists(files)&!grepl("/__pycache__/|\\.pyc$",files)]
  setNames(md5(files),substring(files,nchar(path)+2L))
}
verify_storage <- function(run_root) {
  path <- file.path(run_root,"storage-record.json")
  jobs <- file.path(run_root,"jobs")
  if(!file.exists(path)) return(if(dir.exists(jobs)) normalizePath(jobs) else normalizePath(run_root))
  record <- jsonlite::read_json(path,simplifyVector=TRUE)
  if(!dir.exists(jobs)||!file.exists(record$sentinel_path)) stop("Declared archive volume is unavailable")
  if(!identical(normalizePath(jobs),record$jobs_path) ||
      !identical(readLines(record$sentinel_path,warn=FALSE),record$sentinel_value)) stop("Archive volume identity mismatch")
  normalizePath(jobs)
}
verify_execution_ready <- function(run_root,code=file.path(run_root,"snapshot/study")) {
  config <- readRDS(file.path(run_root,"config.rds"))
  if(isTRUE(config$smoke_only)) return(invisible(TRUE))
  for(name in c("certification","infrastructure-review")) {
    path <- file.path(run_root,paste0(name,".rds"))
    if(!file.exists(path)) stop("Scientific execution requires ",name)
    record <- readRDS(path)
    if(!isTRUE(record$passed)||!length(record$evidence)) stop("Scientific execution has not passed ",name)
    if(!identical(record$source_files,source_manifest(code))) stop("Certified/reviewed source identity mismatch: ",name)
  }
  invisible(TRUE)
}
verify_snapshot <- function(run_root) {
  record <- readRDS(file.path(run_root,"snapshot-record.rds"))
  actual <- source_manifest(file.path(run_root,"snapshot"))
  if(!identical(record$files,actual)) stop("Frozen source hash mismatch")
  if(md5(file.path(run_root,"config.rds"))!=record$config_md5) stop("Frozen config hash mismatch")
  if(length(record$frozen_files)) {
    paths <- file.path(run_root,names(record$frozen_files))
    if(any(!file.exists(paths))||any(md5(paths)!=record$frozen_files)) stop("Frozen protocol/input identity mismatch")
  }
  record
}
verify_package_identity <- function(record) {
  for(r in record) {
    actual <- find.package(r$package,quiet=TRUE)
    if(length(actual)!=1L || !identical(normalizePath(actual),normalizePath(r$path))) stop("Runtime package resolution mismatch: ",r$package)
    if(!identical(as.character(utils::packageVersion(r$package,lib.loc=dirname(actual))),r$version)) stop("Runtime package version mismatch: ",r$package)
    if(r$package %in% loadedNamespaces()) {
      loaded <- getNamespaceInfo(asNamespace(r$package),"path")
      if(!identical(normalizePath(loaded),normalizePath(r$path)) || !identical(as.character(package_version(getNamespaceVersion(r$package))),r$version)) stop("Loaded runtime package mismatch: ",r$package)
    }
    files <- list.files(actual,recursive=TRUE,full.names=TRUE)
    observed <- setNames(md5(files),substring(files,nchar(actual)+2L))
    expected <- setNames(as.character(r$files),substring(names(r$files),nchar(r$path)+2L))
    if(!identical(observed,expected)) stop("Runtime package hash mismatch: ",r$package)
  }
  invisible(TRUE)
}
freeze_sources <- function(code,repo,run_root) {
  target <- file.path(run_root,"snapshot")
  if(file.exists(file.path(run_root,"snapshot-record.rds"))) return(verify_snapshot(run_root))
  if(dir.exists(target)) stop("Incomplete source snapshot preserved; inspect before continuing")
  dir.create(file.path(target,"study"),recursive=TRUE)
  files <- list.files(code,full.names=TRUE); files <- files[basename(files)!="__pycache__"]
  stopifnot(all(file.copy(files,file.path(target,"study"),recursive=TRUE,copy.mode=TRUE)))
  dir.create(file.path(target,"protocol"))
  protocol <- c("dev/superpowers/specs/2026-10-08-spatial-package-comparison-design.md","dev/superpowers/plans/2026-10-08-spatial-package-comparison.md")
  stopifnot(all(file.copy(file.path(repo,protocol),file.path(target,"protocol"))))
  original <- file.path(run_root,"sources/acceptance")
  stopifnot(dir.exists(original),file.copy(original,target,recursive=TRUE))
  frozen <- c("config.rds","comparison-manifest.rds","oracle-manifest.rds","runtime/package-records.rds",
    "input-record.rds","import-manifest.rds","runtime/toolchain-record.rds","runtime/prior-record.rds","certification.rds","infrastructure-review.rds","storage-record.json",
    list.files(file.path(run_root,"training"),full.names=FALSE,pattern="rds$"))
  frozen <- c(frozen[!grepl("^rep",frozen)],
    file.path("training",list.files(file.path(run_root,"training"),pattern="rds$")),
    file.path("bases",list.files(file.path(run_root,"bases"),pattern="rds$")),
    file.path("oracle-inputs",list.files(file.path(run_root,"oracle-inputs"),pattern="rds$")),
    file.path("prediction",list.files(file.path(run_root,"prediction"),pattern="rds$")))
  frozen <- frozen[file.exists(file.path(run_root,frozen))]
  record <- list(schema_version=1L,files=source_manifest(target),frozen_files=setNames(md5(file.path(run_root,frozen)),frozen),config_md5=md5(file.path(run_root,"config.rds")),
    frozen_at=Sys.time(),source_revision=system2("git",c("-C",shQuote(repo),"rev-parse","HEAD"),stdout=TRUE))
  write_once(record,file.path(run_root,"snapshot-record.rds")); record
}
job_identity <- function(run_root,job,basis_path,training_path) {
  c(snapshot=md5(file.path(run_root,"snapshot-record.rds")),config=md5(file.path(run_root,"config.rds")),
    runtime=md5(file.path(run_root,"runtime/package-records.rds")),training=md5(training_path),basis=md5(basis_path),
    job=as.character(job$key))
}
make_jobs <- function(run_root) {
  m <- readRDS(file.path(run_root,"comparison-manifest.rds")); o <- readRDS(file.path(run_root,"oracle-manifest.rds"))
  m$kind <- "comparison"; o$kind <- "oracle"; o$package <- "oracle"; o$arm <- "oracle"
  o$dataset_id <- paste(o$community,o$arrangement,sep="-"); o$key <- paste0("oracle-",o$key)
  cols <- c("key","job_id","kind","package","arm","dataset_id","community","arrangement")
  m$species <- NA_character_; rbind(m[c(cols,"species")],o[c(cols,"species")])
}
prepare_bases <- function(run_root) {
  jobs <- make_jobs(run_root)
  for(dataset in unique(jobs$dataset_id)) {
    tr <- readRDS(file.path(run_root,"training",paste0(dataset,".rds"))); validate_training(tr)
    for(arm in unique(jobs$arm[jobs$dataset_id==dataset & jobs$kind=="comparison"])) {
      dest <- file.path(run_root,"bases",paste0(dataset,"-",arm,".rds"))
      write_once(fit_basis(tr,arm),dest)
    }
  }
}
# This preprocessing function is the only non-oracle fitting route allowed to read full inputs.
prepare_oracles <- function(run_root) {
  jobs <- subset(make_jobs(run_root),kind=="oracle")
  for(comm in unique(jobs$community)) {
    input <- readRDS(file.path(run_root,"inputs",paste0(comm,".rds")))
    for(dataset in unique(jobs$dataset_id[jobs$community==comm])) {
      arr <- jobs$arrangement[match(dataset,jobs$dataset_id)]
      tr <- input$surveys[[arr]]$truth
      minimal <- list(landscape=input$landscape[c("B0","B","range")],surveys=setNames(list(list(truth=tr[c("xy","environment","z") ])),arr))
      write_once(minimal,file.path(run_root,"oracle-inputs",paste0(dataset,".rds")))
    }
  }
}
atomic_json <- function(object,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(tmpdir=dirname(path)); on.exit(unlink(tmp),add=TRUE)
  jsonlite::write_json(object,tmp,auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
  stopifnot(file.rename(tmp,path))
}
read_checked <- function(path,hash) {
  if(!file.exists(path) || !identical(md5(path),unname(hash))) stop("Declared input hash mismatch: ",path)
  readRDS(path)
}
collect_attempts <- function(job_dir,identity) {
  paths <- list.dirs(file.path(job_dir,"attempts"),recursive=FALSE,full.names=TRUE)
  out <- list()
  for(path in paths) {
    state <- resume_state(path,identity)
    if(state$status=="interrupted") {
      request <- readRDS(file.path(path,"request.rds"))
      exit <- file.path(path,"exit.json")
      if(!file.exists(exit)) {
        owner <- file.path(path,"owner.json")
        if(!file.exists(owner)) next
        pid <- jsonlite::read_json(owner,simplifyVector=TRUE)$pid
        command <- suppressWarnings(system2("/bin/ps",c("-p",as.integer(pid),"-o","command="),stdout=TRUE))
        if(length(command) && any(grepl(file.path(path,"request.rds"),command,fixed=TRUE))) stop("Owned worker is still active: ",path)
        atomic_json(list(status="interrupted",exit_status=NA_integer_,elapsed_seconds=NA_real_),exit)
      }
      e <- jsonlite::read_json(exit,simplifyVector=TRUE)
      s <- list(attempt_id=request$settings$attempt_id,status=if(e$status=="resource_interrupted") "resource_interrupted" else "interrupted",
        acceptable=FALSE,objective=NA_real_,elapsed_seconds=e$elapsed_seconds,reason=paste("Worker exit",e$exit_status))
      complete_attempt(path,identity,s); state <- resume_state(path,identity)
    }
    out[[length(out)+1L]] <- state$summary
  }
  out
}
check_attempt_limits <- function(paths,config) {
  paths <- paths[grepl("/attempts/",paths)]
  counts <- c(comparison=sum(!grepl("^oracle-",paths)),oracle=sum(grepl("^oracle-",paths)),total=length(paths))
  exceeded <- names(counts)[counts>config$limits[names(counts)]]
  if(length(exceeded)) stop("Attempt ceiling exceeded: ",paste(exceeded,collapse=", "))
  invisible(TRUE)
}
stage_tasks <- function(tasks,phase) {
  if(phase!="full" || !length(tasks)) return(tasks)
  first <- vapply(tasks,function(t)t$package %in% c("occJSDM","sjSDM","oracle"),logical(1))
  if(any(first)) tasks[first] else tasks
}
plan_queue <- function(run_root,phase="pilot",kind="all",jobs_filter=NULL) {
  snapshot <- verify_snapshot(run_root); config <- readRDS(file.path(run_root,"config.rds"))
  verify_execution_ready(run_root); verify_storage(run_root)
  validate_phase(phase,config)
  jobs <- make_jobs(run_root)
  if(phase=="pilot") jobs <- subset(jobs,community=="rep01" & arrangement %in% c("spread","clustered"))
  if(phase=="full") {
    gate <- file.path(run_root,"runtime-checkpoint.rds")
    if(!file.exists(gate)||!isTRUE(readRDS(gate)$passed)) stop("Runtime checkpoint has not passed; full dispatch forbidden")
  }
  if(!kind %in% c("all","comparison","oracle")) stop("Unknown job kind")
  if(kind!="all") jobs <- jobs[jobs$kind==kind,,drop=FALSE]
  if(!is.null(jobs_filter)) jobs <- jobs[jobs$key %in% jobs_filter,,drop=FALSE]
  tasks <- list(); states <- character()
  for(i in seq_len(nrow(jobs))) {
    job <- as.list(jobs[i,,drop=FALSE]); directory <- file.path(run_root,"jobs",job$key)
    training <- file.path(run_root,"training",paste0(job$dataset_id,".rds"))
    basis <- if(job$kind=="oracle") file.path(run_root,"oracle-inputs",paste0(job$dataset_id,".rds")) else
      file.path(run_root,"bases",paste0(job$dataset_id,"-",job$arm,".rds"))
    identity <- job_identity(run_root,job,basis,training)
    a <- collect_attempts(directory,identity)
    selection_path <- file.path(directory,"selection.rds")
    if(file.exists(selection_path)) {
      selection <- readRDS(selection_path)
      if(!identical(selection$identity,identity)) stop("Immutable selection identity mismatch")
    } else {
      selection <- select_attempts(a,config,if(job$package %in% c("sjSDM","gllvm")) "point" else "bayesian")
      if(selection$status!="pending") {
        selection$identity <- identity; selection$job <- job; selection$created <- Sys.time()
        selection$attempt_hashes <- setNames(md5(file.path(directory,"attempts",sprintf("%02d",selection$attempts),"complete.rds")),as.character(selection$attempts))
        write_once(selection,selection_path)
      }
    }
    states <- c(states,selection$status)
    for(id in selection$next_attempts) {
      path <- file.path(directory,"attempts",sprintf("%02d",id))
      settings <- attempt_settings(job$package,job$job_id,id,config,snapshot$files)
      request <- list(mode="fit",run_root=run_root,job=job,identity=identity,directory=path,training=training,basis=basis,settings=settings)
      request_path <- file.path(path,"request.rds"); write_once(request,request_path)
      tasks[[length(tasks)+1L]] <- list(request=request_path,directory=path,package=job$package,timeout_seconds=config$attempt_hours*3600)
    }
    if(!is.na(selection$selected_attempt) && selection$status!="pending") {
      path <- file.path(directory,"prediction")
      if(!file.exists(file.path(path,"complete.rds"))) {
        if(file.exists(file.path(path,"exit.json"))) {
          # Prediction failures are retained and cannot cause an automatic infinite retry.
          states[length(states)] <- "prediction_failed"; next
        }
        if(file.exists(file.path(path,"owner.json"))) {
          pid <- jsonlite::read_json(file.path(path,"owner.json"),simplifyVector=TRUE)$pid
          command <- suppressWarnings(system2("/bin/ps",c("-p",as.integer(pid),"-o","command="),stdout=TRUE))
          if(length(command) && any(grepl(file.path(path,"request.rds"),command,fixed=TRUE))) stop("Owned prediction worker is still active: ",path)
          atomic_json(list(status="interrupted",exit_status=NA_integer_,elapsed_seconds=NA_real_),file.path(path,"exit.json"))
          states[length(states)] <- "prediction_failed"; next
        }
        attempt_path <- file.path(directory,"attempts",sprintf("%02d",selection$selected_attempt))
        request <- list(mode="predict",run_root=run_root,job=job,identity=identity,directory=path,
          selected=attempt_path,selection_hash=md5(selection_path),prediction=file.path(run_root,"prediction",paste0(job$dataset_id,".rds")),
          prediction_hash=md5(file.path(run_root,"prediction",paste0(job$dataset_id,".rds"))))
        request_path <- file.path(path,"request.rds"); write_once(request,request_path)
        tasks[[length(tasks)+1L]] <- list(request=request_path,directory=path,package=job$package,timeout_seconds=config$attempt_hours*3600)
      } else {
        result <- resume_state(path,identity)
        if(result$status!="complete") states[length(states)] <- "prediction_failed"
      }
    }
  }
  attempts <- list.files(file.path(run_root,"jobs"),pattern="request.rds$",recursive=TRUE)
  check_attempt_limits(attempts,config)
  tasks <- stage_tasks(tasks,phase)
  list(tasks=tasks,status=if(any(states %in% c("budget_amendment_required","prediction_failed"))) "needs_attention" else if(length(tasks)) "running" else "complete",states=states)
}
worker_task <- function(path) {
  request <- readRDS(path); root <- request$run_root; verify_snapshot(root); verify_toolchain(root)
  verify_execution_ready(root); verify_storage(root)
  packages <- readRDS(file.path(root,"runtime/package-records.rds")); verify_package_identity(packages)
  cfg <- readRDS(file.path(root,"config.rds")); started <- Sys.time(); directory <- request$directory
  # Separate process input contracts: fitting has no truth/prediction path; prediction has no outcomes.
  result <- tryCatch({
    if(request$mode=="fit") {
      tr <- read_checked(request$training,request$identity[["training"]])
      b <- read_checked(request$basis,request$identity[["basis"]]); validate_training(tr)
      if(request$job$kind=="oracle") {
        species <- match(request$job$species,tr$species_id)
        tr$species_id <- tr$species_id[species]; tr$y <- tr$y[,species,drop=FALSE]
        fit <- fit_oracle(tr,b,request$job$arrangement,species,request$settings,
          file.path(root,"snapshot/acceptance"),file.path(directory,"cpp-cache"))
        b <- NULL
      } else fit <- get(paste0("fit_",tolower(request$job$package)))(tr,b,request$settings)
      # Save the expensive native result before diagnostics, while keeping completion atomic.
      write_once(list(fit=fit,basis=b),file.path(directory,"fit.rds"))
      diagnostics <- diagnose_fit(fit,tr,b,cfg)
      prediction <- tr; prediction$y <- NULL
      options(spatial.retain_raw_field=fit$package=="oracle" || (fit$package=="occJSDM" && fit$arm=="gp"))
      probabilities <- get(paste0("predict_",tolower(fit$package)))(fit,prediction,b)
      options(spatial.retain_raw_field=FALSE)
      write_once(diagnostics,file.path(directory,"diagnostics.rds"))
      write_once(probabilities,file.path(directory,"training-prediction.rds"))
      point <- fit$package %in% c("sjSDM","gllvm")
      list(attempt_id=request$settings$attempt_id,status="complete",acceptable=!any(diagnostics$failed),
        objective=if(point) diagnostics$objective[1] else NA_real_,
        stationarity=if(point) diagnostics$stationarity[1] else NA_real_,
        training_probability=probabilities$p_full,diagnostics=diagnostics,
        fit_seconds=fit$elapsed_seconds,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")))
    } else {
      selpath <- file.path(dirname(directory),"selection.rds")
      if(md5(selpath)!=request$selection_hash) stop("Selection changed before prediction")
      resume_state(request$selected,request$identity)
      saved <- readRDS(file.path(request$selected,"fit.rds")); f <- saved$fit
      p <- read_checked(request$prediction,request$prediction_hash)
      if(request$job$kind=="oracle") p$species_id <- request$job$species
      # Predictions are opened only after the immutable diagnostic selection exists.
      audit <- p$site_id[cfg$audit_lattice_cells[cfg$audit_lattice_cells<=length(p$site_id)]]
      options(spatial.audit_ids=audit,spatial.draw_sink=function(ids,draws) {
        keep <- which(ids %in% audit)
        if(length(keep)) write_once(list(site_id=ids[keep],block_site_id=ids,rng_before=attr(draws,"rng_before"),species_id=p$species_id,
          draws=lapply(draws,function(x)x[keep,,,,drop=FALSE]),draw_ids=list(iteration=seq_len(dim(draws[[1]])[3]),chain=seq_len(dim(draws[[1]])[4]))),
          file.path(directory,"audit",paste0(ids[keep[1]],".rds")))
      })
      probabilities <- get(paste0("predict_",tolower(f$package)))(f,p,saved$basis)
      write_once(probabilities,file.path(directory,"prediction.rds"))
      list(status="complete",selection_hash=request$selection_hash,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),audit_ids=audit)
    }
  },error=function(e) list(status="failed",attempt_id=if(request$mode=="fit") request$settings$attempt_id else NA_integer_,
    acceptable=FALSE,objective=NA_real_,reason=conditionMessage(e),calls=vapply(sys.calls(),function(x) paste(deparse(x),collapse=" "),character(1)),
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs"))))
  verify_snapshot(root); verify_package_identity(packages); verify_toolchain(root)
  complete_attempt(directory,request$identity,result)
  cat("SPATIAL WORKER COMPLETE",request$job$key,result$status,"\n")
  if(result$status!="complete") cat(result$reason,"\n")
  invisible(result)
}

validate_phase <- function(phase,config) {
  if(!phase %in% c("pilot","full","smoke")) stop("Unknown execution phase")
  if(phase=="smoke" && !isTRUE(config$smoke_only)) stop("Smoke phase requires synthetic inputs")
  phase
}
project_runtime <- function(attempt_times,prediction_times,remaining,config=study_config()) {
  arms <- list(); ratio <- sum(config$mcmc$extended[c("nburn","niter")])/sum(config$mcmc$initial[c("nburn","niter")])
  draw_ratio <- config$mcmc$extended[["niter"]]/config$mcmc$initial[["niter"]]
  for(i in seq_len(nrow(remaining))) {
    row <- remaining[i,,drop=FALSE]; a <- subset(attempt_times,package==row$package & arm==row$arm)
    p <- subset(prediction_times,package==row$package & arm==row$arm)
    bayes <- !row$package %in% c("sjSDM","gllvm")
    maximum <- function(x) if(length(x)&&all(is.finite(x))) max(x) else Inf
    initial <- maximum(a$seconds[if(bayes) a$attempt_id==1 else rep(TRUE,nrow(a))])
    extended <- if(bayes) maximum(a$seconds[a$attempt_id==2]) else 0
    label <- if(bayes && is.finite(extended)) "measured diagnostic-triggered extension" else if(bayes) "iteration-scaled estimate; no measured extension" else "six starts conservatively budgeted"
    if(bayes && !is.finite(extended)) extended <- initial*ratio
    # A selected initial posterior has one quarter as many draws as an extension.
    predict <- maximum(p$seconds*if(bayes) ifelse(p$attempt_id==1,draw_ratio,1) else 1)
    cost <- if(bayes) initial+extended+predict else config$diagnostics$max_starts*initial+predict
    arms[[i]] <- cbind(row,initial_seconds=initial,extended_seconds=extended,prediction_seconds=predict,
      seconds_per_cell=cost,extension_cost=label)
  }
  arms <- do.call(rbind,arms)
  remaining_jobs <- rep(arms$seconds_per_cell,arms$cells)
  # List-scheduling upper bound for four workers, including an indivisible last cell.
  hours <- if(length(remaining_jobs)) (sum(remaining_jobs)/config$workers+(1-1/config$workers)*max(remaining_jobs))/3600 else 0
  c(runtime_gate(hours,config$remaining_wall_hours),list(arms=arms,workers=config$workers,
    convention="maximum observed attempt costs by arm; all extensions/six starts; extended prediction cost; list-scheduling bound"))
}
timing_checkpoint <- function(run_root) {
  verify_snapshot(run_root); config <- readRDS(file.path(run_root,"config.rds"))
  verify_storage(run_root)
  if(isTRUE(config$smoke_only)) stop("Synthetic timing cannot certify the scientific runtime gate")
  jobs <- make_jobs(run_root); pilot <- subset(jobs,community=="rep01" & arrangement %in% c("spread","clustered"))
  attempts <- predictions <- list(); issues <- character(); extension_rates <- list()
  for(i in seq_len(nrow(pilot))) {
    job <- pilot[i,]; dir <- file.path(run_root,"jobs",job$key)
    sp <- file.path(dir,"selection.rds")
    if(!file.exists(sp)) { issues <- c(issues,paste(job$key,"missing immutable selection")); next }
    selection <- readRDS(sp)
    if(selection$status=="budget_amendment_required") issues <- c(issues,paste(job$key,"resource interruption"))
    extension_rates[[i]] <- data.frame(package=job$package,arm=job$arm,key=job$key,attempts=length(selection$attempts),extended=if(job$package %in% c("sjSDM","gllvm")) length(selection$attempts)>3L else 2L %in% selection$attempts,selected_status=selection$status)
    for(id in selection$attempts) {
      ap <- file.path(dir,"attempts",sprintf("%02d",id)); s <- readRDS(file.path(ap,"summary.rds"))
      if(s$status=="complete") {
        e <- jsonlite::read_json(file.path(ap,"exit.json"),simplifyVector=TRUE)
        attempts[[length(attempts)+1L]] <- data.frame(package=job$package,arm=job$arm,attempt_id=id,seconds=e$elapsed_seconds)
      }
    }
    pp <- file.path(dir,"prediction/summary.rds")
    if(file.exists(pp) && readRDS(pp)$status=="complete") {
      e <- jsonlite::read_json(file.path(dirname(pp),"exit.json"),simplifyVector=TRUE)
      predictions[[length(predictions)+1L]] <- data.frame(package=job$package,arm=job$arm,attempt_id=selection$selected_attempt,seconds=e$elapsed_seconds)
    } else if(!is.na(selection$selected_attempt)) issues <- c(issues,paste(job$key,"prediction not complete"))
  }
  empty <- data.frame(package=character(),arm=character(),attempt_id=integer(),seconds=numeric())
  a <- if(length(attempts)) do.call(rbind,attempts) else empty
  p <- if(length(predictions)) do.call(rbind,predictions) else empty
  remaining <- aggregate(key~package+arm,data=jobs[!jobs$key %in% pilot$key,],FUN=length); names(remaining)[3] <- "cells"
  out <- project_runtime(a,p,remaining,config)
  if(length(issues)) out$passed <- FALSE
  out$issues <- issues; out$attempt_times <- a; out$prediction_times <- p; out$pilot_cells <- do.call(rbind,extension_rates)
  out$storage <- storage_projection(run_root,pilot,remaining,config)
  out$time_passed <- out$passed
  out$passed <- out$passed && out$storage$passed
  if(!out$storage$passed) out$issues <- c(out$issues,"Insufficient measured/estimated archive storage; choose another volume before full dispatch")
  out$snapshot_hash <- md5(file.path(run_root,"snapshot-record.rds")); out$truth_scores_opened <- FALSE
  write_once(out,file.path(run_root,"runtime-checkpoint.rds")); out
}
archive_bytes <- function(path) {
  f <- list.files(path,recursive=TRUE,full.names=TRUE,all.files=TRUE)
  sum(file.info(f)$size,na.rm=TRUE)
}
available_bytes <- function(path) {
  out <- system2("/bin/df",c("-Pk",shQuote(path)),stdout=TRUE)
  fields <- strsplit(trimws(tail(out,1)),"[[:space:]]+")[[1]]
  as.numeric(fields[4])*1024
}
storage_projection <- function(run_root,pilot,remaining,config) {
  a <- p <- list()
  for(i in seq_len(nrow(pilot))) {
    job <- pilot[i,]; directory <- file.path(run_root,"jobs",job$key)
    for(ap in list.dirs(file.path(directory,"attempts"),recursive=FALSE,full.names=TRUE)) {
      s <- file.path(ap,"summary.rds")
      if(file.exists(s)&&readRDS(s)$status=="complete") a[[length(a)+1L]] <- data.frame(package=job$package,arm=job$arm,attempt_id=as.integer(basename(ap)),seconds=archive_bytes(ap))
    }
    sp <- file.path(directory,"selection.rds"); pp <- file.path(directory,"prediction/summary.rds")
    if(file.exists(sp)&&file.exists(pp)&&readRDS(pp)$status=="complete") p[[length(p)+1L]] <- data.frame(package=job$package,arm=job$arm,attempt_id=readRDS(sp)$selected_attempt,seconds=archive_bytes(dirname(pp)))
  }
  empty <- data.frame(package=character(),arm=character(),attempt_id=integer(),seconds=numeric())
  cfg <- config; cfg$mcmc$initial[["nburn"]] <- 0L; cfg$mcmc$extended[["nburn"]] <- 0L
  sizes <- project_runtime(if(length(a)) do.call(rbind,a) else empty,if(length(p)) do.call(rbind,p) else empty,remaining,cfg)$arms
  required <- sum(sizes$seconds_per_cell*sizes$cells)
  # Reserve 10 GiB for the OS, transient block arrays, exports and ordinary work.
  storage_path <- verify_storage(run_root)
  free <- available_bytes(storage_path); reserve <- 10*1024^3
  list(passed=is.finite(required)&&required+reserve<free,remaining_bytes=required,available_bytes=free,reserve_bytes=reserve,
    storage_path=storage_path,
    by_arm=transform(sizes,bytes_per_cell=seconds_per_cell)[c("package","arm","cells","bytes_per_cell")],
    convention="Compressed files measured in pilot; maximum by arm; every extension/six starts and fourfold retained-draw scaling when extension is unmeasured")
}
