#!/usr/bin/env Rscript
# Immutable saved inputs; isolated installation; at most four serial samplers.
args <- commandArgs(trailingOnly=TRUE)
option <- function(name,default=NULL) {
  hit <- args[startsWith(args,paste0('--',name,'='))]
  if(length(hit)>1L) stop('Repeated option: ',name)
  if(!length(hit)) {if(is.null(default))stop('Missing --',name);return(default)}
  substring(hit,nchar(name)+4L)
}
known <- c('repo','study','reference','mode','prior','workers','keys')
stopifnot(all(sub('^--([^=]+)=.*$','\\1',args) %in% known))
repo <- normalizePath(option('repo')); study <- normalizePath(option('study'))
reference <- normalizePath(option('reference'))
mode <- option('mode','initial'); prior <- option('prior','half_cauchy')
workers <- as.integer(option('workers','4'))
stopifnot(mode %in% c('pilot','initial','long','starts-initial','starts-long','baseline'),
  prior %in% c('half_cauchy','inverse_gamma'),workers %in% 1:4)
scripts <- file.path(repo,'dev/simstudy/spatial-amplitude-prior')
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/score.R'))
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
source(file.path(scripts,'metrics.R'))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',
  OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package('occJSDM'))==normalizePath(file.path(study,'library/occJSDM')))
RcppParallel::setThreadOptions(numThreads=1)
manifest <- readRDS(file.path(reference,'input-manifest.rds'))
manifest <- manifest[manifest$knots==100L,]
keys <- option('keys','')
if(nzchar(keys)) {
  keys <- strsplit(keys,',',fixed=TRUE)[[1]]
  stopifnot(!anyDuplicated(keys),all(keys %in% manifest$key))
  manifest <- manifest[match(keys,manifest$key),]
} else manifest <- manifest[manifest$arm=='binary',]
starts <- if(startsWith(mode,'starts-')) c(.1,.3,1,3) else NULL
if(!is.null(starts)) stopifnot(prior=='half_cauchy',nrow(manifest)==1L,
  manifest$key=='range6-rep01-binary-k100')
production <- c(file.path(study,'source',c('DESCRIPTION','NAMESPACE')),
  list.files(file.path(study,'source/R'),pattern='\\.R$',full.names=TRUE),
  list.files(file.path(study,'source/src'),pattern='\\.(cpp|h)$|^Makevars',full.names=TRUE),
  file.path(find.package('occJSDM'),'libs/occJSDM.so'),
  file.path(find.package('occJSDM'),c('DESCRIPTION','NAMESPACE')),
  list.files(file.path(find.package('occJSDM'),'R'),full.names=TRUE))
stopifnot(length(production)>10L,all(file.exists(production)))
fit_hashes <- tools::md5sum(production)
script_hashes <- tools::md5sum(c(file.path(scripts,c('run.R','metrics.R','PLAN.md')),
  file.path(repo,'dev/simstudy/spatial-targeted-recheck',c('score.R','analysis.R'))))
priors <- if(prior=='half_cauchy')list(sigma_bs_prior=prior,sigma_bs_scale=1) else list()
mcmc <- switch(mode,pilot=list(nchain=2L,nburn=40L,niter=60L,nthin=1L),
  initial=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L),
  long=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L),
  'starts-initial'=list(nchain=4L,nburn=3000L,niter=5000L,nthin=1L),
  'starts-long'=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L))
fun <- if(is.null(starts)) occJSDM::runOccJSDM else with_spatial_starts(occJSDM::runOccJSDM,starts)
out <- file.path(study,prior,mode);dir.create(out,recursive=TRUE,showWarnings=FALSE)
settings <- list(fit_hashes=fit_hashes,script_hashes=script_hashes,mcmc=mcmc,prior=prior,
  priors=priors,starts=starts,threads_per_fit=1L,session=sessionInfo())
if(!is.null(starts)) {
  bodyfile <- file.path(out,'fitting-function.R')
  bodytext <- deparse(fun,width.cutoff=100)
  if(file.exists(bodyfile)) stopifnot(identical(readLines(bodyfile),bodytext)) else writeLines(bodytext,bodyfile)
  settings$start_function_md5 <- tools::md5sum(bodyfile)
}
settingsfile <- file.path(out,'settings.rds')
if(file.exists(settingsfile)) {
  previous <- readRDS(settingsfile)
  stopifnot(identical(previous[names(previous)!='session'],settings[names(settings)!='session']))
} else saveRDS(settings,settingsfile)
atomic_save <- function(x,path) {
  temp <- paste0(path,'.tmp');saveRDS(x,temp);stopifnot(file.rename(temp,path))
}
score_and_save <- function(saved,input,dest,source_fit) {
  scores <- score_spatial_fit(saved$fit,input,saved$job$arm,saved$job$knots)
  spatial <- score_amplitude_fit(saved$fit,input)
  stopifnot(abs(scores$field$rmse-spatial$score$raw_rmse)<1e-10)
  result <- c(saved[setdiff(names(saved),'fit')],scores,
    list(spatial=spatial,script_hashes=script_hashes,
      source_fit=normalizePath(source_fit),source_fit_md5=unname(tools::md5sum(source_fit))))
  result$reasons <- unique(c(diagnostic_reasons(result),spatial_flags(spatial$diagnostics)))
  atomic_save(result,dest)
  result
}
run_job <- function(row) {
  key <- row$key;inputfile <- file.path(reference,'inputs',paste0(row$community,'.rds'))
  stopifnot(identical(unname(tools::md5sum(inputfile)),row$input_md5))
  input <- readRDS(inputfile)
  # Baseline diagnostics are re-extracted from archived fits, never refitted here.
  if(mode=='baseline') {
    for(phase in c('initial','long')) {
      fitfile <- file.path(reference,phase,paste0(key,'-fit.rds'))
      if(!file.exists(fitfile)) {
        if(phase=='initial') stop('Missing baseline initial fit: ',key)
        next
      }
      originalfile <- file.path(reference,phase,paste0(key,'-result.rds'))
      stopifnot(file.exists(originalfile))
      dest <- file.path(out,paste0(key,'-',phase,'-result.rds'))
      if(file.exists(dest)) {
        old <- readRDS(dest)
        stopifnot(identical(old$script_hashes,script_hashes),
          identical(old$source_fit_md5,unname(tools::md5sum(fitfile))),
          identical(old$source_result_md5,unname(tools::md5sum(originalfile))))
        next
      }
      saved <- readRDS(fitfile)
      validate_result_protocol(saved,row,phase)
      original <- readRDS(originalfile)
      stopifnot(identical(original$job,saved$job),identical(original$fit_hashes,saved$fit_hashes))
      spatial <- score_amplitude_fit(saved$fit,input)
      stopifnot(abs(original$field$rmse-spatial$score$raw_rmse)<1e-10)
      original$spatial <- spatial;original$script_hashes <- script_hashes
      original$source_result <- normalizePath(originalfile)
      original$source_result_md5 <- unname(tools::md5sum(originalfile))
      original$source_fit <- normalizePath(fitfile);original$source_fit_md5 <- unname(tools::md5sum(fitfile))
      original$reasons <- unique(c(diagnostic_reasons(original),spatial_flags(spatial$diagnostics)))
      atomic_save(original,dest)
      cat(key,phase,'baseline diagnostics saved; flags',length(original$reasons),'\n');flush.console()
    }
    return(key)
  }
  dest <- file.path(out,paste0(key,'-result.rds'));fitfile <- file.path(out,paste0(key,'-fit.rds'))
  job <- as.list(row[c('key','community','grid_index','replicate','arm','knots','input_md5')])
  job$input_file <- inputfile
  metadata <- list(job=job,mcmc=mcmc,fit_hashes=fit_hashes,priors=priors,starts=starts)
  verify <- function(saved)stopifnot(identical(saved[names(metadata)],metadata))
  if(file.exists(dest)) {
    old <- readRDS(dest);verify(old);stopifnot(identical(old$script_hashes,script_hashes),
      identical(old$source_fit_md5,unname(tools::md5sum(fitfile))))
    return(key)
  }
  cat(format(Sys.time()),key,prior,mode,'started\n');flush.console()
  if(file.exists(fitfile)) {saved <- readRDS(fitfile);verify(saved)} else {
    assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
    warnings <- character();started <- Sys.time()
    fit <- withCallingHandlers(suppressMessages(fun(input$data[[job$arm]],
      listParams=list(n_factors=0L,n_lattrait=0L,n_supportpoints=job$knots),
      listPriors=priors,threshold=1,occCovariates='environment',
      collCovariates=if(job$arm=='binary')NULL else 'collection',
      spatCovariates=c('longitude','latitude'),MCMCparams=mcmc)),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    finished <- Sys.time()
    stopifnot(identical(fit$infos$spatial_sd_prior$type,prior))
    saved <- c(metadata,list(fit=fit,warnings=warnings,started=started,finished=finished,
      rng_initial=input$fit_rng,rng_final=.Random.seed,threads_per_fit=1L))
    atomic_save(saved,fitfile)
  }
  result <- score_and_save(saved,input,dest,fitfile)
  cat(format(Sys.time()),key,prior,mode,'complete; fit seconds',
    round(as.numeric(difftime(saved$finished,saved$started,units='secs')),1),
    '; diagnostic flags',length(result$reasons),'\n');flush.console()
  key
}
work <- function(i)tryCatch(run_job(manifest[i,,drop=FALSE]),error=function(e){
  cat(manifest$key[i],'ERROR:',conditionMessage(e),'\n');list(key=manifest$key[i],error=conditionMessage(e))})
status <- if(workers==1L)lapply(seq_len(nrow(manifest)),work) else parallel::mclapply(
  seq_len(nrow(manifest)),work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(status,file.path(out,paste0('status-',format(Sys.time(),'%Y%m%d%H%M%S'),'.rds')))
stopifnot(length(status)==nrow(manifest),all(vapply(status,is.character,logical(1))))
cat('Completed',length(status),prior,mode,'jobs.\n')
