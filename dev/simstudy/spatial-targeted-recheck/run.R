#!/usr/bin/env Rscript
# Fits only a separately installed frozen snapshot. Never alters package code.
args <- commandArgs(trailingOnly=TRUE)
option <- function(name,default=NULL) {
  hit <- args[startsWith(args,paste0('--',name,'='))]
  if(length(hit)>1L) stop('Repeated option: ',name)
  if(!length(hit)) { if(is.null(default)) stop('Missing --',name);return(default) }
  substring(hit,nchar(name)+4L)
}
known <- c('repo','study','mode','workers','keys')
stopifnot(all(sub('^--([^=]+)=.*$','\\1',args) %in% known))
repo <- normalizePath(option('repo','.')); study <- normalizePath(option('study'))
mode <- option('mode','pilot'); workers <- as.integer(option('workers','4'))
stopifnot(mode %in% c('prepare','pilot','initial','long'),workers %in% 1:4)
scripts <- file.path(repo,'dev/simstudy/spatial-targeted-recheck')
source(file.path(scripts,'generator.R'));source(file.path(scripts,'score.R'))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',
  OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package('occJSDM'))==normalizePath(file.path(study,'library/occJSDM')))
RcppParallel::setThreadOptions(numThreads=1)
production <- c(file.path(study,'source-main',c('DESCRIPTION','NAMESPACE')),
  list.files(file.path(study,'source-main/R'),pattern='\\.R$',full.names=TRUE),
  list.files(file.path(study,'source-main/src'),pattern='\\.(cpp|h)$|^Makevars',full.names=TRUE),
  file.path(find.package('occJSDM'),'libs/occJSDM.so'))
fit_hashes <- tools::md5sum(c(production,file.path(scripts,c('generator.R','run.R'))))
score_hash <- tools::md5sum(file.path(scripts,'score.R'))
dir.create(file.path(study,'inputs'),showWarnings=FALSE)
jobs <- list(); manifest <- list()
for(grid in c(4L,6L,8L)) for(r in 1:3) {
  community <- sprintf('range%d-rep%02d',grid,r)
  input <- make_spatial_input(grid,r)
  input_file <- file.path(study,'inputs',paste0(community,'.rds'))
  if(file.exists(input_file)) stopifnot(identical(readRDS(input_file),input)) else saveRDS(input,input_file)
  for(arm in c('binary','low','high')) for(knots in c(20L,50L,100L)) {
    key <- sprintf('%s-%s-k%03d',community,arm,knots)
    jobs[[key]] <- list(key=key,community=community,grid_index=grid,replicate=r,
      arm=arm,knots=knots,input_file=input_file,input_md5=unname(tools::md5sum(input_file)))
    manifest[[key]] <- data.frame(key=key,community=community,grid_index=grid,
      replicate=r,arm=arm,knots=knots,input_md5=unname(tools::md5sum(input_file)),
      data_seed=input$settings$data_seed,fit_seed=input$settings$fit_seed)
  }
}
stopifnot(length(jobs)==81L)
manifest <- do.call(rbind,manifest);rownames(manifest) <- NULL
manifest_file <- file.path(study,'input-manifest.rds')
if(file.exists(manifest_file)) stopifnot(identical(readRDS(manifest_file),manifest)) else {
  saveRDS(manifest,manifest_file);write.csv(manifest,file.path(study,'input-manifest.csv'),row.names=FALSE)
}
if(mode=='prepare') {cat('Saved and validated nine paired communities and 81 fit jobs.\n');quit(status=0)}
keys <- option('keys','')
if(nzchar(keys)) {
  keys <- strsplit(keys,',',fixed=TRUE)[[1]]
  stopifnot(all(keys %in% names(jobs)));jobs <- jobs[keys]
} else if(mode=='pilot') jobs <- jobs[grepl('range6-rep01-(binary|high)',names(jobs)) |
  names(jobs)=='range6-rep01-low-k100'] else if(mode=='long') stop('Long runs require explicit selected keys')
mcmc <- switch(mode,pilot=list(nchain=2L,nburn=40L,niter=60L,nthin=1L),
  initial=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L),
  long=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L))
out <- file.path(study,mode);dir.create(out,showWarnings=FALSE)
settings <- list(revision=readLines(file.path(study,'source-revision.txt')),
  fit_hashes=fit_hashes,score_hash=score_hash,mcmc=mcmc,mode=mode,
  workers=workers,threads_per_fit=1L,priors='Unchanged defaults',session=sessionInfo())
settings_file <- file.path(out,'settings.rds')
if(file.exists(settings_file)) {
  old <- readRDS(settings_file)
  stopifnot(identical(old$fit_hashes,fit_hashes),identical(old$mcmc,mcmc),identical(old$score_hash,score_hash))
} else saveRDS(settings,settings_file)
run_job <- function(job) {
  stopifnot(identical(unname(tools::md5sum(job$input_file)),job$input_md5))
  dest <- file.path(out,paste0(job$key,'-result.rds'))
  fitfile <- file.path(out,paste0(job$key,'-fit.rds'))
  if(file.exists(dest)) {
    old <- readRDS(dest)
    stopifnot(identical(old$job,job),identical(old$mcmc,mcmc),
      identical(old$fit_hashes,fit_hashes),identical(old$score_hash,score_hash))
    return(job$key)
  }
  input <- readRDS(job$input_file)
  cat(format(Sys.time()),job$key,'started\n');flush.console()
  if(file.exists(fitfile)) {
    saved <- readRDS(fitfile)
    stopifnot(identical(saved$job,job),identical(saved$mcmc,mcmc),identical(saved$fit_hashes,fit_hashes))
    fit <- saved$fit; warnings <- saved$warnings; started <- saved$started; finished <- saved$finished
  } else {
    assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
    warnings <- character();started <- Sys.time()
    fit <- withCallingHandlers(suppressMessages(occJSDM::runOccJSDM(input$data[[job$arm]],
      listParams=list(n_factors=0L,n_lattrait=0L,n_supportpoints=job$knots),
      listPriors=list(),threshold=1,occCovariates='environment',
      collCovariates=if(job$arm=='binary')NULL else 'collection',
      spatCovariates=c('longitude','latitude'),MCMCparams=mcmc)),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    finished <- Sys.time()
    tmp <- paste0(fitfile,'.tmp')
    saveRDS(list(fit=fit,job=job,mcmc=mcmc,fit_hashes=fit_hashes,
      warnings=warnings,started=started,finished=finished),tmp)
    stopifnot(file.rename(tmp,fitfile))
  }
  scores <- score_spatial_fit(fit,input,job$arm,job$knots)
  result <- c(list(job=job,mcmc=mcmc,fit_hashes=fit_hashes,score_hash=score_hash,
    warnings=warnings,started=started,finished=finished),scores)
  tmp <- paste0(dest,'.tmp');saveRDS(result,tmp);stopifnot(file.rename(tmp,dest))
  cat(format(Sys.time()),job$key,'complete; fit seconds',
    round(as.numeric(difftime(finished,started,units='secs')),1),'; warnings',length(warnings),
    '; max group Rhat',round(max(scores$groups$rhat,na.rm=TRUE),3),'\n');flush.console()
  job$key
}
work <- function(job) tryCatch(run_job(job),error=function(e) {
  cat(job$key,'ERROR:',conditionMessage(e),'\n');list(key=job$key,error=conditionMessage(e))
})
status <- if(workers==1L) lapply(jobs,work) else parallel::mclapply(jobs,work,
  mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(status,file.path(out,paste0('status-',format(Sys.time(),'%Y%m%d%H%M%S'),'.rds')))
stopifnot(length(status)==length(jobs),all(vapply(status,is.character,logical(1))))
cat('Completed',length(status),mode,'fits.\n')
