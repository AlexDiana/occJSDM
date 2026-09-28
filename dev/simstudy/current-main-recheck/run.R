#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
arg <- function(name,default=NULL) {
  x <- grep(paste0('^--',name,'='),args,value=TRUE)
  if(length(x)) sub(paste0('^--',name,'='),'',x[1]) else default
}
repo <- normalizePath(arg('repo','.'))
study <- normalizePath(arg('study'))
archives <- normalizePath(arg('archives'))
mode <- arg('mode','pilot');workers <- as.integer(arg('workers','4'))
stopifnot(mode %in% c('pilot','initial','long'), workers %in% 1:4)
source(file.path(repo,'dev/simstudy/current-main-recheck/helpers.R'))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package('occJSDM'))==normalizePath(file.path(study,'library/occJSDM')))
RcppParallel::setThreadOptions(numThreads=1)
scoring <- load_scoring(repo)
jobs <- make_jobs(archives)
keys <- arg('keys','')
if(nzchar(keys)) {
  keys <- strsplit(keys,',',fixed=TRUE)[[1]]
  stopifnot(all(keys %in% names(jobs)));jobs <- jobs[keys]
} else if(mode=='pilot') {
  jobs <- jobs[c('jsdm-n0100-01','jsdm-n0300-01','jsdm-n1000-01',
    'design-qnear_K6-q20-01','design-qnear_K6-field4-01',
    'design-qnear_K6-sites300-01','design-qnear_K6-knownU-01')]
} else if(mode=='long') stop('Long runs require explicit diagnostic-selected keys.')
mcmc <- switch(mode,pilot=list(nchain=2,nburn=60,niter=80,nthin=1),
  initial=list(nchain=2,nburn=3000,niter=5000,nthin=1),
  long=list(nchain=4,nburn=6000,niter=12000,nthin=1))
out <- file.path(study,mode);dir.create(out,showWarnings=FALSE)
production <- list.files(file.path(study,'source-main'),pattern='\\.(R|cpp|h)$',recursive=TRUE,full.names=TRUE)
helpers <- file.path(repo,c('dev/simstudy/current-main-recheck/helpers.R','dev/simstudy/current-main-recheck/run.R',
  'dev/simstudy/jsdm-sample-size-recheck/helpers.R','dev/simstudy/nonspatial-design-recheck/design_helpers.R',
  'dev/simstudy/nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R'))
source_hashes <- tools::md5sum(c(production,helpers,file.path(find.package('occJSDM'),'libs/occJSDM.so')))
settings <- list(source_revision=readLines(file.path(study,'source-revision.txt')),source_hashes=source_hashes,
  mcmc=mcmc,mode=mode,workers=workers,session=sessionInfo(),jobs=jobs)
sfile <- file.path(out,'settings.rds')
if(file.exists(sfile)) {
 old<-readRDS(sfile);stopifnot(identical(old$source_hashes,source_hashes),identical(old$mcmc,mcmc))
} else saveRDS(settings,sfile)
run_job <- function(job) {
  dest <- file.path(out,paste0(job$key,'-result.rds'))
  fitfile <- file.path(out,paste0(job$key,'-fit.rds'))
  stopifnot(identical(unname(tools::md5sum(job$input_file)),job$input_md5))
  input <- readRDS(job$input_file)
  if(file.exists(dest)) {
    old<-readRDS(dest);stopifnot(identical(old$job,job),identical(old$mcmc,mcmc),identical(old$source_hashes,source_hashes));return(job$key)
  }
  cat(format(Sys.time()),job$key,'started\n');flush.console()
  if(file.exists(fitfile)) {
    saved<-readRDS(fitfile);stopifnot(identical(saved$job,job),identical(saved$mcmc,mcmc),identical(saved$source_hashes,source_hashes))
    fit<-saved$fit;warnings<-saved$warnings;started<-saved$started;finished<-saved$finished
  } else {
    assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
    warnings<-character();started<-Sys.time()
    fit <- withCallingHandlers(suppressMessages({
      if(job$family=='jsdm') {
        occJSDM::runOccJSDM(input$data,listParams=list(n_factors=input$n_factors,n_lattrait=input$n_lattrait),
          occCovariates=input$covariates,collCovariates=NULL,spatCovariates=NULL,MCMCparams=mcmc)
      } else {
        s<-input$scenario
        fitter<-if(job$arm=='knownU') scoring$design$known_u_fitter(input$sim$true_params$jsdmParams_true$U) else occJSDM::runOccJSDM
        fitter(input$sim$data_list,listParams=list(n_factors=s$d,n_lattrait=s$gt),listPriors=job$priors,
          threshold=1,occCovariates=paste0('X_psi.EnvCov.',seq_len(s$ncov_psi)),collCovariates='X_theta',
          spatCovariates=NULL,MCMCparams=mcmc)
      }
    }),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    finished<-Sys.time()
    tmp<-paste0(fitfile,'.tmp')
    saveRDS(list(fit=fit,job=job,mcmc=mcmc,source_hashes=source_hashes,warnings=warnings,started=started,finished=finished),tmp)
    stopifnot(file.rename(tmp,fitfile))
  }
  check_current_fit(fit,input,job)
  control <- if(job$arm=='knownU') scoring$design$check_known_u_output(fit,input$sim$true_params$jsdmParams_true$U) else NULL
  scores <- score_current_fit(fit,input,job,scoring)
  all_rhat <- c(scores$groups$rhat,scores$elements$rhat,scores$additional_diagnostics$rhat)
  diag <- list(max_group_rhat=safe_max(scores$groups$rhat),max_element_rhat=safe_max(c(scores$elements$rhat,scores$additional_diagnostics$rhat)),
    unresolved_rhat=sum(!is.finite(all_rhat) | all_rhat<=0))
  tmp<-paste0(dest,'.tmp')
  saveRDS(c(list(job=job,mcmc=mcmc,source_hashes=source_hashes,warnings=warnings,
    started=started,finished=finished,control_check=control,diagnostics=diag),scores),tmp)
  stopifnot(file.rename(tmp,dest))
  cat(format(Sys.time()),job$key,'complete; fit',round(as.numeric(difftime(finished,started,units='secs')),1),
    's;',length(warnings),'warnings; group Rhat',round(diag$max_group_rhat,3),
    '; element Rhat',round(diag$max_element_rhat,3),'\n');flush.console()
  job$key
}
work <- function(job) tryCatch(run_job(job),error=function(e){cat(job$key,'ERROR:',conditionMessage(e),'\n');list(key=job$key,error=conditionMessage(e))})
status <- if(workers==1L) lapply(jobs,work) else parallel::mclapply(jobs,work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(status,file.path(out,paste0('status-',format(Sys.time(),'%Y%m%d%H%M%S'),'.rds')))
stopifnot(length(status)==length(jobs),all(vapply(status,is.character,logical(1))))
cat('Completed',length(status),mode,'fits.\n')
