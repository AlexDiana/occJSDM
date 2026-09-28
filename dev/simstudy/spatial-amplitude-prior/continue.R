#!/usr/bin/env Rscript
# Finite continuation of this experiment, never an unattended recurring task.
a<-commandArgs(TRUE);stopifnot(length(a) %in% c(3L,4L))
workers<-if(length(a)==4L)as.integer(a[4]) else 4L
stopifnot(length(workers)==1L,is.finite(workers),workers>=4L,workers<=8L)
repo<-normalizePath(a[1]);study<-normalizePath(a[2]);reference<-normalizePath(a[3])
scripts<-file.path(repo,'dev/simstudy/spatial-amplitude-prior')
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
source(file.path(scripts,'metrics.R'));source(file.path(scripts,'analysis.R'))
keys<-unlist(lapply(c(4,6,8),function(g)sprintf('range%d-rep%02d-binary-k100',g,1:3)))
initials<-file.path(study,'half_cauchy/initial',paste0(keys,'-result.rds'))
run_summary<-function(phase) {
  status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(c(file.path(scripts,'summarise.R'),
    paste0('--repo=',repo),paste0('--study=',study),paste0('--phase=',phase))),
    stdout=file.path(study,paste0('summary-',phase,'.log')),stderr=TRUE)
  stopifnot(status==0L)
}
# Prepare shared settings atomically in the parent before separate long processes.
prepare_long<-function(prior) {
  settings<-readRDS(file.path(study,'half_cauchy/initial/settings.rds'))
  settings$mcmc<-list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L)
  settings$prior<-prior
  settings$priors<-if(prior=='half_cauchy')list(sigma_bs_prior=prior,sigma_bs_scale=1) else list()
  folder<-file.path(study,prior,'long');dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  path<-file.path(folder,'settings.rds')
  if(file.exists(path)) {
    old<-readRDS(path);fields<-setdiff(names(settings),'session')
    stopifnot(identical(old[fields],settings[fields]))
  } else {
    temp<-paste0(path,'.tmp');saveRDS(settings,temp);stopifnot(file.rename(temp,path))
  }
}
run_fit<-function(task) {
  log<-file.path(study,paste(task$key,task$prior,task$mode,'log',sep='.'))
  cat(format(Sys.time()),task$key,task$prior,task$mode,'dispatched\n');flush.console()
  status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(c(file.path(scripts,'run.R'),
    paste0('--repo=',repo),paste0('--study=',study),paste0('--reference=',reference),
    paste0('--mode=',task$mode),paste0('--prior=',task$prior),paste0('--keys=',task$key),'--workers=1')),
    stdout=log,stderr=TRUE)
  if(status!=0L)stop('Fit failed, see ',log)
  result<-read_checked_result(file.path(study,task$prior,task$mode,paste0(task$key,'-result.rds')))
  cat(format(Sys.time()),task$key,task$mode,'complete; flags',length(result$reasons),'\n');flush.console()
  if(task$mode=='starts-initial' && length(spatial_flags(result$spatial$diagnostics))) {
    task$mode<-'starts-long';run_fit(task)
  }
  task$key
}
work<-function(task)tryCatch(run_fit(task),error=function(e)list(task=task,error=conditionMessage(e)))
id<-function(task)paste(task$key,task$prior,task$mode,sep=':')
# Both jobs are prespecified. Extra capacity can run them alongside the original
# four-worker initial pool. Only diagnostics determine other early longer jobs.
queue<-list(list(key='range6-rep01-binary-k100',prior='half_cauchy',mode='starts-initial'),
  list(key='range6-rep03-binary-k100',prior='half_cauchy',mode='long'))
active<-list();answers<-list();launched<-character();planned<-FALSE;checked<-character()
new_fit_hashes<-readRDS(file.path(study,'half_cauchy/initial/settings.rds'))$fit_hashes
enqueue<-function(task) {
  if(!id(task) %in% c(launched,vapply(queue,id,character(1))))queue[[length(queue)+1L]]<<-task
}
repeat {
  done<-file.exists(initials)
  statusfiles<-list.files(file.path(study,'half_cauchy/initial'),pattern='^status-',full.names=TRUE)
  if(length(statusfiles)) {
    z<-readRDS(tail(statusfiles,1));if(any(!vapply(z,is.character,logical(1))))stop('Initial fitting failed')
  }
  for(key in setdiff(keys[done],checked)) {
    hc<-read_checked_result(file.path(study,'half_cauchy/initial',paste0(key,'-result.rds')))
    ig<-read_checked_result(file.path(study,'inverse_gamma/baseline',paste0(key,'-initial-result.rds')))
    stopifnot(identical(hc$job$input_md5,ig$job$input_md5),hc$job$key==key,ig$job$key==key)
    validate_choice(hc,hc,'half_cauchy','initial',new_fit_hashes)
    validate_choice(ig,ig,'inverse_gamma','initial',new_fit_hashes)
    if(pair_needs_long(ig,hc)) {
      enqueue(list(key=key,prior='half_cauchy',mode='long'))
      if(!file.exists(file.path(study,'inverse_gamma/baseline',paste0(key,'-long-result.rds'))))
        enqueue(list(key=key,prior='inverse_gamma',mode='long'))
    }
    checked<-c(checked,key)
  }
  if(all(done) && !planned) {
    run_summary('select')
    selection<-read.csv(file.path(study,'summary-binary-select/selection.csv'),stringsAsFactors=FALSE)
    for(prior in c('half_cauchy','inverse_gamma')) {
      needed<-selection$key[selection[[if(prior=='half_cauchy')'run_hc_long' else 'run_ig_long']]]
      for(key in needed) {
        task<-list(key=key,prior=prior,mode='long')
        enqueue(task)
      }
    }
    saveRDS(list(queued=queue,already_dispatched=launched),file.path(study,'continuation-plan.rds'))
    planned<-TRUE
  }
  initial_active_bound<-min(4L,sum(!done))
  slots<-workers-initial_active_bound
  while(length(queue) && length(active)<slots) {
    task<-queue[[1]];queue<-queue[-1L]
    if(task$mode=='long')prepare_long(task$prior)
    child<-parallel::mcparallel(work(task),mc.set.seed=FALSE)
    active[[as.character(child$pid)]]<-child;launched<-c(launched,id(task))
  }
  if(length(active)) {
    completed<-parallel::mccollect(active,wait=FALSE)
    if(length(completed)) {
      answers<-c(answers,completed);active[names(completed)]<-NULL
      saveRDS(answers,file.path(study,'continuation-status.rds'))
    }
  }
  if(planned && !length(active) && !length(queue))break
  Sys.sleep(15)
}
stopifnot(all(vapply(answers,is.character,logical(1))),
  !any(vapply(answers,inherits,logical(1),what='try-error')))
run_summary('initial');run_summary('final')
cat(format(Sys.time()),'Binary comparison and initialization check complete. Review gate and audit before any extension.\n')
