#!/usr/bin/env Rscript
# Finite continuation of this experiment, never an unattended recurring task.
a<-commandArgs(TRUE);stopifnot(length(a)==3L)
repo<-normalizePath(a[1]);study<-normalizePath(a[2]);reference<-normalizePath(a[3])
scripts<-file.path(repo,'dev/simstudy/spatial-amplitude-prior')
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
source(file.path(scripts,'metrics.R'));source(file.path(scripts,'analysis.R'))
keys<-unlist(lapply(c(4,6,8),function(g)sprintf('range%d-rep%02d-binary-k100',g,1:3)))
initials<-file.path(study,'half_cauchy/initial',paste0(keys,'-result.rds'))
while(!all(file.exists(initials))) {
  statusfiles<-list.files(file.path(study,'half_cauchy/initial'),pattern='^status-',full.names=TRUE)
  if(length(statusfiles)) {
    z<-readRDS(tail(statusfiles,1));if(any(!vapply(z,is.character,logical(1))))stop('Initial fitting failed')
  }
  Sys.sleep(30)
}
# The preceding run's pool has no more fits to start once all nine results exist.
run_summary<-function(phase) {
  status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(c(file.path(scripts,'summarise.R'),
    paste0('--repo=',repo),paste0('--study=',study),paste0('--phase=',phase))),
    stdout=file.path(study,paste0('summary-',phase,'.log')),stderr=TRUE)
  stopifnot(status==0L)
}
run_summary('select')
selection<-read.csv(file.path(study,'summary-binary-select/selection.csv'),stringsAsFactors=FALSE)
tasks<-list(list(key='range6-rep01-binary-k100',prior='half_cauchy',mode='starts-initial'))
for(prior in c('half_cauchy','inverse_gamma')) {
  needed<-selection$key[selection[[if(prior=='half_cauchy')'run_hc_long' else 'run_ig_long']]]
  for(key in needed)tasks[[length(tasks)+1L]]<-list(key=key,prior=prior,mode='long')
}
# Prepare shared settings once before separate long-job processes can race.
for(prior in unique(vapply(Filter(function(t)t$mode=='long',tasks),function(t)t$prior,character(1)))) {
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
answers<-parallel::mclapply(tasks,work,mc.cores=4L,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(answers,file.path(study,'continuation-status.rds'))
stopifnot(all(vapply(answers,is.character,logical(1))))
run_summary('initial');run_summary('final')
cat(format(Sys.time()),'Binary comparison and initialization check complete. Review gate and audit before any extension.\n')
