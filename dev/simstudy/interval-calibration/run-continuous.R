#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
arg<-function(k,d){z<-grep(paste0('^--',k,'='),args,value=TRUE);if(length(z))sub(paste0('^--',k,'='),'',z[1])else d}
repo<-normalizePath(arg('repo','.'))
study<-normalizePath(arg('study',file.path(repo,'dev/simstudy/results/interval-calibration-20261008')))
mode<-arg('mode','pilot');reps<-as.integer(arg('reps','100'));workers<-as.integer(arg('workers','4'))
shard<-as.integer(arg('shard','0'));shards<-as.integer(arg('shards','4'))
stopifnot(mode %in% c('pilot','initial','long'),reps>0L,workers %in% 1:4)
lib<-file.path(study,'library');.libPaths(c(lib,.libPaths()))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package('occJSDM'))==normalizePath(file.path(lib,'occJSDM')))
source(file.path(repo,'tests/testthat/helper-simstudy.R'))
source(file.path(repo,'dev/simstudy/interval-calibration/helpers.R'))
RcppParallel::setThreadOptions(numThreads=1)
sourcefiles<-list.files(file.path(study,'source'),pattern='\\.(R|cpp|h)$',recursive=TRUE,full.names=TRUE)
hashes<-tools::md5sum(c(sourcefiles,file.path(lib,'occJSDM/libs/occJSDM.so'),
  file.path(repo,'dev/simstudy/interval-calibration',c('helpers.R','run-continuous.R')),
  file.path(repo,'tests/testthat/helper-simstudy.R')))
mcmc<-switch(mode,pilot=list(nchain=2,nburn=50,niter=100,nthin=1),
  initial=list(nchain=2,nburn=1000,niter=2000,nthin=1),long=list(nchain=4,nburn=3000,niter=6000,nthin=1))
out<-file.path(study,mode);dir.create(out,showWarnings=FALSE)
inputdir<-file.path(study,'inputs');dir.create(inputdir,showWarnings=FALSE)
settings<-list(revision=readLines(file.path(study,'source-revision.txt')),hashes=hashes,mcmc=mcmc,reps=reps,
  mode=mode,session=sessionInfo(),workers=workers)
settingsfile<-file.path(out,'settings.rds')
if(file.exists(settingsfile)) {
  old<-readRDS(settingsfile);stopifnot(identical(old$hashes,hashes),identical(old$mcmc,mcmc),old$reps==reps)
}else {
  pending<-paste0(settingsfile,'.',Sys.getpid(),'.tmp')
  saveRDS(settings,pending);stopifnot(file.rename(pending,settingsfile))
}
snapshot_sources<-file.path(repo,c('dev/simstudy/interval-calibration/run-continuous.R',
  'dev/simstudy/interval-calibration/helpers.R','tests/testthat/helper-simstudy.R'))
snapshot_paths<-file.path(study,c('runner-executed.R','helpers-executed.R','simstudy-helper-executed.R'))
for(i in seq_along(snapshot_paths)) {
  if(file.exists(snapshot_paths[i])) {
    stopifnot(unname(tools::md5sum(snapshot_paths[i]))==unname(hashes[snapshot_sources[i]]))
  }else {
    pending<-paste0(snapshot_paths[i],'.',Sys.getpid(),'.tmp')
    stopifnot(file.copy(snapshot_sources[i],pending),file.rename(pending,snapshot_paths[i]))
  }
}
jobs<-list()
for(r in seq_len(reps)) {
  scenario<-Filter(function(s)s$label=='continuous',simstudy_scenarios())[[1]]
  scenario$useSpatField<-FALSE;scenario$n_supportpoints<-NULL
  inputfile<-file.path(inputdir,sprintf('community-%03d.rds',r))
  if(!file.exists(inputfile)) {
    seed<-simstudy_seed_for(scenario,r);truth<-draw_truth(scenario,seed)
    sim<-simstudy_simulate(truth);fit_rng<-.Random.seed
    zero<-sim;jp<-sim$true_params$jsdmParams_true
    factor_effect<-jp$U%*%jp$L
    zero$data_list$OTU<-sim$data_list$OTU-factor_effect
    zero$true_params$jsdmParams_true$eta<-jp$eta-factor_effect
    zero$true_params$jsdmParams_true$U<-matrix(numeric(),scenario$n,0)
    zero$true_params$jsdmParams_true$L<-matrix(numeric(),0,scenario$S)
    zero$true_params$jsdmParams_true$varPart<-NULL
    zero$true_params$z_true<-zero$data_list$OTU
    stopifnot(max(abs((zero$data_list$OTU-zero$true_params$jsdmParams_true$eta)-
      (sim$data_list$OTU-jp$eta)))<1e-12)
    pending<-paste0(inputfile,'.',Sys.getpid(),'.tmp')
    saveRDS(list(seed=seed,fit_rng=fit_rng,truth=truth,scenario=scenario,sim=list(d2=sim,d0=zero)),pending)
    stopifnot(file.rename(pending,inputfile))
  }
  for(arm in c('d2','d0')) {
    key<-sprintf('%s-%03d',arm,r)
    jobs[[key]]<-list(key=key,arm=arm,replicate=r,input=inputfile,input_md5=unname(tools::md5sum(inputfile)))
  }
}
bad_fit<-function(x)length(x$warnings)>0L || any(!is.finite(x$rows$rhat)) || any(x$rows$rhat>1.05) ||
  any(!is.finite(x$rows$ess_bulk)) || any(x$rows$ess_bulk<400) || any(!is.finite(x$rows$ess_tail)) || any(x$rows$ess_tail<400)
if(mode=='pilot')jobs<-jobs[seq_len(min(4L,length(jobs)))]
if(mode=='long') {
  picked<-vapply(jobs,function(j)bad_fit(readRDS(file.path(study,'initial',paste0(j$key,'-result.rds')))),logical(1))
  jobs<-jobs[picked]
  writeLines(names(jobs),file.path(out,'diagnostic-selected-keys.txt'))
}
if(shard>0L) {
  stopifnot(shard<=shards)
  jobs<-jobs[(seq_along(jobs)-1L)%%shards==shard-1L]
}
run_job<-function(job) {
  resultfile<-file.path(out,paste0(job$key,'-result.rds'))
  if(file.exists(resultfile)) {
    old<-readRDS(resultfile);stopifnot(identical(old$hashes,hashes),identical(old$job,job),identical(old$mcmc,mcmc));return(job$key)
  }
  stopifnot(unname(tools::md5sum(job$input))==job$input_md5)
  input<-readRDS(job$input);sim<-input$sim[[job$arm]];s<-input$scenario
  assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
  warnings<-character();started<-Sys.time()
  fit<-withCallingHandlers(suppressMessages(runOccJSDM(sim$data_list,
    listParams=list(n_factors=if(job$arm=='d2')2L else 0L,n_lattrait=s$gt),
    occCovariates=paste0('X_psi.EnvCov.',seq_len(s$ncov_psi)),
    collCovariates=NULL,spatCovariates=NULL,MCMCparams=mcmc)),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  finished<-Sys.time();jp<-sim$true_params$jsdmParams_true;jo<-fit$results_output$jsdm_output
  stopifnot(fit$infos$model=='continuous',fit$infos$ps==0L,is.null(fit$results_output$beta_theta_output),
    is.null(fit$results_output$p_output),is.null(fit$results_output$q_output),fit$infos$intercept_prior$sd==1)
  eta<-sweep(fit$X_psi%*%jp$B+jp$U%*%jp$L,2,jp$B0,'+')
  stopifnot(max(abs(eta-jp$eta))<1e-10)
  rows<-rbind(interval_rows(jo$B0_output,jp$B0,'B0',job$arm,job$replicate),
    interval_rows(jo$tau_output,jp$tau,'tau',job$arm,job$replicate))
  result<-list(job=job,hashes=hashes,mcmc=mcmc,warnings=warnings,rows=rows,started=started,finished=finished)
  saveRDS(c(result,list(fit=fit)),file.path(out,paste0(job$key,'-fit.rds')))
  saveRDS(result,resultfile)
  cat(format(finished),job$key,'complete',round(as.numeric(difftime(finished,started,units='secs')),1),
    'seconds; max Rhat',round(max(rows$rhat),3),'min ESS',round(min(rows$ess_bulk,rows$ess_tail)),
    'warnings',length(warnings),'\n');flush.console()
  job$key
}
work<-function(j)tryCatch(run_job(j),error=function(e){cat(j$key,'ERROR',conditionMessage(e),'\n');list(error=conditionMessage(e),job=j$key)})
if(workers==1L)status<-lapply(jobs,work) else {
  launcher<-file.path(study,'R-4.5-runtime/bin/Rscript')
  if(!file.exists(launcher))launcher<-file.path(R.home('bin'),'Rscript')
  cl<-parallel::makePSOCKcluster(workers,rscript=launcher,outfile='')
  parallel::clusterExport(cl,c('repo','study','lib','out','hashes','mcmc','run_job','work','interval_rows'))
  parallel::clusterEvalQ(cl,{
    .libPaths(c(lib,.libPaths()));Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    suppressPackageStartupMessages(library(occJSDM));RcppParallel::setThreadOptions(numThreads=1)
    stopifnot(normalizePath(find.package('occJSDM'))==normalizePath(file.path(lib,'occJSDM')))
  })
  status<-parallel::parLapplyLB(cl,jobs,work);parallel::stopCluster(cl)
}
saveRDS(status,file.path(out,if(shard>0L)paste0('status-shard-',shard,'.rds')else 'status.rds'))
stopifnot(all(vapply(status,is.character,logical(1))),length(status)==length(jobs))
cat('Completed',length(jobs),mode,'jobs.\n')
