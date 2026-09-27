#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
arg<-function(key,default){x<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(x))sub(paste0('^--',key,'='),'',x[1]) else default}
here<-normalizePath(arg('study','work/nonspatial-design-recheck-20260919'))
root<-normalizePath(arg('baseline',file.path(here,'../nonspatial-bias-recheck-20260914')))
mode<-arg('mode','run');workers<-as.integer(arg('workers','4'));stopifnot(mode %in% c('prepare','pilot','run'),workers>=1L)
source(file.path(here,'design_helpers.R'));frozen<-load_frozen(root);score_design_fit<-make_design_scorer()
dir.create(file.path(here,'inputs'),showWarnings=FALSE)
mcmc<-if(mode=='pilot')list(nchain=2,nburn=60,niter=80,nthin=1) else list(nchain=2,nburn=3000,niter=5000,nthin=1)
out<-file.path(here,if(mode=='pilot')'pilot' else 'fits');dir.create(out,showWarnings=FALSE)
source_hashes<-c(frozen$hashes,tools::md5sum(file.path(here,c('design_helpers.R','run_design_recheck.R'))))
jobs<-list();inputs_manifest<-list()
for(sc in c('qnear_K6','qfar_K6'))for(r in seq_len(10L)) {
 base_file<-file.path(root,'operational-k6',sprintf('%s-data-%02d.rds',sc,r));a<-readRDS(base_file)
 base_result<-file.path(root,'operational-k6',sprintf('%s-q20-%02d-result.rds',sc,r))
 original<-readRDS(base_result)
 stopifnot(identical(a$fit_rng,readRDS(original$job$input_file)$fit_rng))
 for(arm in c('field4','sites300','knownU')) {
  key<-sprintf('%s-%s-%02d',sc,arm,r);b<-make_design_input(a,arm);check_nested_input(a,b,arm)
  b$design$baseline_input<-base_file;b$design$baseline_result<-base_result
  input_file<-file.path(here,'inputs',paste0(key,'.rds'))
  if(file.exists(input_file))stopifnot(identical(readRDS(input_file),b)) else saveRDS(b,input_file)
  inputs_manifest[[key]]<-data.frame(key=key,scenario=sc,replicate=r,arm=arm,n=b$scenario$n,M=b$scenario$M,K=6,P=2,input_md5=unname(tools::md5sum(input_file)),baseline_input_md5=unname(tools::md5sum(base_file)),baseline_result_md5=unname(tools::md5sum(base_result)),seed=a$seed,original_sites_preserved=100,original_field_samples_preserved=200,original_pcr_rows_preserved=nrow(a$sim$data_list$info))
  if(mode!='pilot'||r==1L)jobs[[key]]<-list(key=key,input_file=input_file,baseline_result=base_result,arm=arm,priors=original$job$priors,mcmc=mcmc,out=out,source_hashes=source_hashes)
 }
}
write.csv(do.call(rbind,inputs_manifest),file.path(here,'input-manifest.csv'),row.names=FALSE)
if(mode=='prepare'){cat('Validated and saved 60 nested/control inputs.\n');quit(status=0)}
settings<-list(source_revision=frozen$source_revision,baseline_root=root,mcmc=mcmc,mode=mode,workers=workers,threads_per_fit=1,source_hashes=source_hashes,keys=names(jobs),session=sessionInfo(),started=Sys.time(),runner=readLines(file.path(here,'run_design_recheck.R')),helpers=readLines(file.path(here,'design_helpers.R')))
settings_file<-file.path(out,'settings.rds')
if(file.exists(settings_file)){
 old<-readRDS(settings_file);stopifnot(identical(old$mcmc,settings$mcmc),identical(old$source_hashes,settings$source_hashes))
}else saveRDS(settings,settings_file)
run_job<-function(job) {
 result_file<-file.path(job$out,paste0(job$key,'-result.rds'))
 if(file.exists(result_file)){old<-readRDS(result_file);stopifnot(identical(old$job,job));return(job$key)}
 input<-readRDS(job$input_file);s<-input$scenario;fit_file<-file.path(job$out,paste0(job$key,'-fit.rds'))
 cat(format(Sys.time()),job$key,'started\n');flush.console()
 if(file.exists(fit_file)){
  saved<-readRDS(fit_file);stopifnot(identical(saved$job,job));fit<-saved$fit;warnings<-saved$warnings;started<-saved$started;finished<-saved$finished
 }else{
  assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
  fitter<-if(job$arm=='knownU')known_u_fitter(input$sim$true_params$jsdmParams_true$U) else occJSDM::runOccJSDM
  warnings<-character();started<-Sys.time()
  fit<-withCallingHandlers(suppressMessages(fitter(input$sim$data_list,listParams=list(n_factors=s$d,n_lattrait=s$gt),listPriors=job$priors,occCovariates=paste0('X_psi.EnvCov.',1:s$ncov_psi),collCovariates='X_theta',spatCovariates=NULL,MCMCparams=job$mcmc)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  finished<-Sys.time();saveRDS(list(fit=fit,job=job,warnings=warnings,started=started,finished=finished),fit_file)
 }
 control_check<-if(job$arm=='knownU')check_known_u_output(fit,input$sim$true_params$jsdmParams_true$U) else NULL
 scores<-score_design_fit(fit,input,job$priors)
 result<-c(list(job=job,started=started,finished=finished,warnings=warnings,control_check=control_check),scores)
 temporary<-paste0(result_file,'.tmp');saveRDS(result,temporary);stopifnot(file.rename(temporary,result_file))
 cat(format(Sys.time()),job$key,'complete;',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'s;',length(warnings),'warnings; max scored Rhat',round(max(scores$elements$rhat),3),'\n');flush.console()
 job$key
}
work<-function(j)tryCatch(run_job(j),error=function(e){cat(j$key,'ERROR:',conditionMessage(e),'\n');list(key=j$key,error=conditionMessage(e))})
status<-if(workers==1L)lapply(jobs,work) else parallel::mclapply(jobs,work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(status,file.path(out,'job-status.rds'))
stopifnot(length(status)==length(jobs),all(vapply(status,is.character,logical(1))))
cat('Completed',length(status),mode,'fits.\n')
