#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
arg<-function(n,d=NULL){x<-grep(paste0('^--',n,'='),args,value=TRUE);if(length(x))sub(paste0('^--',n,'='),'',x[1])else d}
repo<-normalizePath(arg('repo','.'));study<-normalizePath(arg('study'));mode<-arg('mode','pilot')
source(file.path(repo,'dev/simstudy/current-main-recheck/helpers.R'))
settings<-readRDS(file.path(study,if(mode=='selected')'initial' else mode,'settings.rds'))
files<-if(mode=='selected')readRDS(file.path(study,'selection.rds'))$result_file else list.files(file.path(study,mode),pattern='-result\\.rds$',full.names=TRUE)
stopifnot(length(files)>0)
if(mode=='selected') {
  selection<-readRDS(file.path(study,'selection.rds'))
  expected<-names(readRDS(file.path(study,'initial/settings.rds'))$jobs)
  stopifnot(length(files)==110L,!anyDuplicated(files),!anyDuplicated(selection$key),setequal(selection$key,expected))
}
records<-list()
for(file_index in seq_along(files)) {
 f<-files[file_index]
 r<-readRDS(f);job<-r$job
 expected_key<-if(mode=='selected')selection$key[file_index] else sub('-result\\.rds$','',basename(f))
 schedule<-if(mode=='selected')selection$schedule[file_index] else mode
 expected_mcmc<-switch(schedule,pilot=list(nchain=2,nburn=60,niter=80,nthin=1),
   initial=list(nchain=2,nburn=3000,niter=5000,nthin=1),long=list(nchain=4,nburn=6000,niter=12000,nthin=1))
 stopifnot(identical(job,settings$jobs[[expected_key]]),identical(r$source_hashes,settings$source_hashes),
   identical(r$mcmc,expected_mcmc),!job$key %in% names(records))
 stopifnot(identical(unname(tools::md5sum(job$input_file)),job$input_md5))
 input<-readRDS(job$input_file);saved<-readRDS(sub('-result\\.rds$','-fit.rds',f));fit<-saved$fit
 stopifnot(identical(saved$job,job),identical(saved$source_hashes,r$source_hashes),identical(saved$mcmc,expected_mcmc))
 j<-fit$results_output$jsdm_output;ni<-dim(j$B0_output)[2];nc<-dim(j$B0_output)[3]
 stopifnot(ni==expected_mcmc$niter,nc==expected_mcmc$nchain)
 # Independent per-species reconstruction, summing draw-level probabilities.
 n<-nrow(fit$X_psi);S<-dim(j$B0_output)[1];estimate<-matrix(0,n,S)
 for(s in seq_len(S))for(ch in seq_len(nc)) {
   beta<-matrix(j$B_output[,s,,ch],nrow=ncol(fit$X_psi),ncol=ni)
   env<-sweep(fit$X_psi%*%beta,2,j$B0_output[s,,ch],'+')
   hidden<-matrix(0,n,ni)
   for(k in seq_len(dim(j$U_output)[2]))
     hidden<-hidden+sweep(matrix(j$U_output[,k,,ch],n,ni),2,j$L_output[k,s,,ch],'*')
   estimate[,s]<-estimate[,s]+rowSums(plogis(env+hidden))/(ni*nc)
 }
 truth<-if(job$family=='jsdm')plogis(input$truth$eta) else plogis(input$sim$true_params$jsdmParams_true$eta)
 primary<-r$groups[if(job$family=='jsdm') r$groups$scope=='original100' else r$groups$metric=='occupancy_original_sites',]
 stopifnot(nrow(primary)==4L,!anyDuplicated(primary$group),
   setequal(sub('middle','medium',primary$group,fixed=TRUE),c('all','low','medium','high')))
 for(i in seq_len(nrow(primary))) {
   g<-primary[i,];t<-truth[1:100,,drop=FALSE];e<-estimate[1:100,,drop=FALSE]
   mask<-switch(g$group,all=matrix(TRUE,100,S),low=t<.2,medium=t>=.2&t<=.8,middle=t>=.2&t<=.8,high=t>.8)
   stopifnot(!is.null(mask),any(mask),abs(mean(t[mask])-g$truth)<1e-12,
     abs(mean(e[mask])-g$estimate)<1e-12,abs(mean(e[mask]-t[mask])-g$bias)<1e-12,
     abs(mean(abs(e[mask]-t[mask]))-g$mae)<1e-12,
     abs(sqrt(mean((e[mask]-t[mask])^2))-g$rmse)<1e-12)
 }
 difference<-max(abs(estimate-r$estimate))
 stopifnot(difference<1e-12,max(abs(truth-r$truth))<1e-12)
 check_current_fit(fit,input,job)
 records[[job$key]]<-data.frame(key=job$key,input_md5=job$input_md5,
   result_md5=unname(tools::md5sum(f)),fit_md5=unname(tools::md5sum(sub('-result\\.rds$','-fit.rds',f))),
   max_probability_difference=difference,groups_checked=nrow(primary))
 cat(job$key,'verified; max mean difference',format(difference,digits=4),'\n');flush.console()
 rm(saved,fit,j,estimate);gc(FALSE)
}
out<-do.call(rbind,records)
stopifnot(nrow(out)==length(files),setequal(out$key,if(mode=='selected')selection$key else names(settings$jobs)))
write.csv(out,file.path(study,paste0('verification-',mode,'.csv')),row.names=FALSE)
cat('PASS:',nrow(out),'fits;',sum(out$groups_checked),'original-site group scores independently verified.\n')
