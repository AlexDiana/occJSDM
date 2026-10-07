#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
arg<-function(n,d=NULL){x<-grep(paste0('^--',n,'='),args,value=TRUE);if(length(x))sub(paste0('^--',n,'='),'',x[1])else d}
repo<-normalizePath(arg('repo','.'));study<-normalizePath(arg('study'));mode<-arg('mode','plan')
source(file.path(repo,'dev/simstudy/current-main-recheck/helpers.R'))
settings<-readRDS(file.path(study,'initial/settings.rds'));jobs<-settings$jobs
stopifnot(length(jobs)==110L,!anyDuplicated(names(jobs)))
dm<-read.csv(file.path(repo,'dev/simstudy/nonspatial-design-recheck/results/fit-manifest.csv'))
jm<-read.csv(file.path(repo,'dev/simstudy/jsdm-sample-size-recheck/results/fit-manifest.csv'))
historical<-c(paste0('design-',dm$key[dm$longer_fit]),paste0('jsdm-',jm$key[jm$schedule!='initial']))
stopifnot(length(historical)==10L,all(historical %in% names(jobs)))
if(mode=='plan') {
 rows<-lapply(jobs,function(job){
   f<-file.path(study,'initial',paste0(job$key,'-result.rds'));stopifnot(file.exists(f))
   x<-readRDS(f);stopifnot(identical(x$job,job),identical(x$source_hashes,settings$source_hashes),identical(x$mcmc,settings$mcmc))
   d<-x$diagnostics
   flag<-length(x$warnings)>0L || d$unresolved_rhat>0L ||
     !is.finite(d$max_group_rhat) || !is.finite(d$max_element_rhat) ||
     d$max_group_rhat>1.05 || d$max_element_rhat>1.05
   data.frame(key=job$key,historical_long=job$key %in% historical,current_flag=flag,
     warnings=length(x$warnings),max_group_rhat=d$max_group_rhat,max_element_rhat=d$max_element_rhat,
     unresolved_rhat=d$unresolved_rhat,run_long=flag || job$key %in% historical)
 })
 plan<-do.call(rbind,rows)
 output<-file.path(study,'long-selection.csv')
 if(file.exists(output))stopifnot(isTRUE(all.equal(read.csv(output),plan,check.attributes=FALSE))) else write.csv(plan,output,row.names=FALSE)
 writeLines(paste(plan$key[plan$run_long],collapse=','),file.path(study,'long-keys.txt'))
 cat('Selected',sum(plan$run_long),'long runs before comparing ecological errors.\n')
 print(plan[plan$run_long,])
} else {
 stopifnot(mode=='final')
 plan<-read.csv(file.path(study,'long-selection.csv'))
 stopifnot(nrow(plan)==110L,setequal(plan$key,names(jobs)),!anyDuplicated(plan$key))
 rows<-lapply(jobs,function(job){
   p<-plan[plan$key==job$key,];schedule<-if(p$run_long)'long' else 'initial'
   f<-file.path(study,schedule,paste0(job$key,'-result.rds'));stopifnot(file.exists(f))
   r<-readRDS(f)
   expected<-readRDS(file.path(study,schedule,'settings.rds'))$mcmc
   stopifnot(identical(r$job,job),identical(r$source_hashes,settings$source_hashes),identical(r$mcmc,expected))
   old<-job$old_initial_file
   if(job$key %in% historical)old<-sub('/fits/','/long-fits/',old,fixed=TRUE)
   stopifnot(file.exists(old))
   data.frame(key=job$key,family=job$family,arm=job$arm,scenario=job$scenario,replicate=job$replicate,
     schedule=schedule,result_file=f,old_result_file=old,old_initial_file=job$old_initial_file,
     input_file=job$input_file,input_md5=job$input_md5,
     warnings=length(r$warnings),max_group_rhat=r$diagnostics$max_group_rhat,
     max_element_rhat=r$diagnostics$max_element_rhat,unresolved_rhat=r$diagnostics$unresolved_rhat)
 })
 selection<-do.call(rbind,rows);stopifnot(nrow(selection)==110L,!anyDuplicated(selection$key))
 saveRDS(selection,file.path(study,'selection.rds'))
 write.csv(selection,file.path(study,'selected-manifest.csv'),row.names=FALSE)
 print(selection[selection$warnings>0 | selection$max_group_rhat>1.05 | selection$max_element_rhat>1.05 | selection$unresolved_rhat>0,])
 cat('Selected exactly 110 distinct current fits and paired historical results.\n')
}
