#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
arg<-function(key,default){x<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(x))sub(paste0('^--',key,'='),'',x[1]) else default}
here<-normalizePath(arg('study','work/nonspatial-design-recheck-20260919'))
keys<-strsplit(arg('keys',''),',',fixed=TRUE)[[1]];stopifnot(length(keys)>0L,all(nzchar(keys)))
workers<-as.integer(arg('workers','2'));settings<-readRDS(file.path(here,'fits/settings.rds'))
stopifnot(identical(unname(tools::md5sum(names(settings$source_hashes))),unname(settings$source_hashes)))
source(file.path(here,'design_helpers.R'));load_frozen(settings$baseline_root);score_design_fit<-make_design_scorer()
# Use the executed runner's fit-and-score implementation, without executing its jobs.
exprs<-parse(file.path(here,'run_design_recheck.R'))
selected<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('run_job')),as.list(exprs))
stopifnot(length(selected)==1L);eval(selected[[1]],.GlobalEnv)
out<-file.path(here,'long-fits');dir.create(out,showWarnings=FALSE)
jobs<-lapply(keys,function(key){old<-readRDS(file.path(here,'fits',paste0(key,'-result.rds')));j<-old$job;j$mcmc<-list(nchain=4,nburn=6000,niter=12000,nthin=1);j$out<-out;j$source_hashes<-c(j$source_hashes,tools::md5sum(file.path(here,'run_long_design_checks.R')));j})
for(j in jobs)saveRDS(list(job=j,session=sessionInfo(),runner=readLines(file.path(here,'run_long_design_checks.R'))),file.path(out,paste0(j$key,'-settings.rds')))
work<-function(j)tryCatch(run_job(j),error=function(e){cat(j$key,'ERROR:',conditionMessage(e),'\n');list(key=j$key,error=conditionMessage(e))})
status<-parallel::mclapply(jobs,work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(all(vapply(status,is.character,logical(1))));cat('Completed',length(jobs),'longer checks.\n')
