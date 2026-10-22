#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
arg<-function(k,d){z<-grep(paste0('^--',k,'='),args,value=TRUE);if(length(z))sub(paste0('^--',k,'='),'',z[1])else d}
repo<-normalizePath(arg('repo','.'));bs<-as.numeric(strsplit(arg('b','20'),',',fixed=TRUE)[[1]])
shard<-as.integer(arg('shard','1'));shards<-as.integer(arg('shards','1'))
stopifnot(all(bs %in% c(20,30)),shard>=1L,shard<=shards,shards<=4L)
base<-file.path(repo,'dev/simstudy/results/interval-calibration-20261008')
raw<-file.path(repo,'dev/simstudy/results/theta0-prior-mean-20261008')
folder<-file.path(repo,'dev/simstudy/theta0-prior-mean')
.libPaths(c(file.path(base,'library'),.libPaths()));suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package('occJSDM'))==normalizePath(file.path(base,'library/occJSDM')))
RcppParallel::setThreadOptions(numThreads=1)
source(file.path(folder,'helpers.R'))
production<-list.files(file.path(base,'source'),pattern='\\.(R|cpp|h)$',recursive=TRUE,full.names=TRUE)
current<-sub(paste0(base,'/source/'),paste0(repo,'/'),production,fixed=TRUE)
stopifnot(all(unname(tools::md5sum(production))==unname(tools::md5sum(current))))
sources<-c(production,file.path(find.package('occJSDM'),'libs/occJSDM.so'),file.path(folder,c('helpers.R','run-long.R')))
hashes<-tools::md5sum(sources)
for(n in c('helpers.R','run-long.R')) {
 dest<-file.path(raw,paste0('executed-',n));src<-file.path(folder,n)
 if(file.exists(dest))stopifnot(unname(tools::md5sum(dest))==unname(tools::md5sum(src)))else {
   tmp<-paste0(dest,'.',Sys.getpid(),'.tmp');stopifnot(file.copy(src,tmp),file.rename(tmp,dest))
 }
}
mcmc<-list(nchain=4,nburn=6000,niter=12000,nthin=1)
manifest<-read.csv(file.path(repo,'dev/simstudy/interval-calibration/results/selected-fit-provenance.csv'))
manifest<-manifest[manifest$key %in% c('design-qnear_K6-q20-01','design-qnear_K6-field4-01','design-qnear_K6-sites300-01'),]
stopifnot(nrow(manifest)==3L)
jobs<-list()
for(i in seq_len(3L))for(m in c(0,1))for(b in bs) {
 key<-paste0(manifest$key[i],'-mean',m,'-b',b)
 jobs[[key]]<-list(key=key,source_key=manifest$key[i],input=manifest$input[i],input_md5=manifest$input_md5[i],
   design=c('baseline','field4','sites300')[match(manifest$key[i],c('design-qnear_K6-q20-01','design-qnear_K6-field4-01','design-qnear_K6-sites300-01'))],
   collection_mean=m,theta_b=b,replicate=1L)
}
flagged<-read.csv(file.path(folder,'results/full-fit-provenance.csv'))
keys<-flagged$key[flagged$source_flag]
stopifnot(length(keys)==4L,all(keys %in% names(jobs)))
jobs<-jobs[keys]
jobs<-jobs[(seq_along(jobs)-1L)%%shards==shard-1L]
out<-file.path(raw,'long-fits');dir.create(out,showWarnings=FALSE)
status<-list()
for(j in jobs) {
 dest<-file.path(out,paste0(j$key,'.rds'))
 if(file.exists(dest)) {
   old<-readRDS(dest);stopifnot(identical(old$hashes,hashes),identical(old$job,j),identical(old$mcmc,mcmc))
   status[[j$key]]<-'cached';next
 }
 stopifnot(unname(tools::md5sum(j$input))==j$input_md5);x<-readRDS(j$input)
 assign('.Random.seed',x$fit_rng,envir=.GlobalEnv);warnings<-character();started<-Sys.time()
 cat(format(started),j$key,'started\n');flush.console()
 fit<-withCallingHandlers(suppressMessages(fit_input(x,mean_fitter(j$collection_mean),mcmc,j$theta_b)),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 finished<-Sys.time();stopifnot(fit$infos$ps==0L,fit$infos$model=='two_stage',dim(fit$results_output$theta0_output)[2]==12000L,
   dim(fit$results_output$theta0_output)[3]==4L)
 result<-list(job=j,mcmc=mcmc,hashes=hashes,revision=readLines(file.path(base,'source-revision.txt')),
   fit=fit,warnings=warnings,started=started,finished=finished,session=sessionInfo())
 tmp<-paste0(dest,'.',Sys.getpid(),'.tmp');saveRDS(result,tmp);stopifnot(file.rename(tmp,dest))
 status[[j$key]]<-'complete'
 cat(format(finished),j$key,'finished',round(as.numeric(difftime(finished,started,units='secs')),1),'seconds; warnings',length(warnings),'\n');flush.console()
 rm(fit,result);gc(FALSE)
}
saveRDS(status,file.path(raw,paste0('status-long-',paste(bs,collapse='_'),'-shard',shard,'.rds')))
cat('Completed',length(jobs),'validation fits.\n')
