#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
arg <- function(key,default) { x<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(x))sub(paste0('^--',key,'='),'',x[1]) else default }
here <- normalizePath(arg('study','work/jsdm-sample-size-20260919'))
root <- normalizePath(file.path(here,'../nonspatial-bias-recheck-20260914'))
design <- normalizePath(file.path(here,'../nonspatial-design-recheck-20260919'))
source(file.path(here,'helpers.R'));frozen <- load_frozen(root)
mode <- arg('mode','run');workers <- as.integer(arg('workers','3'))
stopifnot(mode %in% c('prepare','pilot','run','long'),workers>=1)
inputdir <- file.path(here,'inputs');dir.create(inputdir,showWarnings=FALSE)
out <- file.path(here,switch(mode,pilot='pilot',long='long-fits','fits'));dir.create(out,showWarnings=FALSE)
mcmc <- switch(mode,pilot=list(nchain=2,nburn=60,niter=80,nthin=1),long=list(nchain=4,nburn=6000,niter=12000,nthin=1),list(nchain=2,nburn=3000,niter=5000,nthin=1))
source_hashes <- c(frozen$hashes,tools::md5sum(file.path(here,c('helpers.R','run_study.R'))))
jobs <- list();manifest <- list()
for(r in 1:10) {
 afile <- file.path(root,'operational-k6',sprintf('qnear_K6-data-%02d.rds',r))
 bfile <- file.path(design,'inputs',sprintf('qnear_K6-sites300-%02d.rds',r))
 a <- readRDS(afile);b <- readRDS(bfile)
 for(n in c(100L,300L,1000L)) {
  input <- make_input(a,b,n);key <- sprintf('n%04d-%02d',n,r)
  file <- file.path(inputdir,paste0(key,'.rds'))
  if(file.exists(file))stopifnot(identical(readRDS(file),input)) else saveRDS(input,file)
  manifest[[key]] <- data.frame(key=key,n=n,replicate=r,base_seed=a$seed,extension_seed=input$extension_seed,input_md5=unname(tools::md5sum(file)),source100_md5=unname(tools::md5sum(afile)),source300_md5=unname(tools::md5sum(bfile)))
  if(mode!='pilot'||r==1L) jobs[[key]] <- list(key=key,n=n,replicate=r,input_file=file,input_md5=unname(tools::md5sum(file)),mcmc=mcmc,source_hashes=source_hashes,out=out)
 }
}
write.csv(do.call(rbind,manifest),file.path(here,'input-manifest.csv'),row.names=FALSE)
if(mode=='prepare') { cat('Prepared and checked 30 paired exact-state inputs.\n');quit(status=0) }
selected <- arg('keys','')
if(nchar(selected)) { keys<-strsplit(selected,',',fixed=TRUE)[[1]];stopifnot(all(keys %in% names(jobs)));jobs<-jobs[keys] }
settings <- list(source_revision=frozen$source_revision,source_hashes=source_hashes,mcmc=mcmc,threads=1,session=sessionInfo(),started=Sys.time(),runner=readLines(file.path(here,'run_study.R')),helpers=readLines(file.path(here,'helpers.R')))
sfile <- file.path(out,'settings.rds')
if(file.exists(sfile)) { old<-readRDS(sfile);stopifnot(identical(old$mcmc,settings$mcmc),identical(old$source_hashes,source_hashes)) } else saveRDS(settings,sfile)
run_job <- function(job) {
 dest <- file.path(job$out,paste0(job$key,'-result.rds'))
 if(file.exists(dest)) { stopifnot(identical(readRDS(dest)$job,job));return(job$key) }
 stopifnot(identical(unname(tools::md5sum(job$input_file)),job$input_md5))
 input <- readRDS(job$input_file);fitfile<-file.path(job$out,paste0(job$key,'-fit.rds'))
 cat(format(Sys.time()),job$key,'started\n');flush.console()
 if(file.exists(fitfile)) { saved<-readRDS(fitfile);stopifnot(identical(saved$job,job));fit<-saved$fit;warnings<-saved$warnings;started<-saved$started;finished<-saved$finished } else {
  assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
  warnings<-character();started<-Sys.time()
  fit <- withCallingHandlers(suppressMessages(occJSDM::runOccJSDM(input$data,listParams=list(n_factors=input$n_factors,n_lattrait=input$n_lattrait),occCovariates=input$covariates,collCovariates=NULL,spatCovariates=NULL,MCMCparams=job$mcmc)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  finished<-Sys.time()
  tmp<-paste0(fitfile,'.tmp');saveRDS(list(fit=fit,job=job,warnings=warnings,started=started,finished=finished),tmp);stopifnot(file.rename(tmp,fitfile))
 }
 scores <- score_fit(fit,input)
 result <- c(list(job=job,warnings=warnings,started=started,finished=finished),scores)
 tmp <- paste0(dest,'.tmp');saveRDS(result,tmp);stopifnot(file.rename(tmp,dest))
 cat(format(Sys.time()),job$key,'complete;',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'s;',length(warnings),'warnings; max group Rhat',round(max(scores$groups$rhat),3),'; max element Rhat',round(max(scores$elements$rhat),3),'\n');flush.console()
 job$key
}
work <- function(j)tryCatch(run_job(j),error=function(e){cat(j$key,'ERROR:',conditionMessage(e),'\n');list(key=j$key,error=conditionMessage(e))})
status <- if(workers==1L)lapply(jobs,work) else parallel::mclapply(jobs,work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(status,file.path(out,paste0('job-status-',format(Sys.time(),'%Y%m%d%H%M%S'),'.rds')))
stopifnot(length(status)==length(jobs),all(vapply(status,is.character,logical(1))))
cat('Completed',length(status),mode,'fits.\n')
