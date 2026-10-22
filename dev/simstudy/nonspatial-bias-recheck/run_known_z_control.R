#!/usr/bin/env Rscript
# Paired diagnostic: observe the exact latent occupancy states from the
# collection-prior datasets, then fit the same non-spatial community model.
args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(if(length(args))args[1] else '.')
.libPaths(c(file.path(root,'library'),.libPaths()))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages(library(occJSDM));RcppParallel::setThreadOptions(numThreads=1)
out<-file.path(root,'known-z');dir.create(out,showWarnings=FALSE)
files<-list.files(file.path(root,'main-results'),pattern='^(traits_isolated|occupancy)-data-[0-9]+.rds$',full.names=TRUE)
stopifnot(length(files)==20)
settings<-list(source_revision=readLines(file.path(root,'source-revision.txt')),mcmc=list(nchain=2,nburn=3000,niter=5000,nthin=1),threads=1,session=sessionInfo(),runner=readLines(file.path(root,'run_known_z_control.R')))
saveRDS(settings,file.path(out,'settings.rds'))
trace_summary<-function(x)c(rhat=posterior::rhat(x),ess_mean=posterior::ess_mean(x),mcse=posterior::mcse_mean(x),chain_gap=abs(mean(x[,1])-mean(x[,2])))
for(file in files) {
 a<-readRDS(file);key<-sprintf('%s-%02d',a$scenario$label,a$replicate)
 dest<-file.path(out,paste0(key,'-result.rds'));if(file.exists(dest))next
 info<-a$sim$data_list$info;info<-info[!duplicated(info$Site),,drop=FALSE];info<-info[order(info$Site),,drop=FALSE]
 info$Sample<-NULL;info$Primer<-NULL
 data<-list(info=info,OTU=a$sim$true_params$z_true,traits=a$sim$data_list$traits)
 colnames(data$OTU)<-colnames(a$sim$data_list$OTU)
 assign('.Random.seed',a$fit_rng,envir=.GlobalEnv)
 warnings<-character();started<-Sys.time()
 fit<-withCallingHandlers(suppressMessages(runOccJSDM(data,listParams=list(n_factors=a$scenario$d,n_lattrait=a$scenario$gt),occCovariates=paste0('X_psi.EnvCov.',1:a$scenario$ncov_psi),collCovariates=NULL,spatCovariates=NULL,MCMCparams=settings$mcmc)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(fit$infos$model=='binary',fit$infos$ps==0)
 jo<-fit$results_output$jsdm_output;jp<-a$sim$true_params$jsdmParams_true
 n<-a$scenario$n;S<-a$scenario$S;ni<-5000;nc<-2
 true_eta<-sweep(fit$X_psi%*%jp$B+jp$U%*%jp$L,2,jp$B0,'+')
 stopifnot(max(abs(true_eta-jp$eta))<1e-10)
 truth<-plogis(true_eta)
 traces<-list();rows<-list()
 add<-function(metric,group,t,x){mu<-rowMeans(matrix(x,nrow=length(t)));trace<-apply(x,c(2,3),mean);e<-mu-as.vector(t);traces[[paste(metric,group,sep=':')]]<<-trace;rows[[length(rows)+1L]]<<-data.frame(metric=metric,group=group,truth=mean(t),estimate=mean(mu),bias=mean(e),mae=mean(abs(e)),rmse=sqrt(mean(e^2)),as.list(trace_summary(trace)))}
 add('B0','all',jp$B0,jo$B0_output)
 pd<-array(NA_real_,c(n*S,ni,nc))
 for(ch in 1:nc)for(it in 1:ni){eta<-sweep(fit$X_psi%*%jo$B_output[,,it,ch]+jo$U_output[,,it,ch]%*%jo$L_output[,,it,ch],2,jo$B0_output[,it,ch],'+');pd[,it,ch]<-as.vector(plogis(eta))}
 reconstructed<-matrix(rowMeans(matrix(pd,nrow=n*S)),n,S)
 # Binary fits store eta-derived occupancy in the public psi running mean.
 if(!is.null(fit$results_output$psi_output))stopifnot(max(abs(reconstructed-fit$results_output$psi_output))<1e-10)
 add('occupancy','all',as.vector(truth),pd)
 for(band in c('low','medium','high')){ix<-which(switch(band,low=truth<.2,medium=truth>=.2 & truth<=.8,high=truth>.8));if(length(ix))add('occupancy',band,as.vector(truth)[ix],pd[ix,,,drop=FALSE])}
 result<-list(key=key,scenario=a$scenario$label,replicate=a$replicate,input=file,groups=do.call(rbind,rows),traces=traces,warnings=warnings,started=started,finished=Sys.time())
 saveRDS(list(fit=fit,data=data,result=result),file.path(out,paste0(key,'-fit.rds')))
 saveRDS(result,dest)
 cat(format(Sys.time()),key,'complete;',length(warnings),'warnings\n');flush.console()
}
