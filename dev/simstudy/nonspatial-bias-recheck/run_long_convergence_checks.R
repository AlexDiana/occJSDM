#!/usr/bin/env Rscript
# Target only the two primary-grid fits with individual scored-parameter
# rank-normalized Rhat above1.1. Original outputs remain unchanged.
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
lib<-file.path(root,'library');source_dir<-file.path(root,'source-main');out<-file.path(root,'convergence-checks');dir.create(out,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()));Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages(library(occJSDM));source(file.path(source_dir,'tests/testthat/helper-simstudy.R'));RcppParallel::setThreadOptions(numThreads=1)
expressions<-parse(file.path(root,'run_nonspatial_recheck_balanced.R'),keep.source=FALSE)
for(name in c('score_fit','run_job')) {
 hit<-Filter(function(e)is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name(name)),as.list(expressions));stopifnot(length(hit)==1);eval(hit[[1]])
}
trace_summary<-function(x){x<-as.matrix(x);c(rhat=posterior::rhat(x),ess_bulk=posterior::ess_bulk(x),ess_mean=posterior::ess_mean(x),mcse=posterior::mcse_mean(x),chain_gap=diff(range(colMeans(x))))}
keys<-c('qfar_K3-q20-08','qfar_K3-q9-09')
jobs<-lapply(keys,function(key){j<-readRDS(file.path(root,'main-results',paste0(key,'-result.rds')))$job;j$out<-out;j$mcmc<-list(nchain=4,nburn=6000,niter=12000,nthin=1);j})
settings<-list(source_revision=readLines(file.path(root,'source-revision.txt')),selection='Individual scored-parameter rank-normalized Rhat above 1.1 in original grid',keys=keys,mcmc=jobs[[1]]$mcmc,threads=1,source_hashes=readRDS(file.path(root,'main-results/settings.rds'))$hashes,runner=readLines(file.path(root,'run_long_convergence_checks.R')),session=sessionInfo())
saveRDS(settings,file.path(out,'settings.rds'))
cl<-parallel::makePSOCKcluster(2,outfile=file.path(out,'workers.log'))
parallel::clusterExport(cl,c('lib','source_dir','score_fit','run_job','trace_summary'))
parallel::clusterEvalQ(cl,{.libPaths(c(lib,.libPaths()));Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');suppressPackageStartupMessages(library(occJSDM));source(file.path(source_dir,'tests/testthat/helper-simstudy.R'));RcppParallel::setThreadOptions(numThreads=1);NULL})
status<-parallel::parLapplyLB(cl,jobs,function(j)tryCatch(run_job(j),error=function(e)list(key=j$key,error=conditionMessage(e))),chunk.size=1L)
parallel::stopCluster(cl);saveRDS(status,file.path(out,'job-status.rds'));stopifnot(all(vapply(status,is.character,logical(1))))
cat('Completed two longer convergence checks.\n')
