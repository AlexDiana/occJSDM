#!/usr/bin/env Rscript
# Rescore completed immutable fits without changing the running fit protocol.
args<-commandArgs(TRUE);stopifnot(length(args) %in% c(2L,3L))
repo<-normalizePath(args[1]);study<-normalizePath(args[2]);workers<-if(length(args)==3L)as.integer(args[3]) else 1L
stopifnot(workers %in% 1:4)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()));library(occJSDM)
scripts<-file.path(repo,'dev/simstudy/spatial-amplitude-prior')
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/score.R'))
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
source(file.path(scripts,'metrics.R'));source(file.path(scripts,'analysis.R'));source(file.path(scripts,'robust.R'))
hashes<-tools::md5sum(c(file.path(scripts,c('robust.R','rescore.R','ESTIMAND-AMENDMENT.md')),
  file.path(repo,'dev/simstudy/spatial-targeted-recheck/score.R')))
paths<-unlist(lapply(c('inverse_gamma','half_cauchy'),function(prior)
  list.files(file.path(study,prior),pattern='-result[.]rds$',recursive=TRUE,full.names=TRUE)))
# The short implementation pilot is excluded from scientific rescoring.
paths<-paths[!grepl('/pilot/',paths,fixed=TRUE)]
rescore<-function(path) {
  relative<-substring(path,nchar(study)+2L);dest<-file.path(study,'robust-v1',relative)
  md5<-unname(tools::md5sum(path))
  if(file.exists(dest)) {
    a<-readRDS(dest)
    stopifnot(identical(a$scoring_hashes,hashes),identical(a$legacy_result_md5,md5),
      identical(a$source_fit_md5,unname(tools::md5sum(a$source_fit))))
    return(data.frame(file=relative,status='verified existing'))
  }
  a<-read_checked_result(path);raw<-readRDS(a$source_fit)
  stopifnot(identical(raw$job,a$job),identical(raw$fit_hashes,a$fit_hashes))
  input<-readRDS(a$job$input_file)
  stopifnot(identical(unname(tools::md5sum(a$job$input_file)),a$job$input_md5))
  js<-raw$fit$results_output$jsdm_output
  spatial<-robust_field_draws(js$Bs_output,js$idx_ls_output,independent_bases(raw$fit),input$truth$field,js$sigmabs_output)
  trace_difference<-max(abs(spatial$traces-a$spatial$traces))
  stopifnot(is.finite(trace_difference),trace_difference<1e-9)
  # Preserve exactly the original trace values for rank-tie consistency.
  spatial$traces<-a$spatial$traces
  spatial$diagnostics<-do.call(rbind,lapply(seq_len(dim(spatial$traces)[1]),function(i)
    data.frame(quantity=dimnames(spatial$traces)[[1]][i],
      as.list(robust_trace_diagnostics(matrix(spatial$traces[i,,],a$mcmc$niter,a$mcmc$nchain))))))
  a$legacy_reasons<-a$reasons;a$spatial<-spatial;a$field<-spatial$score
  # Label mixed estimands explicitly and remove undefined spatial mean diagnostics.
  for(table in c('groups','elements')) {
    a[[table]]$point_estimator<-'posterior mean'
    id<-a[[table]]$metric=='spatial_sd';stopifnot(sum(id)==1L)
    a[[table]]$point_estimator[id]<-'posterior median'
    a[[table]]$estimate[id]<-spatial$amplitude$median
    a[[table]]$bias[id]<-spatial$amplitude$median-input$truth$field_sd
    for(nm in intersect(c('ess_mean','mcse','chain_gap'),names(a[[table]])))a[[table]][id,nm]<-NA_real_
    for(nm in c('rhat','ess_bulk','constant_chain'))a[[table]][id,nm]<-spatial$diagnostics[1,nm]
    if(table=='groups') {
      a[[table]]$mae[id]<-abs(a[[table]]$bias[id]);a[[table]]$rmse[id]<-abs(a[[table]]$bias[id])
    }
  }
  a$reasons<-robust_reasons(a);a$legacy_result<-normalizePath(path);a$legacy_result_md5<-md5
  a$scoring_hashes<-hashes;a$estimand_version<-'robust-v1';a$rescored_at<-Sys.time()
  a$trace_reconstruction_difference<-trace_difference
  validate_robust_spatial(a)
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  tmp<-paste0(dest,'.tmp');saveRDS(a,tmp);stopifnot(file.rename(tmp,dest))
  cat(relative,'median/quantile extraction complete; diagnostic flags',length(a$reasons),'\n');flush.console()
  data.frame(file=relative,status='rescored')
}
work<-function(path)tryCatch(rescore(path),error=function(e)list(file=path,error=conditionMessage(e)))
answer<-if(workers==1L)lapply(paths,work) else parallel::mclapply(paths,work,mc.cores=workers,
  mc.preschedule=FALSE,mc.set.seed=FALSE)
bad<-!vapply(answer,is.data.frame,logical(1))
if(any(bad)){print(answer[bad]);stop('Robust rescoring failed')}
cat('Verified',length(answer),'completed results under robust-v1.\n')
