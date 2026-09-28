#!/usr/bin/env Rscript
# Independent native reconstruction for the amended median/quantile analysis.
args<-commandArgs(TRUE);stopifnot(length(args) %in% c(2L,3L))
study<-normalizePath(args[1]);summary<-normalizePath(args[2]);workers<-if(length(args)==3L)as.integer(args[3]) else 1L
stopifnot(workers %in% 1:4)
script<-normalizePath(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))]))
repo<-dirname(dirname(dirname(dirname(script))))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()));library(occJSDM)
rows<-read.csv(file.path(summary,'fits.csv'),stringsAsFactors=FALSE)
columns<-c('key','prior','result_file','result_md5','fit_file','fit_md5')
rows<-rows[columns];paired_n<-nrow(rows)
probefile<-file.path(summary,'initialization-selection.csv')
if(file.exists(probefile)) {
  probe<-read.csv(probefile,stringsAsFactors=FALSE);stopifnot(nrow(probe)==1L)
  rows<-rbind(rows,data.frame(key='range6-rep01-binary-k100',prior='half_cauchy_startcheck',
    result_file=probe$file,result_md5=probe$md5,fit_file=probe$fit,fit_md5=probe$fit_md5))
}
stopifnot(nrow(rows)>0L,!anyDuplicated(rows[c('key','prior')]))

audit<-function(i) {
  row<-rows[i,];r<-readRDS(row$result_file);raw<-readRDS(row$fit_file);fit<-raw$fit
  stopifnot(identical(r$estimand_version,'robust-v1'),
    identical(unname(tools::md5sum(row$result_file)),row$result_md5),
    identical(unname(tools::md5sum(row$fit_file)),row$fit_md5),
    identical(unname(tools::md5sum(r$legacy_result)),r$legacy_result_md5),
    identical(unname(tools::md5sum(names(r$scoring_hashes))),unname(r$scoring_hashes)),
    identical(raw$job,r$job),identical(raw$fit_hashes,r$fit_hashes))
  original_sources<-names(r$fit_hashes);resolved_sources<-original_sources
  relocated<-!file.exists(original_sources) & grepl('/dev/simstudy/',original_sources,fixed=TRUE)
  resolved_sources[relocated]<-paste0(repo,sub('^.*(/dev/simstudy/.*)$','\\1',original_sources[relocated]))
  stopifnot(all(file.exists(resolved_sources)),
    identical(unname(tools::md5sum(resolved_sources)),unname(r$fit_hashes)))
  input<-readRDS(r$job$input_file)
  stopifnot(identical(unname(tools::md5sum(r$job$input_file)),r$job$input_md5))
  js<-fit$results_output$jsdm_output;n<-nrow(fit$Xs);S<-dim(js$B0_output)[1]
  ni<-dim(js$B0_output)[2];nc<-dim(js$B0_output)[3]
  stopifnot(ni==r$mcmc$niter,nc==r$mcmc$nchain,fit$infos$n_factors==0L,
    identical(fit$infos$speciesNames,colnames(input$truth$z)),
    max(abs(fit$infos$OTU-input$data[[r$job$arm]]$OTU))==0)
  native<-suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid,fit$infos$list_Xs))
  probability<-fields<-array(NA_real_,c(n*S,ni,nc))
  labels<-c('amplitude',paste0('field_mean_',1:S),paste0('field_rms_',1:S),paste0('field_projection_',1:S))
  traces<-array(NA_real_,c(length(labels),ni,nc),dimnames=list(labels,NULL,NULL))
  tc<-sweep(input$truth$field,2,colMeans(input$truth$field),'-')
  for(ch in 1:nc)for(it in 1:ni) {
    field<-occJSDM:::KsBproduct(native$Ks_all[,,js$idx_ls_output[it,ch]],
      matrix(js$Bs_output[,,it,ch],fit$infos$ps,S),fit$infos$list_Xs$Xs_centers)
    fields[,it,ch]<-c(field)
    eta<-outer(rep(1,n),js$B0_output[,it,ch])+fit$X_psi%*%matrix(js$B_output[,,it,ch],ncol(fit$X_psi),S)+field
    probability[,it,ch]<-c(plogis(eta))
    traces[,it,ch]<-c(js$sigmabs_output[it,ch],colMeans(field),sqrt(colMeans(field^2)),
      colSums(sweep(field,2,colMeans(field),'-')*tc)/colSums(tc^2))
  }
  quantiles<-t(vapply(1:(n*S),function(k)unname(quantile(fields[k,,],c(.025,.5,.975))),numeric(3)))
  field<-matrix(quantiles[,2],n,S);fc<-sweep(field,2,colMeans(field),'-')
  field_scores<-c(raw_rmse=sqrt(mean((field-input$truth$field)^2)),raw_mae=mean(abs(field-input$truth$field)),
    centred_rmse=sqrt(mean((fc-tc)^2)),centred_mae=mean(abs(fc-tc)),
    centred_correlation=cor(c(fc),c(tc)),centred_slope=sum(fc*tc)/sum(tc^2),
    median_field_rms=sqrt(mean(field^2)),true_rms=sqrt(mean(input$truth$field^2)),
    raw_bias=mean(field-input$truth$field),
    raw_coverage=mean(quantiles[,1]<=c(input$truth$field) & c(input$truth$field)<=quantiles[,3]),
    raw_interval_width=mean(quantiles[,3]-quantiles[,1]))
  field_error<-max(abs(field-r$spatial$field_median),abs(quantiles[,1]-c(r$spatial$field_lower)),
    abs(quantiles[,3]-c(r$spatial$field_upper)),abs(field_scores-unlist(r$spatial$score[names(field_scores)])))
  amp<-unname(quantile(js$sigmabs_output,c(.025,.5,.975)))
  amplitude_error<-max(abs(amp-unlist(r$spatial$amplitude[c('lower','median','upper')])))
  trace_error<-max(abs(traces-r$spatial$traces))
  # These formulas are written here independently, not imported from robust.R.
  diag<-function(x) {
    if(any(apply(x,2,function(z)all(z==z[1]))))return(c(rhat=NA,ess_bulk=NA,
      ess_median=NA,ess_q025=NA,ess_q975=NA))
    c(rhat=posterior::rhat(x),ess_bulk=posterior::ess_bulk(x),
      ess_median=unname(posterior::ess_quantile(x,.5)),
      ess_q025=unname(posterior::ess_quantile(x,.025)),ess_q975=unname(posterior::ess_quantile(x,.975)))
  }
  crossing<-function(a,b) {
    if(!identical(is.finite(a),is.finite(b)))return(sum(is.finite(a)!=is.finite(b)))
    good<-is.finite(a);limit<-c(1.05,100,100,100,100)
    sum(((a>limit)!=(b>limit))[good & names(a)=='rhat'])+
      sum(((a<limit)!=(b<limit))[good & names(a)!='rhat'])
  }
  diag_error<-0;native_rhat_difference<-0;native_threshold_crossings<-0L;field_diagnostic_difference<-0
  for(j in seq_along(labels)) {
    expected<-diag(matrix(r$spatial$traces[j,,],ni,nc))
    recorded<-unlist(r$spatial$diagnostics[r$spatial$diagnostics$quantity==labels[j],names(expected)])
    stopifnot(identical(is.na(expected),is.na(recorded)))
    diag_error<-max(diag_error,abs(expected-recorded),na.rm=TRUE)
    nd<-diag(matrix(traces[j,,],ni,nc))
    native_rhat_difference<-max(native_rhat_difference,abs(nd['rhat']-expected['rhat']),na.rm=TRUE)
    native_threshold_crossings<-native_threshold_crossings+crossing(nd,expected)
  }
  fd<-r$spatial$field_diagnostics
  stopifnot(identical(fd$quantity,paste0('field_',seq_len(n*S))))
  for(k in 1:(n*S)) {
    nd<-diag(matrix(fields[k,,],ni,nc));stored<-unlist(fd[k,names(nd)])
    stopifnot(identical(is.na(nd),is.na(stored)))
    native_rhat_difference<-max(native_rhat_difference,abs(nd['rhat']-stored['rhat']),na.rm=TRUE)
    field_diagnostic_difference<-max(field_diagnostic_difference,abs(nd-stored),na.rm=TRUE)
    native_threshold_crossings<-native_threshold_crossings+crossing(nd,stored)
  }
  stopifnot(identical(quantiles[,1]<=c(input$truth$field) & c(input$truth$field)<=quantiles[,3],fd$covered))
  mu<-matrix(apply(probability,1,mean),n,S)
  true<-plogis(outer(rep(1,n),input$truth$B0)+tcrossprod(input$truth$X,input$truth$B)+input$truth$field)
  stopifnot(max(abs(true-input$truth$psi))<1e-12)
  interval<-t(vapply(1:(n*S),function(k)unname(quantile(probability[k,,],c(.025,.975))),numeric(2)))
  el<-r$elements[r$elements$metric=='occupancy',]
  covered<-interval[,1]<=c(true) & c(true)<=interval[,2];stopifnot(identical(covered,el$covered))
  probability_error<-max(abs(mu-r$probability),abs(interval-cbind(el$lower,el$upper)))
  masks<-list(all=rep(TRUE,n*S),low=c(true)<.2,medium=c(true)>=.2 & c(true)<=.8,high=c(true)>.8)
  for(p in c(.01,.05,.25,.75))masks[[paste0('prevalence_',p*100,'pct')]]<-rep(abs(colMeans(true)-p)<1e-10,each=n)
  group_error<-0
  for(g in names(masks)) {
    id<-which(masks[[g]]);error<-c(mu-true)[id]
    expected<-c(bias=mean(error),mae=mean(abs(error)),rmse=sqrt(mean(error^2)),coverage=mean(covered[id]))
    stored<-r$groups[r$groups$metric=='occupancy' & r$groups$group==g,names(expected)]
    stopifnot(nrow(stored)==1L);group_error<-max(group_error,abs(expected-unlist(stored)))
  }
  stopifnot(max(field_error,amplitude_error,trace_error,probability_error,group_error)<1e-9,diag_error<1e-12)
  cat(row$key,row$prior,'median/interval audit verified\n');flush.console()
  data.frame(key=row$key,prior=row$prior,relocated_source_count=sum(relocated),field_error=field_error,
    amplitude_error=amplitude_error,trace_error=trace_error,diagnostic_error=diag_error,
    native_rhat_difference=native_rhat_difference,field_diagnostic_difference=field_diagnostic_difference,
    native_threshold_crossings=native_threshold_crossings,probability_error=probability_error,group_error=group_error,
    result_md5=row$result_md5,fit_md5=row$fit_md5)
}
work<-function(i)tryCatch(audit(i),error=function(e)list(key=rows$key[i],prior=rows$prior[i],error=conditionMessage(e)))
answer<-if(workers==1L)lapply(seq_len(nrow(rows)),work) else parallel::mclapply(seq_len(nrow(rows)),work,
  mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
bad<-!vapply(answer,is.data.frame,logical(1))
if(any(bad)){print(answer[bad]);stop('Independent robust audit failed')}
audit<-do.call(rbind,answer)
write.csv(audit[seq_len(paired_n),],file.path(summary,'independent-audit.csv'),row.names=FALSE)
if(nrow(audit)>paired_n)write.csv(audit[nrow(audit),],file.path(summary,'initialization-audit.csv'),row.names=FALSE)
write.csv(data.frame(file=script,md5=unname(tools::md5sum(script))),file.path(summary,'audit-source.csv'),row.names=FALSE)
gatefile<-file.path(summary,'extension-gate.csv')
if(file.exists(gatefile)) {
  stopifnot(paired_n==18L,nrow(audit)==19L)
  gate<-read.csv(gatefile,stringsAsFactors=FALSE)
  crossings<-sum(audit$native_threshold_crossings)
  decision<-rbind(gate,data.frame(criterion='All paired fits and start check audited without unresolved diagnostic threshold crossings',
    value=crossings,pass=crossings==0L))
  write.csv(decision,file.path(summary,'extension-decision.csv'),row.names=FALSE)
  cat('Audited amended extension decision:',if(all(decision$pass))'PASS' else 'FAIL','\n')
}
cat('All',nrow(rows),'selected fits independently verified.\n')
