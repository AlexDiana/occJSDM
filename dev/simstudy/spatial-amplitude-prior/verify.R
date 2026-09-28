#!/usr/bin/env Rscript
# Full-draw audit using native spatial projection, independent of the scorer.
a<-commandArgs(TRUE);stopifnot(length(a)==2L)
study<-normalizePath(a[1]);summary<-normalizePath(a[2])
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()));library(occJSDM)
rows<-read.csv(file.path(summary,'fits.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(rows)>0L,!anyDuplicated(rows[c('key','prior')]))
audit<-function(i) {
  row<-rows[i,];r<-readRDS(row$result_file);raw<-readRDS(row$fit_file);fit<-raw$fit
  stopifnot(identical(unname(tools::md5sum(row$result_file)),row$result_md5),
    identical(unname(tools::md5sum(row$fit_file)),row$fit_md5),
    identical(raw$job,r$job),identical(raw$fit_hashes,r$fit_hashes),
    identical(unname(tools::md5sum(names(r$fit_hashes))),unname(r$fit_hashes)))
  input<-readRDS(r$job$input_file)
  stopifnot(identical(unname(tools::md5sum(r$job$input_file)),r$job$input_md5))
  js<-fit$results_output$jsdm_output;n<-nrow(fit$Xs);S<-dim(js$B0_output)[1]
  ni<-dim(js$B0_output)[2];nc<-dim(js$B0_output)[3]
  stopifnot(ni==r$mcmc$niter,nc==r$mcmc$nchain,fit$infos$n_factors==0L,
    identical(fit$infos$speciesNames,colnames(input$truth$z)),
    max(abs(fit$infos$OTU-input$data[[r$job$arm]]$OTU))==0)
  native<-suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid,fit$infos$list_Xs))
  probability<-array(NA_real_,c(n*S,ni,nc));fieldsum<-matrix(0,n,S)
  labels<-c('amplitude',paste0('field_mean_',1:S),paste0('field_rms_',1:S),paste0('field_projection_',1:S))
  traces<-array(NA_real_,c(length(labels),ni,nc),dimnames=list(labels,NULL,NULL))
  tc<-sweep(input$truth$field,2,colMeans(input$truth$field),'-')
  for(ch in 1:nc)for(it in 1:ni) {
    field<-occJSDM:::KsBproduct(native$Ks_all[,,js$idx_ls_output[it,ch]],
      matrix(js$Bs_output[,,it,ch],fit$infos$ps,S),fit$infos$list_Xs$Xs_centers)
    eta<-outer(rep(1,n),js$B0_output[,it,ch])+fit$X_psi%*%matrix(js$B_output[,,it,ch],ncol(fit$X_psi),S)+field
    probability[,it,ch]<-as.vector(plogis(eta));fieldsum<-fieldsum+field
    traces[,it,ch]<-c(js$sigmabs_output[it,ch],colMeans(field),sqrt(colMeans(field^2)),
      colSums(sweep(field,2,colMeans(field),'-')*tc)/colSums(tc^2))
  }
  field<-fieldsum/(ni*nc);fc<-sweep(field,2,colMeans(field),'-')
  field_scores<-c(raw_rmse=sqrt(mean((field-input$truth$field)^2)),raw_mae=mean(abs(field-input$truth$field)),
    centred_rmse=sqrt(mean((fc-tc)^2)),centred_mae=mean(abs(fc-tc)),
    centred_correlation=cor(c(fc),c(tc)),centred_slope=sum(fc*tc)/sum(tc^2),
    posterior_mean_rms=sqrt(mean(field^2)),true_rms=sqrt(mean(input$truth$field^2)))
  field_error<-max(abs(field-r$spatial$field_mean),abs(field_scores-unlist(r$spatial$score[names(field_scores)])))
  trace_error<-max(abs(traces-r$spatial$traces))
  diag_error<-0;native_rhat_difference<-0;native_threshold_crossings<-0L;folded_rank_changes<-0L
  for(j in seq_along(labels)) {
    native_x<-matrix(traces[j,,],ni,nc)
    x<-matrix(r$spatial$traces[j,,],ni,nc)
    expected<-c(rhat=posterior::rhat(x),ess_mean=posterior::ess_mean(x))
    recorded<-unlist(r$spatial$diagnostics[r$spatial$diagnostics$quantity==labels[j],names(expected)])
    stopifnot(identical(is.na(expected),is.na(recorded)))
    diag_error<-max(diag_error,abs(expected-recorded),na.rm=TRUE)
    native_rhat<-posterior::rhat(native_x)
    native_rhat_difference<-max(native_rhat_difference,abs(native_rhat-expected['rhat']),na.rm=TRUE)
    if(is.finite(native_rhat) && is.finite(expected['rhat'])) {
      native_threshold_crossings<-native_threshold_crossings+as.integer((native_rhat>1.05)!=(expected['rhat']>1.05))
      if(abs(native_rhat-expected['rhat'])>1e-8)
        folded_rank_changes<-folded_rank_changes+sum(rank(abs(native_x-median(native_x)))!=rank(abs(x-median(x))))
    }
  }
  mu<-matrix(apply(probability,1,mean),n,S)
  true<-plogis(outer(rep(1,n),input$truth$B0)+tcrossprod(input$truth$X,input$truth$B)+input$truth$field)
  stopifnot(max(abs(true-input$truth$psi))<1e-12)
  interval<-t(vapply(1:(n*S),function(k)unname(quantile(probability[k,,],c(.025,.975))),numeric(2)))
  el<-r$elements[r$elements$metric=='occupancy',]
  covered<-interval[,1]<=c(true) & c(true)<=interval[,2]
  stopifnot(identical(covered,el$covered))
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
  stopifnot(max(field_error,trace_error,probability_error,group_error)<1e-9,diag_error<1e-6)
  cat(row$key,row$prior,'verified\n');flush.console()
  data.frame(key=row$key,prior=row$prior,field_error=field_error,trace_error=trace_error,
    diagnostic_error=diag_error,native_rhat_difference=native_rhat_difference,
    native_threshold_crossings=native_threshold_crossings,folded_rank_changes=folded_rank_changes,probability_error=probability_error,group_error=group_error,
    result_md5=row$result_md5,fit_md5=row$fit_md5)
}
# Sequential audit keeps memory use bounded alongside any remaining fits.
answer<-lapply(seq_len(nrow(rows)),audit)
write.csv(do.call(rbind,answer),file.path(summary,'independent-audit.csv'),row.names=FALSE)
write.csv(data.frame(file=normalizePath('verify.R'),md5=unname(tools::md5sum('verify.R'))),
  file.path(summary,'audit-source.csv'),row.names=FALSE)
cat('All',nrow(rows),'selected fits independently verified.\n')
