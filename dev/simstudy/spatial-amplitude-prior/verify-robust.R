#!/usr/bin/env Rscript
# Independent native reconstruction for the amended median/quantile analysis.
args<-commandArgs(TRUE);stopifnot(length(args) %in% c(2L,3L))
study<-normalizePath(args[1]);summary<-normalizePath(args[2]);workers<-if(length(args)==3L)as.integer(args[3]) else 1L
stopifnot(workers %in% 1:4)
script<-normalizePath(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))]))
audit_hash<-unname(tools::md5sum(script))
repo<-dirname(dirname(dirname(dirname(script))))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()));library(occJSDM)
rows<-read.csv(file.path(summary,'fits.csv'),stringsAsFactors=FALSE)
columns<-c('key','prior','result_file','result_md5','fit_file','fit_md5')
rows<-rows[columns];paired_n<-nrow(rows)
for(name in c('extension-decision.csv','independent-audit.csv','initialization-audit.csv')) {
  old<-file.path(summary,name);if(file.exists(old))stopifnot(file.remove(old))
}
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
  draw_diag<-function(x) {
    constant<-any(apply(x,2,function(z)all(z==z[1])))
    call<-function(f)if(constant)NA_real_ else suppressWarnings(as.numeric(f(x)))
    c(rhat=call(posterior::rhat),ess_bulk=call(posterior::ess_bulk),
      ess_median=call(function(z)posterior::ess_quantile(z,.5)),
      ess_q025=call(function(z)posterior::ess_quantile(z,.025)),
      ess_q975=call(function(z)posterior::ess_quantile(z,.975)),
      mcse_median=call(function(z)posterior::mcse_quantile(z,.5)),
      mcse_q025=call(function(z)posterior::mcse_quantile(z,.025)),
      mcse_q975=call(function(z)posterior::mcse_quantile(z,.975)),
      chain_median_gap=diff(range(apply(x,2,median))),constant_chain=as.numeric(constant))
  }
  crossing<-function(a,b) {
    check<-c('rhat','ess_bulk','ess_median','ess_q025','ess_q975');a<-a[check];b<-b[check]
    if(!identical(is.finite(a),is.finite(b)))return(sum(is.finite(a)!=is.finite(b)))
    good<-is.finite(a);limit<-c(1.05,100,100,100,100)
    sum(((a>limit)!=(b>limit))[good & names(a)=='rhat'])+
      sum(((a<limit)!=(b<limit))[good & names(a)!='rhat'])
  }
  diag_error<-0;native_rhat_difference<-0;native_threshold_crossings<-0L;field_diagnostic_difference<-0
  for(j in seq_along(labels)) {
    expected<-draw_diag(matrix(r$spatial$traces[j,,],ni,nc))
    recorded<-unlist(r$spatial$diagnostics[r$spatial$diagnostics$quantity==labels[j],names(expected)])
    stopifnot(identical(is.na(expected),is.na(recorded)))
    diag_error<-max(diag_error,abs(expected-recorded),na.rm=TRUE)
    nd<-draw_diag(matrix(traces[j,,],ni,nc))
    native_rhat_difference<-max(native_rhat_difference,abs(nd['rhat']-expected['rhat']),na.rm=TRUE)
    native_threshold_crossings<-native_threshold_crossings+crossing(nd,expected)
  }
  fd<-r$spatial$field_diagnostics
  stopifnot(identical(fd$quantity,paste0('field_',seq_len(n*S))))
  # Reproduce the exact independent-matrix representation for strict recorded
  # diagnostic checks, then separately audit native reconstruction/tie sensitivity.
  xs<-fit$infos$list_Xs
  kernel<-function(a,b,l)exp(-(outer(a[,1],b[,1],'-')^2+outer(a[,2],b[,2],'-')^2)/(2*l^2))
  bases<-lapply(fit$infos$l_s_grid,function(l) {
    lower<-t(chol(kernel(xs$X_tilde,xs$X_tilde,l)+diag(1e-5,fit$infos$ps)))
    t(solve(lower,t(kernel(xs$X_s,xs$X_tilde,l))))[xs$Xs_index,,drop=FALSE]
  })
  field_diagnostic_error<-0;pointwise_reconstruction_error<-0
  for(s in seq_len(S)) {
    exact<-matrix(NA_real_,n,ni*nc)
    for(ch in seq_len(nc))for(g in unique(js$idx_ls_output[,ch])) {
      its<-which(js$idx_ls_output[,ch]==g)
      blocks<-split(its,ceiling(seq_along(its)/100))
      for(block in blocks)exact[,(ch-1L)*ni+block]<-bases[[g]]%*%
        matrix(js$Bs_output[,s,block,ch,drop=FALSE],fit$infos$ps,length(block))
    }
    for(site in seq_len(n)) {
      k<-(s-1L)*n+site;x<-matrix(exact[site,],ni,nc);nx<-matrix(fields[k,,],ni,nc)
      pointwise_reconstruction_error<-max(pointwise_reconstruction_error,abs(x-nx))
      expected<-draw_diag(x);stored<-unlist(fd[k,names(expected)])
      stopifnot(identical(is.na(expected),is.na(stored)))
      field_diagnostic_error<-max(field_diagnostic_error,abs(expected-stored),na.rm=TRUE)
      nd<-draw_diag(nx)
      native_rhat_difference<-max(native_rhat_difference,abs(nd['rhat']-stored['rhat']),na.rm=TRUE)
      field_diagnostic_difference<-max(field_diagnostic_difference,abs(nd-stored),na.rm=TRUE)
      native_threshold_crossings<-native_threshold_crossings+crossing(nd,stored)
    }
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
  range_draws<-fit$infos$l_s_grid[js$idx_ls_output]
  rq<-quantile(range_draws,c(.025,.975),names=FALSE);rb<-mean(range_draws)-input$truth$range
  expected_range<-c(estimate=mean(range_draws),bias=rb,mae=abs(rb),rmse=abs(rb),
    coverage=as.numeric(rq[1]<=input$truth$range && input$truth$range<=rq[2]),interval_width=rq[2]-rq[1])
  stored_range<-r$groups[r$groups$metric=='range' & r$groups$group=='all',names(expected_range)]
  stopifnot(nrow(stored_range)==1L)
  group_error<-max(group_error,abs(expected_range-unlist(stored_range)))
  stopifnot(max(field_error,amplitude_error,trace_error,pointwise_reconstruction_error,probability_error,group_error)<1e-9,
    max(diag_error,field_diagnostic_error)<1e-12)
  cat(row$key,row$prior,'median/interval audit verified\n');flush.console()
  data.frame(key=row$key,prior=row$prior,relocated_source_count=sum(relocated),field_error=field_error,
    amplitude_error=amplitude_error,trace_error=trace_error,diagnostic_error=diag_error,
    native_rhat_difference=native_rhat_difference,field_diagnostic_difference=field_diagnostic_difference,
    field_diagnostic_error=field_diagnostic_error,pointwise_reconstruction_error=pointwise_reconstruction_error,
    native_threshold_crossings=native_threshold_crossings,probability_error=probability_error,group_error=group_error,
    result_md5=row$result_md5,fit_md5=row$fit_md5)
}
work<-function(i)tryCatch(audit(i),error=function(e)list(key=rows$key[i],prior=rows$prior[i],error=conditionMessage(e)))
answer<-if(workers==1L)lapply(seq_len(nrow(rows)),work) else parallel::mclapply(seq_len(nrow(rows)),work,
  mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
bad<-!vapply(answer,is.data.frame,logical(1))
if(any(bad)){print(answer[bad]);stop('Independent robust audit failed')}
audit<-do.call(rbind,answer)
stopifnot(identical(unname(tools::md5sum(script)),audit_hash))
write.csv(audit[seq_len(paired_n),],file.path(summary,'independent-audit.csv'),row.names=FALSE)
if(nrow(audit)>paired_n)write.csv(audit[nrow(audit),],file.path(summary,'initialization-audit.csv'),row.names=FALSE)
write.csv(data.frame(file=script,md5=audit_hash),file.path(summary,'audit-source.csv'),row.names=FALSE)
gatefile<-file.path(summary,'extension-gate.csv')
if(file.exists(gatefile)) {
  stopifnot(paired_n==18L,nrow(audit)==19L)
  scripts<-file.path(repo,'dev/simstudy/spatial-amplitude-prior')
  source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/analysis.R'))
  source(file.path(scripts,'analysis.R'));source(file.path(scripts,'robust.R'))
  provenance<-readRDS(file.path(summary,'provenance.rds'))
  stopifnot(identical(unname(tools::md5sum(names(provenance$script_hashes))),unname(provenance$script_hashes)))
  for(hashes in list(provenance$dependencies$dependency_hashes,provenance$dependencies$control_hashes))
    stopifnot(identical(unname(tools::md5sum(names(hashes))),unname(hashes)))
  stopifnot(all(rows$result_file %in% names(provenance$dependencies$result_hashes)),
    identical(unname(provenance$dependencies$result_hashes[rows$result_file]),rows$result_md5))
  selected<-lapply(rows$result_file[seq_len(paired_n)],readRDS)
  names(selected)<-paste(rows$key[seq_len(paired_n)],rows$prior[seq_len(paired_n)],sep=':')
  stopifnot(identical(unname(tools::md5sum(rows$result_file)),rows$result_md5))
  recovered_rows<-do.call(rbind,lapply(seq_len(paired_n),function(i)
    robust_result_row(selected[[i]],rows$prior[i],if(selected[[i]]$mcmc$nchain==4L)'long' else 'initial')))
  audited_probe<-readRDS(rows$result_file[nrow(rows)])
  gate<-robust_extension_gate(pair_rows(recovered_rows),selected,audited_probe)
  stored_gate<-read.csv(gatefile,stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(gate,stored_gate,tolerance=1e-12,check.attributes=FALSE)))
  crossings<-sum(audit$native_threshold_crossings)
  decision<-rbind(gate,data.frame(criterion='All paired fits and start check audited without unresolved diagnostic threshold crossings',
    value=crossings,pass=crossings==0L))
  write.csv(decision,file.path(summary,'extension-decision.csv'),row.names=FALSE)
  write.csv(rows,file.path(summary,'extension-inputs.csv'),row.names=FALSE)
  cat('Audited amended extension decision:',if(all(decision$pass))'PASS' else 'FAIL','\n')
}
cat('All',nrow(rows),'selected fits independently verified.\n')
