#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
repo<-normalizePath(if(length(args))args[1]else'.')
study<-normalizePath(if(length(args)>1L)args[2]else file.path(repo,'dev/simstudy/results/interval-calibration-20261008'))
out<-file.path(repo,'dev/simstudy/interval-calibration/results')
source(file.path(repo,'dev/simstudy/interval-calibration/helpers.R'))
checks<-list();check<-function(name,value){stopifnot(isTRUE(value));checks[[name]]<<-data.frame(check=name,passed=TRUE)}
oldfile<-file.path(repo,'dev/simstudy/results/simstudy-20260821-135619.rds')
old<-subset(readRDS(oldfile)$rows,block=='q')
new<-read.csv(file.path(out,'historical-q-rows.csv'))
# Independent integration of the log-intensity density, not the scorer's CDF.
retained<-integrate(function(v)dnorm(v,1.5,1),log(1.5),Inf,rel.tol=1e-12)$value
check('historical interval rows and order preserved',identical(old$scenario,new$scenario)&&identical(old$element,new$element)&&identical(old$replicate,new$replicate))
check('historical q positive-read truth from density integration',max(abs(new$truth-old$truth*retained))<1e-14)
check('historical posterior interval endpoints preserved',max(abs(new$lower-old$lower),abs(new$upper-old$upper),abs(new$post_mean-old$post_mean))<1e-14)
check('historical corrected coverage decisions',all(new$covered==as.integer(old$lower<=old$truth*retained & old$upper>=old$truth*retained)))

rows<-read.csv(file.path(out,'selected-nonspatial-rows.csv'))
manifest<-read.csv(file.path(out,'selected-fit-provenance.csv'))
check('all 60 selected ordinary fits present',nrow(manifest)==60L&&!anyDuplicated(manifest$key))
for(i in seq_len(nrow(manifest))) {
  m<-manifest[i,];x<-readRDS(m$input)
  check(paste('input identity',m$key),unname(tools::md5sum(m$input))==m$input_md5)
  f<-readRDS(m$fit)$fit
  check(paste('fit identity',m$key),unname(tools::md5sum(m$fit))==m$fit_md5)
  info<-x$sim$data_list$info
  samples<-info[!duplicated(info[c('Site','Sample')]),,drop=FALSE]
  samples<-samples[order(samples$Site,samples$Sample),,drop=FALSE]
  raw<-samples$X_theta;bt<-x$sim$true_params$beta_theta_true
  direct<-rbind(bt[1,]+mean(raw)*bt[2,],sd(raw)*bt[2,])
  check(paste('collection truth preserves predictor',m$key),max(abs(cbind(1,raw)%*%bt-f$X_theta%*%direct))<1e-10)
  scenario<-sub('^design-','',m$key);scenario<-sub('-[0-9]+$','',scenario)
  arm<-if(grepl('-q20$',scenario))'baseline'else sub('^.*-','',scenario)
  sc<-paste(sub('-[^-]+$','',scenario),arm,sep=':')
  rr<-subset(rows,scenario==sc & replicate==x$replicate & group=='all')
  for(block in c('collection_slope','collection_intercept','q')) {
    z<-rr[rr$block==block,];target<-switch(block,collection_slope=direct[2,],collection_intercept=direct[1,],q=as.vector(x$sim$true_params$q_true)*retained)
    check(paste('direct truth',m$key,block),length(target)==nrow(z)&&max(abs(target-z$truth))<1e-12)
    values<-if(block=='q') {
      unlist(lapply(seq_len(ncol(x$sim$true_params$q_true)),function(s)
        lapply(seq_len(nrow(x$sim$true_params$q_true)),function(p)as.vector(f$results_output$q_output[p,s,,]))),recursive=FALSE)
    }else lapply(seq_len(ncol(bt)),function(s)as.vector(f$results_output$beta_theta_output[if(block=='collection_slope')2L else 1L,s,,]))
    bounds<-vapply(values,quantile,numeric(2),probs=c(.025,.975),names=FALSE)
    centers<-vapply(values,mean,numeric(1))
    check(paste('direct interval endpoints',m$key,block),max(abs(bounds[1,]-z$lower),abs(bounds[2,]-z$upper),abs(centers-z$post_mean))<1e-11)
  }
  rm(f);gc(FALSE)
}
selectionfile<-file.path(out,'continuous-selection.csv')
if(file.exists(selectionfile)) {
  selected<-read.csv(selectionfile);scored<-read.csv(file.path(out,'continuous-selected-rows.csv'))
  check('complete continuous selection',nrow(selected)==200L&&!anyDuplicated(selected$key)&&all(table(selected$arm)==100L))
  provenance<-readRDS(selected$result[1])$hashes
  retained<-grepl(paste0(study,'/(source|library)/'),names(provenance),fixed=FALSE)
  check('frozen production source and installed library identity',any(retained)&&identical(unname(tools::md5sum(names(provenance)[retained])),unname(provenance[retained])))
  for(i in seq_len(nrow(selected))) {
    s<-selected[i,];savedfit<-readRDS(s$fit);f<-savedfit$fit
    input<-readRDS(savedfit$job$input);sim<-input$sim[[s$arm]];jp<-sim$true_params$jsdmParams_true
    check(paste('continuous fit identity',s$key),unname(tools::md5sum(s$fit))==s$fit_md5)
    check(paste('continuous input identity',s$key),unname(tools::md5sum(savedfit$job$input))==savedfit$job$input_md5)
    check(paste('continuous zero spatial and correct factors',s$key),f$infos$ps==0L&&dim(f$results_output$jsdm_output$L_output)[1]==if(s$arm=='d2')2L else 0L)
    y2<-input$sim$d2;y0<-input$sim$d0
    check(paste('paired Gaussian residuals',s$key),max(abs((y2$data_list$OTU-y2$true_params$jsdmParams_true$eta)-
      (y0$data_list$OTU-y0$true_params$jsdmParams_true$eta)))<1e-12)
    original<-readRDS(file.path(study,'initial',paste0(s$key,'-result.rds')))
    snapshots<-c('runner-executed.R','helpers-executed.R','simstudy-helper-executed.R')
    originalpaths<-file.path(repo,c('dev/simstudy/interval-calibration/run-continuous.R',
      'dev/simstudy/interval-calibration/helpers.R','tests/testthat/helper-simstudy.R'))
    check(paste('executed analysis source snapshots',s$key),identical(unname(tools::md5sum(file.path(study,snapshots))),
      unname(original$hashes[originalpaths])))
    direct_flag<-function(z,warnings)any(!is.finite(z$rhat)|z$rhat>1.05|!is.finite(z$ess_bulk)|z$ess_bulk<400|
      !is.finite(z$ess_tail)|z$ess_tail<400)||length(warnings)>0L
    initial_flag<-direct_flag(original$rows,original$warnings)
    check(paste('diagnostic-only fit selection',s$key),identical(s$initial_flagged,initial_flag)&&
      identical(s$schedule,if(initial_flag)'long'else'initial')&&
      identical(s$selected_flagged,direct_flag(savedfit$rows,savedfit$warnings)))
    for(block in c('B0','tau')) {
      z<-scored[scored$scenario==s$arm & scored$replicate==s$replicate & scored$block==block,]
      arr<-f$results_output$jsdm_output[[if(block=='B0')'B0_output'else'tau_output']]
      truth<-jp[[block]]
      bounds<-vapply(seq_along(truth),function(k)quantile(as.vector(arr[k,,]),c(.025,.975),names=FALSE),numeric(2))
      centers<-vapply(seq_along(truth),function(k)mean(arr[k,,]),numeric(1))
      check(paste('continuous direct intervals and truth',s$key,block),nrow(z)==length(truth)&&
        max(abs(z$truth-truth),abs(z$lower-bounds[1,]),abs(z$upper-bounds[2,]),abs(z$post_mean-centers))<1e-11)
    }
    rm(savedfit,f);gc(FALSE)
  }
  known_nuisance<-read.csv(file.path(out,'continuous-known-nuisance-rows.csv'))
  known_covariance<-read.csv(file.path(out,'continuous-known-covariance-rows.csv'))
  for(r in 1:100) {
    input<-readRDS(file.path(study,'inputs',sprintf('community-%03d.rds',r)))
    for(arm in c('d0','d2')) {
      sim<-input$sim[[arm]];jp<-sim$true_params$jsdmParams_true
      design<-scale(as.matrix(sim$data_list$info[paste0('X_psi.EnvCov.',1:2)]))
      residual<-sim$data_list$OTU-design%*%jp$B-jp$U%*%jp$L
      n<-nrow(residual)
      variance<-jp$tau^2/(n+jp$tau^2)
      center<-colMeans(residual)*n/(n+jp$tau^2)
      z<-known_nuisance[known_nuisance$scenario==arm & known_nuisance$replicate==r,]
      check(paste('independent known-nuisance posterior',arm,r),max(abs(z$truth-jp$B0),abs(z$post_mean-center),
        abs(z$lower-center+qnorm(.975)*sqrt(variance)),abs(z$upper-center-qnorm(.975)*sqrt(variance)))<1e-12)
      sigma<-diag(jp$tau^2)+input$truth$jsdmParams$sigma_h^2*crossprod(jp$L)
      precision<-chol2inv(chol(sigma+diag(n,ncol(residual))))
      covariance<-diag(ncol(residual))-n*precision
      center<-as.vector(precision%*%colSums(sim$data_list$OTU-design%*%jp$B))
      z<-known_covariance[known_covariance$scenario==arm & known_covariance$replicate==r,]
      check(paste('independent known-covariance posterior',arm,r),max(abs(z$truth-jp$B0),abs(z$post_mean-center),
        abs(z$lower-center+qnorm(.975)*sqrt(diag(covariance))),abs(z$upper-center-qnorm(.975)*sqrt(diag(covariance))))<1e-12)
    }
  }
}
# Cluster-based Monte Carlo uncertainty can be reproduced without helper code.
prefixes<-c('historical-q','selected-nonspatial')
if(file.exists(selectionfile))prefixes<-c(prefixes,'continuous-initial','continuous-selected','continuous-known-nuisance','continuous-known-covariance')
for(prefix in prefixes) {
  rr<-read.csv(file.path(out,paste0(prefix,'-replicates.csv')))
  ss<-read.csv(file.path(out,paste0(prefix,'-summary.csv')))
  for(i in seq_len(nrow(ss))) {
    z<-subset(rr,scenario==ss$scenario[i]&block==ss$block[i]&group==ss$group[i])
    check(paste('community coverage uncertainty',prefix,i),abs(mean(z$covered)-ss$coverage[i])<1e-12 && abs(sd(z$covered)/sqrt(nrow(z))-ss$coverage_mcse[i])<1e-12)
  }
}
# The screened sensitivity must remain available even when all fits fail.
expressions<-as.list(parse(file.path(repo,'dev/simstudy/interval-calibration/summarise.R')))
definition<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('write_screened')),expressions)
stopifnot(length(definition)==1L);eval(definition[[1]])
scratch<-tempfile('interval-empty-');dir.create(scratch)
screened<-write_screened(new[0,,drop=FALSE],new,'empty',scratch)
check('empty screened sensitivity records every missing cell',nrow(screened$summary)==4L&&all(screened$summary$communities==0L)&&all(is.na(screened$summary$coverage)))
unlink(scratch,recursive=TRUE)
write.csv(do.call(rbind,checks),file.path(out,'verification.csv'),row.names=FALSE)
cat(length(checks),'independent checks passed.\n')
