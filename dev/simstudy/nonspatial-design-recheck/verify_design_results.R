#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/nonspatial-design-recheck-20260919')
use_long<-length(args)>1L && args[2]=='long';out<-file.path(here,if(use_long)'summary-long' else 'summary')
d<-read.csv(file.path(out,'dataset-groups.csv'));s<-read.csv(file.path(out,'summary.csv'));p<-read.csv(file.path(out,'paired-comparisons.csv'));species<-read.csv(file.path(out,'species-occupancy-original-sites.csv'))
files<-list.files(file.path(here,'fits'),'-result.rds$',full.names=TRUE);stopifnot(length(files)==60L)
if(use_long){lf<-file.path(here,'long-fits',basename(files));files[file.exists(lf)]<-lf[file.exists(lf)]}
manifest<-read.csv(file.path(here,'input-manifest.csv'));stopifnot(nrow(manifest)==60L)
near<-function(x,y){stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),max(abs(x-y))<1e-10)}
seen<-character();checked<-0L
check_fit<-function(result_file,wanted_arm){
 r<-readRDS(result_file);a<-readRDS(r$job$input_file);fit<-readRDS(sub('-result.rds$','-fit.rds',result_file))$fit
 sc<-a$scenario$label;rep<-a$replicate
 truth<-plogis(a$sim$true_params$jsdmParams_true$eta);est<-fit$results_output$psi_output
 for(wanted_scope in c('original_sites','all_fitted_sites')){
  sites<-if(wanted_scope=='original_sites')1:100 else seq_len(nrow(truth));t<-truth[sites,,drop=FALSE];e<-est[sites,,drop=FALSE]
  for(g in c('all','low','medium','high')){
   sel<-switch(g,all=matrix(TRUE,nrow(t),ncol(t)),low=t<.2,medium=t>=.2&t<=.8,high=t>.8)
   z<-subset(d,scenario==sc & replicate==rep & d$arm==wanted_arm & d$scope==wanted_scope & metric=='occupancy' & group==g)
   stopifnot(nrow(z)==1L,z$n_elements==sum(sel));err<-e[sel]-t[sel]
   near(c(z$truth,z$estimate,z$bias,z$mae,z$rmse),c(mean(t[sel]),mean(e[sel]),mean(err),mean(abs(err)),sqrt(mean(err^2))))
  }
 }
 si<-a$sim$data_list$info;si<-si[!duplicated(si$Sample),];si<-si[order(si$Sample),]
 ix<-if(wanted_arm=='baseline')1:200 else a$design$original_samples
 tc<-plogis(cbind(1,si$X_theta[ix])%*%a$sim$true_params$beta_theta_true);ec<-fit$results_output$theta_output[ix,,drop=FALSE];err<-ec-tc
 z<-subset(d,scenario==sc & replicate==rep & d$arm==wanted_arm & metric=='collection_probability')
 stopifnot(nrow(z)==1L);near(c(z$truth,z$estimate,z$bias,z$mae,z$rmse),c(mean(tc),mean(ec),mean(err),mean(abs(err)),sqrt(mean(err^2))))
 sp<-species[species$scenario==sc & species$replicate==rep & species$arm==wanted_arm,];stopifnot(nrow(sp)==10L)
 near(sp$prevalence,colMeans(truth[1:100,]));near(sp$mae,colMeans(abs(est[1:100,]-truth[1:100,])))
 checked<<-checked+1L
}
for(f in files){
 r<-readRDS(f);check_fit(f,r$job$arm)
 a<-readRDS(r$job$input_file);mi<-manifest[manifest$key==r$job$key,];stopifnot(nrow(mi)==1L,identical(unname(tools::md5sum(r$job$input_file)),mi$input_md5),identical(unname(tools::md5sum(r$job$baseline_result)),mi$baseline_result_md5),identical(unname(tools::md5sum(a$design$baseline_input)),mi$baseline_input_md5))
 if(!r$job$baseline_result %in% seen){check_fit(r$job$baseline_result,'baseline');seen<-c(seen,r$job$baseline_result)}
}
for(i in seq_len(nrow(s))){
 z<-s[i,];ix<-rep(TRUE,nrow(d));for(k in c('scenario','arm','scope','metric','group'))ix<-ix & d[[k]]==z[[k]];x<-d[ix,];stopifnot(nrow(x)==10L)
 near(c(z$bias,z$mae,z$rmse,z$bias_se,z$mae_se),c(mean(x$bias),mean(x$mae),sqrt(mean(x$rmse^2)),sd(x$bias)/sqrt(10),sd(x$mae)/sqrt(10)))
}
for(i in seq_len(nrow(p))){
 z<-p[i,];wanted_scope<-if(z$metric=='occupancy')'original_sites' else if(z$metric=='collection_probability')'original_samples' else 'all_fitted_sites'
 a<-subset(d,scenario==z$scenario & arm==z$arm & metric==z$metric & group==z$group & d$scope==wanted_scope);b<-subset(d,scenario==z$scenario & arm=='baseline' & metric==z$metric & group==z$group & d$scope==wanted_scope)
 a<-a[order(a$replicate),];b<-b[order(b$replicate),];stopifnot(nrow(a)==10L,nrow(b)==10L,identical(a$replicate,b$replicate));dm<-a$mae-b$mae;db<-a$bias-b$bias
 near(c(z$change_mae,z$change_bias,z$change_mae_se,z$change_bias_se,z$change_mse),c(mean(dm),mean(db),sd(dm)/sqrt(10),sd(db)/sqrt(10),mean(a$rmse^2-b$rmse^2)))
 stopifnot(z$mae_improved_datasets==sum(dm<0))
}
stopifnot(checked==80L,nrow(species)==800L)
cat('PASS: 80 saved-fit probability means, all occupancy bands, original-sample collection errors, 800 species rows, aggregate and paired errors, and unchanged baseline/input hashes.\n')
