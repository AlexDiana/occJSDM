#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/nonspatial-design-recheck-20260919')
files<-list.files(file.path(here,'fits'),'-result.rds$',full.names=TRUE)
stopifnot(length(files)==60L)
use_long<-length(args)>1L && args[2]=='long'
if(use_long){long_files<-file.path(here,'long-fits',basename(files));files[file.exists(long_files)]<-long_files[file.exists(long_files)]}
out<-file.path(here,if(use_long)'summary-long' else 'summary');dir.create(out,showWarnings=FALSE)
settings<-readRDS(file.path(here,'fits/settings.rds'));records<-lapply(files,readRDS)
stopifnot(identical(unname(tools::md5sum(names(settings$source_hashes))),unname(settings$source_hashes)))
levels<-list();paired<-list();manifest<-list();flags<-list();warns<-list();seen<-character();species_rows<-list()
add_rows<-function(rec,scenario,replicate,arm,scope,source_metric) {
 d<-rec$groups[rec$groups$metric==source_metric & rec$groups$group %in% c('all','low','medium','high'),]
 d$metric<-if(source_metric=='occupancy_original_sites')'occupancy' else source_metric
 data.frame(scenario=scenario,replicate=replicate,arm=arm,scope=scope,d)
}
for(i in seq_along(records)) {
 r<-records[[i]];stopifnot(identical(unname(tools::md5sum(names(r$job$source_hashes))),unname(r$job$source_hashes)));input<-readRDS(r$job$input_file);base<-readRDS(r$job$baseline_result)
 fit<-readRDS(sub('-result.rds$','-fit.rds',files[i]))$fit
 tp<-input$sim$true_params;psi<-fit$results_output$psi_output[1:100,,drop=FALSE];truth_psi<-plogis(tp$jsdmParams_true$eta[1:100,,drop=FALSE])
 info<-input$sim$data_list$info;si<-info[!duplicated(info$Sample),];si<-si[order(si$Sample),]
 ix<-input$design$original_samples;theta<-fit$results_output$theta_output[ix,,drop=FALSE]
 truth_theta<-plogis(cbind(1,si$X_theta[ix])%*%tp$beta_theta_true)
 collection<-data.frame(metric='collection_probability',group='all',n_elements=length(theta),truth=mean(truth_theta),estimate=mean(theta),bias=mean(theta-truth_theta),mae=mean(abs(theta-truth_theta)),rmse=sqrt(mean((theta-truth_theta)^2)),rhat=NA_real_,ess_bulk=NA_real_,ess_mean=NA_real_,mcse=NA_real_,chain_gap=NA_real_)
 sc<-input$scenario$label;rep<-input$replicate;arm<-r$job$arm;key<-r$job$key;bk<-paste(sc,rep)
 if(!bk %in% seen) {
  for(scope in c('original_sites','all_fitted_sites'))levels[[length(levels)+1L]]<-add_rows(base,sc,rep,'baseline',scope,'occupancy')
  for(m in c('p','q'))levels[[length(levels)+1L]]<-add_rows(base,sc,rep,'baseline','all_fitted_sites',m)
  levels[[length(levels)+1L]]<-add_rows(base,sc,rep,'baseline','original_samples','collection_probability')
  species_rows[[length(species_rows)+1L]]<-data.frame(scenario=sc,replicate=rep,arm='baseline',base$species)
  seen<-c(seen,bk)
 }
 levels[[length(levels)+1L]]<-add_rows(r,sc,rep,arm,'original_sites','occupancy_original_sites')
 levels[[length(levels)+1L]]<-add_rows(r,sc,rep,arm,'all_fitted_sites','occupancy')
 for(m in c('p','q'))levels[[length(levels)+1L]]<-add_rows(r,sc,rep,arm,'all_fitted_sites',m)
 levels[[length(levels)+1L]]<-data.frame(scenario=sc,replicate=rep,arm=arm,scope='original_samples',collection)
 species_rows[[length(species_rows)+1L]]<-data.frame(scenario=sc,replicate=rep,arm=arm,species=colnames(input$sim$data_list$OTU),prevalence=colMeans(truth_psi),occupancy_estimate=colMeans(psi),bias=colMeans(psi-truth_psi),mae=colMeans(abs(psi-truth_psi)),rmse=sqrt(colMeans((psi-truth_psi)^2)))
 for(m in c('occupancy','p','q','collection_probability')) {
  original<-base$groups[base$groups$metric==m & base$groups$group %in% c('all','low','medium','high'),]
  updated<-if(m=='collection_probability')collection else r$groups[r$groups$metric==if(m=='occupancy')'occupancy_original_sites' else m,]
  p<-merge(original,updated,by='group',suffixes=c('_baseline','_new'))
  stopifnot(nrow(p)==nrow(original),all(p$n_elements_baseline==p$n_elements_new),max(abs(p$truth_baseline-p$truth_new))<1e-10)
  p$scenario<-sc;p$replicate<-rep;p$arm<-arm;p$metric<-m
  p$change_bias<-p$bias_new-p$bias_baseline;p$change_mae<-p$mae_new-p$mae_baseline;p$change_mse<-p$rmse_new^2-p$rmse_baseline^2
  p$mean_probability_mcse_difference<-sqrt(p$mcse_baseline^2+p$mcse_new^2)
  paired[[length(paired)+1L]]<-p
 }
 manifest[[i]]<-data.frame(key=key,scenario=sc,replicate=rep,arm=arm,n=input$scenario$n,M=input$scenario$M,P=input$scenario$P,K=input$scenario$K,
   longer_fit=r$job$mcmc$nchain==4,nchain=r$job$mcmc$nchain,nburn=r$job$mcmc$nburn,niter=r$job$mcmc$niter,nthin=r$job$mcmc$nthin,
   started=as.character(r$started),finished=as.character(r$finished),warnings=length(r$warnings),max_element_rhat=max(r$elements$rhat),max_group_rhat=max(r$groups$rhat),
   max_occupancy_group_mcse=max(r$groups$mcse[r$groups$metric=='occupancy_original_sites']),input_md5=unname(tools::md5sum(r$job$input_file)),result_md5=unname(tools::md5sum(files[i])),
   known_u_projection_error=if(is.null(r$control_check))NA_real_ else unname(r$control_check['projection']))
 ef<-r$elements[r$elements$rhat>1.1,];if(nrow(ef))flags[[length(flags)+1L]]<-data.frame(key=key,type='parameter',metric=ef$metric,item=as.character(ef$element),rhat=ef$rhat,mcse=ef$mcse)
 gf<-r$groups[r$groups$rhat>1.05,];if(nrow(gf))flags[[length(flags)+1L]]<-data.frame(key=key,type='group',metric=gf$metric,item=gf$group,rhat=gf$rhat,mcse=gf$mcse)
 warns[[i]]<-data.frame(key=key,warning=if(length(r$warnings))r$warnings else NA_character_)
}
d<-do.call(rbind,levels);p<-do.call(rbind,paired);m<-do.call(rbind,manifest)
se<-function(x)sd(x)/sqrt(length(x))
s<-do.call(rbind,lapply(split(d,interaction(d[c('scenario','arm','scope','metric','group')],drop=TRUE)),function(x){stopifnot(nrow(x)==10L);data.frame(x[1,c('scenario','arm','scope','metric','group')],datasets=nrow(x),truth=mean(x$truth),estimate=mean(x$estimate),bias=mean(x$bias),bias_se=se(x$bias),mae=mean(x$mae),mae_se=se(x$mae),rmse=sqrt(mean(x$rmse^2)),mean_probability_mcse=sqrt(sum(x$mcse^2))/nrow(x),max_rhat=if(all(is.na(x$rhat)))NA_real_ else max(x$rhat,na.rm=TRUE))}))
c<-do.call(rbind,lapply(split(p,interaction(p[c('scenario','arm','metric','group')],drop=TRUE)),function(x){stopifnot(nrow(x)==10L);data.frame(x[1,c('scenario','arm','metric','group')],datasets=nrow(x),change_bias=mean(x$change_bias),change_bias_se=se(x$change_bias),change_mae=mean(x$change_mae),change_mae_se=se(x$change_mae),mae_improved_datasets=sum(x$change_mae<0),change_mse=mean(x$change_mse),approx_mean_probability_mcse=sqrt(sum(x$mean_probability_mcse_difference^2))/nrow(x))}))
write.csv(d,file.path(out,'dataset-groups.csv'),row.names=FALSE);write.csv(s,file.path(out,'summary.csv'),row.names=FALSE)
write.csv(p,file.path(out,'paired-dataset-comparisons.csv'),row.names=FALSE);write.csv(c,file.path(out,'paired-comparisons.csv'),row.names=FALSE)
write.csv(m,file.path(out,'fit-manifest.csv'),row.names=FALSE);write.csv(do.call(rbind,warns),file.path(out,'warnings.csv'),row.names=FALSE)
fl<-if(length(flags))do.call(rbind,flags) else data.frame(key=character(),type=character(),metric=character(),item=character(),rhat=numeric(),mcse=numeric())
write.csv(fl,file.path(out,'diagnostic-flags.csv'),row.names=FALSE)
write.csv(do.call(rbind,species_rows),file.path(out,'species-occupancy-original-sites.csv'),row.names=FALSE)
options(width=180)
x<-subset(s,scope=='original_sites' & metric=='occupancy');x[c('truth','estimate','bias','bias_se','mae','mae_se','rmse','mean_probability_mcse')]<-100*x[c('truth','estimate','bias','bias_se','mae','mae_se','rmse','mean_probability_mcse')]
cat('Occupancy levels at the same original 100 sites, percentages / percentage points:\n');print(x,row.names=FALSE,digits=4)
x<-subset(c,metric=='occupancy');x[c('change_bias','change_bias_se','change_mae','change_mae_se','approx_mean_probability_mcse')]<-100*x[c('change_bias','change_bias_se','change_mae','change_mae_se','approx_mean_probability_mcse')]
cat('\nPaired changes, percentage points (negative MAE change = improvement):\n');print(x,row.names=FALSE,digits=4)
cat('\nDiagnostic flags:\n');print(fl,row.names=FALSE,digits=4)
cat('\nCompleted fits:',nrow(m),'; fits with warnings:',sum(m$warnings>0),'; warning messages:',sum(m$warnings),'\n')
