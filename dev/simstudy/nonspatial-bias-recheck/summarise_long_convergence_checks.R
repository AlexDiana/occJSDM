#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
out<-file.path(root,'main-results/summary');settings<-readRDS(file.path(root,'convergence-checks/settings.rds'));keys<-settings$keys
status<-readRDS(file.path(root,'convergence-checks/job-status.rds'));stopifnot(length(status)==length(keys),all(vapply(status,is.character,logical(1))),setequal(unlist(status),keys))
summary<-read.csv(file.path(out,'summary.csv'));rows<-list();elements<-list();warnings<-list();manifest<-list()
for(key in keys) {
 a<-readRDS(file.path(root,'main-results',paste0(key,'-result.rds')));b<-readRDS(file.path(root,'convergence-checks',paste0(key,'-result.rds')))
 stopifnot(identical(a$job$input_file,b$job$input_file),identical(a$job$priors,b$job$priors),identical(b$job$mcmc,settings$mcmc),max(abs(b$reconstruction))<1e-10)
 p<-merge(a$groups,b$groups,by=c('metric','group'),suffixes=c('_original','_long'));stopifnot(max(abs(p$truth_original-p$truth_long))<1e-10)
 p$key<-key;p$change_bias<-p$bias_long-p$bias_original;p$change_in_ten_dataset_average<-p$change_bias/10
 prior<-a$job$prior_name;sc<-readRDS(a$job$input_file)$scenario$label
 original_mean<-summary[summary$scenario==sc & summary$prior==prior,c('metric','group','datasets','bias','mae','rmse')];names(original_mean)[4:6]<-c('original_ten_dataset_average','original_ten_dataset_mae','original_ten_dataset_rmse')
 p<-merge(p,original_mean,by=c('metric','group'));stopifnot(nrow(p)==nrow(a$groups),all(p$datasets==10L))
 p$long_run_replacement_average<-p$original_ten_dataset_average+p$change_bias/p$datasets
 p$change_mae<-p$mae_long-p$mae_original
 p$change_in_ten_dataset_mae<-p$change_mae/p$datasets
 p$long_run_replacement_mae<-p$original_ten_dataset_mae+p$change_in_ten_dataset_mae
 # The primary summary averages squared errors before taking the square root.
 p$change_mse<-p$rmse_long^2-p$rmse_original^2
 p$long_run_replacement_rmse<-sqrt(p$original_ten_dataset_rmse^2+p$change_mse/p$datasets)
 p$change_in_ten_dataset_rmse<-p$long_run_replacement_rmse-p$original_ten_dataset_rmse
 rows[[key]]<-p
 e<-merge(a$elements,b$elements,by=c('metric','element'),suffixes=c('_original','_long'));e$key<-key;elements[[key]]<-e
 warnings[[key]]<-data.frame(key=key,warning=if(length(b$warnings))b$warnings else NA_character_)
 manifest[[key]]<-data.frame(key=key,source_revision=settings$source_revision,input_md5=unname(tools::md5sum(b$job$input_file)),result_md5=unname(tools::md5sum(file.path(root,'convergence-checks',paste0(key,'-result.rds')))),nchain=b$job$mcmc$nchain,nburn=b$job$mcmc$nburn,niter=b$job$mcmc$niter,nthin=b$job$mcmc$nthin,started=as.character(b$started),finished=as.character(b$finished),warnings=length(b$warnings),original_max_element_rhat=max(a$elements$rhat),long_max_element_rhat=max(b$elements$rhat),original_max_group_rhat=max(a$groups$rhat),long_max_group_rhat=max(b$groups$rhat))
}
write.csv(do.call(rbind,rows),file.path(out,'long-run-group-comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,elements),file.path(out,'long-run-parameter-comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,warnings),file.path(out,'long-run-warnings.csv'),row.names=FALSE)
write.csv(do.call(rbind,manifest),file.path(out,'long-run-fit-manifest.csv'),row.names=FALSE)
x<-do.call(rbind,rows);options(width=170);print(x[x$metric %in% c('p','q','occupancy','B0') & x$group %in% c('all','low','medium','high'),c('key','metric','group','bias_original','bias_long','rhat_original','rhat_long','change_in_ten_dataset_average')],row.names=FALSE,digits=4)
cat('Worst individual scored Rhat before:',max(vapply(elements,function(e)max(e$rhat_original),numeric(1))),'after:',max(vapply(elements,function(e)max(e$rhat_long),numeric(1))),'\n')
