#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
out<-file.path(root,'main-results/summary');keys<-readRDS(file.path(root,'convergence-checks/settings.rds'))$keys
summary<-read.csv(file.path(out,'summary.csv'));rows<-list();elements<-list();warnings<-list()
for(key in keys) {
 a<-readRDS(file.path(root,'main-results',paste0(key,'-result.rds')));b<-readRDS(file.path(root,'convergence-checks',paste0(key,'-result.rds')))
 p<-merge(a$groups,b$groups,by=c('metric','group'),suffixes=c('_original','_long'));stopifnot(max(abs(p$truth_original-p$truth_long))<1e-10)
 p$key<-key;p$change_bias<-p$bias_long-p$bias_original;p$change_in_ten_dataset_average<-p$change_bias/10
 prior<-a$job$prior_name;sc<-readRDS(a$job$input_file)$scenario$label
 original_mean<-summary[summary$scenario==sc & summary$prior==prior,c('metric','group','datasets','bias')];names(original_mean)[4]<-'original_ten_dataset_average'
 p<-merge(p,original_mean,by=c('metric','group'));p$long_run_replacement_average<-p$original_ten_dataset_average+p$change_bias/p$datasets;rows[[key]]<-p
 e<-merge(a$elements,b$elements,by=c('metric','element'),suffixes=c('_original','_long'));e$key<-key;elements[[key]]<-e
 warnings[[key]]<-data.frame(key=key,warning=if(length(b$warnings))b$warnings else NA_character_)
}
write.csv(do.call(rbind,rows),file.path(out,'long-run-group-comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,elements),file.path(out,'long-run-parameter-comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,warnings),file.path(out,'long-run-warnings.csv'),row.names=FALSE)
x<-do.call(rbind,rows);options(width=170);print(x[x$metric %in% c('p','q','occupancy','B0') & x$group %in% c('all','low','medium','high'),c('key','metric','group','bias_original','bias_long','rhat_original','rhat_long','change_in_ten_dataset_average')],row.names=FALSE,digits=4)
cat('Worst individual scored Rhat before:',max(vapply(elements,function(e)max(e$rhat_original),numeric(1))),'after:',max(vapply(elements,function(e)max(e$rhat_long),numeric(1))),'\n')
