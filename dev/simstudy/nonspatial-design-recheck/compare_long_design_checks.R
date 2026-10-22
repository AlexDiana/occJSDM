#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/nonspatial-design-recheck-20260919')
a<-read.csv(file.path(here,'summary/summary.csv'));b<-read.csv(file.path(here,'summary-long/summary.csv'))
z<-merge(a,b,by=c('scenario','arm','scope','metric','group'),suffixes=c('_original','_with_longer_fits'))
z$change_bias<-z$bias_with_longer_fits-z$bias_original;z$change_mae<-z$mae_with_longer_fits-z$mae_original;z$change_rmse<-z$rmse_with_longer_fits-z$rmse_original
out<-file.path(here,'summary-long');write.csv(z,file.path(out,'longer-fit-sensitivity.csv'),row.names=FALSE)
files<-list.files(file.path(here,'long-fits'),'-result.rds$',full.names=TRUE);rows<-list();parameters<-list()
for(f in files){
 b<-readRDS(f);a<-readRDS(file.path(here,'fits',basename(f)))
 x<-merge(a$groups,b$groups,by=c('metric','group'),suffixes=c('_original','_longer'));x$key<-b$job$key;x$change_bias<-x$bias_longer-x$bias_original;x$change_mae<-x$mae_longer-x$mae_original;rows[[b$job$key]]<-x
 x<-merge(a$elements,b$elements,by=c('metric','element'),suffixes=c('_original','_longer'));x$key<-b$job$key;parameters[[b$job$key]]<-x
}
write.csv(do.call(rbind,rows),file.path(out,'longer-fit-group-comparison.csv'),row.names=FALSE);write.csv(do.call(rbind,parameters),file.path(out,'longer-fit-parameter-comparison.csv'),row.names=FALSE)
sel<-subset(z,scope=='original_sites' & metric=='occupancy');cat('Longer fits:',length(files),'\nLargest shift in a ten-dataset signed occupancy error:',100*max(abs(sel$change_bias)),'percentage points\nLargest shift in a ten-dataset occupancy MAE:',100*max(abs(sel$change_mae)),'percentage points\n')
