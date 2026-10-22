#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
out<-file.path(root,'main-results/summary');x<-read.csv(file.path(out,'dataset-groups.csv'))
keys<-c('prior','replicate','metric','group')
results<-list();dataset_rows<-list()
for(level in c('qnear','qfar')) {
 a<-x[x$scenario==paste0(level,'_K3'),];b<-x[x$scenario==paste0(level,'_K6'),]
 p<-merge(a,b,by=keys,suffixes=c('_K3','_K6'))
 stopifnot(length(unique(p$replicate))==10L,all(abs(p$truth_K3-p$truth_K6)<1e-10))
 p$contamination<-level;p$change_bias<-p$bias_K6-p$bias_K3;p$change_mae<-p$mae_K6-p$mae_K3;p$change_mse<-p$rmse_K6^2-p$rmse_K3^2
 dataset_rows[[level]]<-p
 for(g in split(p,interaction(p[c('prior','metric','group')],drop=TRUE)))results[[length(results)+1L]]<-data.frame(contamination=level,prior=g$prior[1],metric=g$metric[1],group=g$group[1],datasets=nrow(g),K3_bias=mean(g$bias_K3),K6_bias=mean(g$bias_K6),change_bias=mean(g$change_bias),paired_se=if(nrow(g)>1)sd(g$change_bias)/sqrt(nrow(g)) else NA_real_,K3_mae=mean(g$mae_K3),K6_mae=mean(g$mae_K6),change_mae=mean(g$change_mae),change_mse=mean(g$change_mse))
}
write.csv(do.call(rbind,results),file.path(out,'pcr-three-vs-six-summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,dataset_rows),file.path(out,'pcr-three-vs-six-dataset-comparison.csv'),row.names=FALSE)
cat('Saved paired K3 versus K6 comparisons; generating truths match.\n')
