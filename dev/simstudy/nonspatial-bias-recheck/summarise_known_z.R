args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(if(length(args))args[1] else '.')
files<-list.files(file.path(root,'known-z'),pattern='-result.rds$',full.names=TRUE)
z<-do.call(rbind,lapply(files,function(f){a<-readRDS(f);data.frame(scenario=a$scenario,replicate=a$replicate,a$groups)}))
out<-file.path(root,'main-results','summary');dir.create(out,showWarnings=FALSE,recursive=TRUE)
se<-function(x)if(length(x)>1)sd(x)/sqrt(length(x)) else NA_real_
summary<-do.call(rbind,lapply(split(z,interaction(z[c('scenario','metric','group')],drop=TRUE)),function(x)data.frame(x[1,c('scenario','metric','group')],datasets=nrow(x),truth=mean(x$truth),estimate=mean(x$estimate),bias=mean(x$bias),between_dataset_se=se(x$bias),mae=mean(x$mae),rmse=sqrt(mean(x$rmse^2)),combined_mcse=sqrt(sum(x$mcse^2))/nrow(x),max_rhat=max(x$rhat))))
x<-read.csv(file.path(out,'dataset-groups.csv'))
x<-x[x$prior=='var2',]
paired<-merge(x,z,by=c('scenario','replicate','metric','group'),suffixes=c('_full','_known_z'))
paired$change_bias<-paired$bias_full-paired$bias_known_z
paired$change_mae<-paired$mae_full-paired$mae_known_z
paired_summary<-do.call(rbind,lapply(split(paired,interaction(paired[c('scenario','metric','group')],drop=TRUE)),function(x)data.frame(x[1,c('scenario','metric','group')],datasets=nrow(x),full_model_bias=mean(x$bias_full),known_occupancy_bias=mean(x$bias_known_z),difference=mean(x$change_bias),paired_se=se(x$change_bias),change_mae=mean(x$change_mae))))
write.csv(z,file.path(out,'known-occupancy-dataset-groups.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'known-occupancy-summary.csv'),row.names=FALSE)
write.csv(paired_summary,file.path(out,'known-occupancy-paired-comparison.csv'),row.names=FALSE)
cat('Summarised',length(files),'known-occupancy controls\n')
