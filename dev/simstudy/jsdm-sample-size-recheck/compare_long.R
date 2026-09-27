args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/jsdm-sample-size-20260919')
files<-list.files(file.path(here,'long-fits'),'-result.rds$',full.names=TRUE);stopifnot(length(files)>0)
rows<-list()
for(f in files) {
 b<-readRDS(f);a<-readRDS(file.path(here,'fits',basename(f)))
 stopifnot(identical(a$job$input_md5,b$job$input_md5),identical(a$groups[c('scope','group')],b$groups[c('scope','group')]))
 rows[[length(rows)+1L]]<-data.frame(key=b$job$key,n=b$job$n,replicate=b$job$replicate,b$groups[c('scope','group')],initial_bias=a$groups$bias,long_bias=b$groups$bias,initial_mae=a$groups$mae,long_mae=b$groups$mae,change_bias=b$groups$bias-a$groups$bias,change_mae=b$groups$mae-a$groups$mae,initial_rhat=a$groups$rhat,long_rhat=b$groups$rhat)
}
x<-do.call(rbind,rows);write.csv(x,file.path(here,'summary','long-check-comparisons.csv'),row.names=FALSE)
changes<-aggregate(x[c('change_bias','change_mae')],x[c('n','scope','group')],sum);changes$change_bias<-changes$change_bias/10;changes$change_mae<-changes$change_mae/10
write.csv(changes,file.path(here,'summary','long-check-summary-changes.csv'),row.names=FALSE)
cat('Largest change in ten-community average bias:',100*max(abs(changes$change_bias)),'points; MAE:',100*max(abs(changes$change_mae)),'points.\n')
