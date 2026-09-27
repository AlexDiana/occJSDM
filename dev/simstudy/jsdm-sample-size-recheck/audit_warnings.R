args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/jsdm-sample-size-20260919')
rows<-list()
for(folder in c('fits','long-fits'))for(file in list.files(file.path(here,folder),'-result.rds$',full.names=TRUE)) {
 r<-readRDS(file)
 if(length(r$warnings)) {
  fit<-readRDS(sub('-result.rds','-fit.rds',file,fixed=TRUE))$fit
  for(w in r$warnings) {
   nm<-sub(".*'([^']+)'.*",'\\1',w);x<-fit$results_output$jsdm_output[[nm]]
   if(is.null(x))x<-fit$results_output$jsdm_output[[paste0(nm,'_output')]]
   rh<-if(is.matrix(x))posterior::rhat(x) else NA_real_
   rows[[length(rows)+1L]]<-data.frame(key=r$job$key,schedule=folder,warning=w,rank_rhat_if_scalar=rh)
  }
 }
}
result<-if(length(rows))do.call(rbind,rows) else data.frame(key=character(),schedule=character(),warning=character(),rank_rhat_if_scalar=numeric())
write.csv(result,file.path(here,'summary','warnings.csv'),row.names=FALSE);print(result,row.names=FALSE)
