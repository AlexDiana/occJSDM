#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/jsdm-sample-size-20260919')
out<-file.path(here,'summary');dir.create(out,showWarnings=FALSE)
manifest<-read.csv(file.path(here,'input-manifest.csv'));groups<-species<-elements<-points<-fits<-list()
for(k in seq_len(nrow(manifest))) {
 m<-manifest[k,];base<-file.path(here,'fits',paste0(m$key,'-result.rds'));long<-file.path(here,'long-fits',paste0(m$key,'-result.rds'))
 path<-if(file.exists(long))long else base;stopifnot(file.exists(path));r<-readRDS(path)
 stopifnot(r$job$input_md5==m$input_md5,all(is.finite(r$estimate)))
 labels<-data.frame(key=m$key,n=m$n,replicate=m$replicate)
 groups[[k]]<-cbind(labels,r$groups);species[[k]]<-cbind(labels,r$species);elements[[k]]<-cbind(labels,r$elements)
 points[[k]]<-cbind(labels,expand.grid(site=1:100,species=r$species$species),truth=as.vector(r$truth[1:100,]),estimate=as.vector(r$estimate[1:100,]))
 fits[[k]]<-cbind(labels,schedule=if(file.exists(long))'long' else 'initial',chains=r$job$mcmc$nchain,burn=r$job$mcmc$nburn,retained=r$job$mcmc$niter,warnings=length(r$warnings),max_group_rhat=max(r$groups$rhat),max_element_rhat=max(r$elements$rhat),min_element_ess=min(r$elements$ess_mean),max_group_mcse=max(r$groups$mcse),input_md5=r$job$input_md5,result_md5=unname(tools::md5sum(path)),fit_md5=unname(tools::md5sum(sub('-result.rds','-fit.rds',path,fixed=TRUE))),seconds=as.numeric(difftime(r$finished,r$started,units='secs')))
}
groups<-do.call(rbind,groups);species<-do.call(rbind,species);elements<-do.call(rbind,elements);points<-do.call(rbind,points);fits<-do.call(rbind,fits)
write<-function(x,name)write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE)
write(groups,'dataset-groups');write(species,'dataset-species');write(elements,'element-diagnostics');write(points,'original-site-probabilities');write(fits,'fit-manifest')
se<-function(x)sd(x)/sqrt(length(x))
summaries<-lapply(split(groups,interaction(groups$n,groups$scope,groups$group,drop=TRUE)),function(x) {
 stopifnot(nrow(x)==10L,length(unique(x$replicate))==10L)
 vals<-list();for(v in c('truth','estimate','bias','mae','rmse')){vals[[v]]<-mean(x[[v]]);vals[[paste0(v,'_se')]]<-se(x[[v]])}
 cbind(x[1,c('n','scope','group')],communities=nrow(x),min_cells=min(x$cells),max_cells=max(x$cells),as.data.frame(vals),combined_mean_mcse=sqrt(sum(x$mcse^2))/nrow(x),max_group_rhat=max(x$rhat))
})
summary<-do.call(rbind,summaries);rownames(summary)<-NULL;write(summary,'summary')
paired<-list()
for(pair in list(c(100,300),c(100,1000),c(300,1000)))for(g in c('all','low','middle','high')) {
 a<-groups[groups$n==pair[1]&groups$scope=='original100'&groups$group==g,];b<-groups[groups$n==pair[2]&groups$scope=='original100'&groups$group==g,]
 a<-a[order(a$replicate),];b<-b[order(b$replicate),];stopifnot(identical(a$replicate,b$replicate),identical(a$cells,b$cells),max(abs(a$truth-b$truth))<1e-12)
 paired[[length(paired)+1L]]<-data.frame(from=pair[1],to=pair[2],group=g,replicate=a$replicate,reduction_mae=a$mae-b$mae,change_bias=b$bias-a$bias)
}
paired<-do.call(rbind,paired);write(paired,'paired-datasets')
pairs<-do.call(rbind,lapply(split(paired,interaction(paired$from,paired$to,paired$group,drop=TRUE)),function(x) {
 avg<-mean(x$reduction_mae);err<-se(x$reduction_mae)
 cbind(x[1,c('from','to','group')],communities=nrow(x),improved=sum(x$reduction_mae>0),reduction_mae=avg,se=err,lower=avg-qt(.975,9)*err,upper=avg+qt(.975,9)*err)
}));rownames(pairs)<-NULL;write(pairs,'paired-summary')
# Each plotted bin first averages within a community, then across communities.
points$bin<-pmin(9L,floor(points$truth*10))
bins<-do.call(rbind,lapply(split(points,interaction(points$n,points$replicate,points$bin,drop=TRUE)),function(x)cbind(x[1,c('n','replicate','bin')],cells=nrow(x),truth=mean(x$truth),estimate=mean(x$estimate))))
write(bins,'dataset-calibration')
calibration<-do.call(rbind,lapply(split(bins,interaction(bins$n,bins$bin,drop=TRUE)),function(x)cbind(x[1,c('n','bin')],communities=nrow(x),truth=mean(x$truth),estimate=mean(x$estimate),estimate_se=se(x$estimate))))
write(calibration,'calibration')
cat('Summarized',nrow(fits),'fits,',nrow(groups),'groups,',nrow(points),'original-site probabilities.\n')
print(summary[summary$scope=='original100' & summary$group=='all',c('n','bias','mae','mae_se','max_group_rhat')],row.names=FALSE)
print(pairs[pairs$group=='all',],row.names=FALSE)
cat('Maximum element Rhat:',max(fits$max_element_rhat),'; warnings:',sum(fits$warnings),'\n')
