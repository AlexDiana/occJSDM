#!/usr/bin/env Rscript
args<-commandArgs(TRUE);repo<-normalizePath(if(length(args))args[1]else '.')
folder<-file.path(repo,'dev/simstudy/theta0-prior-mean');out<-file.path(folder,'results')
raw<-file.path(repo,'dev/simstudy/results/theta0-prior-mean-20261008')
phase<-if(length(args)>1L)args[2]else 'initial';stopifnot(phase %in% c('initial','selected'))
source(file.path(folder,'helpers.R'))
rr<-read.csv(file.path(out,'importance-rows.csv'));cc<-summaries(rr)$communities
contrasts<-list();add<-function(label,a,b,filter=rep(TRUE,nrow(cc))) {
 ids<-c('contamination','design','replicate')
 from<-cc[filter & cc$collection_mean==a[1] & cc$theta_b==a[2],]
 to<-cc[filter & cc$collection_mean==b[1] & cc$theta_b==b[2],]
 z<-merge(from,to,by=ids,suffixes=c('_from','_to'))
 if(!nrow(z))return(invisible(NULL))
 # Both arms must pass on the same input for the screened comparison.
 okay<-aggregate(as.integer(!rr$source_flag & !rr$weight_flag),
   rr[c(ids,'collection_mean','theta_b')],min);names(okay)[6]<-'okay'
 aa<-okay[okay$collection_mean==a[1] & okay$theta_b==a[2],c(ids,'okay')]
 bb<-okay[okay$collection_mean==b[1] & okay$theta_b==b[2],c(ids,'okay')]
 ok<-merge(aa,bb,by=ids,suffixes=c('_from','_to'));ok$screened<-ok$okay_from==1 & ok$okay_to==1
 z<-merge(z,ok[c(ids,'screened')],by=ids)
 metrics<-c('covered','bias','abs_error','width','lower_miss','upper_miss')
 for(v in metrics)z[[paste0('change_',v)]]<-z[[paste0(v,'_to')]]-z[[paste0(v,'_from')]]
 z$change_width_percent<-100*(z$width_to/z$width_from-1)
 z$contrast<-label;contrasts[[length(contrasts)+1L]]<<-z
}
for(b in sort(unique(cc$theta_b)))add(paste0('collection mean 1 to 0, theta b=',b),c(1,b),c(0,b))
if(30 %in% cc$theta_b)for(m in c(0,1))add(paste0('theta b 30 to 20, collection mean=',m),c(m,30),c(m,20))
paired<-do.call(rbind,contrasts);write.csv(paired,file.path(out,'paired-communities.csv'),row.names=FALSE)
paired_summary<-function(p,scope) {
 gg<-split(p,interaction(p$contrast,p$contamination,p$design,drop=TRUE))
 do.call(rbind,lapply(gg,function(z){
   n<-nrow(z);crit<-if(n>1)qt(.975,n-1)else NA_real_
   do.call(rbind,lapply(grep('^change_',names(z),value=TRUE),function(v){
     se<-sd(z[[v]])/sqrt(n);est<-mean(z[[v]])
     data.frame(scope=scope,contrast=z$contrast[1],contamination=z$contamination[1],design=z$design[1],
       communities=n,metric=sub('change_','',v),change=est,mcse=se,lower=est-crit*se,upper=est+crit*se)
   }))
 }))
}
ss<-paired_summary(paired,'all')
if(any(paired$screened))ss<-rbind(ss,paired_summary(subset(paired,screened),'common screened inputs'))
write.csv(ss,file.path(out,'paired-summary.csv'),row.names=FALSE)
if(30 %in% cc$theta_b) {
  ids<-c('contamination','design','replicate')
  a<-subset(paired,contrast=='theta b 30 to 20, collection mean=0')
  b<-subset(paired,contrast=='theta b 30 to 20, collection mean=1')
  ab<-merge(a,b,by=ids,suffixes=c('_mean0','_mean1'))
  interaction_rows<-do.call(rbind,lapply(c('covered','bias','abs_error','width'),function(v)
    data.frame(ab[ids],metric=v,interaction=ab[[paste0('change_',v,'_mean0')]]-
      ab[[paste0('change_',v,'_mean1')]],screened=ab$screened_mean0 & ab$screened_mean1)))
  write.csv(interaction_rows,file.path(out,'factorial-interaction-communities.csv'),row.names=FALSE)
  interaction_summary<-function(z,scope) {
    gg<-split(z,interaction(z$contamination,z$design,z$metric,drop=TRUE))
    do.call(rbind,lapply(gg,function(x){
      n<-nrow(x);est<-mean(x$interaction);se<-sd(x$interaction)/sqrt(n)
      crit<-if(n>1)qt(.975,n-1)else NA_real_
      data.frame(scope=scope,x[1,c('contamination','design','metric')],communities=n,
        interaction=est,mcse=se,lower=est-crit*se,upper=est+crit*se)
    }))
  }
  zz<-interaction_summary(interaction_rows,'all')
  if(any(interaction_rows$screened))zz<-rbind(zz,interaction_summary(subset(interaction_rows,screened),'common screened inputs'))
  write.csv(zz,file.path(out,'factorial-interaction-summary.csv'),row.names=FALSE)
}


files<-list.files(file.path(raw,'fits'),pattern='[.]rds$',full.names=TRUE)
stopifnot(length(files)==3L*2L*length(unique(rr$theta_b)))
if(phase=='selected') {
 replacements<-file.path(raw,'long-fits',basename(files));have<-file.exists(replacements)
 initial_manifest<-read.csv(file.path(out,'factorial-full-fit-provenance.csv'))
 expected<-paste0(initial_manifest$key[initial_manifest$source_flag],'.rds')
 stopifnot(setequal(basename(files[have]),expected));files[have]<-replacements[have]
}
prefix<-if(phase=='initial')'factorial-full-fit'else 'factorial-full-fit-selected'
fullrows<-fulldiag<-comparison<-fm<-list()
old<-read.csv(file.path(out,'selected-fit-provenance.csv'))
for(i in seq_along(files)) {
 saved<-readRDS(files[i]);j<-saved$job;x<-readRDS(j$input);ro<-saved$fit$results_output
 meta<-data.frame(key=j$key,source_key=j$source_key,contamination='low',design=j$design,replicate=j$replicate,
   collection_mean=j$collection_mean,theta_b=j$theta_b)
 diag<-source_diagnostics(saved$fit,saved$warnings);fulldiag[[i]]<-cbind(meta,diag$elements)
 rows<-score_theta(ro$theta0_output,x$truth$params$theta0,rep(1/prod(dim(ro$theta0_output)[-1L]),prod(dim(ro$theta0_output)[-1L])))
 fullrows[[i]]<-cbind(meta,source_flag=diag$flag,weight_flag=FALSE,rows)
 fm[[i]]<-cbind(meta,file=files[i],md5=unname(tools::md5sum(files[i])),warnings=length(saved$warnings),
   max_rhat=diag$max_rhat,min_ess=diag$min_ess,source_flag=diag$flag,
   elapsed_seconds=as.numeric(difftime(saved$finished,saved$started,units='secs')))
 orig<-old[old$key==j$source_key,];prop<-readRDS(orig$fit);p<-prop$fit$results_output
 ni<-dim(p$theta0_output)[2L];nc<-dim(p$theta0_output)[3L]
 wi<-weight_info(joint_log_ratio(p$beta_theta_output,p$theta0_output,j$collection_mean,j$theta_b),ni,nc)
 means<-draw_matrix(p$theta0_output)%*%wi$weights
 approx<-subset(rr,key==j$source_key & collection_mean==j$collection_mean & theta_b==j$theta_b)
 stopifnot(nrow(approx)==10L,max(abs(approx$post_mean-as.vector(means)))<1e-10)
 for(k in 1:10) {
   values<-matrix(p$theta0_output[k,,],ni,nc)
   influence<-matrix(wi$weights*length(wi$weights),ni,nc)*(values-as.numeric(means[k]))
   imp_mcse<-posterior::mcse_mean(influence)
   full_mcse<-posterior::mcse_mean(matrix(ro$theta0_output[k,,],dim(ro$theta0_output)[2L],dim(ro$theta0_output)[3L]))
   comparison[[length(comparison)+1L]]<-cbind(meta,species=k,importance_mean=approx$post_mean[k],
     full_mean=rows$post_mean[k],difference=rows$post_mean[k]-approx$post_mean[k],
     importance_mcse=imp_mcse,full_mcse=full_mcse,
     difference_mcse_upper_bound=imp_mcse+full_mcse,
     importance_width=approx$upper[k]-approx$lower[k],full_width=rows$upper[k]-rows$lower[k],
     lower_difference=rows$lower[k]-approx$lower[k],upper_difference=rows$upper[k]-approx$upper[k],
     source_flag=orig$source_flag,full_flag=diag$flag,weight_flag=wi$weight_flag)
 }
 cat(i,'/',length(files),j$key,'compared\n');flush.console()
}
fr<-do.call(rbind,fullrows);save_tables(fr,prefix,out)
write.csv(do.call(rbind,fulldiag),file.path(out,paste0(prefix,'-element-diagnostics.csv')),row.names=FALSE)
write.csv(do.call(rbind,fm),file.path(out,paste0(prefix,'-provenance.csv')),row.names=FALSE)
cp<-do.call(rbind,comparison)
# Source and refit share starting RNG states. Summed MCSE is a conservative
# uncertainty scale under unknown covariance; there is no independence Z score.
write.csv(cp,file.path(out,paste0('factorial-importance-versus-full-',phase,'.csv')),row.names=FALSE)
cat('Paired community contrasts and full-fit validation saved.\n')
