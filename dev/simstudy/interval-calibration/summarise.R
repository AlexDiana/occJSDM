#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
repo<-normalizePath(if(length(args))args[1]else'.')
study<-normalizePath(if(length(args)>1L)args[2]else file.path(repo,'dev/simstudy/results/interval-calibration-20261008'))
out<-file.path(repo,'dev/simstudy/interval-calibration/results')
source(file.path(repo,'dev/simstudy/interval-calibration/helpers.R'))
flag_rows<-function(z)any(!is.finite(z$rhat)|z$rhat>1.05|!is.finite(z$ess_bulk)|z$ess_bulk<400|!is.finite(z$ess_tail)|z$ess_tail<400)
write_screened<-function(rows,primary_rows,prefix,out) {
  reference<-summarise_intervals(primary_rows)
  result<-if(nrow(rows))summarise_intervals(rows)else lapply(reference,function(x)x[0,,drop=FALSE])
  key<-function(x)paste(x$scenario,x$block,x$group,sep=':')
  missing<-reference$summary[!key(reference$summary)%in%key(result$summary),,drop=FALSE]
  if(nrow(missing)) {
    missing[,setdiff(names(missing),c('scenario','block','group'))]<-NA_real_
    missing$communities<-0L
    result$summary<-rbind(result$summary,missing)
  }
  for(nm in names(result))write.csv(result[[nm]],file.path(out,paste0(prefix,'-',nm,'.csv')),row.names=FALSE)
  invisible(result)
}
saved<-read.csv(file.path(out,'selected-nonspatial-rows.csv'))
manifest<-read.csv(file.path(out,'selected-fit-provenance.csv'))
groups<-split(subset(saved,group=='all'),interaction(subset(saved,group=='all')$scenario,subset(saved,group=='all')$replicate,drop=TRUE))
diag<-do.call(rbind,lapply(groups,function(z)data.frame(scenario=z$scenario[1],replicate=z$replicate[1],
  flagged=flag_rows(z),max_rhat=max(z$rhat),min_bulk_ess=min(z$ess_bulk),min_tail_ess=min(z$ess_tail))))
write.csv(diag,file.path(out,'selected-diagnostics.csv'),row.names=FALSE)
keep<-merge(saved,diag[c('scenario','replicate','flagged')])
write_screened(subset(keep,!flagged),saved,'selected-screened-sensitivity',out)
historical<-read.csv(file.path(out,'historical-q-replicates.csv'))
paired<-do.call(rbind,lapply(c('qnear','qfar'),function(family) {
  a<-subset(historical,scenario==paste0(family,'_K3'))
  b<-subset(historical,scenario==paste0(family,'_K30'))
  z<-merge(a,b,by=c('replicate','block','group'),suffixes=c('_K3','_K30'))
  delta<-z$covered_K30-z$covered_K3
  se<-sd(delta)/sqrt(nrow(z));critical<-qt(.975,nrow(z)-1L)
  data.frame(family,communities=nrow(z),coverage_delta_K30_minus_K3=mean(delta),mcse=se,
    mc_lower=mean(delta)-critical*se,mc_upper=mean(delta)+critical*se)
}))
write.csv(paired,file.path(out,'historical-q-paired.csv'),row.names=FALSE)

initialfiles<-list.files(file.path(study,'initial'),pattern='-result.rds$',full.names=TRUE)
stopifnot(length(initialfiles)==200L)
initial<-lapply(initialfiles,readRDS)
initialrows<-do.call(rbind,lapply(initial,`[[`,'rows'))
write_interval_tables(initialrows,'continuous-initial',out)
selection<-lapply(seq_along(initial),function(i) {
  x<-initial[[i]];j<-x$job
  flagged<-flag_rows(x$rows)||length(x$warnings)>0L
  path<-initialfiles[i];schedule<-'initial'
  if(flagged) {
    path<-file.path(study,'long',paste0(j$key,'-result.rds'))
    stopifnot(file.exists(path));schedule<-'long'
  }
  y<-readRDS(path)
  stopifnot(identical(x$job,y$job),identical(x$hashes,y$hashes))
  list(rows=y$rows,manifest=data.frame(key=j$key,arm=j$arm,replicate=j$replicate,schedule=schedule,
    initial_flagged=flagged,selected_flagged=flag_rows(y$rows)||length(y$warnings)>0L,
    result=path,result_md5=unname(tools::md5sum(path)),
    fit=sub('-result.rds$','-fit.rds',path),fit_md5=unname(tools::md5sum(sub('-result.rds$','-fit.rds',path))),
    warnings=length(y$warnings),max_rhat=max(y$rows$rhat),min_bulk_ess=min(y$rows$ess_bulk),min_tail_ess=min(y$rows$ess_tail)))
})
rows<-do.call(rbind,lapply(selection,`[[`,'rows'))
selected<-do.call(rbind,lapply(selection,`[[`,'manifest'))
stopifnot(!anyDuplicated(selected$key),all(table(selected$arm)==100L))
write.csv(selected,file.path(out,'continuous-selection.csv'),row.names=FALSE)
write_interval_tables(rows,'continuous-selected',out)
initial_summary<-read.csv(file.path(out,'continuous-initial-summary.csv'))
selected_summary<-read.csv(file.path(out,'continuous-selected-summary.csv'))
substitution<-merge(initial_summary,selected_summary,by=c('scenario','block','group'),suffixes=c('_initial','_selected'))
for(metric in c('coverage','bias','width','rmse'))substitution[[paste0(metric,'_change')]]<-substitution[[paste0(metric,'_selected')]]-substitution[[paste0(metric,'_initial')]]
write.csv(substitution,file.path(out,'continuous-substitution.csv'),row.names=FALSE)
keep<-merge(rows,selected[c('arm','replicate','selected_flagged')],by.x=c('scenario','replicate'),by.y=c('arm','replicate'))
write_screened(subset(keep,!selected_flagged),rows,'continuous-screened-sensitivity',out)

# A conditional Gaussian reference removes true slopes and hidden effects,
# uses the known noise SD and the actual Normal(0,1) intercept prior.
oracle<-list();marginal<-list()
for(r in 1:100) {
  x<-readRDS(file.path(study,'inputs',sprintf('community-%03d.rds',r)))
  for(arm in c('d2','d0')) {
    sim<-x$sim[[arm]];jp<-sim$true_params$jsdmParams_true
    rawx<-as.matrix(sim$data_list$info[paste0('X_psi.EnvCov.',1:2)])
    nuisance<-scale(rawx)%*%jp$B+jp$U%*%jp$L
    residual<-sim$data_list$OTU-nuisance
    variance<-1/(1+nrow(residual)/jp$tau^2)
    center<-colSums(residual)/jp$tau^2*variance
    oracle[[paste(arm,r)]]<-data.frame(scenario=arm,replicate=r,block='B0',element=1:10,truth=jp$B0,
      post_mean=center,lower=center-qnorm(.975)*sqrt(variance),upper=center+qnorm(.975)*sqrt(variance))
    # Integrate iid Gaussian site factors with their known generating covariance.
    sigma<-diag(jp$tau^2)+x$truth$jsdmParams$sigma_h^2*crossprod(jp$L)
    precision<-solve(sigma)
    covariance<-solve(diag(10)+nrow(residual)*precision)
    factor_integrated<-sim$data_list$OTU-scale(rawx)%*%jp$B
    center<-as.vector(covariance%*%precision%*%colSums(factor_integrated))
    marginal[[paste(arm,r)]]<-data.frame(scenario=arm,replicate=r,block='B0',element=1:10,truth=jp$B0,
      post_mean=center,lower=center-qnorm(.975)*sqrt(diag(covariance)),upper=center+qnorm(.975)*sqrt(diag(covariance)))
  }
}
write_interval_tables(do.call(rbind,oracle),'continuous-known-nuisance',out)
write_interval_tables(do.call(rbind,marginal),'continuous-known-covariance',out)
cat('Summarized 200 initial fits and diagnostic-selected longer fits.\n')
