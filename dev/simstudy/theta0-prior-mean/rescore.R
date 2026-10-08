#!/usr/bin/env Rscript
args<-commandArgs(TRUE);repo<-normalizePath(if(length(args))args[1]else '.')
bs<-if(length(args)>1L)as.numeric(strsplit(args[2],',',fixed=TRUE)[[1]])else 20
stopifnot(all(bs %in% c(20,30)))
folder<-file.path(repo,'dev/simstudy/theta0-prior-mean');out<-file.path(folder,'results')
raw<-file.path(repo,'dev/simstudy/results/theta0-prior-mean-20261008')
source(file.path(folder,'helpers.R'))
p<-read.csv(file.path(repo,'dev/simstudy/interval-calibration/results/selected-fit-provenance.csv'))
selection<-readRDS(file.path(repo,'dev/simstudy/results/pr11-current-20260927/selection.rds'))
rows<-diagnostics<-weights<-references<-manifest<-list();count<-0L
for(i in seq_len(nrow(p))) {
 job<-selection[match(p$key[i],selection$key),];stopifnot(nrow(job)==1L)
 input<-readRDS(p$input[i]);saved<-readRDS(p$fit[i]);fit<-saved$fit;ro<-fit$results_output
 stopifnot(unname(tools::md5sum(p$input[i]))==p$input_md5[i],unname(tools::md5sum(p$fit[i]))==p$fit_md5[i],
   fit$infos$ps==0L,fit$infos$model=='two_stage',identical(as.numeric(saved$job$priors$b_q),20),
   identical(dim(ro$beta_theta_output)[-1L],dim(ro$theta0_output)))
 ni<-dim(ro$theta0_output)[2];nc<-dim(ro$theta0_output)[3];truth<-input$truth$params$theta0
 sd<-source_diagnostics(fit,saved$warnings)
 meta<-data.frame(key=job$key,contamination=if(job$scenario=='qnear_K6')'low'else'high',
   design=job$arm,replicate=job$replicate)
 diagnostics[[i]]<-cbind(meta,sd$elements)
 manifest[[i]]<-cbind(meta,p[i,c('schedule','input','input_md5','fit','fit_md5')],
   source_flag=sd$flag,warnings=sd$warnings,max_rhat=sd$max_rhat,min_ess=sd$min_ess)
 for(m in c(0,1))for(b in bs) {
   count<-count+1L;lw<-joint_log_ratio(ro$beta_theta_output,ro$theta0_output,m,b)
   wi<-weight_info(lw,ni,nc)
   rows[[count]]<-cbind(meta,collection_mean=m,theta_b=b,source_flag=sd$flag,
     weight_flag=wi$weight_flag,score_theta(ro$theta0_output,truth,wi$weights))
   weights[[count]]<-cbind(meta,collection_mean=m,theta_b=b,draws=ni*nc,
     pareto_k=wi$pareto_k,weight_ess=wi$weight_ess,psis_ess=wi$psis_ess,r_eff=wi$r_eff,
     min_chain_weight=min(wi$chain_weight),max_chain_weight=max(wi$chain_weight),weight_flag=wi$weight_flag)
 }
 # The simulator stores samples in site-major order, independently verified later.
 z<-input$sim$true_params$z_true;w<-input$sim$true_params$w_true
 stopifnot(nrow(w)==nrow(z)*input$scenario$M,ncol(w)==ncol(z))
 zz<-z[rep(seq_len(nrow(z)),each=input$scenario$M),,drop=FALSE]
 for(b in bs) {
   n1<-colSums(zz==0 & w==1);n0<-colSums(zz==0 & w==0)
   q<-qbeta(c(.025,.975),rep(1+n1,each=2L),rep(b+n0,each=2L));q<-matrix(q,nrow=2L)
   references[[length(references)+1L]]<-cbind(meta,theta_b=b,species=seq_along(truth),truth=truth,
     absent_collected=n1,absent_uncollected=n0,post_mean=(1+n1)/(1+b+n1+n0),lower=q[1,],upper=q[2,])
 }
 cat(i,'/60',job$key,'complete; source flag',sd$flag,'\n');flush.console()
 rm(saved,fit);gc(FALSE)
}
rr<-do.call(rbind,rows);save_tables(rr,'importance',out)
if(any(!rr$source_flag & !rr$weight_flag))save_tables(subset(rr,!source_flag & !weight_flag),'importance-screened',out)
write.csv(do.call(rbind,diagnostics),file.path(out,'source-element-diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,weights),file.path(out,'importance-diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,manifest),file.path(out,'selected-fit-provenance.csv'),row.names=FALSE)
ref<-do.call(rbind,references);ref$collection_mean<-0;ref$covered<-as.integer(ref$truth>=ref$lower & ref$truth<=ref$upper)
save_tables(ref,'known-state-reference',out)
saveRDS(list(bs=bs,hashes=tools::md5sum(file.path(folder,c('helpers.R','rescore.R'))),session=sessionInfo()),file.path(raw,'rescore-settings.rds'))
cat('All 60 selected sources scored.\n')
