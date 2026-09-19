#!/usr/bin/env Rscript
# Independent validation: no study scoring helpers are sourced.
args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else 'work/jsdm-sample-size-20260919')
out<-file.path(here,'summary');manifest<-read.csv(file.path(out,'fit-manifest.csv'));g<-read.csv(file.path(out,'dataset-groups.csv'));summary<-read.csv(file.path(out,'summary.csv'));pp<-read.csv(file.path(out,'paired-summary.csv'));points<-read.csv(file.path(out,'original-site-probabilities.csv'))
stopifnot(nrow(manifest)==30L,nrow(points)==30000L)
maximum_difference<-0
for(k in seq_len(nrow(manifest))) {
 m<-manifest[k,];folder<-if(m$schedule=='long')'long-fits' else 'fits';file<-file.path(here,folder,paste0(m$key,'-fit.rds'));inputfile<-file.path(here,'inputs',paste0(m$key,'.rds'))
 stopifnot(unname(tools::md5sum(file))==m$fit_md5,unname(tools::md5sum(inputfile))==m$input_md5)
 saved<-readRDS(file);f<-saved$fit;input<-readRDS(inputfile);r<-readRDS(sub('-fit.rds','-result.rds',file,fixed=TRUE))
 stopifnot(identical(saved$job,r$job),identical(unname(tools::md5sum(names(saved$job$source_hashes))),unname(saved$job$source_hashes)),f$infos$model=='binary',f$infos$ps==0,is.null(f$results_output$p_output),is.null(f$results_output$q_output),is.na(r$public_mean_difference))
 X<-scale(as.matrix(input$data$info[,input$covariates]));t<-input$truth;truth<-plogis(X%*%t$B+t$U%*%t$L+matrix(t$B0,input$n,10,byrow=TRUE))
 stopifnot(max(abs(truth-r$truth))<1e-12,max(abs(X-f$X_psi))<1e-12)
 j<-f$results_output$jsdm_output;ni<-dim(j$B0_output)[2];nc<-dim(j$B0_output)[3]
 # Reconstruct all draws independently, one species at a time, using vectorized products.
 meanp<-matrix(0,input$n,10L)
 for(s in 1:10) {
  for(ch in seq_len(nc)) {
   env<-X%*%j$B_output[,s,,ch]+matrix(j$B0_output[s,,ch],input$n,ni,byrow=TRUE)
   hidden<-matrix(0,input$n,ni)
   for(d in seq_len(input$n_factors))hidden<-hidden+sweep(j$U_output[,d,,ch],2,j$L_output[d,s,,ch],'*')
   meanp[,s]<-meanp[,s]+rowMeans(plogis(env+hidden))/nc
  }
 }
 delta<-max(abs(meanp-r$estimate));maximum_difference<-max(maximum_difference,delta);stopifnot(delta<1e-11)
 for(scope in c('original100','allsites'))for(band in c('all','low','middle','high')) {
  ix<-matrix(FALSE,input$n,10);ix[if(scope=='original100')1:100 else seq_len(input$n),]<-TRUE
  if(band=='low')ix<-ix & truth<.2
  if(band=='middle')ix<-ix & truth>=.2 & truth<=.8
  if(band=='high')ix<-ix & truth>.8
  e<-meanp[ix]-truth[ix];row<-g[g$key==m$key & g$scope==scope & g$group==band,]
  stopifnot(nrow(row)==1,row$cells==sum(ix),abs(row$bias-mean(e))<1e-11,abs(row$mae-mean(abs(e)))<1e-11,abs(row$rmse-sqrt(mean(e^2)))<1e-11)
 }
 p<-points[points$key==m$key,];stopifnot(max(abs(p$truth-as.vector(truth[1:100,])))<1e-12,max(abs(p$estimate-as.vector(meanp[1:100,])))<1e-11)
 cat(m$key,'verified\n');flush.console()
}
for(k in seq_len(nrow(summary))) {
 s<-summary[k,];x<-g[g$n==s$n & g$scope==s$scope & g$group==s$group,];stopifnot(nrow(x)==10)
 for(v in c('truth','estimate','bias','mae','rmse'))stopifnot(abs(s[[v]]-mean(x[[v]]))<1e-12,abs(s[[paste0(v,'_se')]]-sd(x[[v]])/sqrt(10))<1e-12)
}
for(k in seq_len(nrow(pp))) {
 s<-pp[k,];a<-g[g$n==s$from & g$scope=='original100' & g$group==s$group,];b<-g[g$n==s$to & g$scope=='original100' & g$group==s$group,];a<-a[order(a$replicate),];b<-b[order(b$replicate),];d<-a$mae-b$mae
 stopifnot(abs(s$reduction_mae-mean(d))<1e-12,abs(s$se-sd(d)/sqrt(10))<1e-12,abs(s$lower-(mean(d)-qt(.975,9)*sd(d)/sqrt(10)))<1e-12,abs(s$upper-(mean(d)+qt(.975,9)*sd(d)/sqrt(10)))<1e-12)
}
cat('PASS: all 30 raw fits, 240 group scores, 30,000 original-site estimates, summary means and paired intervals independently verified. Maximum posterior-mean discrepancy:',format(maximum_difference,digits=8),'\n')
