#!/usr/bin/env Rscript
args<-commandArgs(TRUE);repo<-normalizePath(if(length(args))args[1]else '.')
folder<-file.path(repo,'dev/simstudy/theta0-prior-mean');out<-file.path(folder,'results')
raw<-file.path(repo,'dev/simstudy/results/theta0-prior-mean-20261008')
source(file.path(folder,'helpers.R'))
checks<-0L;chain_weights<-list();check<-function(x){stopifnot(isTRUE(x));checks<<-checks+1L}
rr<-read.csv(file.path(out,'importance-rows.csv'));mf<-read.csv(file.path(out,'selected-fit-provenance.csv'))
refs<-read.csv(file.path(out,'known-state-reference-rows.csv'));bs<-sort(unique(rr$theta_b))
check(nrow(mf)==60L);check(nrow(rr)==60L*2L*length(bs)*10L);check(!anyDuplicated(rr[c('key','collection_mean','theta_b','species')]))
for(i in seq_len(nrow(mf))) {
 saved<-readRDS(mf$fit[i]);x<-readRDS(mf$input[i]);ro<-saved$fit$results_output
 check(unname(tools::md5sum(mf$input[i]))==mf$input_md5[i]);check(unname(tools::md5sum(mf$fit[i]))==mf$fit_md5[i])
 info<-x$sim$data_list$info
 samples<-info[!duplicated(info$Sample),c('Site','Sample')]
 # Sample IDs are sequential and independently recover the simulator's field-state order.
 check(identical(as.integer(samples$Sample),seq_len(nrow(x$sim$true_params$w_true))))
 zz<-x$sim$true_params$z_true[as.integer(samples$Site),,drop=FALSE];ww<-x$sim$true_params$w_true
 t<-matrix(ro$theta0_output,nrow=10L);btheta<-matrix(ro$beta_theta_output[2,,,,drop=FALSE],nrow=10L)
 for(m in c(0,1))for(b in bs) {
  # Independent density evaluation, rather than the simplified scoring formula.
  lw<-colSums(dnorm(btheta,m,sqrt(2),log=TRUE)-dnorm(btheta,0,sqrt(2),log=TRUE))+
    colSums(dbeta(t,1,b,log=TRUE)-dbeta(t,1,20,log=TRUE))
  actual<-joint_log_ratio(ro$beta_theta_output,ro$theta0_output,m,b)
  check(max(abs(actual-lw))<1e-10)
  wt<-exp(lw-max(lw));wt<-wt/sum(wt)
  ni<-dim(ro$theta0_output)[2L];nc<-dim(ro$theta0_output)[3L]
  chain_weights[[length(chain_weights)+1L]]<-data.frame(key=mf$key[i],collection_mean=m,theta_b=b,chain=seq_len(nc),weight_fraction=colSums(matrix(wt,ni,nc)))
  rows<-subset(rr,key==mf$key[i] & collection_mean==m & theta_b==b)
  for(s in 1:10) {
   r<-rows[rows$species==s,];check(abs(r$truth-x$truth$params$theta0[s])<1e-12)
   check(abs(sum(wt*t[s,])-r$post_mean)<1e-10)
   # An inverse-CDF quantile brackets the target mass with and without its atom.
   for(j in 1:2) {
     v<-if(j==1)r$lower else r$upper;p<-c(.025,.975)[j]
     check(min(abs(t[s,]-v))<1e-12)
     # CSV endpoints round-trip to about 15 digits; allow only value rounding.
     check(sum(wt[t[s,]<v-1e-12])<p+1e-10);check(sum(wt[t[s,]<=v+1e-12])>=p-1e-10)
   }
   check(r$covered==as.integer(r$lower<=r$truth & r$truth<=r$upper))
  }
 }
 for(b in bs) {
  r<-subset(refs,key==mf$key[i] & theta_b==b);r<-r[order(r$species),]
  n1<-colSums(zz==0 & ww==1);n0<-colSums(zz==0 & ww==0)
  check(identical(as.numeric(r$absent_collected),as.numeric(n1)))
  check(identical(as.numeric(r$absent_uncollected),as.numeric(n0)))
  check(max(abs(r$post_mean-(1+n1)/(1+b+n1+n0)))<1e-10)
  check(max(abs(pbeta(r$lower,1+n1,b+n0)-.025))<1e-10)
  check(max(abs(pbeta(r$upper,1+n1,b+n0)-.975))<1e-10)
 }
}
compare<-file.path(out,'factorial-full-fit-provenance.csv')
if(file.exists(compare)) {
 fm<-read.csv(compare);check(nrow(fm)==6L*length(bs))
 selected<-file.path(out,'factorial-full-fit-selected-provenance.csv')
 if(file.exists(selected)) {
   fm<-rbind(fm,read.csv(selected));fm<-fm[!duplicated(fm$file),]
   initial_flags<-read.csv(compare)
   check(nrow(fm)==6L*length(bs)+sum(initial_flags$source_flag))
 }
 for(i in seq_len(nrow(fm))) {
  x<-readRDS(fm$file[i]);check(unname(tools::md5sum(fm$file[i]))==fm$md5[i])
  check(x$job$key==fm$key[i] && x$job$theta_b==fm$theta_b[i] && x$job$collection_mean==fm$collection_mean[i])
  check(unname(tools::md5sum(x$job$input))==x$job$input_md5)
  long<-basename(dirname(fm$file[i]))=='long-fits'
  ni<-if(long)12000L else 6000L;burn<-if(long)6000 else 3000
  check(identical(x$mcmc,list(nchain=4,nburn=burn,niter=as.numeric(ni),nthin=1)))
  check(identical(dim(x$fit$results_output$theta0_output),c(10L,ni,4L)))
  check(all(unname(tools::md5sum(names(x$hashes)))==unname(x$hashes)))
  runner<-if(!long)'run.R'else if(x$job$theta_b==30)'run-historical-long.R'else 'run-long.R'
  for(n in c('helpers.R',runner))check(unname(tools::md5sum(file.path(raw,paste0('executed-',n))))==
    unname(x$hashes[file.path(folder,n)]))
 }
}
# Reconstruct paired community summaries from element rows without using summaries().
ps<-read.csv(file.path(out,'paired-summary.csv'))
for(i in seq_len(nrow(ps))) {
 r<-ps[i,];z<-subset(rr,contamination==r$contamination & design==r$design)
 if(grepl('collection mean 1 to 0',r$contrast,fixed=TRUE)) {
   b<-as.numeric(sub('.*theta b=','',r$contrast));from_mean<-1;to_mean<-0;from_b<-to_b<-b
 }else {
   check(grepl('theta b 30 to 20',r$contrast,fixed=TRUE))
   m<-as.numeric(sub('.*collection mean=','',r$contrast));from_mean<-to_mean<-m;from_b<-30;to_b<-20
 }
 metric<-as.character(r$metric);values<-numeric()
 for(k in sort(unique(z$replicate))) {
   from<-subset(z,replicate==k & collection_mean==from_mean & theta_b==from_b)
   to<-subset(z,replicate==k & collection_mean==to_mean & theta_b==to_b)
   check(nrow(from)==10L && nrow(to)==10L)
   if(r$scope=='common screened inputs' && any(c(from$source_flag,from$weight_flag,to$source_flag,to$weight_flag)))next
   scalar<-function(x) switch(metric,covered=mean(x$covered),bias=mean(x$post_mean-x$truth),
     abs_error=mean(abs(x$post_mean-x$truth)),width=mean(x$upper-x$lower),
     lower_miss=mean(x$truth<x$lower),upper_miss=mean(x$truth>x$upper))
   v<-if(metric=='width_percent')100*(mean(to$upper-to$lower)/mean(from$upper-from$lower)-1)else scalar(to)-scalar(from)
   values<-c(values,v)
 }
 check(length(values)==r$communities);check(abs(mean(values)-r$change)<1e-10)
 check(abs(sd(values)/sqrt(length(values))-r$mcse)<1e-10)
 if(length(values)>1) {
   crit<-qt(.975,length(values)-1);check(abs(r$lower-(mean(values)-crit*r$mcse))<1e-10)
   check(abs(r$upper-(mean(values)+crit*r$mcse))<1e-10)
 }
}
# Direct four-arm reconstruction of the interaction, independent of contrast merges.
if(length(bs)==2L) {
  ii<-read.csv(file.path(out,'factorial-interaction-communities.csv'))
  for(i in seq_len(nrow(ii))) {
    r<-ii[i,];z<-subset(rr,contamination==r$contamination & design==r$design & replicate==r$replicate)
    scalar<-function(x)switch(as.character(r$metric),covered=mean(x$covered),bias=mean(x$post_mean-x$truth),
      abs_error=mean(abs(x$post_mean-x$truth)),width=mean(x$upper-x$lower))
    a<-subset(z,collection_mean==0 & theta_b==20);b<-subset(z,collection_mean==0 & theta_b==30)
    c<-subset(z,collection_mean==1 & theta_b==20);d<-subset(z,collection_mean==1 & theta_b==30)
    check(abs((scalar(a)-scalar(b))-(scalar(c)-scalar(d))-r$interaction)<1e-10)
    check(r$screened==!any(c(z$source_flag,z$weight_flag)))
  }
  ss<-read.csv(file.path(out,'factorial-interaction-summary.csv'))
  for(i in seq_len(nrow(ss))) {
    r<-ss[i,];z<-subset(ii,contamination==r$contamination & design==r$design & metric==r$metric)
    if(r$scope=='common screened inputs')z<-subset(z,screened)
    check(nrow(z)==r$communities);check(abs(mean(z$interaction)-r$interaction)<1e-10)
    check(abs(sd(z$interaction)/sqrt(nrow(z))-r$mcse)<1e-10)
    if(nrow(z)>1) {
      crit<-qt(.975,nrow(z)-1);check(abs(r$lower-(r$interaction-crit*r$mcse))<1e-10)
      check(abs(r$upper-(r$interaction+crit*r$mcse))<1e-10)
    }
  }
}
write.csv(do.call(rbind,chain_weights),file.path(out,'chain-weight-fractions.csv'),row.names=FALSE)
history<-read.csv(file.path(out,'history-prior-blocks.csv'))
check(nrow(history)==4L)
check(all(unname(tools::md5sum(file.path(repo,history$source)))==history$md5))
rescore<-readRDS(file.path(raw,'rescore-settings.rds'))
check(all(unname(tools::md5sum(names(rescore$hashes)))==unname(rescore$hashes)))
check(unname(tools::md5sum(file.path(raw,'executed-rescore.R')))==unname(rescore$hashes[file.path(folder,'rescore.R')]))
result<-paste('PASS:',checks,'independent density, weighted-CDF, truth, known-state reference, paired-community and provenance checks.')
writeLines(result,file.path(out,'verification.txt'));cat(result,'\n')
