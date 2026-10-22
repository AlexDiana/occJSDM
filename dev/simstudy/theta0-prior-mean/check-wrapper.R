#!/usr/bin/env Rscript
repo<-normalizePath(if(length(commandArgs(TRUE)))commandArgs(TRUE)[1]else '.')
base<-file.path(repo,'dev/simstudy/results/interval-calibration-20261008')
.libPaths(c(file.path(base,'library'),.libPaths()))
suppressPackageStartupMessages(library(occJSDM));RcppParallel::setThreadOptions(numThreads=1)
source(file.path(repo,'dev/simstudy/theta0-prior-mean/helpers.R'))
f0<-mean_fitter(0);f1<-mean_fitter(1)
stopifnot(identical(body(f0),body(occJSDM::runOccJSDM)),identical(formals(f0),formals(f1)))
a<-read.csv(file.path(repo,'dev/simstudy/interval-calibration/results/selected-fit-provenance.csv'))
x<-readRDS(a$input[1]);m<-list(nchain=2,nburn=30,niter=60,nthin=1)
fit<-function(f){assign('.Random.seed',x$fit_rng,envir=.GlobalEnv);suppressMessages(suppressWarnings(fit_input(x,f,m)))}
f<-fit(occJSDM::runOccJSDM);g<-fit(f0);h<-fit(f1)
stopifnot(identical(f$results_output,g$results_output),identical(f$X_theta,g$X_theta),
  !identical(f$results_output$beta_theta_output,h$results_output$beta_theta_output),f$infos$ps==0L)
# Density-ratio identities use arbitrary joint draws independent of a fitted posterior.
beta<-array(seq(-2,2,length.out=2*3*5*2),c(2,3,5,2));theta<-array(seq(.01,.7,length.out=3*5*2),c(3,5,2))
for(mean in c(0,1))for(b in c(20,30)) {
 actual<-joint_log_ratio(beta,theta,mean,b)
 independent<-colSums(dnorm(draw_matrix(beta[2,,,,drop=FALSE]),mean,sqrt(2),log=TRUE)-
  dnorm(draw_matrix(beta[2,,,,drop=FALSE]),0,sqrt(2),log=TRUE))+
  colSums(dbeta(draw_matrix(theta),1,b,log=TRUE)-dbeta(draw_matrix(theta),1,20,log=TRUE))
 stopifnot(max(abs(actual-independent))<1e-12)
}
stopifnot(identical(weighted_quantile(c(5,1,9),c(.2,.7,.1),c(.025,.5,.975)),c(1,1,9)))
cat('PASS: mean-zero clone equals exported seeded fit; mean-one affects collection output; joint prior-density identities and weighted quantiles verified.\n')
