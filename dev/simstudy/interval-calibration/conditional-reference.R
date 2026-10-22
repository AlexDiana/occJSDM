#!/usr/bin/env Rscript
# A distribution check of the actual Gaussian coefficient update. This is
# separate from full-model calibration and conditions on the noise scale.
args<-commandArgs(trailingOnly=TRUE)
repo<-normalizePath(if(length(args))args[1]else'.')
study<-normalizePath(if(length(args)>1L)args[2]else file.path(repo,'dev/simstudy/results/interval-calibration-20261008'))
.libPaths(c(file.path(study,'library'),.libPaths()))
suppressPackageStartupMessages(library(occJSDM))
draw<-getFromNamespace('sampleB_SoR','occJSDM')
seed<-getFromNamespace('setOccJSDMSeed','occJSDM')
results<-list()
for(noise_sd in c(.5,1,2)) {
  set.seed(20261008);seed(20261008L)
  # A noncentered design produces correlated coefficient draws, detecting
  # a wrong Cholesky orientation as well as incorrect marginal variance.
  X<-cbind(1,seq(-1,1,length.out=31)+.4)
  precision_prior<-diag(c(1,1/2))
  prior_mean<-c(0,.2)
  response<-sin(seq_len(nrow(X)))+.4+.8*X[,2]
  omega<-rep(1/noise_sd^2,nrow(X))
  k<-response*omega
  precision<-crossprod(X,X*omega)+precision_prior
  covariance<-solve(precision)
  center<-as.vector(covariance%*%(crossprod(X,k)+precision_prior%*%prior_mean))
  samples<-replicate(20000,draw(X,precision_prior,prior_mean,k,omega,matrix(numeric(),0,0),matrix(numeric(),0,0),0L))
  samples<-matrix(samples,nrow=ncol(X))
  estimated<-rowMeans(samples);observed_covariance<-cov(t(samples))
  mean_se<-sqrt(diag(covariance)/ncol(samples))
  standardized_mean_error<-(estimated-center)/mean_se
  variance_ratio<-diag(observed_covariance)/diag(covariance)
  stopifnot(max(abs(standardized_mean_error))<5,max(abs(variance_ratio-1))<.05,
    max(abs(observed_covariance-covariance))/max(diag(covariance))<.05)
  results[[as.character(noise_sd)]]<-data.frame(noise_sd,element=1:2,expected_mean=center,observed_mean=estimated,
    mean_error_in_mcse=standardized_mean_error,expected_variance=diag(covariance),observed_variance=diag(observed_covariance),variance_ratio,
    expected_covariance12=covariance[1,2],observed_covariance12=observed_covariance[1,2])
}
write.csv(do.call(rbind,results),file.path(repo,'dev/simstudy/interval-calibration/results/gaussian-update-reference.csv'),row.names=FALSE)
cat('Gaussian joint coefficient mean/covariance checked at three noise scales.\n')
