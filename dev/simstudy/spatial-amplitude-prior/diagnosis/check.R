args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
repo <- normalizePath(args[1]); archive <- args[2]
dir.create(archive,recursive=TRUE,showWarnings=FALSE)
folder <- file.path(repo,'dev/simstudy/spatial-amplitude-prior/diagnosis')
Rcpp::sourceCpp(file.path(folder,'ellipse.cpp'),cacheDir=file.path(archive,'cpp-cache'))
source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/robust.R'))

# Independent likelihood check, including extreme linear predictors and counts.
for(n in c(1L,10L))for(b in c(FALSE,TRUE)) {
  set.seed(2911+n+b);f<-rnorm(9+b);offset<-seq(-100,100,length.out=9)
  eta<-f[1:9]+offset+if(b)f[10] else 0
  y<-rbinom(9,n,plogis(eta))
  observed<-ellipse_loglik_cpp(f,y,n,offset,b)
  expected<-sum(dbinom(y,n,plogis(eta),log=TRUE)-lchoose(n,y))
  stopifnot(abs(observed-expected)<1e-10)
}
eta<-c(-1000,-100,-10,0,10,100,1000);y<-c(10,5,0,5,10,5,0)
expected<-sum(y*plogis(eta,log.p=TRUE)+(10-y)*plogis(eta,lower.tail=FALSE,log.p=TRUE))
stopifnot(abs(ellipse_loglik_cpp(rep(0,7),y,10L,eta,FALSE)-expected)<1e-10,
  inherits(try(ellipse_loglik_cpp(0,.5,1L,0,FALSE),silent=TRUE),'try-error'))

quadrature <- function(y,n,offset,variance) {
  prior_variance<-variance
  logdensity<-function(x)dnorm(x,sd=sqrt(prior_variance),log=TRUE)+
    y*(offset+x)-n*(pmax(offset+x,0)+log1p(exp(-abs(offset+x))))
  peak<-optimize(function(x)-logdensity(x),c(-30,30))$minimum
  density<-function(x)exp(logdensity(x)-logdensity(peak))
  norm<-integrate(density,-Inf,Inf,rel.tol=1e-11)$value
  moment<-function(k)integrate(function(x)x^k*density(x),-Inf,Inf,rel.tol=1e-10)$value/norm
  mean<-moment(1);posterior_variance<-moment(2)-mean^2
  cdf<-function(q)if(q<=peak)integrate(density,-Inf,q,rel.tol=1e-10)$value/norm else
    1-integrate(density,q,Inf,rel.tol=1e-10)$value/norm
  quantile<-function(p)uniroot(function(q)cdf(q)-p,
    c(-30,30),tol=1e-9)$root
  list(mean=mean,variance=posterior_variance,quantiles=vapply(c(.025,.5,.975),quantile,numeric(1)))
}

rows<-list()
for(case in 1:3) {
  joint<-case==3L; y<-if(case==1L)0 else 3; n<-if(case==1L)1L else 8L
  offset<-if(case==1L)-5 else .4; variance<-if(joint)2 else 1
  expected<-quadrature(y,n,offset,variance)
  initial<-matrix(c(-2,0,1,3),if(joint)2 else 1,4,byrow=TRUE)
  set.seed(29120+case)
  result<-ellipse_draws_cpp(diag(nrow(initial)),y,n,offset,joint,initial,1000L,20000L)
  x<-if(joint)apply(result$draws,c(2,3),sum) else result$draws[1,,]
  diagnostics<-robust_trace_diagnostics(x)
  observed<-mean(x);mcse<-as.numeric(posterior::mcse_mean(x))
  q<-quantile(x,c(.025,.5,.975),names=FALSE)
  stopifnot(abs(observed-expected$mean)<6*mcse+.002,diagnostics['rhat']<1.01,
    diagnostics['ess_bulk']>1000,abs(var(as.vector(x))-expected$variance)<.04,
    all(abs(q-expected$quantiles)<.08))
  if(joint)for(j in 1:2)stopifnot(abs(mean(result$draws[j,,])-expected$mean/2)<.025,
    abs(var(as.vector(result$draws[j,,]))-(.5+expected$variance/4))<.035)
  rows[[case]]<-data.frame(case=case,mean=observed,expected_mean=expected$mean,mcse=mcse,
    variance=var(as.vector(x)),expected_variance=expected$variance,
    max_quantile_error=max(abs(q-expected$quantiles)),rhat=diagnostics['rhat'],bulk_ess=diagnostics['ess_bulk'])
}

# Independent Gaussian contrast field - intercept must retain N(0,2).
contrast<-result$draws[1,,]-result$draws[2,,]
stopifnot(abs(mean(contrast))<.035,abs(var(as.vector(contrast))-2)<.07)
# Offset -1000 makes the zero-success likelihood numerically constant. This
# independently checks correlated prior draws, beyond the joint contrast.
covariance<-matrix(c(1,.6,.6,2),2,2)
set.seed(7612)
prior<-ellipse_draws_cpp(t(chol(covariance)),c(0,0),1L,c(-1000,-1000),FALSE,
  matrix(0,2,4),1000L,10000L)$draws
stopifnot(max(abs(rowMeans(matrix(prior,2))))<.04,
  max(abs(cov(t(matrix(prior,2)))-covariance))<.07)
set.seed(781);a<-ellipse_draws_cpp(diag(2),3,8L,.4,TRUE,matrix(0,2,4),10L,20L)
set.seed(781);b<-ellipse_draws_cpp(diag(2),3,8L,.4,TRUE,matrix(0,2,4),10L,20L)
stopifnot(identical(a,b),identical(dim(a$draws),c(2L,20L,4L)))
write.csv(do.call(rbind,rows),file.path(archive,'quadrature-checks.csv'),row.names=FALSE)
files<-c(file.path(folder,c('ellipse.cpp','check.R')),file.path(repo,'dev/simstudy/spatial-amplitude-prior/robust.R'))
write.csv(data.frame(file=files,md5=unname(tools::md5sum(files))),
  file.path(archive,'validation-source-md5.csv'),row.names=FALSE)
cat('Independent likelihood, quadrature, Gaussian contrast and reproducibility checks passed.\n')
