# Standalone data generator: no occJSDM simulation or fitting internals.
# Uses the existing simstudy_seed(label, replicate) scheme.
spatial_seed <- function(label, replicate) {
  as.integer(sum(as.integer(charToRaw(label))) * 1000L + replicate)
}

make_spatial_input <- function(grid_index, replicate) {
  stopifnot(length(grid_index)==1L, grid_index %in% 1:10,
            length(replicate)==1L, replicate>=1L, replicate==as.integer(replicate))
  data_seed <- spatial_seed(paste0('spatial-targeted-data-range',grid_index),replicate)
  fit_seed <- spatial_seed(paste0('spatial-targeted-fit-range',grid_index),replicate)
  set.seed(data_seed)
  n <- 100L; S <- 8L; M <- 2L; P <- 2L; K <- 6L
  coordinates <- matrix(runif(n*2L),n,2L)
  Xs <- unname(scale(coordinates))
  environment <- rnorm(n); X <- as.numeric(scale(environment))
  range <- seq(.01,.3,length.out=10L)[grid_index]
  covariance <- exp(-as.matrix(dist(Xs))^2/(2*range^2))+diag(1e-10,n)
  innovations <- matrix(rnorm(n*S),n,S)
  field <- t(chol(covariance)) %*% innovations
  B <- rep(c(.4,-.4),4L)
  prevalence <- c(.01,.01,.05,.05,.25,.25,.75,.75)
  offset <- X %o% B+field
  B0 <- vapply(seq_len(S),function(s) uniroot(function(b)
    mean(plogis(b+offset[,s]))-prevalence[s],c(-30,30),tol=1e-12)$root,numeric(1))
  psi <- plogis(sweep(offset,2,B0,'+'))
  species <- paste0('species',seq_len(S))
  z <- matrix(as.integer(runif(n*S)<psi),n,S,dimnames=list(NULL,species))
  site <- data.frame(Site=seq_len(n),environment=environment,
    longitude=coordinates[,1],latitude=coordinates[,2])
  sample_site <- rep(seq_len(n),each=M)
  collection <- rnorm(n*M)
  Xt <- cbind(1,as.numeric(scale(collection)))
  beta_theta <- rbind(qlogis(runif(S,.2,.5)),rep(c(1,-1),4L))
  theta <- plogis(Xt %*% beta_theta)
  theta0 <- runif(S,.02,.1)
  p <- matrix(runif(P*S,.3,.6),P,S)
  q_uniform <- matrix(runif(P*S),P,S)
  q <- list(low=.01+.04*q_uniform,high=.15+.15*q_uniform)
  wprob <- ifelse(z[sample_site,]==1L,theta,
    matrix(theta0,n*M,S,byrow=TRUE))
  w <- matrix(as.integer(runif(n*M*S)<wprob),n*M,S)
  sample <- rep(seq_len(n*M),each=P*K)
  primer <- rep(rep(seq_len(P),each=K),times=n*M)
  row_site <- sample_site[sample]
  info <- data.frame(Site=row_site,Sample=sample,Primer=primer,
    environment=environment[row_site],longitude=coordinates[row_site,1],
    latitude=coordinates[row_site,2],collection=collection[sample])
  u <- matrix(runif(length(sample)*S),length(sample),S)
  observations <- lapply(q,function(qi) {
    probability <- ifelse(w[sample,]==1L,p[primer,],qi[primer,])
    y <- matrix(as.integer(u<probability),length(sample),S,dimnames=list(NULL,species))
    list(info=info,OTU=y,traits=NULL)
  })
  set.seed(fit_seed)
  list(data=c(list(binary=list(info=site,OTU=z,traits=NULL)),observations),
    truth=list(psi=psi,z=z,w=w,field=field,field_sd=1,range=range,
      grid_index=grid_index,X=X,Xs=Xs,B0=B0,B=B,Xt=Xt,
      beta_theta=beta_theta,theta=theta,theta0=theta0,p=p,q=q,
      p_positive=p,q_positive=q,target_prevalence=prevalence,
      innovations=innovations),
    settings=list(n=n,S=S,M=M,P=P,K=K,replicate=replicate,
      data_seed=data_seed,fit_seed=fit_seed,
      observation_model='Direct Bernoulli detections; p/q are positive-observation probabilities'),
    fit_rng=.Random.seed)
}
