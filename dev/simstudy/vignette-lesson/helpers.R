# Helpers for the teaching experiment, not part of the package API.
effective_detection_rate <- function(event_rate, mu, sigma, threshold=1) {
  stopifnot(threshold >= 1, threshold == as.integer(threshold))
  event_rate * stats::pnorm(log(threshold + .5), mu, sigma, lower.tail=FALSE)
}

make_lesson <- function() {
  n <- 100L; S <- 10L; M <- 2L; P <- 2L; K <- 6L
  settings <- list(n=n, S=S, g=2L, M=rep(M,n), P=P,
                   K=rep(K,n*M*P), ncov_psi=2L, ncov_theta=1L)
  p_base <- seq(.35,.85,length.out=S)
  params <- list(p=rbind(p_base,pmin(.95,p_base+.10)),
                 q=rbind(seq(.025,.065,length.out=S),
                         seq(.065,.025,length.out=S)),
                 theta0=seq(.02,.08,length.out=S),
                 theta_baseline=seq(.35,.75,length.out=S),
                 mu1=5, sigma1=1, mu0=1.5, sigma0=1)
  jsdm <- list(gt=1L,d=2L,ds=0L,sigma_b=.5,sigma_bs=.5,
               sigma_ts=.5,sigma_h=1,sigma_s=.5,l_s=.3,
               tau=rep(1,S),useSpatField=FALSE)
  seed <- 20260919L
  set.seed(seed)
  sim <- occJSDM::simulateOccJSDMData(settings,params,jsdm,model="two_stage")
  sp <- colnames(sim$data_list$OTU)
  dimnames(sim$true_params$z_true) <- list(as.character(seq_len(n)),sp)
  dimnames(sim$true_params$w_true) <- list(as.character(seq_len(n*M)),sp)
  dimnames(sim$true_params$jsdmParams_true$eta) <- list(as.character(seq_len(n)),sp)
  out <- list(sim=sim,settings=settings,params=params,jsdm=jsdm,
              seed=seed,threshold=1L)
  obs <- lesson_observations(out)
  out$cases <- select_lesson_cases(obs)
  out
}

lesson_observations <- function(bundle) {
  dat <- bundle$sim$data_list
  info <- dat$info
  z <- bundle$sim$true_params$z_true
  w <- bundle$sim$true_params$w_true
  zi <- match(as.character(info$Site),rownames(z))
  wi <- match(as.character(info$Sample),rownames(w))
  sz <- match(colnames(dat$OTU),colnames(z))
  sw <- match(colnames(dat$OTU),colnames(w))
  if (anyNA(c(zi,wi,sz,sw))) stop("Unmatched truth identity")
  pcr <- ave(seq_len(nrow(info)), interaction(info$Site,info$Sample,info$Primer),
             FUN=seq_along)
  do.call(rbind,lapply(seq_len(ncol(dat$OTU)), function(s) {
    reads <- dat$OTU[,s]; positive <- as.integer(reads >= bundle$threshold)
    zs <- z[zi,sz[s]]; ws <- w[wi,sw[s]]
    source <- ifelse(is.na(positive),"Missing",
                     ifelse(positive==0,"No detection",
                            ifelse(ws==0,"Laboratory false positive",
                                   ifelse(zs==0,"Field-stage false positive","True detection"))))
    data.frame(row=seq_len(nrow(info)),species=colnames(dat$OTU)[s],
               Site=info$Site,Sample=info$Sample,Primer=info$Primer,
               PCR=as.integer(pcr),reads=reads,positive=positive,z=zs,w=ws,source=source)
  }))
}

lesson_samples <- function(obs) {
  keys <- interaction(obs$species,obs$Site,obs$Sample,drop=TRUE)
  out <- do.call(rbind,lapply(split(obs,keys),function(x) {
    data.frame(species=x$species[1],Site=x$Site[1],Sample=x$Sample[1],
               z=x$z[1],w=x$w[1],positives=sum(x$positive,na.rm=TRUE),
               observed=sum(!is.na(x$positive)))
  }))
  rownames(out) <- NULL
  out <- out[order(out$species,out$Site,out$Sample),]
  out$positive_true_samples <- ave(as.integer(out$w==1 & out$positives>=3),
                                  interaction(out$species,out$Site),FUN=sum)
  out
}

select_lesson_cases <- function(obs) {
  x <- lesson_samples(obs)
  eligible <- list(
    "Weak true detection"=x$z==1 & x$w==1 & x$positives %in% 1:2,
    "Laboratory false positive"=x$w==0 & x$positives>0,
    "Strong true detection"=x$z==1 & x$w==1 & x$positives>=6 & x$positive_true_samples>=2,
    "Field-stage false positive"=x$z==0 & x$w==1 & x$positives>0)
  do.call(rbind,lapply(names(eligible),function(label) {
    ix <- which(eligible[[label]])
    if(!length(ix)) stop("No teaching case found for: ",label)
    cbind(case=label,eligible=length(ix),x[ix[1],,drop=FALSE],row.names=NULL)
  }))
}

lesson_binary_data <- function(bundle) {
  dat <- bundle$sim$data_list
  info <- dat$info[!duplicated(dat$info$Site),,drop=FALSE]
  info <- info[order(info$Site),c("Site",grep("^X_psi",names(info),value=TRUE)),drop=FALSE]
  list(info=info,OTU=bundle$sim$true_params$z_true[as.character(info$Site),,drop=FALSE],
       traits=dat$traits)
}

lesson_source_hashes <- function(root=".") {
  files <- c("DESCRIPTION","NAMESPACE",list.files("R",full.names=TRUE),
             list.files("src",pattern="\\.(cpp|h|hpp)$",full.names=TRUE))
  tools::md5sum(files)
}

lesson_fit <- function(bundle,arm,mcmc=list(nchain=4L,nburn=3000L,niter=6000L,nthin=1L)) {
  stopifnot(arm %in% c("perfect","default","alternative"))
  dat <- if(arm=="perfect") lesson_binary_data(bundle) else bundle$sim$data_list
  seed <- c(perfect=20260920L,default=20260921L,alternative=20260922L)[[arm]]
  priors <- if(arm=="alternative") list(a_q=1,b_q=4,a_theta0=1,b_theta0=4) else list()
  args <- list(data=dat,listParams=list(n_factors=2L,n_lattrait=1L),
               threshold=1,occCovariates=grep("^X_psi",names(dat$info),value=TRUE),
               collCovariates=if(arm=="perfect") NULL else "X_theta",
               spatCovariates=NULL,MCMCparams=mcmc,listPriors=priors,
               summarisedLatentPresences=TRUE)
  warnings <- character()
  set.seed(seed)
  started <- Sys.time()
  fit <- withCallingHandlers(do.call(occJSDM::runOccJSDM,args),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w))
  })
  list(arm=arm,fit=fit,seed=seed,mcmc=mcmc,priors=priors,warnings=warnings,
       started=started,finished=Sys.time(),session=sessionInfo(),
       source_hashes=lesson_source_hashes())
}
