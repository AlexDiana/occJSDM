# Research-only helpers. The production package and its priors are unchanged.
load_frozen <- function(root) {
 .libPaths(c(file.path(root,'library'),.libPaths()))
 Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 suppressPackageStartupMessages(library(occJSDM))
 RcppParallel::setThreadOptions(numThreads=1)
 s <- readRDS(file.path(root,'operational-k6/settings.rds'))
 stopifnot(identical(unname(tools::md5sum(names(s$hashes))),unname(s$hashes)))
 s
}

make_input <- function(a,b,n) {
 stopifnot(n %in% c(100L,300L,1000L),a$scenario$n==100L,b$scenario$n==300L)
 src <- if(n==100L) a else b
 info <- src$sim$data_list$info
 info <- info[!duplicated(info$Site),,drop=FALSE]
 info <- info[order(info$Site),,drop=FALSE]
 covariates <- paste0('X_psi.EnvCov.',seq_len(a$scenario$ncov_psi))
 info <- info[,c('Site',covariates),drop=FALSE];rownames(info)<-NULL
 jp <- src$sim$true_params$jsdmParams_true
 z <- src$sim$true_params$z_true
 extension_seed <- a$seed+3010000L
 if(n==1000L) {
  set.seed(extension_seed)
  xold <- as.matrix(info[,covariates]); nadd <- n-nrow(xold)
  x <- rbind(xold,matrix(rnorm(nadd*ncol(xold),sd=10),nadd,ncol(xold)))
  jp$U <- rbind(jp$U,matrix(rnorm(nadd*ncol(jp$U),sd=a$truth$jsdmParams$sigma_h),nadd,ncol(jp$U)))
  ratio <- apply(x,2,sd)/apply(xold,2,sd)
  jp$B0 <- jp$B0+as.vector(((colMeans(x)-colMeans(xold))/apply(xold,2,sd))%*%jp$B)
  jp$B <- sweep(jp$B,1,ratio,'*')
  for(nm in c('G','C','Bt'))jp[[nm]]<-sweep(jp[[nm]],2,ratio,'*')
  jp$eta <- sweep(scale(x)%*%jp$B+jp$U%*%jp$L,2,jp$B0,'+')
  z <- rbind(z,matrix(rbinom(nadd*ncol(z),1,plogis(jp$eta[301:n,])),nadd,ncol(z)))
  info <- data.frame(Site=seq_len(n));info[,covariates]<-x
 }
 colnames(z) <- colnames(a$sim$data_list$OTU)
 # Retain only relevant truth; do not carry stale spatial or observation arrays.
 truth <- jp[c('B0','B','G','C','Bt','A','U','L','eta')]
 stopifnot(max(abs(sweep(scale(as.matrix(info[,covariates]))%*%truth$B+truth$U%*%truth$L,2,truth$B0,'+')-truth$eta))<1e-12)
 list(n=n,replicate=a$replicate,seed=a$seed,extension_seed=extension_seed,
      data=list(info=info,OTU=z,traits=a$sim$data_list$traits),truth=truth,
      covariates=covariates,n_factors=a$scenario$d,n_lattrait=a$scenario$gt,fit_rng=a$fit_rng)
}

error_metrics <- function(truth,estimate) {
 e <- estimate-truth
 c(cells=length(e),truth=mean(truth),estimate=mean(estimate),bias=mean(e),mae=mean(abs(e)),rmse=sqrt(mean(e^2)))
}
trace_stats <- function(x) {
 c(rhat=posterior::rhat(x),ess_mean=posterior::ess_mean(x),mcse=posterior::mcse_mean(x),chain_gap=diff(range(colMeans(x))))
}

score_fit <- function(fit,input) {
 stopifnot(fit$infos$model=='binary',fit$infos$ps==0,
           is.null(fit$results_output$p_output),is.null(fit$results_output$q_output))
 j <- fit$results_output$jsdm_output
 n <- input$n;S <- ncol(input$data$OTU);ni <- dim(j$B0_output)[2];nc <- dim(j$B0_output)[3]
 X <- fit$X_psi;tp <- input$truth
 true_eta <- sweep(X%*%tp$B+tp$U%*%tp$L,2,tp$B0,'+')
 stopifnot(max(abs(true_eta-tp$eta))<1e-10)
 truth <- plogis(true_eta)
 original <- unlist(lapply(seq_len(S),function(s)(s-1L)*n+seq_len(100L)))
 masks <- list()
 for(scope in c('original100','allsites')) {
  base <- if(scope=='original100')original else seq_len(n*S)
  for(band in c('all','low','middle','high')) {
   p <- truth[base]
   pick <- switch(band,all=rep(TRUE,length(p)),low=p<.2,middle=p>=.2&p<=.8,high=p>.8)
   masks[[paste(scope,band,sep=':')]] <- base[pick]
  }
 }
 traces <- array(0,c(length(masks),ni,nc),dimnames=list(names(masks),NULL,NULL))
 # Store original-site draws only; all-site means and group traces are accumulated.
 original_draws <- array(0,c(100L*S,ni,nc))
 sum_psi <- matrix(0,n,S)
 for(ch in seq_len(nc))for(it in seq_len(ni)) {
  eta <- sweep(X%*%j$B_output[,,it,ch]+j$U_output[,,it,ch]%*%j$L_output[,,it,ch],2,j$B0_output[,it,ch],'+')
  p <- plogis(eta);sum_psi <- sum_psi+p
  original_draws[,it,ch] <- p[original]
  traces[,it,ch] <- vapply(masks,function(ix)mean(p[ix]),numeric(1))
 }
 mean_psi <- sum_psi/(ni*nc)
 stopifnot(all(is.finite(mean_psi)))
 public_error <- NA_real_
 if(!is.null(fit$results_output$psi_output)) {
  public_error <- max(abs(mean_psi-fit$results_output$psi_output))
  stopifnot(is.finite(public_error),public_error<1e-10)
 }
 groups <- do.call(rbind,lapply(seq_along(masks),function(k) {
  labels <- strsplit(names(masks)[k],':',fixed=TRUE)[[1]];ix<-masks[[k]]
  data.frame(scope=labels[1],group=labels[2],as.list(error_metrics(truth[ix],mean_psi[ix])),as.list(trace_stats(traces[k,,])))
 }))
 species <- do.call(rbind,lapply(seq_len(S),function(s)data.frame(species=colnames(input$data$OTU)[s],as.list(error_metrics(truth[1:100,s],mean_psi[1:100,s])))))
 diagnostic_block <- function(arr,metric) {
  dd <- dim(arr);np <- prod(head(dd,-2));flat<-array(arr,c(np,ni,nc))
  do.call(rbind,lapply(seq_len(np),function(k)data.frame(metric=metric,element=k,as.list(trace_stats(flat[k,,])))))
 }
 elements <- rbind(diagnostic_block(j$B0_output,'intercept'),diagnostic_block(j$B_output,'environment_slope'),diagnostic_block(original_draws,'original_probability'))
 list(groups=groups,species=species,elements=elements,traces=traces,truth=truth,estimate=mean_psi,public_mean_difference=public_error)
}
