args <- commandArgs(trailingOnly=TRUE); code <- args[1]; root <- args[2]; package <- args[3]
source(file.path(code,"..","config.R")); for(f in c("basis.R","adapter-contract.R","adapter-occjsdm.R","adapter-hmsc.R","adapter-gllvm.R")) source(file.path(code,f))
set.seed(7); n <- 20L; S <- 8L
tr <- list(schema_version=1L,dataset_id="native-synthetic",site_id=sprintf("s%03d",1:n),species_id=paste0("sp",1:S),
  environment=matrix(rnorm(n),ncol=1,dimnames=list(NULL,"environment")),xy=matrix(runif(2*n),n,2),input_hashes=c(synthetic="fixture"))
tr$y <- matrix(rbinom(n*S,1,.4),n,S,dimnames=list(NULL,tr$species_id))
pred <- tr; pred$y <- NULL
settings <- list(seed=33L,prediction_seed=733L,job_id=1L,attempt_id=1L,source_hashes=c(fixture="synthetic"),
  mcmc=c(nchain=2L,nburn=10L,niter=20L,nthin=1L),gllvm=list(maxit=100L))
if(length(args)>3 && args[4]=="predict") {
  for(arm in c("env","trend","mem","spatial2")) {
    saved <- readRDS(file.path(root,paste0(package,"-",arm,".rds")))
    result <- get(paste0("predict_",tolower(package)))(saved$fit,saved$input,saved$basis)
    stopifnot(max(abs(result$p_full-saved$prediction$p_full))<1e-10)
  }
  cat(package,"fresh process prediction PASS\n"); quit(status=0)
}
for(arm in c("env","trend","mem","spatial2")) {
  basis <- fit_basis(tr,arm)
  f <- get(paste0("fit_",tolower(package)))(tr,basis,settings)
  p <- get(paste0("predict_",tolower(package)))(f,pred,basis)
  stopifnot(identical(dim(p$p_full),c(n,S)),identical(f$species_id,tr$species_id))
  if(arm=="env") stopifnot(max(abs(p$p_full-p$p_space_zero))==0)
  if(package=="Hmsc") {
    native_fun <- getS3method("predict","Hmsc",envir=asNamespace("Hmsc"))
    env <- new.env(parent=environment(native_fun)); env$detectCores <- function(...) 1L; environment(native_fun)<-env
    native <- native_fun(f$native,expected=TRUE,nParallel=1)
    stopifnot(max(abs(p$p_full-Reduce(`+`,native)/length(native)))<1e-10)
    if(arm=="spatial2") {
      post <- f$native$postList[[1]][[1]]; level <- f$native$ranLevels$site
      # Force a positive saved grid value to test native conditional variation.
      q <- which(level$alphapw[,1]>0)[3]; a <- level$alphapw[q,1]
      far <- subset_prediction(pred,1:2); far$site_id <- c("newA","newB"); far$xy <- rbind(c(2,2),c(3,3))
      newlevel <- hmsc_prediction_level(f,far)
      Ktt <- exp(-xy_distances(tr$xy,tr$xy)/a); Knt <- exp(-xy_distances(far$xy,tr$xy)/a)
      Knn <- exp(-xy_distances(far$xy,far$xy)/a)
      truth <- conditional_spatial_moments(post$Eta[[1]][,1],matrix(0,n,n),Ktt,Knt,Knn)
      set.seed(444L)
      native_mean <- Hmsc:::predictLatentFactor(newlevel$ids,tr$site_id,list(post$Eta[[1]]),list(c(q,q)),newlevel$level,predictMean=TRUE)[[1]][,1]
      stopifnot(max(abs(native_mean-truth$mean))<1e-10)
      draws <- Hmsc:::predictLatentFactor(newlevel$ids,tr$site_id,rep(list(post$Eta[[1]]),2000),rep(list(c(q,q)),2000),newlevel$level)
      set.seed(444L)
      cached <- hmsc_conditional_draws(newlevel$ids,tr$site_id,rep(list(post$Eta[[1]]),2000),rep(list(c(q,q)),2000),newlevel$level)
      stopifnot(max(abs(unlist(cached)-unlist(draws)))<1e-9)
      vals <- vapply(draws,function(x)x[,1],numeric(2))
      stopifnot(max(abs(rowMeans(vals)-truth$mean)) < 5*sqrt(max(diag(truth$covariance))/2000)+1e-6,
        max(abs(cov(t(vals))-truth$covariance))<.12)
      np <- predict_hmsc(f,far,basis); stopifnot(all(is.finite(np$p_full)),all(np$intervals$field_upper>np$intervals$field_lower))
    }
  } else {
    fit <- f$native$fit
    X <- as.matrix(occ_design(tr,basis)); raw <- sweep(X %*% t(fit$params$Xcoef),2,fit$params$beta0,"+")
    if(arm=="spatial2") raw <- raw+sweep(fit$lvs,2,fit$params$sigma.lv,"*") %*% t(fit$params$theta)
    native_raw <- predict(fit,type="link")
    stopifnot(max(abs(raw-native_raw))<1e-8)
    if(arm=="spatial2") {
      moments <- gllvm_field_moments(f,pred)
      stopifnot(max(abs(moments$mean-(sweep(fit$lvs,2,fit$params$sigma.lv,"*") %*% t(fit$params$theta))))<1e-8)
      far <- subset_prediction(pred,1L); far$xy[,] <- 20; far$site_id <- "remote"
      fm <- gllvm_field_moments(f,far)
      stopifnot(all(fm$variance>0),max(abs(fm$mean))<1e-6)
      np <- predict_gllvm(f,far,basis); stopifnot(identical(dim(np$p_full),c(1L,S)))
    } else stopifnot(max(abs(p$p_full-plogis(native_raw)))<1e-10)
  }
  path <- file.path(root,paste0(package,"-",arm,".rds")); saveRDS(list(fit=f,basis=basis,prediction=p,input=pred),path)
  saved <- readRDS(path); p2 <- get(paste0("predict_",tolower(package)))(saved$fit,pred,saved$basis)
  stopifnot(max(abs(p2$p_full-p$p_full))<1e-10)
  cat(package,arm,"PASS\n")
}
