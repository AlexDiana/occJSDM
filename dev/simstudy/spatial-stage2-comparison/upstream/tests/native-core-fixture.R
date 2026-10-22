args <- commandArgs(trailingOnly=TRUE); code <- args[1]; root <- args[2]; package <- args[3]; mode <- args[4]
source(file.path(code,"..","config.R")); for(f in c("basis.R","adapter-contract.R","adapter-occjsdm.R")) source(file.path(code,f))
tr <- readRDS(file.path(root,"training.rds")); pred <- tr; pred$y <- NULL
cfg <- study_config(); arms <- if(package=="occJSDM") c("env","trend","mem","gp") else c("env","trend","mem","dnn")
for(arm in arms) {
  basis <- fit_basis(tr,arm)
  path <- file.path(root,paste0(package,"-",arm,".rds"))
  if(mode=="fit") {
    settings <- list(seed=289L,job_id=1L,attempt_id=1L,source_hashes=c(fixture="synthetic"),
      mcmc=c(nchain=2L,nburn=10L,niter=20L,nthin=1L),support_points=10L)
    f <- get(paste0("fit_",tolower(package)))(tr,basis,settings)
    p <- get(paste0("predict_",tolower(package)))(f,pred,basis)
    stopifnot(identical(dim(p$p_full),c(nrow(tr$y),2L)),identical(f$species_id,tr$species_id),
      identical(f$input_hashes,tr$input_hashes),identical(f$source_hashes,settings$source_hashes))
    if(arm=="env") stopifnot(max(abs(p$p_full-p$p_space_zero))==0,all(p$field_mean==0))
    if(package=="occJSDM") {
      stored <- f$native$results_output$psi_output
      if(length(dim(stored))==4L) stored <- apply(stored,c(1,2),mean)
      stopifnot(max(abs(p$p_full-stored))<1e-10)
      # Check block reconstruction against the native new-site median target.
      if(arm=="gp") {
        X <- occ_design(pred,basis)
        native <- suppressMessages(occJSDM::predictNewSites(f$native,X_psi=X,
          X_s=data.frame(longitude=pred$xy[,1],latitude=pred$xy[,2]),
          useEnvCov=TRUE,useSpatial=TRUE,useBiotic=FALSE,verbose=FALSE))
        d <- occ_draw_extractor(f,pred,basis)(seq_len(nrow(tr$y)))
        stopifnot(max(abs(apply(d$p_full,c(1,2),median)-native[2,,]))<1e-8)
      }
    } else {
      stopifnot(f$native$diagonal_verified,all(f$native$calibration$passed),
        max(abs(p$p_full-f$native$objective$probabilities))<1e-6)
      if(arm!="env") {
        flatten <- function(net) unlist(lapply(net,function(x)c(x$weight,x$bias)))
        stopifnot(max(abs(flatten(f$native$network$spatial)-flatten(f$native$first_phase$spatial)))>0,
          f$native$before$spatial_gradient>0)
      }
    }
    saveRDS(list(fit=f,basis=basis,prediction=p),path)
  } else {
    saved <- readRDS(path); f <- saved$fit; basis <- saved$basis
    p <- get(paste0("predict_",tolower(package)))(f,pred,basis)
    stopifnot(max(abs(p$p_full-saved$prediction$p_full))<1e-12)
    idx <- rev(seq_len(nrow(tr$y))); rev <- subset_prediction(pred,idx)
    p2 <- get(paste0("predict_",tolower(package)))(f,rev,basis)
    stopifnot(max(abs(p2$p_full-p$p_full[idx,,drop=FALSE]))<1e-10)
    one <- get(paste0("predict_",tolower(package)))(f,subset_prediction(pred,1L),basis)
    stopifnot(identical(dim(one$p_full),c(1L,2L)),max(abs(one$p_full-p$p_full[1,,drop=FALSE]))<1e-10)
    wrong <- pred; wrong$species_id <- rev(pred$species_id)
    stopifnot(inherits(try(get(paste0("predict_",tolower(package)))(f,wrong,basis),silent=TRUE),"try-error"))
  }
  cat(package,arm,mode,"PASS\n")
}
