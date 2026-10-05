# Fresh verifier oracles call upstream directly, without importing production dispatch.
a <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(a[1]); out <- file.path(root,'reports/verification/quickstart')
Sys.setenv(P2A_R_PROJECT=file.path(root,'r-runtime'))
source(file.path(root,'r-runtime','activate.R'))
suppressPackageStartupMessages(library(occJSDM))
suppressPackageStartupMessages(library(jsonlite))
make_summary <- function(f) {
  r <- returnOccupancyRates(f)
  data.frame(species=colnames(r),mean=apply(r,2,mean),median=apply(r,2,quantile,.5),
             q2.5=apply(r,2,quantile,.025),q97.5=apply(r,2,quantile,.975),row.names=NULL)
}
specs <- list(
 binary_single=list(base='binary',S=1L,seed=2401L,oc='habitat_changed',cc=character(),nf=0L,threshold=1L,thin=2L),
 occupancy_changed=list(base='occupancy',S=3L,seed=2402L,oc=c('habitat_changed','habitat_factor'),cc='collection_changed',nf=1L,threshold=2L,thin=2L),
 two_stage_single=list(base='two_stage',S=1L,seed=2403L,oc='habitat_changed',cc='collection_changed',nf=0L,threshold=3L,thin=3L),
 two_stage_changed=list(base='two_stage',S=3L,seed=2404L,oc=c('habitat_changed','habitat_factor'),cc='collection_changed',nf=1L,threshold=3L,thin=3L))
result <- list()
for(nm in names(specs)) {
 s <- specs[[nm]]
 d <- readRDS(file.path(root,'notebooks/quickstart/data',paste0(s$base,'_data.rds')))
 d$OTU <- d$OTU[,seq_len(s$S),drop=FALSE]
 colnames(d$OTU) <- paste0('verifier_species_',seq_len(s$S))
 d$traits <- NULL
 d$info$habitat_changed <- d$info$X_psi.EnvCov.1 * 1.7 + .2
 unit <- if(is.null(d$info$Site))seq_len(nrow(d$info)) else d$info$Site
 d$info$habitat_factor <- factor(ifelse(as.integer(factor(unit))%%2L==0L,'oak','pine'),levels=c('pine','oak'))
 if(!is.null(d$info$X_theta)) d$info$collection_changed <- -d$info$X_theta + .7
 rownames(d$info) <- rownames(d$OTU) <- paste0('vrow_',seq_len(nrow(d$info)))
 input <- file.path(out,paste0(nm,'_data.rds'));saveRDS(d,input)
 opts <- list(occCovariates=as.list(s$oc),collCovariates=as.list(s$cc),threshold=s$threshold,
              listParams=list(n_factors=s$nf),MCMCparams=list(nchain=2L,nburn=20L,niter=20L,nthin=s$thin))
 info <- list(data_path=input,options=opts,seed=s$seed,model=occJSDM:::inferDataModel(d))
 fit <- tryCatch({set.seed(s$seed);runOccJSDM(data=d,occCovariates=s$oc,collCovariates=s$cc,
                     spatCovariates=character(),threshold=s$threshold,listParams=opts$listParams,
                     MCMCparams=opts$MCMCparams,summarisedLatentPresences=TRUE)},error=function(e)e)
 if(inherits(fit,'error')) info$error <- conditionMessage(fit) else {
   saveRDS(fit,file.path(out,paste0(nm,'_oracle_fit.rds')))
   write.csv(returnConvergenceDiagnostics(fit),file.path(out,paste0(nm,'_oracle_diagnostics.csv')),row.names=FALSE,na='NA')
   su <- tryCatch(make_summary(fit),error=function(e)e)
   if(inherits(su,'error')) info$summary_error <- conditionMessage(su) else write.csv(su,file.path(out,paste0(nm,'_oracle_summary.csv')),row.names=FALSE,na='NA')
   info$dimensions <- as.list(dim(d$OTU));info$species <- as.list(fit$infos$speciesNames)
 }
 result[[nm]] <- info
}
# Source-backed invalid input fixtures from a valid fresh verifier input.
d <- readRDS(result$occupancy_changed$data_path)
bad <- d; names(bad)[names(bad)=='info'] <- 'info_extra';saveRDS(bad,file.path(out,'missing_info.rds'))
bad <- d; bad$traitsMatrix <- matrix(1,ncol(d$OTU),1);saveRDS(bad,file.path(out,'traits_alias.rds'))
bad <- d; rownames(bad$OTU) <- rev(rownames(bad$OTU));saveRDS(bad,file.path(out,'row_mismatch.rds'))
write_json(result,file.path(out,'independent_inputs.json'),auto_unbox=TRUE,pretty=TRUE,na='null',digits=17)
cat('Independent direct upstream oracles generated\n')
