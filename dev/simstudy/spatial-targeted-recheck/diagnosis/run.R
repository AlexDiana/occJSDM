#!/usr/bin/env Rscript
# Diagnosis only: leaves the fitted package, saved fits and original study unchanged.
args <- commandArgs(trailingOnly=TRUE)
option <- function(name) {
  hit <- args[startsWith(args,paste0('--',name,'='))]
  if(length(hit)!=1L) stop('Supply exactly one --',name,'=PATH')
  substring(hit,nchar(name)+4L)
}
stopifnot(length(args)==2L,all(sub('^--([^=]+)=.*$','\\1',args)%in%c('study','out')))
study <- normalizePath(option('study'));out <- option('out')
if(dir.exists(out)) stop('Use a fresh output directory; existing results are preserved.')
dir.create(out,recursive=TRUE);out <- normalizePath(out)
scripts <- dirname(normalizePath(sub('^--file=','',commandArgs()[startsWith(commandArgs(),'--file=')][1])))
source(file.path(scripts,'oracle.R'))
Rcpp::sourceCpp(file.path(scripts,'quadrature.cpp'))
source(file.path(scripts,'check.R'),local=TRUE)
source_files <- file.path(scripts,c('run.R','oracle.R','quadrature.cpp','check.R'))
write.csv(data.frame(file=basename(source_files),md5=unname(tools::md5sum(source_files))),file.path(out,'source-hashes.csv'),row.names=FALSE)
sel <- read.csv(file.path(study,'summary-final/selected-fits.csv'))
sel <- sel[sel$knots==20,]
sel$selected_file <- file.path(study,sel$phase,paste0(sel$key,'-result.rds'))
stopifnot(nrow(sel)==27L,!anyDuplicated(sel$key))
provenance <- list()
rows <- list()
for(i in seq_len(nrow(sel))) {
  job <- sel[i,];input <- readRDS(file.path(study,'inputs',paste0(job$community,'.rds')))
  stopifnot(unname(tools::md5sum(job$selected_file))==job$selected_md5)
  r <- readRDS(job$selected_file); t <- input$truth
  input_file <- file.path(study,'inputs',paste0(job$community,'.rds'))
  stopifnot(unname(tools::md5sum(input_file))==r$job$input_md5,input$settings$n==100L,input$settings$M==2L)
  fit_file <- sub('-result.rds','-fit.rds',job$selected_file)
  saved <- if(job$arm!='binary') readRDS(fit_file) else NULL
  if(!is.null(saved)) stopifnot(identical(saved$job,r$job),identical(saved$mcmc,r$mcmc),identical(saved$fit_hashes,r$fit_hashes))
  fit <- if(is.null(saved))NULL else saved$fit
  provenance[[i]] <- data.frame(key=job$key,phase=job$phase,result_md5=job$selected_md5,input_md5=r$job$input_md5,
    fit_md5=if(is.null(saved))NA_character_ else unname(tools::md5sum(fit_file)))
  getblock <- function(nm) r$elements[r$elements$metric==nm,]
  for(s in 1:8) {
    th <- getblock('collection_coefficient');e <- r$species[s,]
    z_sample <- rep(t$z[,s],each=2)
    d <- data.frame(community=job$community,arm=job$arm,species=s,target=t$target_prevalence[s],
      occupied=e$occupied,psi_estimate=e$estimate,B0_true=t$B0[s],B0_estimate=getblock('intercept')$estimate[s],
      B0_lower=getblock('intercept')$lower[s],B0_upper=getblock('intercept')$upper[s],
      w_true_occupied=sum(t$w[,s]*z_sample),w_true_unoccupied=sum(t$w[,s]*(1-z_sample)),
      theta_true=mean(t$theta[,s]),theta0_true=t$theta0[s],
      theta_estimate=if(is.null(fit))NA else mean(fit$results_output$theta_output[,s]),
      theta0_estimate=if(is.null(fit))NA else getblock('theta0')$estimate[s],
      collection_intercept_estimate=if(is.null(fit))NA else th$estimate[2*s-1],
      collection_slope_estimate=if(is.null(fit))NA else th$estimate[2*s],
      p_true=mean(t$p[,s]),p_estimate=if(is.null(fit))NA else mean(getblock('p')$estimate[(2*s-1):(2*s)]),
      q_true=if(job$arm=='binary')NA else mean(t$q[[job$arm]][,s]),
      q_estimate=if(is.null(fit))NA else mean(getblock('q')$estimate[(2*s-1):(2*s)]),
      z_estimate=if(is.null(fit))NA else mean(fit$results_output$z_output[,s]),
      w_estimate=if(is.null(fit))NA else mean(fit$results_output$w_output[,s]))
    if(!is.null(fit)) {
      ro <- fit$results_output
      d$correlation_B0_collection <- cor(as.vector(ro$jsdm_output$B0_output[s,,]),as.vector(ro$beta_theta_output[1,s,,]))
      d$correlation_B0_theta0 <- cor(as.vector(ro$jsdm_output$B0_output[s,,]),as.vector(ro$theta0_output[s,,]))
      d$signed_collection_slope <- mean(ro$beta_theta_output[2,s,,])*t$beta_theta[2,s]
    } else {
      d$correlation_B0_collection <- d$correlation_B0_theta0 <- d$signed_collection_slope <- NA_real_
    }
    pair <- colSums(matrix(t$w[,s],2,100));occupied <- t$z[,s]==1
    d$double_positive_occupied <- sum(pair[occupied]==2)
    d$double_positive_unoccupied <- sum(pair[!occupied]==2)
    rows[[length(rows)+1]] <- d
  }
  cat(job$key,'extracted\n')
}
write.csv(do.call(rbind,provenance),file.path(out,'input-identities.csv'),row.names=FALSE)
d <- do.call(rbind,rows);write.csv(d,file.path(out,'saved-fit-diagnostic.csv'),row.names=FALSE)
num <- names(d)[vapply(d,is.numeric,logical(1))];num <- setdiff(num,c('species','target'))
a <- aggregate(d[num],d[c('arm','target')],function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE))
write.csv(a,file.path(out,'saved-fit-diagnostic-aggregate.csv'),row.names=FALSE)
print(a[a$target<=.05,],row.names=FALSE)

inputs <- sort(list.files(file.path(study,'inputs'),pattern='rds$',full.names=TRUE))
rows <- list();resolution <- list()
for(path in inputs) {
 input <- readRDS(path);community <- sub('.rds$','',basename(path));t <- input$truth
 for(s in 1:4) {
  offset <- t$X*t$B[s]+t$field[,s]
  for(arm in c('binary','field','low','high')) {
   L <- site_likelihood(input,s,arm)
   for(sd in c(1,2.5)) {
    x <- intercept_posterior(L,offset,sd)
    fine <- intercept_posterior(L,offset,sd,step=.0125)
    rows[[length(rows)+1]] <- data.frame(community,species=s,target=t$target_prevalence[s],arm,
      nuisance='known',prior_sd=sd,estimate=x$mean,collection_mean=mean(t$theta[,s]),boundary=x$boundary)
    resolution[[length(resolution)+1]] <- data.frame(community,species=s,arm,nuisance='known',prior_sd=sd,
      difference=abs(x$mean-fine$mean),boundary=max(x$boundary,fine$boundary))
   }
  }
  for(arm in c('field','low','high')) {
   x <- collection_marginal(input,s,arm,step=.1)
   coarse <- collection_marginal(input,s,arm,step=.2)
   rows[[length(rows)+1]] <- cbind(data.frame(community,species=s,target=t$target_prevalence[s],arm,
      nuisance='collection_intercept_unknown'),x)
   resolution[[length(resolution)+1]] <- data.frame(community,species=s,arm,
      nuisance='collection_intercept_unknown',prior_sd=x$prior_sd,
      difference=abs(x$estimate-coarse$estimate),boundary=pmax(x$boundary,coarse$boundary))
  }
  cat(format(Sys.time()),community,s,'complete\n');flush.console()
 }
}
d <- do.call(rbind,rows);write.csv(d,file.path(out,'oracle-results.csv'),row.names=FALSE)
v <- do.call(rbind,resolution);write.csv(v,file.path(out,'oracle-resolution.csv'),row.names=FALSE)
a <- aggregate(d[c('estimate','collection_mean')],d[c('target','arm','nuisance','prior_sd')],mean)
write.csv(a,file.path(out,'oracle-aggregate.csv'),row.names=FALSE);print(a,row.names=FALSE)
print(apply(v[c('difference','boundary')],2,max))
stopifnot(max(v$difference)<1e-6,max(v$boundary)<1e-6)
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
