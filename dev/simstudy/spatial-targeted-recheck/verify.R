#!/usr/bin/env Rscript
# Independent full-draw audit. Does not source the generator or scorer.
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(name,default=NULL) {
  hit <- grep(paste0('^--',name,'='),args,value=TRUE)
  if(length(hit)>1L)stop('Repeated option')
  if(!length(hit)){if(is.null(default))stop('Missing --',name);return(default)}
  substring(hit,nchar(name)+4L)
}
study <- normalizePath(getarg('study'));phase <- getarg('phase','pilot')
workers <- as.integer(getarg('workers','2'));stopifnot(workers %in% 1:4)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()));suppressPackageStartupMessages(library(occJSDM))
if(phase=='selected') {
  selection <- read.csv(getarg('selection'),stringsAsFactors=FALSE)
  files <- selection$selected_file
  stopifnot(length(files)==81L,!anyDuplicated(selection$key))
} else {
  stopifnot(phase %in% c('pilot','initial','long'))
  files <- list.files(file.path(study,phase),pattern='-result.rds$',full.names=TRUE)
  stopifnot(length(files)>0L)
}
check <- function(path) {
  result <- readRDS(path)
  fitfile <- sub('-result.rds$','-fit.rds',path)
  raw <- readRDS(fitfile);input <- readRDS(result$job$input_file);fit <- raw$fit
  stopifnot(identical(raw$job,result$job),identical(raw$fit_hashes,result$fit_hashes),
    identical(unname(tools::md5sum(result$job$input_file)),result$job$input_md5),
    identical(unname(tools::md5sum(names(result$fit_hashes))),unname(result$fit_hashes)))
  n <- nrow(input$truth$z);s <- ncol(input$truth$z)
  jo <- fit$results_output$jsdm_output
  ni <- dim(jo$B0_output)[2];nc <- dim(jo$B0_output)[3]
  stopifnot(ni==result$mcmc$niter,nc==result$mcmc$nchain)
  stopifnot(identical(fit$infos$speciesNames,colnames(input$truth$z)),
    max(abs(fit$infos$OTU-input$data[[result$job$arm]]$OTU))==0)
  native <- suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid,fit$infos$list_Xs))
  predictions <- array(NA_real_,c(n*s,ni,nc))
  for(ch in seq_len(nc)) for(i in seq_len(ni)) {
    field <- occJSDM:::KsBproduct(native$Ks_all[,,jo$idx_ls_output[i,ch]],
      matrix(jo$Bs_output[,,i,ch],result$job$knots,s),fit$infos$list_Xs$Xs_centers)
    environmental <- tcrossprod(fit$X_psi[,1],jo$B_output[1,,i,ch])
    intercept <- outer(rep(1,n),jo$B0_output[,i,ch])
    predictions[,i,ch] <- as.vector(plogis(intercept+environmental+field))
  }
  mu <- matrix(apply(predictions,1,mean),n,s)
  mean_difference <- max(abs(mu-result$probability))
  true_psi <- plogis(outer(rep(1,n),input$truth$B0)+
    tcrossprod(input$truth$X,input$truth$B)+input$truth$field)
  stopifnot(max(abs(true_psi-input$truth$psi))<1e-12)
  interval <- t(vapply(seq_len(n*s),function(i)
    as.numeric(quantile(predictions[i,,],c(.025,.975))),numeric(2)))
  cells <- result$elements[result$elements$metric=='occupancy',]
  coverage <- interval[,1]<=as.vector(true_psi) & as.vector(true_psi)<=interval[,2]
  interval_difference <- max(abs(interval-cbind(cells$lower,cells$upper)))
  stopifnot(identical(as.logical(cells$covered),coverage))
  masks <- list(all=rep(TRUE,n*s),low=as.vector(true_psi)<.2,
    medium=as.vector(true_psi)>=.2 & as.vector(true_psi)<=.8,high=as.vector(true_psi)>.8)
  for(p in c(.01,.05,.25,.75)) masks[[paste0('prevalence_',p*100,'pct')]] <-
    rep(abs(colMeans(true_psi)-p)<1e-10,each=n)
  group_difference <- 0
  for(g in names(masks)) {
    idx <- which(masks[[g]])
    row <- result$groups[result$groups$metric=='occupancy' & result$groups$group==g,]
    stopifnot(nrow(row)==1L,row$n==length(idx))
    error <- as.vector(mu-true_psi)[idx]
    expected <- c(bias=mean(error),mae=mean(abs(error)),rmse=sqrt(mean(error^2)),coverage=mean(coverage[idx]))
    group_difference <- max(group_difference,max(abs(as.numeric(row[names(expected)])-expected)))
  }
  stopifnot(all(is.finite(c(mean_difference,interval_difference,group_difference))),
    max(mean_difference,interval_difference,group_difference)<1e-10)
  cat(result$job$key,'verified\n');flush.console()
  data.frame(key=result$job$key,probability_difference=mean_difference,
    interval_difference=interval_difference,group_difference=group_difference,
    cells=n*s,groups=length(masks),file=path,result_md5=unname(tools::md5sum(path)),
    fit_md5=unname(tools::md5sum(fitfile)))
}
work <- function(path)tryCatch(check(path),error=function(e)list(file=path,error=conditionMessage(e)))
ans <- if(workers==1L)lapply(files,work) else parallel::mclapply(files,work,mc.cores=workers,mc.preschedule=FALSE)
errors <- ans[!vapply(ans,is.data.frame,logical(1))]
if(length(errors)){print(errors);stop('Independent verification failed')}
output <- getarg('out',file.path(study,paste0('verification-',phase,'-independent.csv')))
if(file.exists(output))stop('Refusing to replace verification: ',output)
write.csv(do.call(rbind,ans),output,row.names=FALSE)
cat('Verified',length(ans),'fits using native full-draw reconstruction.\n')
