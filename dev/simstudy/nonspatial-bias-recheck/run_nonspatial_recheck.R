#!/usr/bin/env Rscript
# A targeted diagnostic. It does not modify the fitted package or choose priors.
options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
arg <- function(k, d) { x <- grep(paste0('^--',k,'='),args,value=TRUE); if(length(x)) sub(paste0('^--',k,'='),'',x[1]) else d }
root <- normalizePath(arg('root','.'),mustWork=TRUE)
source_dir <- normalizePath(arg('source',file.path(root,'source-main')),mustWork=TRUE)
lib <- normalizePath(arg('library',file.path(root,'library')),mustWork=TRUE)
out <- arg('out',file.path(root,'results'))
reps <- as.integer(arg('replicates','2')); workers <- as.integer(arg('workers','4'))
burn <- as.integer(arg('burn','1000')); iter <- as.integer(arg('iter','1000'))
mode <- arg('mode','all')
stopifnot(reps>0,workers>0,burn>=0,iter>=20, mode %in% c('all','collection','q'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out)
.libPaths(c(lib,.libPaths()))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages(library(occJSDM))
source(file.path(source_dir,'tests/testthat/helper-simstudy.R'))
RcppParallel::setThreadOptions(numThreads=1)
all_scenarios <- simstudy_scenarios()
find_scenario <- function(label) Filter(function(s) s$label==label,all_scenarios)[[1]]
labels <- c(if(mode %in% c('all','collection')) c('traits_isolated','occupancy'), if(mode %in% c('all','q')) c('qnear_K3','qnear_K30','qfar_K3','qfar_K30'))
settings <- list(source=source_dir,library=lib,source_revision=readLines(file.path(root,'source-revision.txt')),mcmc=list(nchain=2,nburn=burn,niter=iter,nthin=1),threads=1,spatial=FALSE,session=sessionInfo())
paths <- c(file.path(source_dir,c('R/runOccJSDM.R','R/simulateData.R','R/jsdmfun.R','src/functions.cpp','src/jsdm.cpp','src/rng.h','tests/testthat/helper-simstudy.R')),file.path(lib,'occJSDM/libs/occJSDM.so'))
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_path <- if(length(script_arg))normalizePath(sub('^--file=','',script_arg[1])) else NA_character_
settings$runner_path <- script_path
settings$runner_source <- if(!is.na(script_path))readLines(script_path,warn=FALSE) else NULL
settings$hashes <- tools::md5sum(c(paths,if(!is.na(script_path))script_path))
saveRDS(settings,file.path(out,'settings.rds'))
jobs <- list()
for(label in labels) {
  scenario <- find_scenario(label)
  scenario$useSpatField <- FALSE; scenario$n_supportpoints <- NULL
  for(r in seq_len(reps)) {
    input_file <- file.path(out,sprintf('%s-data-%02d.rds',label,r))
    if(!file.exists(input_file)) {
      seed <- simstudy_seed_for(scenario,r)
      truth <- draw_truth(scenario,seed)
      sim <- simstudy_simulate(truth)
      input <- list(scenario=scenario,replicate=r,seed=seed,truth=truth,sim=sim,fit_rng=.Random.seed)
      stopifnot(!truth$jsdmParams$useSpatField)
      saveRDS(input,input_file)
    }
    priors <- if(label %in% c('traits_isolated','occupancy')) list(var2=list(b_betatheta_slope_var=2),var05=list(b_betatheta_slope_var=.5),var01=list(b_betatheta_slope_var=.1)) else list(q20=list(a_q=1,b_q=20),q9=list(a_q=1,b_q=9))
    if(label %in% c('traits_isolated','qnear_K3','qfar_K3')) priors$p32 <- list(b_betatheta_slope_var=2,a_p=3,b_p=2,a_q=1,b_q=20)
    for(prior_name in names(priors)) {
      key <- sprintf('%s-%s-%02d',label,prior_name,r)
      jobs[[key]] <- list(key=key,input_file=input_file,prior_name=prior_name,priors=priors[[prior_name]],mcmc=settings$mcmc,out=out)
    }
  }
}

trace_summary <- function(x) {
  x <- as.matrix(x)
  c(rhat=posterior::rhat(x),ess_bulk=posterior::ess_bulk(x),ess_mean=posterior::ess_mean(x),mcse=posterior::mcse_mean(x),chain_gap=abs(mean(x[,1])-mean(x[,2])))
}
score_fit <- function(fit,input,prior) {
  sim <- input$sim; scenario <- input$scenario; tp <- sim$true_params; jp <- tp$jsdmParams_true
  ro <- fit$results_output; jo <- ro$jsdm_output
  S <- scenario$S; n <- scenario$n; ni <- dim(jo$B0_output)[2]; nc <- dim(jo$B0_output)[3]
  stopifnot(fit$infos$ps==0,fit$infos$model==scenario$model,identical(fit$infos$speciesNames,colnames(sim$data_list$OTU)))
  info <- sim$data_list$info
  site_info <- info[!duplicated(info$Site),,drop=FALSE]
  site_info <- site_info[order(site_info$Site),,drop=FALSE]
  rawx <- as.matrix(site_info[paste0('X_psi.EnvCov.',1:scenario$ncov_psi)])
  stopifnot(max(abs(unname(scale(rawx))-unname(fit$X_psi)))<1e-10)
  # The generating environmental design is already standardised; B0 is on
  # the same scale. Check its reconstructed true linear predictor explicitly.
  true_eta <- sweep(fit$X_psi %*% jp$B + jp$U %*% jp$L,2,jp$B0,'+')
  stopifnot(max(abs(true_eta-jp$eta))<1e-10)
  true_psi <- plogis(true_eta)
  if(scenario$model=='two_stage') {
    sample_info <- info[!duplicated(info[c('Site','Sample')]),,drop=FALSE]
    sample_info <- sample_info[order(sample_info$Site,sample_info$Sample),,drop=FALSE]
  } else sample_info <- info[order(info$Site,seq_len(nrow(info))),,drop=FALSE]
  rawt <- sample_info$X_theta
  stopifnot(max(abs(unname(cbind(1,scale(rawt)))-unname(fit$X_theta)))<1e-10)
  true_bt <- tp$beta_theta_true
  true_bt[1,] <- true_bt[1,]+mean(rawt)*true_bt[2,]
  true_bt[2,] <- sd(rawt)*true_bt[2,]
  true_theta <- plogis(fit$X_theta %*% true_bt)
  stopifnot(max(abs(true_theta-plogis(cbind(1,rawt)%*%tp$beta_theta_true)))<1e-10)
  blocks <- list(B0=list(draws=jo$B0_output,truth=jp$B0),collection_intercept=list(draws=ro$beta_theta_output[1,,,,drop=FALSE],truth=true_bt[1,]),collection_slope=list(draws=ro$beta_theta_output[2,,,,drop=FALSE],truth=true_bt[2,]),theta0=list(draws=ro$theta0_output,truth=input$truth$params$theta0))
  retention <- c(p=pnorm(log(1.5),5,1,lower.tail=FALSE),q=pnorm(log(1.5),1.5,1,lower.tail=FALSE))
  if(scenario$model=='two_stage') {
    blocks$p <- list(draws=ro$p_output,truth=tp$p_true*retention['p'])
    blocks$q <- list(draws=ro$q_output,truth=tp$q_true*retention['q'])
    blocks$q_nominal <- list(draws=ro$q_output,truth=tp$q_true)
  }
  element_rows <- list(); group_rows <- list(); group_traces <- list()
  add_group <- function(metric,group,truth,x) {
    dx <- dim(x); stopifnot(length(dx)==3,length(truth)==dx[1])
    means <- rowMeans(matrix(x,nrow=dx[1]))
    trace <- apply(x,c(2,3),mean)
    errors <- means-as.vector(truth)
    group_rows[[length(group_rows)+1L]] <<- data.frame(metric=metric,group=group,n_elements=length(truth),truth=mean(truth),estimate=mean(means),bias=mean(errors),mae=mean(abs(errors)),rmse=sqrt(mean(errors^2)),as.list(trace_summary(trace)))
    group_traces[[paste(metric,group,sep=':')]] <<- trace
  }
  for(nm in names(blocks)) {
    b <- blocks[[nm]]; truth <- as.vector(b$truth); x <- array(b$draws,c(length(truth),ni,nc))
    estimate <- rowMeans(matrix(x,nrow=length(truth)))
    diagnostics <- t(vapply(seq_along(truth),function(i) trace_summary(x[i,,]),numeric(5)))
    element_rows[[nm]] <- data.frame(metric=nm,element=seq_along(truth),truth=truth,estimate=estimate,bias=estimate-truth,mae=abs(estimate-truth),diagnostics)
    add_group(nm,'all',truth,x)
    if(nm=='collection_slope') for(sg in c(-1,0,1)) {
      ix <- which(sign(truth)==sg)
      if(length(ix)) add_group(nm,c('negative','zero','positive')[sg+2],truth[ix],x[ix,,,drop=FALSE])
    }
  }
  psi_draws <- array(NA_real_,c(n*S,ni,nc)); theta_draws <- array(NA_real_,c(nrow(fit$X_theta)*S,ni,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    eta <- sweep(fit$X_psi%*%jo$B_output[,,it,ch]+jo$U_output[,,it,ch]%*%jo$L_output[,,it,ch],2,jo$B0_output[,it,ch],'+')
    psi_draws[,it,ch] <- as.vector(plogis(eta))
    theta_draws[,it,ch] <- as.vector(plogis(fit$X_theta%*%ro$beta_theta_output[,,it,ch]))
  }
  reconstructed <- matrix(rowMeans(matrix(psi_draws,nrow=n*S)),n,S)
  psi_difference <- max(abs(reconstructed-ro$psi_output))
  theta_difference <- max(abs(rowMeans(matrix(theta_draws,nrow=length(true_theta)))-as.vector(ro$theta_output)))
  stopifnot(psi_difference<1e-10,theta_difference<1e-10)
  add_group('occupancy','all',as.vector(true_psi),psi_draws)
  add_group('collection_probability','all',as.vector(true_theta),theta_draws)
  bands <- list(low=true_psi<.2,medium=true_psi>=.2 & true_psi<=.8,high=true_psi>.8)
  for(band in names(bands)) {
    ix <- which(bands[[band]])
    if(length(ix)) add_group('occupancy',band,as.vector(true_psi)[ix],psi_draws[ix,,,drop=FALSE])
  }
  prevalence <- colMeans(true_psi)
  for(group in c('rare_below_20pct','common_20pct_or_more')) {
    species <- if(group=='rare_below_20pct') which(prevalence<.2) else which(prevalence>=.2)
    ix <- unlist(lapply(species,function(s) ((s-1L)*n+1L):(s*n)))
    if(length(ix)) add_group('occupancy',group,as.vector(true_psi)[ix],psi_draws[ix,,,drop=FALSE])
  }
  species_rows <- data.frame(species=colnames(sim$data_list$OTU),prevalence=prevalence,occupancy_estimate=colMeans(reconstructed),bias=colMeans(reconstructed-true_psi),mae=colMeans(abs(reconstructed-true_psi)),rmse=sqrt(colMeans((reconstructed-true_psi)^2)))
  oracle <- NULL
  if(scenario$model=='two_stage') {
    y <- (sim$data_list$OTU>=1)*1
    rows <- list()
    for(p in seq_len(scenario$P)) for(s in seq_len(S)) for(state in 0:1) {
      ix <- which(info$Primer==p & tp$w_true[info$Sample,s]==state)
      a <- if(state==0) (prior$a_q %||% 1) else (prior$a_p %||% 5)
      b <- if(state==0) (prior$b_q %||% 20) else (prior$b_p %||% 1)
      rows[[length(rows)+1L]] <- data.frame(metric=if(state==0)'q' else 'p',primer=p,species=s,n=length(ix),successes=sum(y[ix,s]),estimate=(a+sum(y[ix,s]))/(a+b+length(ix)),truth=if(state==0)tp$q_true[p,s]*retention['q'] else tp$p_true[p,s]*retention['p'])
    }
    oracle <- do.call(rbind,rows)
  }
  list(groups=do.call(rbind,group_rows),elements=do.call(rbind,element_rows),species=species_rows,traces=group_traces,known_w=oracle,retention=retention,reconstruction=c(psi=psi_difference,theta=theta_difference))
}
run_job <- function(job) {
  result_file <- file.path(job$out,paste0(job$key,'-result.rds'))
  if(file.exists(result_file)) return(job$key)
  input <- readRDS(job$input_file); s <- input$scenario
  assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
  warnings <- character(); started <- Sys.time()
  fit <- withCallingHandlers(suppressMessages(runOccJSDM(input$sim$data_list,listParams=list(n_factors=s$d,n_lattrait=s$gt),listPriors=job$priors,occCovariates=paste0('X_psi.EnvCov.',1:s$ncov_psi),collCovariates='X_theta',spatCovariates=NULL,MCMCparams=job$mcmc)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  saveRDS(list(fit=fit,job=job,warnings=warnings,started=started,finished=Sys.time()),file.path(job$out,paste0(job$key,'-fit.rds')))
  scores <- score_fit(fit,input,job$priors)
  result <- c(list(job=job,warnings=warnings,started=started,finished=Sys.time()),scores)
  temporary <- paste0(result_file,'.tmp')
  saveRDS(result,temporary); stopifnot(file.rename(temporary,result_file))
  cat(format(Sys.time()),job$key,'complete;',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'s; warnings',length(warnings),'\n')
  flush.console()
  job$key
}
if(workers==1L) {
  result <- lapply(jobs,function(j) tryCatch(run_job(j),error=function(e)list(key=j$key,error=conditionMessage(e))))
} else {
  cl <- parallel::makePSOCKcluster(workers,outfile=file.path(out,'workers.log'))
  parallel::clusterExport(cl,c('lib','source_dir'))
  parallel::clusterEvalQ(cl,{
    .libPaths(c(lib,.libPaths()));Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    suppressPackageStartupMessages(library(occJSDM));source(file.path(source_dir,'tests/testthat/helper-simstudy.R'));RcppParallel::setThreadOptions(numThreads=1)
    NULL
  })
  parallel::clusterExport(cl,c('run_job','score_fit','trace_summary'))
  result <- parallel::parLapplyLB(cl,jobs,function(j)tryCatch(run_job(j),error=function(e)list(key=j$key,error=conditionMessage(e))))
  parallel::stopCluster(cl)
}
saveRDS(result,file.path(out,'job-status.rds'))
errors <- Filter(is.list,result)
if(length(errors)) {print(errors);quit(status=1)}
cat('Completed',length(result),'jobs in',out,'\n')
