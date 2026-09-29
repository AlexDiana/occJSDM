#!/usr/bin/env Rscript
# Occupancy-intercept prior study runner: one complete fit per community.
#
#   Rscript run.R --repo=REPO --study=STUDY --phase=A1|A2|B --sd=1|2|3|5
#     --schedule=initial|long|pilot --workers=N [--keys=k1,k2] [--inputs-root=DIR]
#     [--archives=DIR]
#
# Each fit evaluates the control runner's own fitting statement (parsed from the
# tracked script whose md5 was recorded with the control fits), with the fitting
# function wrapped only to add listPriors$sigma_b0. Inputs, input RNG states,
# MCMC settings and priors come from the control archives; A1/B inputs are read
# from the copy under --inputs-root and must match the recorded md5. The
# installed library must match the committed fingerprint before it is loaded.
# One sampler thread per worker, at most eight workers, one lock per key.
# Exits non-zero if any job fails. Launch long batches detached.
args <- commandArgs(trailingOnly=TRUE)
repo_arg <- sub('^--repo=','',grep('^--repo=',args,value=TRUE))
if(length(repo_arg)!=1L) stop('Missing or repeated --repo')
source(file.path(normalizePath(repo_arg),'dev/simstudy/occupancy-intercept-prior/jobs.R'))
opt <- parse_run_args(args)
repo <- normalizePath(opt$repo);study <- normalizePath(opt$study)
archives <- normalizePath(opt$archives);inputs_root <- normalizePath(opt$inputs_root)
phase <- opt$phase;sd <- opt$sd;schedule <- opt$schedule

# Refuse any build other than the tested one before anything is loaded from it.
source_revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))

Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',
  OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(file.path(study,'library'),.libPaths()))
suppressPackageStartupMessages(library(occJSDM))
pkg <- normalizePath(find.package('occJSDM'))
stopifnot(pkg==normalizePath(file.path(study,'library/occJSDM')))
RcppParallel::setThreadOptions(numThreads=1)

fit_hashes <- hash_files(production_files(study,pkg),study)
runner <- archived_runner(phase_runner_name(phase),repo,archives)
expr <- with_intercept_prior(runner$expr)
fit_expression <- deparse(expr,width.cutoff=500L)
mcmc <- schedule_mcmc(phase,schedule,archives)
priors <- phase_priors(phase,archives)
jobs <- phase_jobs(phase,archives,inputs_root,opt$keys)
# pr11's own post-fit invariant (trait preprocessing, species order, no spatial support).
pr11 <- NULL
if(phase %in% c('A1','B')) {
  helpers <- file.path(repo,PR11_HELPERS)
  recorded <- recorded_hash(readRDS(runner$settings_path)$source_hashes,PR11_HELPERS)
  stopifnot(identical(unname(tools::md5sum(helpers)),recorded))
  pr11 <- new.env(parent=globalenv());sys.source(helpers,envir=pr11)
}
script_hashes <- c(hash_files(file.path(repo,runner_script_files(phase)),repo),
  hash_files(unlist(mcmc_sources(phase,archives)),archives,'archives/'))

run_job <- function(spec) {
  dest <- fit_path(study,phase,sd,schedule,spec$key)
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  tag <- paste(spec$key,phase,paste0('sd',sd),schedule)
  metadata <- fit_metadata(phase=phase,key=spec$key,sd=sd,schedule=schedule,mcmc=mcmc,job=spec$job,
    input_md5=spec$input_md5,runner=runner[c('name','file','md5','settings')],fit_expression=fit_expression,
    phase_priors=priors,source_revision=source_revision,fit_hashes=fit_hashes,script_hashes=script_hashes)
  # Informational only; never part of the resume comparison.
  info <- list(job=spec$job,input_file=spec$input_file,runner_file_absolute=runner$path,
    repo=repo,study=study,archives=archives,library=pkg)
  lock <- acquire_fit_lock(dest);on.exit(release_fit_lock(lock),add=TRUE)
  if(existing_fit_status(dest,metadata)=='resume') {
    cat(format(Sys.time()),tag,'already complete with identical settings; resumed\n');flush.console()
    return(spec$key)
  }
  checked_input(spec$input_file,spec$input_md5)
  input <- readRDS(spec$input_file)
  cat(format(Sys.time()),tag,'started\n');flush.console()
  result <- evaluate_fit(expr,input,spec$job,mcmc,bindings=list(.fitter=make_fitter(sd),priors=priors))
  saved <- c(metadata,info,result,list(elapsed_seconds=as.numeric(difftime(result$finished,result$started,units='secs')),
    input_md5_after=unname(tools::md5sum(spec$input_file)),loaded_library=normalizePath(find.package('occJSDM')),
    session=sessionInfo()))
  post_fit <- function(s) {
    check_fit_prior(s$fit,sd)
    if(phase=='A2') {
      if(!identical(s$fit$infos$spatial_sd_prior$type,'inverse_gamma')) stop('Fit does not record the inverse-gamma amplitude prior')
    } else pr11$check_current_fit(s$fit,input,spec$job)
    if(!identical(s$input_md5_after,spec$input_md5)) stop('Input changed during fitting: ',spec$input_file)
    if(!identical(s$loaded_library,pkg)) stop('Loaded library changed during fitting')
  }
  finalise_fit(saved,dest,post_fit)
  cat(format(Sys.time()),tag,'complete; fit seconds',round(saved$elapsed_seconds,1),
    '; warnings',length(saved$warnings),'\n');flush.console()
  spec$key
}
work <- function(spec) tryCatch(run_job(spec),error=function(e) {
  cat(spec$key,'ERROR:',conditionMessage(e),'\n');list(key=spec$key,error=conditionMessage(e))
})
status <- if(opt$workers==1L) lapply(jobs,work) else parallel::mclapply(jobs,work,
  mc.cores=opt$workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
out <- dirname(fit_path(study,phase,sd,schedule,'x'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
saveRDS(status,file.path(out,paste0('status-',format(Sys.time(),'%Y%m%d%H%M%S'),'-',Sys.getpid(),'.rds')))
finish_jobs(status,names(jobs),paste(phase,paste0('sd',sd),schedule))
