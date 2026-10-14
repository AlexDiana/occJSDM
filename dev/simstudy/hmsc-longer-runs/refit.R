#!/usr/bin/env Rscript
# One Hmsc refit per process (README.md in this directory): the frozen worker's
# model for one selected job, sampled four times as long as its longer attempt.
#
#   Rscript refit.R --repo=REPO --root=EXTENSION_ARCHIVE --study=ARCHIVE --job-number=N
#
# Reads the job's input from the extension archive (md5-checked, read only),
# fits with fit_hmsc() (common.R; the frozen worker's statements, Hmsc 3.3-7
# from ARCHIVE/library, seed 26100000 + 100 * job number + 3), extracts the
# draws and runs the extension's own parameter_diagnostics(). Writes
# ARCHIVE/jobs/<job>/attempt-3-parameters.rds, attempt-3-diagnostics.csv and
# result.rds under a per-job lock. Never opens truth/. A job already saved with
# identical settings is reported as resumed; one with different settings is
# never overwritten. Exit 0 on success or resume, 1 on failure.
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(here,'common.R'))

settings_record <- function(job) list(attempt=ATTEMPT,seed=job_seed(job$job_number),warmup=SCHEDULE$warmup,draws=SCHEDULE$draws,
  thin=SCHEDULE$thin,chains=SCHEDULE$chains,factors=2L,latent_traits=0L,link='probit',hmsc_version=HMSC_VERSION,
  fit_job_md5=FIT_JOB_MD5,code_md5=CODE_MD5,input_md5=job$input_md5,job=job$job,job_number=job$job_number)

refit_main <- function(args) {
  o <- parse_options(args,known=c('repo','root','study','job-number'),required=c('repo','root','study','job-number'))
  repo <- normalizePath(o$repo,mustWork=TRUE);root <- normalizePath(o$root,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this refit.R')
  if(!grepl('^[0-9]+$',o$`job-number`)) stop('--job-number must be a whole number')
  jobs <- read_jobs(repo);job <- jobs[jobs$job_number==as.integer(o$`job-number`),,drop=FALSE]
  if(nrow(job)!=1L) stop('Job number ',o$`job-number`,' is not in ',JOBS_FILE)
  manifest <- extension_manifest(root);m <- manifest[manifest$job_number==job$job_number,,drop=FALSE]
  if(nrow(m)!=1L || !identical(m$job,job$job) || !identical(m$input_md5,job$input_md5) || m$package!='Hmsc')
    stop('Manifest row differs from the selection for job ',job$job)
  check_fit_job(root);code <- load_code(repo)
  pkg <- load_hmsc(study)
  out <- study_job_dir(study,job$job);dir.create(out,recursive=TRUE,showWarnings=FALSE)
  lock <- file.path(out,'job.lock')
  if(!dir.create(lock,showWarnings=FALSE)) stop('Another process holds ',lock)
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  settings <- settings_record(job);tag <- paste('job',job$job_number,job$job)
  dest <- file.path(out,'result.rds')
  if(file.exists(dest)) {
    saved <- readRDS(dest)
    if(!identical(saved$settings,settings)) stop('Refusing to overwrite ',dest,': settings differ')
    cat(format(Sys.time()),tag,RESUMED_TEXT,'\n');return(0L)
  }
  input <- job_input(root,job)
  warning_messages <- character()
  cat(format(Sys.time()),tag,'started; seed',settings$seed,'; warmup',SCHEDULE$warmup,'draws',SCHEDULE$draws,'thin',SCHEDULE$thin,'chains',SCHEDULE$chains,'\n');flush.console()
  started <- Sys.time()
  fit <- withCallingHandlers(fit_hmsc(input,SCHEDULE,settings$seed),
    warning=function(w) {warning_messages <<- c(warning_messages,conditionMessage(w));invokeRestart('muffleWarning')})
  fitted <- Sys.time()
  p <- code$extract_bayesian(fit,'Hmsc');rm(fit);invisible(gc(FALSE))
  stopifnot(identical(dim(p$beta)[3:4],c(SCHEDULE$draws,SCHEDULE$chains)),p$link=='probit')
  diagnostics <- code$parameter_diagnostics(p,input)
  finished <- Sys.time()
  saveRDS(p,file.path(out,sprintf('attempt-%d-parameters.rds',ATTEMPT)))
  utils::write.csv(diagnostics$table,file.path(out,sprintf('attempt-%d-diagnostics.csv',ATTEMPT)),row.names=FALSE)
  result <- list(ok=TRUE,diagnostic_pass=diagnostics$passed,integration_error=diagnostics$integration_error,
    max_rhat=diagnostics$max_rhat,min_bulk_ess=diagnostics$min_bulk_ess,min_tail_ess=diagnostics$min_tail_ess,
    settings=settings,job=m,warnings=unique(warning_messages),started=started,fitted=fitted,finished=finished,
    fit_seconds=as.numeric(difftime(fitted,started,units='secs')),elapsed_seconds=as.numeric(difftime(finished,started,units='secs')),
    iterations_per_chain=iterations_per_chain(SCHEDULE),input_md5=job$input_md5,truth_used=FALSE,library=pkg,
    script_md5=hash_files(file.path(repo,STUDY_REL,c('common.R','refit.R')),repo),session=sessionInfo())
  tmp <- paste0(dest,'.partial');saveRDS(result,tmp);stopifnot(file.rename(tmp,dest))
  cat(format(Sys.time()),tag,'complete; fit seconds',round(result$fit_seconds,1),'; diagnostics',if(diagnostics$passed) 'passed' else 'flagged',
    '; failing quantities',sum(!diagnostics$table$pass),'; warnings',length(result$warnings),'\n');flush.console()
  0L
}

if(sys.nframe()==0L) {
  status <- tryCatch(refit_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
