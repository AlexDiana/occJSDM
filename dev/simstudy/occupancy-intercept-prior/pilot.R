#!/usr/bin/env Rscript
# Plumbing pilots: one community per phase at sd 3 on the pilot schedule
# (2 chains, 200 burn-in, 200 retained). Never scored.
#
#   Rscript pilot.R --repo=REPO --study=STUDY [--archives=DIR] [--inputs-root=DIR]
#
# Runs run.R once per phase in separate processes, then checks the installed
# library, the recorded intercept prior, input hashes before and after, resume
# and overwrite refusal on a real fit, and that each community's reused control
# fit is readable with the same posterior output layout as the new fit. Through
# the real run.R it also checks that an identical rerun resumes from the same
# checkout and from a copy of the scripts at another absolute path, that a held
# lock refuses the key without being broken, and that a library or revision
# differing from the committed fingerprint is refused before anything is
# fitted. Exits non-zero unless every check passes.
args <- commandArgs(trailingOnly=TRUE)
repo_arg <- sub('^--repo=','',grep('^--repo=',args,value=TRUE))
if(length(repo_arg)!=1L) stop('Missing or repeated --repo')
scripts <- file.path(normalizePath(repo_arg),'dev/simstudy/occupancy-intercept-prior')
source(file.path(scripts,'jobs.R'));source(file.path(scripts,'verify-helpers.R'))
o <- parse_options(args,known=c('repo','study','archives','inputs-root'),required=c('repo','study'))
repo <- normalizePath(o$repo);study <- normalizePath(o$study)
archives <- normalizePath(o$archives %||% dirname(study))
inputs_root <- normalizePath(o$`inputs-root` %||% file.path(archives,'intercept-prior-inputs'))
source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/execution.R'))
PILOTS <- data.frame(phase=c('A1','A2','B'),
  key=c('jsdm-n0100-01','range6-rep01-binary-k100','design-qnear_K6-sites300-01'),stringsAsFactors=FALSE)
SD <- 3
pdir <- file.path(study,'pilot');dir.create(pdir,showWarnings=FALSE)
library_dir <- normalizePath(file.path(study,'library/occJSDM'))

specs <- unlist(lapply(PHASES,function(p) phase_jobs(p,archives,inputs_root)),recursive=FALSE)
inputs <- data.frame(key=names(specs),input_file=vapply(specs,`[[`,'','input_file'),
  recorded_md5=vapply(specs,`[[`,'','input_md5'),stringsAsFactors=FALSE)
inputs$md5_before <- unname(tools::md5sum(inputs$input_file))
stopifnot(identical(inputs$md5_before,inputs$recorded_md5))

run_args <- function(p,key,repo_dir=repo,study_dir=study) c(paste0('--repo=',repo_dir),paste0('--study=',study_dir),paste0('--phase=',p),
  paste0('--sd=',SD),'--schedule=pilot',paste0('--keys=',key),'--workers=1',
  paste0('--archives=',archives),paste0('--inputs-root=',inputs_root))
cat(format(Sys.time()),'launching',nrow(PILOTS),'pilot fits\n');flush.console()
status <- parallel::mclapply(seq_len(nrow(PILOTS)),function(i)
  run_logged_r(file.path(scripts,'run.R'),run_args(PILOTS$phase[i],PILOTS$key[i]),
    file.path(pdir,paste0(PILOTS$phase[i],'-run.log'))),mc.cores=nrow(PILOTS),mc.preschedule=FALSE)
if(!all(unlist(status)==0L)) stop('A pilot fit failed; see ',pdir)
cat(format(Sys.time()),'pilot fits complete\n');flush.console()

# Resume: an identical rerun must leave the saved fit untouched.
a1 <- fit_path(study,'A1',SD,'pilot',PILOTS$key[1]);a1_md5 <- unname(tools::md5sum(a1))
unchanged_a1 <- function() identical(unname(tools::md5sum(a1)),a1_md5)
log_has <- function(log,pattern) any(grepl(pattern,readLines(log),fixed=TRUE))
resume_log <- file.path(pdir,'A1-resume.log')
resume_status <- run_logged_r(file.path(scripts,'run.R'),run_args('A1',PILOTS$key[1]),resume_log)
resume_ok <- resume_status==0L && unchanged_a1() && log_has(resume_log,'resumed')

hardening <- list()
record <- function(check,pass,detail) hardening[[length(hardening)+1L]] <<-
  data.frame(check=check,pass=isTRUE(pass),detail=detail,stringsAsFactors=FALSE)
record('identical rerun resumes',resume_ok,paste('exit',resume_status,'; fit md5 unchanged',unchanged_a1()))

# The same scripts at another absolute path must resume the same fit.
scratch <- tempfile('pilot-hardening-');dir.create(scratch)
relocated <- file.path(scratch,'relocated-repo')
for(f in runner_script_files('A1')) {
  dir.create(dirname(file.path(relocated,f)),recursive=TRUE,showWarnings=FALSE)
  stopifnot(file.copy(file.path(repo,f),file.path(relocated,f)))
}
relocated_log <- file.path(pdir,'A1-resume-relocated.log')
s1 <- run_logged_r(file.path(relocated,RUN_SCRIPT),run_args('A1',PILOTS$key[1],repo_dir=relocated),relocated_log)
record('rerun from a copy of the scripts at another path resumes',
  s1==0L && unchanged_a1() && log_has(relocated_log,'resumed'),
  paste('repo copy',relocated,'; exit',s1,'; fit md5 unchanged',unchanged_a1()))

# A held (or stale) lock refuses the key and is left in place.
lock <- fit_lock_path(a1);stopifnot(dir.create(lock))
writeLines('planted by pilot.R to test lock refusal',file.path(lock,'owner'))
lock_log <- file.path(pdir,'A1-lock-refusal.log')
s2 <- run_logged_r(file.path(scripts,'run.R'),run_args('A1',PILOTS$key[1]),lock_log)
lock_kept <- dir.exists(lock);release_fit_lock(lock)
record('held lock refuses the key without breaking it',
  s2!=0L && lock_kept && unchanged_a1() && log_has(lock_log,'already exists'),
  paste('exit',s2,'; lock kept',lock_kept,'; fit md5 unchanged',unchanged_a1()))

# A study whose library or revision differs from the committed fingerprint is
# refused before any fit (scratch copy; the real study library is untouched).
tampered <- file.path(scratch,'tampered-study');dir.create(tampered)
stopifnot(all(file.copy(file.path(study,c('source','library','source-revision.txt')),tampered,recursive=TRUE)))
description <- file.path(tampered,'library/occJSDM/DESCRIPTION')
cat('Tampered: yes\n',file=description,append=TRUE)
library_log <- file.path(pdir,'A1-fingerprint-library-refusal.log')
s3 <- run_logged_r(file.path(scripts,'run.R'),run_args('A1',PILOTS$key[1],study_dir=tampered),library_log)
record('changed installed library file is refused',
  s3!=0L && log_has(library_log,'fingerprint mismatch') && log_has(library_log,'changed library/occJSDM/DESCRIPTION') &&
    !dir.exists(file.path(tampered,'fits')),paste('exit',s3,'; no fits directory',!dir.exists(file.path(tampered,'fits'))))
stopifnot(file.copy(file.path(study,'library/occJSDM/DESCRIPTION'),description,overwrite=TRUE))
writeLines(strrep('0',40),file.path(tampered,'source-revision.txt'))
revision_log <- file.path(pdir,'A1-fingerprint-revision-refusal.log')
s4 <- run_logged_r(file.path(scripts,'run.R'),run_args('A1',PILOTS$key[1],study_dir=tampered),revision_log)
record('different source revision is refused',
  s4!=0L && log_has(revision_log,'fingerprint mismatch') && log_has(revision_log,'is not the recorded revision') &&
    !dir.exists(file.path(tampered,'fits')),paste('exit',s4,'; no fits directory',!dir.exists(file.path(tampered,'fits'))))
unlink(scratch,recursive=TRUE)
hardening <- do.call(rbind,hardening)

layout <- function(results,mcmc) {
  x <- posterior_layout(results);n <- mcmc$niter/mcmc$nthin;c <- mcmc$nchain
  x$dim <- sub(paste0('(^|x)',n,'x',c,'$'),'\\1ITERxCHAIN',x$dim);x
}
checks <- lapply(seq_len(nrow(PILOTS)),function(i) {
  p <- PILOTS$phase[i];key <- PILOTS$key[i];spec <- specs[[key]]
  dest <- fit_path(study,p,SD,'pilot',key);saved <- readRDS(dest)
  metadata <- saved_fit_metadata(saved)
  changed <- metadata;changed$sd <- 5
  before <- unname(tools::md5sum(dest))
  refused <- tryCatch({existing_fit_status(dest,changed);FALSE},error=function(e) grepl('Refusing to overwrite',conditionMessage(e)))
  refused <- refused && identical(unname(tools::md5sum(dest)),before) && existing_fit_status(dest,metadata)=='resume'
  input <- readRDS(spec$input_file)
  b0 <- saved$fit$results_output$jsdm_output$B0_output
  sel <- control_schedules(p,archives);sel <- sel[sel$key==key,]
  ctl <- readRDS(sel$control_fit);ctl_md5 <- unname(tools::md5sum(sel$control_fit))
  cb0 <- ctl$fit$results_output$jsdm_output$B0_output
  ctl_job <- if(p=='A2') setequal(names(ctl$job),names(spec$job)) && identical(ctl$job[names(spec$job)],spec$job) else identical(ctl$job,spec$job)
  new_layout <- layout(saved$fit$results_output,saved$mcmc);old_layout <- layout(ctl$fit$results_output,ctl$mcmc)
  layout_same <- identical(new_layout,old_layout)
  if(!layout_same) {
    both <- merge(new_layout,old_layout,by='path',all=TRUE,suffixes=c('_pilot','_control'))
    write.csv(both,file.path(pdir,paste0(p,'-layout-difference.csv')),row.names=FALSE)
  }
  data.frame(phase=p,key=key,sd=SD,schedule='pilot',fit_file=dest,
    library_is_study=identical(saved$library,library_dir) && identical(saved$loaded_library,library_dir),
    fit_hashes_current=identical(saved$fit_hashes,hash_files(production_files(study,library_dir),study)),
    intercept_prior_sd=saved$fit$infos$intercept_prior$sd,
    intercept_prior_ok=identical(saved$fit$infos$intercept_prior,list(mean=0,sd=SD)),
    listPriors_added=paste(names(saved$listPriors_added),unlist(saved$listPriors_added),sep='=',collapse=';'),
    mcmc=paste(unlist(saved$mcmc),collapse='/'),
    input_md5_ok=identical(saved$input_md5,spec$input_md5) && identical(saved$input_md5_after,spec$input_md5),
    rng_initial_is_input_state=identical(saved$rng_initial,input$fit_rng),
    b0_dim=paste(dim(b0),collapse='x'),posterior_finite=all(is.finite(b0)),
    warnings=length(saved$warnings),elapsed_seconds=round(saved$elapsed_seconds,2),
    overwrite_refused=refused,resume_ok=if(p=='A1') resume_ok else NA,
    control_schedule=sel$schedule,control_fit=sel$control_fit,control_fit_md5=ctl_md5,
    control_md5_matches_record=is.na(sel$control_fit_md5) || identical(sel$control_fit_md5,ctl_md5),
    control_job_matches=ctl_job,control_mcmc=paste(unlist(ctl$mcmc),collapse='/'),
    control_b0_dim=paste(dim(cb0),collapse='x'),control_posterior_finite=all(is.finite(cb0)),
    control_outputs=nrow(old_layout),layout_matches_control=layout_same,stringsAsFactors=FALSE)
})
table <- do.call(rbind,checks)
inputs$md5_after <- unname(tools::md5sum(inputs$input_file))
table$all_inputs_unchanged <- identical(inputs$md5_after,inputs$md5_before)
write.csv(inputs,file.path(pdir,'input-hashes.csv'),row.names=FALSE)
table$pass <- table$library_is_study & table$fit_hashes_current & table$intercept_prior_ok & table$input_md5_ok &
  table$rng_initial_is_input_state & table$posterior_finite & table$overwrite_refused &
  (is.na(table$resume_ok) | table$resume_ok) & table$control_md5_matches_record & table$control_job_matches &
  table$control_posterior_finite & table$layout_matches_control & table$all_inputs_unchanged
write.csv(table,file.path(pdir,'pilot-check.csv'),row.names=FALSE)
compact <- table;compact$fit_file <- sub(paste0('^',study,'/'),'',compact$fit_file)
compact$control_fit <- sub(paste0('^',archives,'/'),'',compact$control_fit)
results <- file.path(scripts,'results');dir.create(results,showWarnings=FALSE)
write.csv(compact,file.path(results,'pilot.csv'),row.names=FALSE)
write.csv(hardening,file.path(pdir,'hardening-check.csv'),row.names=FALSE)
write.csv(hardening,file.path(results,'pilot-hardening.csv'),row.names=FALSE)
print(t(table[,setdiff(names(table),c('fit_file','control_fit'))]))
print(hardening)
finish_gate(c(stats::setNames(table$pass,paste('pilot',table$phase)),stats::setNames(hardening$pass,hardening$check)),'pilot')
