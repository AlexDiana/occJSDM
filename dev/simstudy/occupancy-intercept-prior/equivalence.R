#!/usr/bin/env Rscript
# Cross-installation default equivalence of the study library against every
# installed library that produced a reused control fit.
#
#   Rscript equivalence.R --repo=REPO --study=STUDY --mode=launch [--workers=4]
#     [--archives=DIR] [--inputs-root=DIR]
#
# launch runs each fit below as its own R process (so each loads exactly one
# library), then compares. A control run evaluates the archived control
# runner's unmodified fitting statement with the control library; a study run
# evaluates run.R's statement with STUDY/library, default priors (no sigma_b0)
# or explicit sigma_b0 = 1. Same saved input and RNG state; short schedule
# (2 chains, 50 burn-in, 50 retained). The two full-length runs refit the saved
# A1 and B initial controls with STUDY/library at default priors.
args <- commandArgs(trailingOnly=TRUE)
repo_arg <- sub('^--repo=','',grep('^--repo=',args,value=TRUE))
if(length(repo_arg)!=1L) stop('Missing or repeated --repo')
scripts <- file.path(normalizePath(repo_arg),'dev/simstudy/occupancy-intercept-prior')
source(file.path(scripts,'jobs.R'))
o <- parse_options(args,known=c('repo','study','mode','run','archives','inputs-root','workers'),
  required=c('repo','study','mode'))
repo <- normalizePath(o$repo);study <- normalizePath(o$study)
archives <- normalizePath(o$archives %||% dirname(study))
inputs_root <- normalizePath(o$`inputs-root` %||% file.path(archives,'intercept-prior-inputs'))
eqdir <- file.path(study,'equivalence');dir.create(eqdir,showWarnings=FALSE)
TOLERANCE <- 1e-12

CHECKS <- data.frame(check=c('A1-pr11','B-pr11','A2-targeted','A2-amplitude'),phase=c('A1','B','A2','A2'),
  key=c('jsdm-n0100-01','design-qnear_K6-sites300-01','range6-rep01-binary-k100','range6-rep01-binary-k100'),
  control=c('pr11','pr11','targeted','amplitude'),stringsAsFactors=FALSE)
RUNS <- rbind(
  data.frame(run=paste0('control-',CHECKS$check),phase=CHECKS$phase,key=CHECKS$key,control=CHECKS$control,
    variant='control',stringsAsFactors=FALSE),
  do.call(rbind,lapply(c('A1','B','A2'),function(p) data.frame(run=paste0('study-',p,c('-default','-explicit1')),
    phase=p,key=CHECKS$key[match(p,CHECKS$phase)],control=NA_character_,variant=c('default','explicit1'),
    stringsAsFactors=FALSE))),
  data.frame(run=c('study-A1-full-default','study-B-full-default'),phase=c('A1','B'),
    key=c('jsdm-n0100-01','design-qnear_K6-sites300-01'),control=NA_character_,variant='full-default',
    stringsAsFactors=FALSE))
run_file <- function(run) file.path(eqdir,paste0(run,'.rds'))

fit_one <- function(r) {
  control <- r$variant=='control'
  lib <- if(control) file.path(archive_dir(archives,r$control),'library') else file.path(study,'library')
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  .libPaths(c(lib,.libPaths()))
  suppressPackageStartupMessages(library(occJSDM))
  pkg <- normalizePath(find.package('occJSDM'))
  stopifnot(pkg==normalizePath(file.path(lib,'occJSDM')))
  RcppParallel::setThreadOptions(numThreads=1)
  if(control) {
    runner <- archived_runner(r$control,repo,archives);expr <- runner$expr
    # The spatial-amplitude runner calls `fun` (runOccJSDM without start
    # overrides) with its inverse-gamma `priors`, list().
    bindings <- if(r$control=='amplitude') list(fun=occJSDM::runOccJSDM,priors=phase_priors('A2',archives)) else list()
  } else {
    runner <- archived_runner(phase_runner_name(r$phase),repo,archives);expr <- with_intercept_prior(runner$expr)
    bindings <- list(.fitter=make_fitter(if(r$variant=='explicit1') 1 else NULL),priors=phase_priors(r$phase,archives))
  }
  mcmc <- schedule_mcmc(r$phase,'initial',archives)
  if(r$variant!='full-default') mcmc <- short_mcmc(mcmc,50,50)
  spec <- phase_jobs(r$phase,archives,inputs_root,r$key)[[1]]
  checked_input(spec$input_file,spec$input_md5)
  input <- readRDS(spec$input_file)
  result <- evaluate_fit(expr,input,spec$job,mcmc,bindings=bindings)
  revision <- readLines(file.path(if(control) archive_dir(archives,r$control) else study,'source-revision.txt'))
  saved <- c(as.list(r),result,list(library=pkg,revision=revision,library_hashes=tools::md5sum(installed_library_files(pkg)),
    runner=runner[c('name','file','md5')],expression=deparse(expr,width.cutoff=500L),mcmc=mcmc,
    input_file=spec$input_file,input_md5_before=spec$input_md5,input_md5_after=unname(tools::md5sum(spec$input_file)),
    elapsed_seconds=as.numeric(difftime(result$finished,result$started,units='secs'))))
  atomic_save(saved,run_file(r$run))
  cat(r$run,'saved;',round(saved$elapsed_seconds,1),'seconds;',length(saved$warnings),'warnings\n')
}

compare <- function() {
  load <- function(run) readRDS(run_file(run))
  # Design matrices and every other non-posterior component must be identical.
  other_same <- function(a,b) {
    rest <- setdiff(union(names(a),names(b)),c('results_output','infos'))
    identical(a[rest],b[rest])
  }
  # Metadata fields in fit$infos that differ, marked as added in the study fit
  # or changed. Only the documented prior records may be additions.
  infos_diff <- function(a,b) {
    fields <- union(names(a$infos),names(b$infos))
    differ <- fields[!vapply(fields,function(f) identical(a$infos[[f]],b$infos[[f]]),logical(1))]
    paste0(differ,ifelse(differ %in% names(a$infos),'(changed)','(added)'),collapse=';')
  }
  allowed_infos <- function(x) !nzchar(x) || all(strsplit(x,';',fixed=TRUE)[[1]] %in%
    c('intercept_prior(added)','spatial_sd_prior(added)'))
  rows <- lapply(seq_len(nrow(CHECKS)),function(i) {
    ck <- CHECKS[i,];ctl <- load(paste0('control-',ck$check))
    def <- load(paste0('study-',ck$phase,'-default'));one <- load(paste0('study-',ck$phase,'-explicit1'))
    d <- max_abs_difference(ctl$fit$results_output,def$fit$results_output)
    data.frame(check=ck$check,phase=ck$phase,key=ck$key,comparison='control library vs study library default (short)',
      control_library=ctl$library,control_revision=ctl$revision,control_so_md5=unname(ctl$library_hashes[1]),
      study_revision=def$revision,study_so_md5=unname(def$library_hashes[1]),
      mcmc=paste(unlist(def$mcmc),collapse='/'),numeric_values=count_numeric(def$fit$results_output),
      max_abs_diff=d,bitwise_identical=identical(ctl$fit$results_output,def$fit$results_output),
      other_components_identical=other_same(ctl$fit,def$fit),infos_differences=infos_diff(ctl$fit,def$fit),
      warnings_control=length(ctl$warnings),warnings_study=length(def$warnings),
      warnings_identical=identical(ctl$warnings,def$warnings),rng_final_identical=identical(ctl$rng_final,def$rng_final),
      control_prior_recorded=!is.null(ctl$fit$infos$intercept_prior),
      study_default_prior_sd1=identical(def$fit$infos$intercept_prior,list(mean=0,sd=1)),
      explicit1_identical_to_default=identical(one$fit,def$fit) && identical(one$warnings,def$warnings) &&
        identical(one$rng_final,def$rng_final),
      explicit1_prior_sd1=identical(one$fit$infos$intercept_prior,list(mean=0,sd=1)),
      inputs_unchanged=all(c(ctl$input_md5_after,def$input_md5_after,one$input_md5_after)==def$input_md5_before),
      stringsAsFactors=FALSE)
  })
  full <- lapply(c('A1','B'),function(p) {
    f <- load(paste0('study-',p,'-full-default'))
    sel <- control_schedules(p,archives);sel <- sel[sel$key==f$key,]
    stopifnot(sel$schedule=='initial')
    ctl <- readRDS(sel$control_fit)
    stopifnot(identical(ctl$mcmc,f$mcmc))
    d <- max_abs_difference(ctl$fit$results_output,f$fit$results_output)
    data.frame(check=paste0(p,'-pr11-saved-control'),phase=p,key=f$key,
      comparison='saved control fit vs study library default (full initial schedule)',
      control_library=file.path(archive_dir(archives,'pr11'),'library'),
      control_revision=readLines(file.path(archive_dir(archives,'pr11'),'source-revision.txt')),
      control_so_md5=unname(ctl$source_hashes[endsWith(names(ctl$source_hashes),'/libs/occJSDM.so')]),
      study_revision=f$revision,study_so_md5=unname(f$library_hashes[1]),
      mcmc=paste(unlist(f$mcmc),collapse='/'),numeric_values=count_numeric(f$fit$results_output),
      max_abs_diff=d,bitwise_identical=identical(ctl$fit$results_output,f$fit$results_output),
      other_components_identical=other_same(ctl$fit,f$fit),infos_differences=infos_diff(ctl$fit,f$fit),
      warnings_control=length(ctl$warnings),warnings_study=length(f$warnings),
      warnings_identical=identical(ctl$warnings,f$warnings),rng_final_identical=NA,
      control_prior_recorded=!is.null(ctl$fit$infos$intercept_prior),
      study_default_prior_sd1=identical(f$fit$infos$intercept_prior,list(mean=0,sd=1)),
      explicit1_identical_to_default=NA,explicit1_prior_sd1=NA,
      inputs_unchanged=f$input_md5_after==f$input_md5_before,stringsAsFactors=FALSE)
  })
  table <- do.call(rbind,c(rows,full))
  table$pass <- table$max_abs_diff<=TOLERANCE & table$other_components_identical &
    vapply(table$infos_differences,allowed_infos,logical(1)) & table$warnings_identical &
    (is.na(table$rng_final_identical) | table$rng_final_identical) & !table$control_prior_recorded &
    table$study_default_prior_sd1 & (is.na(table$explicit1_identical_to_default) | table$explicit1_identical_to_default) &
    (is.na(table$explicit1_prior_sd1) | table$explicit1_prior_sd1) & table$inputs_unchanged
  write.csv(table,file.path(eqdir,'equivalence.csv'),row.names=FALSE)
  compact <- table[,c('check','phase','key','comparison','control_revision','control_so_md5','study_revision',
    'study_so_md5','mcmc','numeric_values','max_abs_diff','bitwise_identical','other_components_identical',
    'infos_differences','warnings_identical','rng_final_identical','explicit1_identical_to_default','inputs_unchanged','pass')]
  results <- file.path(scripts,'results');dir.create(results,showWarnings=FALSE)
  write.csv(compact,file.path(results,'equivalence.csv'),row.names=FALSE)
  print(table[,c('check','max_abs_diff','bitwise_identical','warnings_identical','rng_final_identical',
    'explicit1_identical_to_default','pass')])
  cat(if(all(table$pass)) 'All equivalence checks pass.\n' else 'EQUIVALENCE FAILED for: ',
    paste(table$check[!table$pass],collapse=', '),'\n')
  invisible(table)
}

if(o$mode=='fit') {
  r <- RUNS[RUNS$run==o$run,];if(nrow(r)!=1L) stop('Unknown --run: ',o$run)
  fit_one(r)
} else if(o$mode=='compare') {
  compare()
} else if(o$mode=='launch') {
  source(file.path(repo,'dev/simstudy/spatial-amplitude-prior/execution.R'))
  workers <- parse_workers(o$workers %||% '4')
  pending <- RUNS$run[!file.exists(run_file(RUNS$run))]
  cat(format(Sys.time()),'launching',length(pending),'equivalence fits\n');flush.console()
  extra <- c(paste0('--repo=',repo),paste0('--study=',study),paste0('--archives=',archives),
    paste0('--inputs-root=',inputs_root))
  status <- parallel::mclapply(pending,function(run) {
    s <- run_logged_r(file.path(scripts,'equivalence.R'),c(extra,'--mode=fit',paste0('--run=',run)),
      file.path(eqdir,paste0(run,'.log')))
    cat(format(Sys.time()),run,'exit',s,'\n');s
  },mc.cores=workers,mc.preschedule=FALSE)
  if(!all(unlist(status)==0L)) stop('An equivalence fit failed; see logs in ',eqdir)
  compare()
} else stop('Unknown --mode: ',o$mode)
