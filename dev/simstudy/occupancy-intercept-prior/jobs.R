# Fitting-path helpers for the occupancy-intercept prior study. run.R sources
# this file and records its md5 in every fit, so it holds only what a fit
# needs; verification helpers live in verify-helpers.R, which run.R neither
# sources nor hashes. No package is loaded here. Archived inputs, settings,
# selections and fits are read, never written.

`%||%` <- function(x,y) if(is.null(x)) y else x

# The pr11 job specs record inputs under this root, which this study cannot
# read; Doug copied them to <archives>/intercept-prior-inputs (same layout).
ORIGINAL_INPUT_ROOT <- '/Users/douglasyu/Documents/Codex/2026-09-10/fam/work'
PHASES <- c('A1','A2','B')
PRIOR_SDS <- c(1,2,3,5)
SCHEDULES <- c('initial','long','pilot')
MAX_WORKERS <- 8L
PROTOCOL_MCMC <- list(initial=c(nchain=2,nburn=3000,niter=5000,nthin=1),
  long=c(nchain=4,nburn=6000,niter=12000,nthin=1),
  pilot=c(nchain=2,nburn=200,niter=200,nthin=1))
ARCHIVES <- c(pr11='pr11-current-20260927',targeted='spatial-targeted-20260927',
  amplitude='spatial-amplitude-20260928')
# Archived runners whose fitting expression is reused, with the settings file
# and hash field that recorded their fingerprint when the controls were fitted.
RUNNERS <- list(
  pr11=list(file='dev/simstudy/current-main-recheck/run.R',settings='initial/settings.rds',field='source_hashes'),
  amplitude=list(file='dev/simstudy/spatial-amplitude-prior/run.R',settings='inverse_gamma/long/settings.rds',field='script_hashes'),
  targeted=list(file='dev/simstudy/spatial-targeted-recheck/run.R',settings='initial/settings.rds',field='fit_hashes'))

# Repo-relative paths of the files run.R reads to fit, all hashed into each fit.
RUN_SCRIPT <- 'dev/simstudy/occupancy-intercept-prior/run.R'
JOBS_SCRIPT <- 'dev/simstudy/occupancy-intercept-prior/jobs.R'
FINGERPRINT_FILE <- 'dev/simstudy/occupancy-intercept-prior/results/library-fingerprint.csv'
PR11_HELPERS <- 'dev/simstudy/current-main-recheck/helpers.R'

archive_dir <- function(archives,name) file.path(archives,ARCHIVES[[name]])
runner_script_files <- function(phase) c(RUN_SCRIPT,JOBS_SCRIPT,FINGERPRINT_FILE,
  RUNNERS[[phase_runner_name(phase)]]$file,if(phase %in% c('A1','B')) PR11_HELPERS)
phase_runner_name <- function(phase) switch(phase,A1=,B='pr11',A2='amplitude',stop('Unknown phase: ',phase))

# Generic --name=value parser: named character list of the options given.
parse_options <- function(args,known,required=character()) {
  bad <- !grepl('^--[A-Za-z][A-Za-z-]*=',args)
  if(any(bad)) stop('Malformed argument: ',args[bad][1])
  opt <- sub('^--([^=]+)=.*$','\\1',args);value <- sub('^--[^=]+=','',args)
  unknown <- setdiff(opt,known)
  if(length(unknown)) stop('Unknown option: --',unknown[1])
  if(anyDuplicated(opt)) stop('Repeated option: --',opt[duplicated(opt)][1])
  missing_opt <- setdiff(required,opt)
  if(length(missing_opt)) stop('Missing --',missing_opt[1])
  as.list(stats::setNames(value,opt))
}

parse_workers <- function(text) {
  if(!grepl('^[0-9]+$',text) || !as.integer(text) %in% seq_len(MAX_WORKERS))
    stop('--workers must be an integer from 1 to ',MAX_WORKERS)
  as.integer(text)
}

# run.R interface. --archives defaults to the directory holding the normalised
# --study, and --inputs-root to <archives>/intercept-prior-inputs.
parse_run_args <- function(args) {
  o <- parse_options(args,
    known=c('repo','study','phase','sd','schedule','workers','keys','inputs-root','archives'),
    required=c('repo','study','phase','sd','schedule','workers'))
  if(!o$phase %in% PHASES) stop('--phase must be one of ',paste(PHASES,collapse=', '))
  if(!o$sd %in% as.character(PRIOR_SDS)) stop('--sd must be one of ',paste(PRIOR_SDS,collapse=', '))
  if(!o$schedule %in% SCHEDULES) stop('--schedule must be one of ',paste(SCHEDULES,collapse=', '))
  keys <- o$keys
  if(!is.null(keys)) {
    if(!nzchar(keys)) stop('Empty --keys')
    keys <- strsplit(keys,',',fixed=TRUE)[[1]]
    if(any(!nzchar(keys))) stop('Empty key in --keys')
    if(anyDuplicated(keys)) stop('Duplicated key in --keys: ',keys[duplicated(keys)][1])
  }
  archives <- o$archives %||% dirname(normalizePath(o$study,mustWork=FALSE))
  list(repo=o$repo,study=o$study,phase=o$phase,sd=as.numeric(o$sd),schedule=o$schedule,
    workers=parse_workers(o$workers),keys=keys,inputs_root=o$`inputs-root` %||% file.path(archives,'intercept-prior-inputs'),
    archives=archives)
}

phase_keys <- function(phase) switch(phase,
  A1=c(sprintf('jsdm-n0100-%02d',1:10),sprintf('jsdm-n0300-%02d',1:10)),
  B=c(sprintf('design-qnear_K6-sites300-%02d',1:10),sprintf('design-qfar_K6-sites300-%02d',1:10)),
  A2=sprintf('range%d-rep%02d-binary-k100',rep(c(4L,6L,8L),each=3L),1:3),
  stop('Unknown phase: ',phase))

remap_input <- function(path,to,from=ORIGINAL_INPUT_ROOT) {
  prefix <- paste0(from,'/')
  if(!startsWith(path,prefix) || nchar(path)<=nchar(prefix))
    stop('Input path is not under the original input root: ',path)
  file.path(to,substring(path,nchar(prefix)+1L))
}

checked_input <- function(path,md5) {
  if(!file.exists(path)) stop('Missing input: ',path)
  found <- unname(tools::md5sum(path))
  if(!identical(found,md5)) stop('Input md5 mismatch for ',path,': recorded ',md5,', found ',found)
  path
}

# One spec per community: the archived job record (exactly as the control runner
# used it), the path actually read and the md5 that path must have.
phase_jobs <- function(phase,archives,inputs_root,keys=NULL) {
  all_keys <- phase_keys(phase);keys <- keys %||% all_keys
  bad <- setdiff(keys,all_keys)
  if(length(bad)) stop('Key not in phase ',phase,': ',bad[1])
  if(phase=='A2') {
    targeted <- archive_dir(archives,'targeted')
    manifest <- readRDS(file.path(targeted,'input-manifest.rds'))
    manifest <- manifest[manifest$knots==100L,]
    rows <- manifest[match(keys,manifest$key),,drop=FALSE]
    if(!identical(rows$key,keys)) stop('Spatial input manifest lacks requested keys')
    specs <- lapply(seq_len(nrow(rows)),function(i) {
      row <- rows[i,,drop=FALSE]
      # Same record the spatial-amplitude runner built for its fits.
      job <- as.list(row[c('key','community','grid_index','replicate','arm','knots','input_md5')])
      job$input_file <- file.path(targeted,'inputs',paste0(row$community,'.rds'))
      if(!identical(job$arm,'binary') || !identical(job$knots,100L)) stop('Not a binary k100 job: ',job$key)
      list(key=job$key,phase=phase,job=job,input_file=job$input_file,input_md5=job$input_md5)
    })
  } else {
    jobs <- readRDS(file.path(archive_dir(archives,'pr11'),'initial/settings.rds'))$jobs
    family <- if(phase=='A1') 'jsdm' else 'design'
    specs <- lapply(keys,function(key) {
      job <- jobs[[key]]
      if(is.null(job) || !identical(job$key,key)) stop('No archived pr11 job for ',key)
      if(!identical(job$family,family) || (phase=='B' && !identical(job$arm,'sites300')))
        stop('Archived job ',key,' is not a phase ',phase,' job')
      list(key=key,phase=phase,job=job,input_file=remap_input(job$input_file,inputs_root),
        input_md5=job$input_md5)
    })
  }
  names(specs) <- keys
  specs
}

mcmc_sources <- function(phase,archives) {
  if(phase %in% c('A1','B')) {
    p <- archive_dir(archives,'pr11')
    return(list(initial=file.path(p,'initial/settings.rds'),long=file.path(p,'long/settings.rds')))
  }
  if(phase!='A2') stop('Unknown phase: ',phase)
  a <- archive_dir(archives,'amplitude');t <- archive_dir(archives,'targeted')
  # First file is authoritative (the spatial-amplitude runner); the second is the
  # spatial-targeted runner that fitted the initial and one longer control.
  list(initial=c(file.path(a,'half_cauchy/initial/settings.rds'),file.path(t,'initial/settings.rds')),
    long=c(file.path(a,'inverse_gamma/long/settings.rds'),file.path(t,'long/settings.rds')))
}

# Shorten a schedule while keeping each field's storage type.
short_mcmc <- function(base,nburn,niter) {
  base$nburn <- as.vector(nburn,mode=typeof(base$nburn))
  base$niter <- as.vector(niter,mode=typeof(base$niter))
  base
}

check_protocol_mcmc <- function(mcmc,schedule,source) {
  expected <- PROTOCOL_MCMC[[schedule]]
  if(!identical(names(mcmc),names(expected)) || !isTRUE(all(unlist(mcmc)==expected)))
    stop('Archived MCMC settings in ',source,' differ from the protocol ',schedule,' schedule')
  invisible(TRUE)
}

schedule_mcmc <- function(phase,schedule,archives) {
  if(!schedule %in% SCHEDULES) stop('Unknown schedule: ',schedule)
  from <- if(schedule=='pilot') 'initial' else schedule
  files <- mcmc_sources(phase,archives)[[from]]
  settings <- lapply(files,function(f) readRDS(f)$mcmc)
  for(i in seq_along(files)) check_protocol_mcmc(settings[[i]],from,files[i])
  mcmc <- settings[[1]]
  if(schedule=='pilot') {
    mcmc <- short_mcmc(mcmc,200,200);check_protocol_mcmc(mcmc,'pilot',files[1])
  }
  mcmc
}

fit_path <- function(study,phase,sd,schedule,key)
  file.path(study,'fits',phase,paste0('sd',format(sd)),schedule,paste0(key,'-fit.rds'))

atomic_save <- function(x,path) {
  temp <- paste0(path,'.tmp');saveRDS(x,temp);stopifnot(file.rename(temp,path))
}

# The fields that identify a fit and must be identical for it to be resumed.
# Only content and settings enter it: absolute paths (repo, study, archives,
# input copy, library) are saved beside it as information and never compared,
# so an unchanged fit resumes from another checkout. The job record drops its
# input path for the same reason; input identity is the recorded md5.
fit_metadata <- function(phase,key,sd,schedule,mcmc,job,input_md5,runner,fit_expression,
  phase_priors,source_revision,fit_hashes,script_hashes,threads_per_fit=1L)
  list(phase=phase,key=key,sd=sd,schedule=schedule,mcmc=mcmc,
    job_record=job[setdiff(names(job),'input_file')],input_md5=input_md5,
    listPriors_added=list(sigma_b0=sd),runner=runner,fit_expression=fit_expression,
    phase_priors=phase_priors,source_revision=source_revision,fit_hashes=fit_hashes,
    script_hashes=script_hashes,threads_per_fit=threads_per_fit)

# md5s named by path relative to root (plus an optional label), so the record
# depends on content and layout but not on where the tree lives.
hash_files <- function(files,root,label='') {
  root <- normalizePath(root,mustWork=TRUE);full <- normalizePath(files,mustWork=TRUE)
  prefix <- paste0(root,'/');inside <- startsWith(full,prefix)
  if(!all(inside)) stop('File outside ',root,': ',full[!inside][1])
  stats::setNames(unname(tools::md5sum(full)),paste0(label,substring(full,nchar(prefix)+1L)))
}

# Atomic per-key lock held while a key is checked, fitted and saved. A lock
# left by a killed run is reported, never broken automatically.
fit_lock_path <- function(dest) sub('-fit\\.rds$','.lock',dest)
acquire_fit_lock <- function(dest) {
  lock <- fit_lock_path(dest)
  if(identical(lock,dest)) stop('Not a fit path: ',dest)
  if(!dir.create(lock,showWarnings=FALSE)) {
    if(!dir.exists(lock)) stop('Could not create fit lock ',lock)
    owner <- file.path(lock,'owner')
    who <- if(file.exists(owner)) paste(readLines(owner,warn=FALSE),collapse=', ') else 'owner unknown'
    stop('Fit lock ',lock,' already exists (',who,'): another run.R process may be fitting this key, ',
      'or the lock is stale from an interrupted run. Check that no process holds it, then remove it by hand. ',
      'Refusing this key.')
  }
  writeLines(c(paste('pid',Sys.getpid()),paste('host',Sys.info()[['nodename']]),
    paste('since',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'))),file.path(lock,'owner'))
  lock
}
release_fit_lock <- function(lock) {unlink(lock,recursive=TRUE);invisible(!dir.exists(lock))}

quarantine_path <- function(dest,time=Sys.time()) {
  q <- sub('-fit\\.rds$',paste0('-fit.QUARANTINE-',format(time,'%Y%m%d%H%M%S'),'.rds'),dest)
  if(identical(q,dest)) stop('Not a fit path: ',dest)
  q
}

# Save a completed fit under its normal name only if every post-fit invariant
# holds; otherwise keep the evidence under a quarantine name and fail the job.
finalise_fit <- function(saved,dest,check) {
  problem <- tryCatch({check(saved);NULL},error=function(e) conditionMessage(e))
  if(is.null(problem)) {atomic_save(saved,dest);return(dest)}
  q <- quarantine_path(dest)
  if(file.exists(q)) q <- sub('\\.rds$',paste0('-',Sys.getpid(),'.rds'),q)
  saved$quarantine_reason <- problem
  atomic_save(saved,q)
  stop('Post-fit invariant failed for ',basename(dest),': ',problem,'; fit kept as ',q)
}

# Keys whose worker did not return exactly its own key (errors, try-error
# objects from mclapply, NULL from killed children, missing results).
job_failures <- function(status,keys) {
  ok <- vapply(seq_along(keys),function(i) {
    s <- if(i<=length(status)) status[[i]] else NULL
    is.character(s) && !inherits(s,'try-error') && length(s)==1L && identical(unname(s),keys[i])
  },logical(1))
  keys[!ok]
}
finish_jobs <- function(status,keys,label) {
  failed <- job_failures(status,keys)
  if(length(failed)) {
    cat('FAILED',length(failed),'of',length(keys),label,'jobs:',paste(failed,collapse=', '),'\n')
    quit(save='no',status=1)
  }
  cat('Completed',length(keys),label,'jobs.\n')
  invisible(TRUE)
}

# 'new' when absent, 'resume' when every metadata field is identical; else stop.
existing_fit_status <- function(path,metadata) {
  if(!file.exists(path)) return('new')
  saved <- readRDS(path)
  same <- vapply(names(metadata),function(f) f %in% names(saved) && identical(saved[[f]],metadata[[f]]),logical(1))
  if(!all(same)) stop('Refusing to overwrite existing fit ',path,'; differing fields: ',
    paste(names(metadata)[!same],collapse=', '))
  'resume'
}

is_missing_arg <- function(x,i) is.symbol(x[[i]]) && !nzchar(as.character(x[[i]]))

# Replace every occurrence of one language object by another, leaving formals,
# empty arguments and constants untouched.
replace_in_call <- function(expr,from,to) {
  if(identical(expr,from)) return(to)
  if(!is.call(expr)) return(expr)
  parts <- as.list(expr)
  # parts[i] <- list(...) keeps NULL arguments such as collCovariates=NULL.
  for(i in seq_along(parts)) if(!is_missing_arg(parts,i)) parts[i] <- list(replace_in_call(parts[[i]],from,to))
  as.call(parts)
}

recorded_hash <- function(hashes,suffix) {
  hit <- hashes[endsWith(names(hashes),paste0('/',suffix))]
  if(length(hit)!=1L) stop('Expected one recorded hash for ',suffix,', found ',length(hit))
  unname(hit)
}

# The unique `fit <- withCallingHandlers(...)` statement of an archived runner,
# after checking that the runner is byte-identical to the recorded one.
archived_fit_expression <- function(file,md5) {
  found <- unname(tools::md5sum(file))
  if(!identical(found,md5)) stop('Archived runner md5 mismatch for ',file,': recorded ',md5,', found ',found)
  hits <- list()
  visit <- function(e) {
    if(!is.call(e)) return(invisible())
    if(identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('fit')) && is.call(e[[3]]) &&
       identical(e[[3]][[1]],as.name('withCallingHandlers'))) {
      hits[[length(hits)+1L]] <<- e;return(invisible())
    }
    for(i in seq_along(e)) if(!is_missing_arg(e,i)) visit(e[[i]])
  }
  for(e in parse(file,keep.source=FALSE)) visit(e)
  if(length(hits)!=1L) stop('Expected one fitting expression in ',file,', found ',length(hits))
  hits[[1]]
}

archived_runner <- function(name,repo,archives) {
  spec <- RUNNERS[[name]]
  if(is.null(spec)) stop('Unknown archived runner: ',name)
  settings_path <- file.path(archive_dir(archives,name),spec$settings)
  md5 <- recorded_hash(readRDS(settings_path)[[spec$field]],spec$file)
  path <- file.path(repo,spec$file)
  # file and settings are relative (repo, archives); path and settings_path absolute.
  list(name=name,file=spec$file,md5=md5,settings=file.path(ARCHIVES[[name]],spec$settings),
    path=path,settings_path=settings_path,expr=archived_fit_expression(path,md5))
}

# The archived expression with its fitting function replaced by `.fitter`,
# which differs from the control call only by adding listPriors$sigma_b0.
with_intercept_prior <- function(expr) {
  out <- replace_in_call(expr,quote(occJSDM::runOccJSDM),as.name('.fitter'))
  out <- replace_in_call(out,as.name('fun'),as.name('.fitter'))
  if(identical(out,expr)) stop('The archived expression contains no fitting call')
  if(any(grepl('runOccJSDM',deparse(out),fixed=TRUE))) stop('An archived fitting call was not replaced')
  out
}

make_fitter <- function(sd,target=NULL) {
  force(sd);force(target)
  function(...,listPriors=list()) {
    if(!is.null(sd)) listPriors$sigma_b0 <- sd
    fitter <- if(is.null(target)) occJSDM::runOccJSDM else target
    fitter(...,listPriors=listPriors)
  }
}

# Evaluate a fitting statement exactly as the archived runner did: saved input
# RNG state assigned immediately before, warnings collected by its handler.
# `expr` is only ever a statement parsed from a tracked research runner whose
# md5 matched the fingerprint recorded with the control fits (see
# archived_fit_expression), never text from inputs or other untrusted data.
evaluate_fit <- function(expr,input,job,mcmc,bindings=list()) {
  env <- new.env(parent=globalenv())
  env$input <- input;env$job <- job;env$mcmc <- mcmc;env$warnings <- character()
  for(name in names(bindings)) assign(name,bindings[[name]],envir=env)
  assign('.Random.seed',input$fit_rng,envir=.GlobalEnv)
  started <- Sys.time()
  eval(expr,env)
  finished <- Sys.time()
  list(fit=env$fit,warnings=env$warnings,started=started,finished=finished,
    rng_initial=input$fit_rng,rng_final=get('.Random.seed',envir=.GlobalEnv))
}

check_fit_prior <- function(fit,sd) {
  expected <- list(mean=0,sd=sd %||% 1)
  if(!identical(fit$infos$intercept_prior,expected))
    stop('Fitted object records intercept prior ',paste(deparse(fit$infos$intercept_prior),collapse=''),
      ', expected ',paste(deparse(expected),collapse=''))
  TRUE
}

# Everything that determines fitting behaviour: exported source and the
# installed package (compiled library, metadata and byte-compiled R code).
production_files <- function(study,pkg,must_exist=TRUE) {
  files <- c(file.path(study,'source',c('DESCRIPTION','NAMESPACE')),
    list.files(file.path(study,'source/R'),pattern='\\.R$',full.names=TRUE),
    list.files(file.path(study,'source/src'),pattern='\\.(cpp|h)$|^Makevars',full.names=TRUE),
    file.path(pkg,'libs/occJSDM.so'),file.path(pkg,c('DESCRIPTION','NAMESPACE')),
    list.files(file.path(pkg,'R'),full.names=TRUE))
  if(must_exist) stopifnot(length(files)>10L,all(file.exists(files)))
  files
}

library_fingerprint_files <- function(study)
  c(production_files(study,file.path(study,'library/occJSDM'),must_exist=FALSE),file.path(study,'source-revision.txt'))

# The installed build must be exactly the tested one: same revision and the
# same md5 for every exported source file and installed library file, with no
# file missing or added. Returns the revision.
check_library_fingerprint <- function(study,csv) {
  if(!file.exists(csv)) stop('Missing library fingerprint: ',csv)
  fp <- utils::read.csv(csv,colClasses='character')
  if(!identical(names(fp),c('file','md5','revision')) || !nrow(fp) || anyDuplicated(fp$file))
    stop('Malformed library fingerprint: ',csv)
  revision <- unique(fp$revision)
  if(length(revision)!=1L || !grepl('^[0-9a-f]{40}$',revision)) stop('Malformed library fingerprint revision in ',csv)
  rfile <- file.path(study,'source-revision.txt')
  actual_revision <- if(file.exists(rfile)) readLines(rfile,warn=FALSE) else character()
  if(!identical(actual_revision,revision)) stop('Library fingerprint mismatch in ',study,': source revision ',
    paste(actual_revision,collapse=' '),' is not the recorded revision ',revision)
  files <- library_fingerprint_files(study);files <- files[file.exists(files)]
  actual <- hash_files(files,study);expected <- stats::setNames(fp$md5,fp$file)
  missing_files <- setdiff(names(expected),names(actual));unexpected <- setdiff(names(actual),names(expected))
  common <- intersect(names(expected),names(actual));changed <- common[actual[common]!=expected[common]]
  if(length(c(missing_files,unexpected,changed))) stop('Library fingerprint mismatch in ',study,':',
    if(length(changed)) paste0(' changed ',paste(changed,collapse=', ')),
    if(length(missing_files)) paste0(' missing ',paste(missing_files,collapse=', ')),
    if(length(unexpected)) paste0(' unexpected ',paste(unexpected,collapse=', ')))
  revision
}

phase_priors <- function(phase,archives) {
  if(phase!='A2') return(NULL)
  priors <- readRDS(file.path(archive_dir(archives,'amplitude'),'inverse_gamma/long/settings.rds'))$priors
  if(!identical(priors,list())) stop('Inverse-gamma control priors are not the defaults')
  priors
}
