#!/usr/bin/env Rscript
# Mirror-frequency check (README.md in this directory): one repeat of the pr11
# long fit of design-qfar_K6-sites300-05 (community 5) per process. The four
# chains run one after another inside this process from one seed, as the
# package runs them on main.
#
#   Rscript run.R --repo=REPO --study=ARCHIVE --inputs-root=DIR
#     --stage=a|b|pilot --repeat=K [--archives=DIR]
#
# --archives defaults to the directory holding ARCHIVE; its
# pr11-current-20260927 archive supplies the job record, the long schedule and
# the md5s recorded for the control runner.
#
# Each repeat evaluates the pr11 control runner's own fitting statement (parsed
# from dev/simstudy/current-main-recheck/run.R after its md5 is checked against
# the one recorded with the pr11 fits), with the fitting function replaced by
# `.fitter`, which calls occJSDM::runOccJSDM with the statement's arguments
# unchanged and records the listPriors passed. Stage a uses the pr11 long
# schedule read from the archive (4 chains, 6,000 burn-in, 12,000 retained,
# no thinning); stage b the same with 24,000 burn-in; the pilot 4 chains of
# 50 burn-in and 50 retained. Seed: set.seed(SEED_BASE[stage] + repeat) with
# the kind of the saved input RNG state, immediately before the statement.
# Output: ARCHIVE/fits/stage-<stage>/repeat-KK-fit.rds under an atomic
# per-fit lock (repeat-KK.lock). A fit already saved with identical metadata
# is left alone and reported as resumed; one with any differing field is never
# overwritten. The installed library must match results/library-fingerprint.csv
# before it is loaded. Exit status 0 on success or resume, 1 on any failure.
# Launch batches with launch.R; analyse a stage with analyse.R.

`%||%` <- function(x,y) if(is.null(x)) y else x

STUDY_DIR <- local({
  here <- NULL
  for(i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if(!is.null(f)) {here <- dirname(normalizePath(f));break}
  }
  if(is.null(here)) {
    a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
    if(length(a)) here <- dirname(normalizePath(sub('^--file=','',a[1])))
  }
  here %||% getwd()
})
SCRIPT_REPO <- normalizePath(file.path(STUDY_DIR,'../../..'))

STUDY_REL <- 'dev/simstudy/mirror-frequency'
RUN_SCRIPT <- file.path(STUDY_REL,'run.R')
FINGERPRINT_FILE <- file.path(STUDY_REL,'results/library-fingerprint.csv')
JOBS_SCRIPT <- 'dev/simstudy/occupancy-intercept-prior/jobs.R'
PR11_RUNNER <- 'dev/simstudy/current-main-recheck/run.R'
PR11_HELPERS <- 'dev/simstudy/current-main-recheck/helpers.R'
PR11_ARCHIVE <- 'pr11-current-20260927'
PR11_SETTINGS <- file.path(PR11_ARCHIVE,c('initial/settings.rds','long/settings.rds'))
# Reused files, frozen at these md5s (jobs.R as committed in c28eedc; helpers.R
# as recorded with the pr11 fits), as the convergence-flag diagnosis froze them.
REUSED_MD5 <- c(jobs='d6dc59ff47a25de08fad60b8bde827d1',helpers='50d5c6bac89f9739509a86c59380dfbf')

KEY <- 'design-qfar_K6-sites300-05'
TARGET_SPECIES <- 6L
TARGET_SPECIES_NAME <- 'OTU_6'
STAGES <- c('a','b','pilot')
REPEATS <- list(a=1:25,b=1:25,pilot=1L)
SEED_BASE <- c(a=20261100L,b=20261200L,pilot=20261050L)
SEED_RULE <- 'set.seed(SEED_BASE[stage] + repeat, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection") immediately before the fitting statement; the four chains then run one after another in this process from that one stream'
RNG_KIND <- c(kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
# runOccJSDM runs nburn + niter * nthin iterations and keeps niter draws, so
# niter counts retained draws.
PR11_LONG <- c(nchain=4,nburn=6000,niter=12000,nthin=1)
STAGE_B_NBURN <- 24000
PILOT_MCMC <- list(nchain=4,nburn=50,niter=50,nthin=1)
RESUMED_TEXT <- 'already complete with identical settings; resumed'

# ---- Reused code ------------------------------------------------------------------

check_reused_md5 <- function(repo) {
  files <- file.path(repo,c(jobs=JOBS_SCRIPT,helpers=PR11_HELPERS))
  found <- unname(tools::md5sum(files))
  bad <- names(REUSED_MD5)[found!=REUSED_MD5]
  if(length(bad)) stop('Reused file changed: ',paste(files[match(bad,names(REUSED_MD5))],collapse=', '))
  invisible(TRUE)
}
check_reused_md5(SCRIPT_REPO)
PR11 <- new.env(parent=globalenv())
sys.source(file.path(SCRIPT_REPO,PR11_HELPERS),envir=PR11)
PR11$extract_functions(file.path(SCRIPT_REPO,JOBS_SCRIPT),
  c('ORIGINAL_INPUT_ROOT','remap_input','checked_input','parse_options','hash_files','atomic_save',
    'fit_lock_path','acquire_fit_lock','release_fit_lock','quarantine_path','finalise_fit',
    'existing_fit_status','is_missing_arg','replace_in_call','recorded_hash','archived_fit_expression',
    'with_intercept_prior','production_files','library_fingerprint_files','check_library_fingerprint'),
  environment())

# ---- Stages ----------------------------------------------------------------------

repeat_seed <- function(stage,k) unname(SEED_BASE[[stage]]+as.integer(k))

# The pr11 long schedule as the archive records it, checked against PR11_LONG.
pr11_long_mcmc <- function(archives) {
  m <- readRDS(file.path(archives,PR11_SETTINGS[2]))$mcmc
  if(!identical(names(m),names(PR11_LONG)) || !all(unlist(m)==PR11_LONG))
    stop('The pr11 long schedule is not ',paste(names(PR11_LONG),PR11_LONG,collapse=' '))
  m
}

stage_mcmc <- function(stage,archives) {
  if(stage=='pilot') return(PILOT_MCMC)
  m <- pr11_long_mcmc(archives)
  if(stage=='b') m$nburn <- as.vector(STAGE_B_NBURN,mode=typeof(m$nburn))
  m
}

run_spec <- function(stage,k,archives) {
  if(length(stage)!=1L || !stage %in% STAGES) stop('--stage must be one of ',paste(STAGES,collapse=', '))
  k <- as.integer(k)
  if(length(k)!=1L || is.na(k) || !k %in% REPEATS[[stage]])
    stop('--repeat for stage ',stage,' must be one of ',min(REPEATS[[stage]]),' to ',max(REPEATS[[stage]]))
  list(stage=stage,repeat_index=k,seed=repeat_seed(stage,k),mcmc=stage_mcmc(stage,archives))
}

parse_run_args <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','stage','repeat'),
    required=c('repo','study','inputs-root','stage','repeat'))
  if(!grepl('^[0-9]+$',o$`repeat`)) stop('--repeat must be a whole number')
  archives <- o$archives %||% dirname(normalizePath(o$study,mustWork=FALSE))
  spec <- run_spec(o$stage,as.integer(o$`repeat`),archives)
  list(repo=o$repo,study=o$study,inputs_root=o$`inputs-root`,archives=archives,spec=spec)
}

fit_file <- function(study,spec)
  file.path(study,'fits',paste0('stage-',spec$stage),sprintf('repeat-%02d-fit.rds',spec$repeat_index))

iterations_run <- function(mcmc) mcmc$nchain*(mcmc$nburn+mcmc$niter*mcmc$nthin)

# ---- Archived job and runner -----------------------------------------------------------

# The community-5 job record as the pr11 runner used it (identical in its
# initial and long settings), with the input path remapped to the copy.
pr11_job <- function(archives,inputs_root) {
  jobs <- lapply(file.path(archives,PR11_SETTINGS),function(f) readRDS(f)$jobs[[KEY]])
  job <- jobs[[1]]
  if(is.null(job) || !identical(job$key,KEY) || !identical(jobs[[1]],jobs[[2]]))
    stop('The pr11 settings do not hold one identical job record for ',KEY)
  if(!identical(job$family,'design') || !identical(job$arm,'sites300')) stop('Unexpected pr11 job record for ',KEY)
  list(job=job,input_file=remap_input(job$input_file,inputs_root),input_md5=job$input_md5)
}

# The pr11 runner's fitting statement, after checking run.R and helpers.R
# against the md5s recorded in both pr11 settings files.
pr11_runner <- function(repo,archives) {
  hashes <- lapply(file.path(archives,PR11_SETTINGS),function(f) readRDS(f)$source_hashes)
  one <- function(suffix) {
    h <- unique(vapply(hashes,recorded_hash,character(1),suffix=suffix))
    if(length(h)!=1L) stop('pr11 settings record different md5s for ',suffix)
    h
  }
  runner_md5 <- one(PR11_RUNNER);helpers_md5 <- one(PR11_HELPERS)
  if(!identical(helpers_md5,REUSED_MD5[['helpers']])) stop('pr11 helpers.R md5 differs from the frozen one')
  found <- unname(tools::md5sum(file.path(repo,PR11_HELPERS)))
  if(!identical(found,helpers_md5)) stop('pr11 helpers.R changed: recorded ',helpers_md5,', found ',found)
  list(file=PR11_RUNNER,md5=runner_md5,helpers_md5=helpers_md5,settings=PR11_SETTINGS,
    expr=archived_fit_expression(file.path(repo,PR11_RUNNER),runner_md5))
}

# ---- Fitting --------------------------------------------------------------------------

# The function bound to `.fitter` in the fitting statement: the package's
# fitter with the statement's arguments unchanged.
make_fitter <- function(recorder,target=occJSDM::runOccJSDM) {
  force(target)
  function(...,listPriors=list()) {
    recorder$listPriors <- listPriors
    target(...,listPriors=listPriors)
  }
}

# Evaluate the fitting statement as the pr11 runner did, after set.seed(seed)
# with RNG_KIND. `expr` is only ever the statement parsed from the md5-checked
# pr11 runner.
evaluate_seeded <- function(expr,input,job,mcmc,bindings=list(),seed) {
  env <- new.env(parent=globalenv())
  env$input <- input;env$job <- job;env$mcmc <- mcmc;env$warnings <- character()
  for(name in names(bindings)) assign(name,bindings[[name]],envir=env)
  if(!is.numeric(seed) || length(seed)!=1L || !is.finite(seed)) stop('seed must be one number')
  set.seed(seed,kind=RNG_KIND[['kind']],normal.kind=RNG_KIND[['normal.kind']],sample.kind=RNG_KIND[['sample.kind']])
  rng_initial <- get('.Random.seed',envir=.GlobalEnv)
  if(!identical(rng_initial[1],input$fit_rng[1])) stop('RNG kind differs from the saved input state')
  started <- Sys.time()
  eval(expr,env)
  finished <- Sys.time()
  list(fit=env$fit,warnings=env$warnings,started=started,finished=finished,
    rng_initial=rng_initial,rng_final=get('.Random.seed',envir=.GlobalEnv))
}

# The fields that identify a fit; all must be identical for a saved fit to be
# accepted as this one. Absolute paths are saved beside them, never compared.
fit_metadata <- function(spec,job,input_md5,fit_expression,runner,source_revision,fit_hashes,script_hashes)
  list(key=KEY,stage=spec$stage,repeat_index=spec$repeat_index,seed=spec$seed,seed_rule=SEED_RULE,
    rng_kind=RNG_KIND,mcmc=spec$mcmc,chains_in_process=as.integer(spec$mcmc$nchain),
    job_record=job[setdiff(names(job),'input_file')],input_md5=input_md5,fit_expression=fit_expression,
    runner=runner[c('file','md5','helpers_md5','settings')],source_revision=source_revision,
    fit_hashes=fit_hashes,script_hashes=script_hashes,threads_per_fit=1L)

fit_script_files <- function() c(RUN_SCRIPT,FINGERPRINT_FILE,JOBS_SCRIPT,PR11_RUNNER,PR11_HELPERS)

load_study_library <- function(study) {
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  .libPaths(c(file.path(study,'library'),.libPaths()))
  suppressPackageStartupMessages(library(occJSDM))
  pkg <- normalizePath(find.package('occJSDM'))
  if(!identical(pkg,normalizePath(file.path(study,'library/occJSDM')))) stop('occJSDM loaded from ',pkg,', not the study library')
  RcppParallel::setThreadOptions(numThreads=1)
  pkg
}

# Post-fit invariants; a failure keeps the fit under a quarantine name.
check_saved_fit <- function(s,spec,prepared,job_priors,pkg) {
  fit <- s$fit;ro <- fit$results_output;mcmc <- spec$mcmc
  PR11$check_current_fit(fit,prepared$input,prepared$job)
  dd <- c(10L,as.integer(mcmc$niter),as.integer(mcmc$nchain))
  if(!identical(dim(ro$jsdm_output$B0_output),dd)) stop('Unexpected B0_output dimensions')
  if(!identical(dim(ro$theta0_output),dd)) stop('Unexpected theta0_output dimensions')
  if(!identical(s$listPriors_used,job_priors)) stop('listPriors passed to runOccJSDM differ from the job priors')
  species <- fit$infos$speciesNames
  if(!identical(species[TARGET_SPECIES],TARGET_SPECIES_NAME)) stop('Species ',TARGET_SPECIES,' is not ',TARGET_SPECIES_NAME)
  varies <- apply(ro$theta0_output,1,function(v) length(unique(v))>1L)
  if(!all(varies)) stop('A theta0 is constant')
  if(!identical(s$input_md5_after,prepared$input_md5)) stop('Input changed during fitting: ',prepared$input_file)
  if(!identical(s$loaded_library,pkg)) stop('Loaded library changed during fitting')
  invisible(TRUE)
}

prepare_input <- function(source) {
  checked_input(source$input_file,source$input_md5)
  list(input=readRDS(source$input_file),job=source$job,input_file=source$input_file,input_md5=source$input_md5)
}

run_main <- function(args) {
  opt <- parse_run_args(args);spec <- opt$spec
  repo <- normalizePath(opt$repo,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo (',repo,') is not the checkout holding this run.R (',SCRIPT_REPO,')')
  study <- normalizePath(opt$study,mustWork=TRUE);archives <- normalizePath(opt$archives,mustWork=TRUE)
  inputs_root <- normalizePath(opt$inputs_root,mustWork=TRUE)
  tag <- paste0(KEY,' stage ',spec$stage,' repeat ',spec$repeat_index)
  # Refuse any build other than the fingerprinted one before loading it.
  source_revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
  pkg <- load_study_library(study)
  fit_hashes <- hash_files(production_files(study,pkg),study)
  runner <- pr11_runner(repo,archives)
  expr <- with_intercept_prior(runner$expr);fit_expression <- deparse(expr,width.cutoff=500L)
  source <- pr11_job(archives,inputs_root)
  script_hashes <- c(hash_files(file.path(repo,fit_script_files()),repo),
    hash_files(file.path(archives,PR11_SETTINGS),archives,'archives/'))
  dest <- fit_file(study,spec);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  lock <- acquire_fit_lock(dest);on.exit(release_fit_lock(lock),add=TRUE)
  prepared <- prepare_input(source)
  metadata <- fit_metadata(spec,prepared$job,prepared$input_md5,fit_expression,runner,source_revision,fit_hashes,script_hashes)
  if(existing_fit_status(dest,metadata)=='resume') {
    cat(format(Sys.time()),tag,RESUMED_TEXT,'\n');flush.console()
    return(0L)
  }
  recorder <- new.env()
  fitter <- make_fitter(recorder)
  cat(format(Sys.time()),tag,'started; seed',spec$seed,'; schedule',paste(names(spec$mcmc),unlist(spec$mcmc),collapse=' '),'\n')
  flush.console()
  result <- evaluate_seeded(expr,prepared$input,prepared$job,spec$mcmc,bindings=list(.fitter=fitter),seed=spec$seed)
  elapsed <- as.numeric(difftime(result$finished,result$started,units='secs'))
  saved <- c(metadata,list(job=prepared$job,input_file=prepared$input_file,repo=repo,study=study,archives=archives,library=pkg),
    result,
    list(listPriors_used=recorder$listPriors,iterations_run=iterations_run(spec$mcmc),elapsed_seconds=elapsed,
      seconds_per_1000_iterations=1000*elapsed/iterations_run(spec$mcmc),
      input_md5_after=unname(tools::md5sum(prepared$input_file)),
      loaded_library=normalizePath(find.package('occJSDM')),session=sessionInfo()))
  finalise_fit(saved,dest,function(s) check_saved_fit(s,spec,prepared,source$job$priors,pkg))
  cat(format(Sys.time()),tag,'complete; fit seconds',round(elapsed,1),'; warnings',length(saved$warnings),'\n')
  flush.console()
  0L
}

if(sys.nframe()==0L) {
  status <- tryCatch(run_main(commandArgs(trailingOnly=TRUE)),
    error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
