#!/usr/bin/env Rscript
# Task 3 diagnostic-fit runner (PLAN.md and README.md in this directory): one
# single-chain fit of design-qfar_K6-sites300-05 (community 5) per process.
#
#   Rscript run.R --repo=REPO --study=ARCHIVE --inputs-root=DIR
#     --run=extended|a|b|c|pilot [--variant=extended|a|b|c] --chain=K
#     [--archives=DIR]
#
# --variant is required with --run=pilot (which fits that variant at the pilot
# schedule) and refused otherwise. --archives defaults to the directory
# holding ARCHIVE; its pr11-current-20260927 archive supplies the job record
# and the md5s recorded for the control runner.
#
# Each fit evaluates the pr11 control runner's own fitting statement (parsed
# from dev/simstudy/current-main-recheck/run.R after its md5 is checked against
# the one recorded with the pr11 fits), with the fitting function replaced by
# `.fitter`, which calls occJSDM::runOccJSDM with the statement's arguments and
# differs only as the run states (VARIANTS below):
#   extended  no change;
#   a         listPriors gains a_theta0 = 1 and b_theta0 = 100;
#   b         the fixed-theta0 clone of fixed-theta0.R holds species 6's theta0
#             at its generating value, input$truth$params$theta0[6];
#   c         no change to the fit, but the input is the community-5 input with
#             species 6 removed (drop_species()), saved once under
#             ARCHIVE/inputs/ and checked on every later use.
# Seeds (ruling R9): set.seed(20261001 + K) with the kind of the saved input
# RNG states (Mersenne-Twister, Inversion, Rejection) immediately before the
# statement; the package then seeds its sampler from R's stream as in every fit.
# Output: ARCHIVE/fits/<run>/chain-KK-fit.rds (pilots: fits/pilot/<variant>/),
# under an atomic per-fit lock (chain-KK.lock). A fit already saved with
# identical metadata is left alone and reported as resumed; one with any
# differing field is never overwritten. The installed library must match
# results/library-fingerprint.csv before it is loaded. Exit status 0 on
# success or resume, 1 on any failure. Launch batches with launch.R.

`%||%` <- function(x,y) if(is.null(x)) y else x

DIAG_DIR <- local({
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
SCRIPT_REPO <- normalizePath(file.path(DIAG_DIR,'../../..'))

STUDY_REL <- 'dev/simstudy/convergence-flag-diagnosis'
RUN_SCRIPT <- file.path(STUDY_REL,'run.R')
FIXED_THETA0_SCRIPT <- file.path(STUDY_REL,'fixed-theta0.R')
FINGERPRINT_FILE <- file.path(STUDY_REL,'results/library-fingerprint.csv')
CLONE_ARCHIVE_DIR <- file.path(STUDY_REL,'results/fixed-theta0')
JOBS_SCRIPT <- 'dev/simstudy/occupancy-intercept-prior/jobs.R'
PR11_RUNNER <- 'dev/simstudy/current-main-recheck/run.R'
PR11_HELPERS <- 'dev/simstudy/current-main-recheck/helpers.R'
PR11_ARCHIVE <- 'pr11-current-20260927'
PR11_SETTINGS <- file.path(PR11_ARCHIVE,c('initial/settings.rds','long/settings.rds'))
# Reused files, frozen at these md5s (jobs.R as committed in c28eedc; helpers.R
# as recorded with the pr11 fits).
REUSED_MD5 <- c(jobs='d6dc59ff47a25de08fad60b8bde827d1',helpers='50d5c6bac89f9739509a86c59380dfbf')

KEY <- 'design-qfar_K6-sites300-05'
TARGET_SPECIES <- 6L
TARGET_SPECIES_NAME <- 'OTU_6'
VARIANTS <- c('extended','a','b','c')
RUN_CHAINS <- list(extended=1:16,a=1:8,b=1:8,c=1:8,pilot=1:4)
SEED_BASE <- 20261001L
SEED_RULE <- 'set.seed(20261001 + chain, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection") immediately before the fitting statement'
RNG_KIND <- c(kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
# runOccJSDM runs nburn + niter * nthin iterations and keeps niter draws
# (R/runOccJSDM.R line 1182), so niter counts retained draws (ruling R11).
SCHEDULES <- list(diagnostic=list(nchain=1,nburn=10000,niter=10000,nthin=4),
  pilot=list(nchain=1,nburn=200,niter=200,nthin=4))
VARIANT_A_PRIORS <- list(a_theta0=1,b_theta0=100)
DERIVED_INPUT <- 'inputs/qfar_K6-sites300-05-without-species6.rds'
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
source(file.path(DIAG_DIR,'fixed-theta0.R'))

# ---- Runs --------------------------------------------------------------------------

chain_seed <- function(chain) SEED_BASE+as.integer(chain)

run_spec <- function(run,variant=NULL,chain) {
  if(!run %in% c(VARIANTS,'pilot')) stop('--run must be one of ',paste(c(VARIANTS,'pilot'),collapse=', '))
  if(run=='pilot') {
    if(is.null(variant) || !variant %in% VARIANTS) stop('--run=pilot needs --variant= one of ',paste(VARIANTS,collapse=', '))
  } else {
    if(!is.null(variant)) stop('--variant is only for --run=pilot')
    variant <- run
  }
  chain <- as.integer(chain)
  if(length(chain)!=1L || is.na(chain) || !chain %in% RUN_CHAINS[[run]])
    stop('--chain for run ',run,' must be one of ',min(RUN_CHAINS[[run]]),' to ',max(RUN_CHAINS[[run]]))
  schedule <- if(run=='pilot') 'pilot' else 'diagnostic'
  list(run=run,variant=variant,chain=chain,seed=chain_seed(chain),schedule=schedule,mcmc=SCHEDULES[[schedule]])
}

parse_run_args <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','run','variant','chain'),
    required=c('repo','study','inputs-root','run','chain'))
  if(!grepl('^[0-9]+$',o$chain)) stop('--chain must be a whole number')
  spec <- run_spec(o$run,o$variant,as.integer(o$chain))
  list(repo=o$repo,study=o$study,inputs_root=o$`inputs-root`,
    archives=o$archives %||% dirname(normalizePath(o$study,mustWork=FALSE)),spec=spec)
}

fit_file <- function(study,spec) {
  dir <- if(spec$run=='pilot') file.path(study,'fits','pilot',spec$variant) else file.path(study,'fits',spec$run)
  file.path(dir,sprintf('chain-%02d-fit.rds',spec$chain))
}

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

# ---- Variant (c): community 5 without species 6 ------------------------------------------

# Every species-indexed element of a design-family input, with the margin that
# indexes species (0 for a vector).
SPECIES_ELEMENTS <- list(
  list(path=c('sim','data_list','OTU'),margin=2L),
  list(path=c('sim','data_list','traits'),margin=1L),
  list(path=c('sim','true_params','jsdmParams_true','B0'),margin=0L),
  list(path=c('sim','true_params','jsdmParams_true','B'),margin=2L),
  list(path=c('sim','true_params','jsdmParams_true','Bt'),margin=1L),
  list(path=c('sim','true_params','jsdmParams_true','A'),margin=1L),
  list(path=c('sim','true_params','jsdmParams_true','Bs'),margin=2L),
  list(path=c('sim','true_params','jsdmParams_true','As'),margin=1L),
  list(path=c('sim','true_params','jsdmParams_true','Bst'),margin=1L),
  list(path=c('sim','true_params','jsdmParams_true','SE'),margin=2L),
  list(path=c('sim','true_params','jsdmParams_true','spatField'),margin=2L),
  list(path=c('sim','true_params','jsdmParams_true','L'),margin=2L),
  list(path=c('sim','true_params','jsdmParams_true','eta'),margin=2L),
  list(path=c('sim','true_params','jsdmParams_true','tau'),margin=0L),
  list(path=c('sim','true_params','beta_theta_true'),margin=2L),
  list(path=c('sim','true_params','z_true'),margin=2L),
  list(path=c('sim','true_params','w_true'),margin=2L),
  list(path=c('sim','true_params','p_true'),margin=2L),
  list(path=c('sim','true_params','q_true'),margin=2L),
  list(path=c('truth','jsdmParams','tau'),margin=0L),
  list(path=c('truth','params','p'),margin=2L),
  list(path=c('truth','params','q'),margin=2L),
  list(path=c('truth','params','theta0'),margin=0L),
  list(path=c('truth','params','theta_baseline'),margin=0L))
SPECIES_COUNTS <- list(c('scenario','S'),c('truth','datasettings','S'))

element_path <- function(p) paste(p,collapse='/')

# Paths of the leaves (non-list elements) with a dimension, or a length when
# they have no dimensions, equal to S.
species_indexed_paths <- function(x,S,prefix=character()) {
  if(is.list(x) && !is.data.frame(x)) {
    out <- character();nm <- names(x) %||% rep('',length(x))
    for(i in seq_along(x)) out <- c(out,species_indexed_paths(x[[i]],S,
      c(prefix,if(nzchar(nm[i])) nm[i] else paste0('[[',i,']]'))))
    return(out)
  }
  if(is.data.frame(x)) return(if(nrow(x)==S) element_path(prefix) else character())
  d <- dim(x)
  hit <- if(is.null(d)) length(x)==S else any(d==S)
  if(hit) element_path(prefix) else character()
}

drop_along <- function(x,margin,j) switch(as.character(margin),
  '0'=x[-j],'1'=x[-j,,drop=FALSE],'2'=x[,-j,drop=FALSE],stop('Unknown margin ',margin))

# The input with species j removed from the observation array and every
# species-indexed element, and the two species counts reduced by one. Refuses
# an input with a species-indexed element it does not know.
drop_species <- function(input,j=TARGET_SPECIES,name=TARGET_SPECIES_NAME) {
  S <- input$scenario$S
  if(!identical(S,input$truth$datasettings$S) || !identical(S,ncol(input$sim$data_list$OTU)))
    stop('Inconsistent species counts in the input')
  known <- vapply(SPECIES_ELEMENTS,function(e) element_path(e$path),character(1))
  found <- species_indexed_paths(input,S)
  if(!setequal(found,known)) stop('Species-indexed elements differ from the known list: ',
    paste(c(setdiff(found,known),setdiff(known,found)),collapse=', '))
  if(!identical(colnames(input$sim$data_list$OTU)[j],name) || !identical(rownames(input$sim$data_list$traits)[j],name))
    stop('Species ',j,' is not ',name)
  out <- input
  for(e in SPECIES_ELEMENTS) {
    x <- input[[e$path]]
    extra <- setdiff(names(attributes(x)),c('dim','dimnames','names'))
    if(length(extra)) stop('Element ',element_path(e$path),' has attributes ',paste(extra,collapse=', '))
    out[[e$path]] <- drop_along(x,e$margin,j)
  }
  for(p in SPECIES_COUNTS) out[[p]] <- S-1L
  left <- species_indexed_paths(out,S)
  if(length(left)) stop('Species-indexed elements left at ',S,': ',paste(left,collapse=', '))
  out
}

# Save the derived input once (under a lock shared by all variant (c)
# processes) or check an existing copy; returns its md5.
save_or_check_derived <- function(derived,file,wait_seconds=600) {
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
  lock <- paste0(file,'.lock');waited <- 0
  while(!dir.create(lock,showWarnings=FALSE)) {
    if(waited>=wait_seconds) stop('Timed out waiting for ',lock,'; remove it by hand if no run.R process holds it')
    Sys.sleep(.5);waited <- waited+.5
  }
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  if(file.exists(file)) {
    if(!identical(readRDS(file),derived)) stop('Existing derived input ',file,' differs from drop_species() of the input')
  } else atomic_save(derived,file)
  unname(tools::md5sum(file))
}

# ---- Fitting --------------------------------------------------------------------------

fixed_theta0_value <- function(input) {
  v <- input$truth$params$theta0[TARGET_SPECIES]
  if(length(v)!=1L || !is.finite(v) || v<=0 || v>=1) stop('No usable generating theta0 for species ',TARGET_SPECIES)
  v
}

# The function bound to `.fitter` in the fitting statement.
make_run_fitter <- function(variant,input,recorder,target=NULL) {
  if(is.null(target)) target <- if(variant=='b') make_fixed_theta0_fitter(TARGET_SPECIES,fixed_theta0_value(input)) else
    occJSDM::runOccJSDM
  added <- if(variant=='a') VARIANT_A_PRIORS else list()
  force(target);force(added)
  function(...,listPriors=list()) {
    clash <- intersect(names(added),names(listPriors))
    if(length(clash)) stop('The archived listPriors already set ',paste(clash,collapse=', '))
    listPriors <- c(listPriors,added)
    recorder$listPriors <- listPriors
    target(...,listPriors=listPriors)
  }
}

# Evaluate the fitting statement as the pr11 runner did, except for the RNG
# state: seed = 'input' assigns the saved input state (equivalence check);
# a number calls set.seed() with RNG_KIND. `expr` is only ever the statement
# parsed from the md5-checked pr11 runner.
evaluate_seeded <- function(expr,input,job,mcmc,bindings=list(),seed) {
  env <- new.env(parent=globalenv())
  env$input <- input;env$job <- job;env$mcmc <- mcmc;env$warnings <- character()
  for(name in names(bindings)) assign(name,bindings[[name]],envir=env)
  if(identical(seed,'input')) assign('.Random.seed',input$fit_rng,envir=.GlobalEnv) else
    set.seed(seed,kind=RNG_KIND[['kind']],normal.kind=RNG_KIND[['normal.kind']],sample.kind=RNG_KIND[['sample.kind']])
  rng_initial <- get('.Random.seed',envir=.GlobalEnv)
  if(!identical(rng_initial[1],input$fit_rng[1])) stop('RNG kind differs from the saved input state')
  started <- Sys.time()
  eval(expr,env)
  finished <- Sys.time()
  list(fit=env$fit,warnings=env$warnings,started=started,finished=finished,
    rng_initial=rng_initial,rng_final=get('.Random.seed',envir=.GlobalEnv))
}

variant_definition <- function(variant,input=NULL,clone_md5=NULL,source_input_md5=NULL,derived_md5=NULL) switch(variant,
  extended=list(change='none'),
  a=list(change='listPriors added',listPriors_added=VARIANT_A_PRIORS),
  b=list(change='theta0 of one species held fixed by the fixed-theta0.R clone of sample_theta0',
    species=TARGET_SPECIES,species_name=TARGET_SPECIES_NAME,fixed_value=fixed_theta0_value(input),
    fixed_value_source='input$truth$params$theta0[6]',clone_md5=clone_md5),
  c=list(change='species removed from the input',species=TARGET_SPECIES,species_name=TARGET_SPECIES_NAME,
    source_input_md5=source_input_md5,derived_input_md5=derived_md5,derived_input=DERIVED_INPUT),
  stop('Unknown variant ',variant))

# The fields that identify a fit; all must be identical for a saved fit to be
# accepted as this one. Absolute paths are saved beside them, never compared.
fit_metadata <- function(spec,job,input_md5,variant_def,fit_expression,runner,source_revision,fit_hashes,script_hashes)
  list(key=KEY,run=spec$run,variant=spec$variant,chain=spec$chain,seed=spec$seed,seed_rule=SEED_RULE,
    rng_kind=RNG_KIND,schedule=spec$schedule,mcmc=spec$mcmc,job_record=job[setdiff(names(job),'input_file')],
    input_md5=input_md5,variant_definition=variant_def,fit_expression=fit_expression,
    runner=runner[c('file','md5','helpers_md5','settings')],source_revision=source_revision,
    fit_hashes=fit_hashes,script_hashes=script_hashes,threads_per_fit=1L)

fit_script_files <- function() c(RUN_SCRIPT,FIXED_THETA0_SCRIPT,FINGERPRINT_FILE,file.path(CLONE_ARCHIVE_DIR,'hashes.csv'),
  JOBS_SCRIPT,PR11_RUNNER,PR11_HELPERS)

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
  S <- if(spec$variant=='c') 9L else 10L
  if(!identical(dim(ro$jsdm_output$B0_output),c(S,as.integer(mcmc$niter),1L))) stop('Unexpected B0_output dimensions')
  if(!identical(dim(ro$theta0_output),c(S,as.integer(mcmc$niter),1L))) stop('Unexpected theta0_output dimensions')
  expected <- c(job_priors,if(spec$variant=='a') VARIANT_A_PRIORS else list())
  if(!identical(s$listPriors_used,expected)) stop('listPriors passed to runOccJSDM differ from the variant definition')
  species <- fit$infos$speciesNames
  varies <- apply(ro$theta0_output[,,1,drop=FALSE],1,function(v) length(unique(v))>1L)
  if(spec$variant=='b') {
    if(!all(ro$theta0_output[TARGET_SPECIES,,1]==fixed_theta0_value(prepared$input))) stop('Species 6 theta0 draws are not the fixed value')
    if(!all(varies[-TARGET_SPECIES])) stop('A theta0 other than the fixed one is constant')
  } else if(!all(varies)) stop('A theta0 is constant')
  if(spec$variant=='c') {
    if(TARGET_SPECIES_NAME %in% species || length(species)!=9L) stop('Variant (c) fit still holds ',TARGET_SPECIES_NAME)
  } else if(!identical(species[TARGET_SPECIES],TARGET_SPECIES_NAME)) stop('Species 6 is not ',TARGET_SPECIES_NAME)
  if(!identical(s$input_md5_after,prepared$input_md5)) stop('Input changed during fitting: ',prepared$input_file)
  if(!identical(s$loaded_library,pkg)) stop('Loaded library changed during fitting')
  invisible(TRUE)
}

prepare_input <- function(spec,source,study) {
  checked_input(source$input_file,source$input_md5)
  original <- readRDS(source$input_file)
  if(spec$variant!='c') return(list(input=original,job=source$job,input_file=source$input_file,input_md5=source$input_md5))
  derived <- drop_species(original)
  file <- file.path(study,DERIVED_INPUT)
  md5 <- save_or_check_derived(derived,file)
  job <- source$job
  job$input_file <- file;job$input_md5 <- md5;job$source_input_md5 <- source$input_md5
  job$dropped_species <- TARGET_SPECIES
  list(input=derived,job=job,input_file=file,input_md5=md5)
}

run_main <- function(args) {
  opt <- parse_run_args(args);spec <- opt$spec
  repo <- normalizePath(opt$repo,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo (',repo,') is not the checkout holding this run.R (',SCRIPT_REPO,')')
  study <- normalizePath(opt$study,mustWork=TRUE);archives <- normalizePath(opt$archives,mustWork=TRUE)
  inputs_root <- normalizePath(opt$inputs_root,mustWork=TRUE)
  tag <- paste0(KEY,' ',spec$run,if(spec$run=='pilot') paste0('/',spec$variant),' chain ',spec$chain)
  # Refuse any build other than the fingerprinted one before loading it.
  source_revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
  pkg <- load_study_library(study)
  fit_hashes <- hash_files(production_files(study,pkg),study)
  runner <- pr11_runner(repo,archives)
  expr <- with_intercept_prior(runner$expr);fit_expression <- deparse(expr,width.cutoff=500L)
  source <- pr11_job(archives,inputs_root)
  clone_md5 <- if(spec$variant=='b') check_fixed_theta0_archive(file.path(repo,CLONE_ARCHIVE_DIR))[['clone_sample_theta0_fixed']] else NULL
  script_hashes <- c(hash_files(file.path(repo,fit_script_files()),repo),
    hash_files(file.path(archives,PR11_SETTINGS),archives,'archives/'))
  dest <- fit_file(study,spec);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  lock <- acquire_fit_lock(dest);on.exit(release_fit_lock(lock),add=TRUE)
  prepared <- prepare_input(spec,source,study)
  def <- variant_definition(spec$variant,prepared$input,clone_md5,source$input_md5,
    if(spec$variant=='c') prepared$input_md5 else NULL)
  metadata <- fit_metadata(spec,prepared$job,prepared$input_md5,def,fit_expression,runner,source_revision,fit_hashes,script_hashes)
  if(existing_fit_status(dest,metadata)=='resume') {
    cat(format(Sys.time()),tag,RESUMED_TEXT,'\n');flush.console()
    return(0L)
  }
  recorder <- new.env()
  fitter <- make_run_fitter(spec$variant,prepared$input,recorder)
  cat(format(Sys.time()),tag,'started; seed',spec$seed,'; schedule',paste(names(spec$mcmc),unlist(spec$mcmc),collapse=' '),'\n')
  flush.console()
  result <- evaluate_seeded(expr,prepared$input,prepared$job,spec$mcmc,bindings=list(.fitter=fitter),seed=spec$seed)
  elapsed <- as.numeric(difftime(result$finished,result$started,units='secs'))
  saved <- c(metadata,list(job=prepared$job,input_file=prepared$input_file,source_input_file=source$input_file,
      repo=repo,study=study,archives=archives,library=pkg),result,
    list(listPriors_used=recorder$listPriors,iterations_run=iterations_run(spec$mcmc),elapsed_seconds=elapsed,
      seconds_per_1000_iterations=1000*elapsed/iterations_run(spec$mcmc),
      input_md5_after=unname(tools::md5sum(prepared$input_file)),
      source_input_md5_after=unname(tools::md5sum(source$input_file)),
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
