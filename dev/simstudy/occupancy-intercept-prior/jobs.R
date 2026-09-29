# Helpers for the occupancy-intercept prior study. No package is loaded here;
# run.R, equivalence.R, controls.R, pilot.R and test-run.R source this file.
# Archived inputs, settings, selections and fits are read, never written.

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

archive_dir <- function(archives,name) file.path(archives,ARCHIVES[[name]])
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

# run.R interface. --archives defaults to the directory holding --study, and
# --inputs-root to <archives>/intercept-prior-inputs.
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
  archives <- o$archives %||% dirname(o$study)
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

# Relocate a recorded control-fit path into the given archives root.
relocate_archive_path <- function(path,archives,default) {
  for(name in names(ARCHIVES)) {
    tag <- paste0('/',ARCHIVES[[name]],'/');at <- regexpr(tag,path,fixed=TRUE)
    if(at>0) return(file.path(archive_dir(archives,name),substring(path,at+nchar(tag))))
  }
  if(startsWith(path,'/')) stop('Control fit is outside the known archives: ',path)
  file.path(default,path)
}

# For each community, the schedule of the control arm's selected fit.
control_schedules <- function(phase,archives) {
  keys <- phase_keys(phase)
  if(phase %in% c('A1','B')) {
    pr11 <- archive_dir(archives,'pr11')
    selection <- readRDS(file.path(pr11,'selection.rds'))
    plan <- read.csv(file.path(pr11,'long-selection.csv'),stringsAsFactors=FALSE)
    i <- match(keys,selection$key);j <- match(keys,plan$key)
    if(anyNA(i) || anyNA(j)) stop('Recorded pr11 selection is missing phase ',phase,' keys')
    schedule <- selection$schedule[i]
    if(!all(schedule %in% c('initial','long')) || !identical(schedule=='long',as.logical(plan$run_long[j])))
      stop('pr11 selection.rds and long-selection.csv disagree for phase ',phase)
    result <- selection$result_file[i]
    if(!identical(basename(result),paste0(keys,'-result.rds')) || !identical(basename(dirname(result)),schedule))
      stop('pr11 selected result files disagree with the selected schedules')
    fit <- file.path(pr11,schedule,paste0(keys,'-fit.rds'));md5 <- rep(NA_character_,length(keys))
  } else if(phase=='A2') {
    amplitude <- archive_dir(archives,'amplitude')
    final <- file.path(amplitude,'robust-v1/summary-binary-final')
    fits <- read.csv(file.path(final,'fits.csv'),stringsAsFactors=FALSE)
    fits <- fits[fits$prior=='inverse_gamma',]
    selection <- read.csv(file.path(final,'selection.csv'),stringsAsFactors=FALSE)
    i <- match(keys,fits$key);j <- match(keys,selection$key)
    if(anyNA(i) || anyNA(j) || anyDuplicated(fits$key)) stop('Recorded spatial-amplitude selection is missing A2 keys')
    schedule <- fits$phase[i]
    if(!all(schedule %in% c('initial','long')) || !identical(schedule=='long',as.logical(selection$needs_long[j])))
      stop('spatial-amplitude fits.csv and selection.csv disagree')
    fit <- vapply(fits$fit_file[i],relocate_archive_path,character(1),archives=archives,default=amplitude,USE.NAMES=FALSE)
    if(!identical(basename(fit),paste0(keys,'-fit.rds')) || !identical(basename(dirname(fit)),schedule))
      stop('spatial-amplitude selected fit files disagree with the selected schedules')
    md5 <- fits$fit_md5[i]
  } else stop('Unknown phase: ',phase)
  data.frame(phase=phase,key=keys,schedule=schedule,control_fit=fit,control_fit_md5=md5,stringsAsFactors=FALSE)
}

fit_path <- function(study,phase,sd,schedule,key)
  file.path(study,'fits',phase,paste0('sd',format(sd)),schedule,paste0(key,'-fit.rds'))

atomic_save <- function(x,path) {
  temp <- paste0(path,'.tmp');saveRDS(x,temp);stopifnot(file.rename(temp,path))
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
  settings <- file.path(archive_dir(archives,name),spec$settings)
  md5 <- recorded_hash(readRDS(settings)[[spec$field]],spec$file)
  file <- file.path(repo,spec$file)
  list(name=name,file=file,md5=md5,settings=settings,expr=archived_fit_expression(file,md5))
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
production_files <- function(study,pkg) {
  files <- c(file.path(study,'source',c('DESCRIPTION','NAMESPACE')),
    list.files(file.path(study,'source/R'),pattern='\\.R$',full.names=TRUE),
    list.files(file.path(study,'source/src'),pattern='\\.(cpp|h)$|^Makevars',full.names=TRUE),
    file.path(pkg,'libs/occJSDM.so'),file.path(pkg,c('DESCRIPTION','NAMESPACE')),
    list.files(file.path(pkg,'R'),full.names=TRUE))
  stopifnot(length(files)>10L,all(file.exists(files)))
  files
}

installed_library_files <- function(pkg)
  c(file.path(pkg,'libs/occJSDM.so'),file.path(pkg,c('DESCRIPTION','NAMESPACE')),list.files(file.path(pkg,'R'),full.names=TRUE))

phase_priors <- function(phase,archives) {
  if(phase!='A2') return(NULL)
  priors <- readRDS(file.path(archive_dir(archives,'amplitude'),'inverse_gamma/long/settings.rds'))$priors
  if(!identical(priors,list())) stop('Inverse-gamma control priors are not the defaults')
  priors
}

# Names and dimensions of every posterior output array, for structure checks.
posterior_layout <- function(results) {
  walk <- function(x,prefix) {
    if(is.list(x) && !is.data.frame(x)) return(do.call(rbind,c(list(NULL),lapply(names(x),function(n) walk(x[[n]],paste0(prefix,n,'/'))))))
    data.frame(path=sub('/$','',prefix),class=class(x)[1],dim=paste(dim(x) %||% length(x),collapse='x'),stringsAsFactors=FALSE)
  }
  walk(results,'')
}

# Largest absolute difference over all numeric leaves; Inf when structures differ.
max_abs_difference <- function(a,b) {
  if(is.list(a) || is.list(b)) {
    if(!is.list(a) || !is.list(b) || !identical(names(a),names(b)) || length(a)!=length(b)) return(Inf)
    if(!length(a)) return(0)
    return(max(vapply(seq_along(a),function(i) max_abs_difference(a[[i]],b[[i]]),numeric(1))))
  }
  if(!identical(dim(a),dim(b)) || length(a)!=length(b)) return(Inf)
  if(is.numeric(a) && is.numeric(b)) {
    if(!identical(is.na(a),is.na(b))) return(Inf)
    ok <- !is.na(a);if(!any(ok)) return(0)
    inf <- is.infinite(a[ok]) | is.infinite(b[ok])
    if(any(inf) && !identical(a[ok][inf],b[ok][inf])) return(Inf)
    d <- abs(a[ok][!inf]-b[ok][!inf]);return(if(length(d)) max(d) else 0)
  }
  if(identical(a,b)) 0 else Inf
}

count_numeric <- function(x) if(is.list(x)) sum(vapply(x,count_numeric,numeric(1))) else if(is.numeric(x)) length(x) else 0
