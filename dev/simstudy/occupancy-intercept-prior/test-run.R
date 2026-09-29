library(testthat)
source('jobs.R')

# Real archives are optional: tests that need them are skipped when absent.
repo <- normalizePath('../../..')
archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
inputs_root <- file.path(archives,'intercept-prior-inputs')
have_archives <- all(dir.exists(file.path(archives,c('pr11-current-20260927',
  'spatial-targeted-20260927','spatial-amplitude-20260928','intercept-prior-inputs'))))
need_archives <- function() if(!have_archives) skip('saved study archives not available')

# A minimal archive tree holding only the settings files the schedule code reads.
fake_archive <- function(pr11_initial=list(nchain=2,nburn=3000,niter=5000,nthin=1),
  targeted_initial=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L)) {
  root <- tempfile('archives');dir.create(root)
  put <- function(path,x) {dir.create(dirname(file.path(root,path)),recursive=TRUE,showWarnings=FALSE)
    saveRDS(x,file.path(root,path))}
  put('pr11-current-20260927/initial/settings.rds',list(mcmc=pr11_initial))
  put('pr11-current-20260927/long/settings.rds',list(mcmc=list(nchain=4,nburn=6000,niter=12000,nthin=1)))
  put('spatial-amplitude-20260928/half_cauchy/initial/settings.rds',
    list(mcmc=list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L)))
  put('spatial-amplitude-20260928/inverse_gamma/long/settings.rds',
    list(mcmc=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L),priors=list()))
  put('spatial-targeted-20260927/initial/settings.rds',list(mcmc=targeted_initial))
  put('spatial-targeted-20260927/long/settings.rds',list(mcmc=list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L)))
  root
}

base_args <- c('--repo=/r','--study=/x/results/study','--phase=A1','--sd=3','--schedule=pilot','--workers=2')
with_arg <- function(name,value) c(base_args[!startsWith(base_args,paste0('--',name,'='))],
  if(!is.null(value)) paste0('--',name,'=',value))

test_that('arguments parse, default and validate', {
  a <- parse_run_args(base_args)
  expect_identical(a$repo,'/r');expect_identical(a$study,'/x/results/study')
  expect_identical(a$phase,'A1');expect_identical(a$sd,3);expect_identical(a$schedule,'pilot')
  expect_identical(a$workers,2L);expect_null(a$keys)
  expect_identical(a$archives,'/x/results')
  expect_identical(a$inputs_root,'/x/results/intercept-prior-inputs')
  b <- parse_run_args(c(base_args,'--keys=k1,k2','--inputs-root=/in','--archives=/arch'))
  expect_identical(b$keys,c('k1','k2'));expect_identical(b$inputs_root,'/in');expect_identical(b$archives,'/arch')
  for(s in c('1','2','3','5')) expect_identical(parse_run_args(with_arg('sd',s))$sd,as.numeric(s))
  for(p in c('A1','A2','B')) expect_identical(parse_run_args(with_arg('phase',p))$phase,p)
  for(s in c('initial','long','pilot')) expect_identical(parse_run_args(with_arg('schedule',s))$schedule,s)
  expect_identical(parse_run_args(with_arg('workers','8'))$workers,8L)
  expect_error(parse_run_args(with_arg('sd',NULL)),'Missing --sd')
  expect_error(parse_run_args(with_arg('repo',NULL)),'Missing --repo')
  expect_error(parse_run_args(c(base_args,'--sd=2')),'Repeated option')
  expect_error(parse_run_args(c(base_args,'--prior=half_cauchy')),'Unknown option')
  expect_error(parse_run_args(c(base_args,'stray')),'Malformed argument')
  for(bad in c('4','0.5','abc','','1e0x')) expect_error(parse_run_args(with_arg('sd',bad)),'--sd')
  expect_error(parse_run_args(with_arg('phase','C')),'--phase')
  expect_error(parse_run_args(with_arg('schedule','short')),'--schedule')
  for(bad in c('0','9','x','2.5')) expect_error(parse_run_args(with_arg('workers',bad)),'--workers')
  expect_error(parse_run_args(c(base_args,'--keys=k1,k1')),'Duplicated')
  expect_error(parse_run_args(c(base_args,'--keys=')),'Empty --keys')
})

test_that('input paths are remapped under the copy root and never read from the original root', {
  orig <- ORIGINAL_INPUT_ROOT
  expect_identical(remap_input(file.path(orig,'jsdm-sample-size-20260919/inputs/n0100-01.rds'),'/copy'),
    '/copy/jsdm-sample-size-20260919/inputs/n0100-01.rds')
  expect_error(remap_input('/elsewhere/a.rds','/copy'),'not under the original input root')
  expect_error(remap_input(paste0(orig,'-other/a.rds'),'/copy'),'not under the original input root')
  expect_error(remap_input(orig,'/copy'),'not under the original input root')
})

test_that('inputs whose md5 differs from the recorded job md5 are refused', {
  f <- tempfile(fileext='.rds');saveRDS(list(a=1),f);md5 <- unname(tools::md5sum(f))
  expect_identical(checked_input(f,md5),f)
  expect_error(checked_input(f,'0123456789abcdef0123456789abcdef'),'md5 mismatch')
  expect_error(checked_input(tempfile(),md5),'Missing input')
  saveRDS(list(a=2),f)
  expect_error(checked_input(f,md5),'md5 mismatch')
})

test_that('every A1 and B job resolves to a copied input with the recorded md5, and A2 reads in place', {
  need_archives()
  for(phase in c('A1','B')) {
    jobs <- phase_jobs(phase,archives,inputs_root)
    expect_identical(names(jobs),phase_keys(phase))
    for(spec in jobs) {
      expect_true(startsWith(spec$job$input_file,paste0(ORIGINAL_INPUT_ROOT,'/')))
      expect_true(startsWith(spec$input_file,paste0(inputs_root,'/')))
      expect_identical(spec$input_md5,spec$job$input_md5)
      expect_identical(checked_input(spec$input_file,spec$input_md5),spec$input_file)
    }
  }
  a2 <- phase_jobs('A2',archives,inputs_root)
  expect_identical(names(a2),phase_keys('A2'))
  expect_identical(names(a2[[1]]$job),c('key','community','grid_index','replicate','arm','knots','input_md5','input_file'))
  for(spec in a2) {
    expect_identical(spec$job$arm,'binary');expect_identical(spec$job$knots,100L)
    expect_identical(spec$input_file,file.path(archives,'spatial-targeted-20260927/inputs',paste0(spec$job$community,'.rds')))
    expect_identical(checked_input(spec$input_file,spec$input_md5),spec$input_file)
  }
  expect_error(phase_jobs('A1',archives,inputs_root,keys='design-qnear_K6-sites300-01'),'not in phase A1')
  expect_identical(names(phase_jobs('B',archives,inputs_root,keys='design-qfar_K6-sites300-03')),
    'design-qfar_K6-sites300-03')
  # A copy root lacking the file is refused before any read.
  empty <- tempfile('copy');dir.create(empty)
  spec <- phase_jobs('A1',archives,empty,keys='jsdm-n0100-01')[[1]]
  expect_error(checked_input(spec$input_file,spec$input_md5),'Missing input')
})

test_that('phase keys are the prespecified communities', {
  expect_identical(phase_keys('A1'),c(sprintf('jsdm-n0100-%02d',1:10),sprintf('jsdm-n0300-%02d',1:10)))
  expect_identical(phase_keys('B'),c(sprintf('design-qnear_K6-sites300-%02d',1:10),
    sprintf('design-qfar_K6-sites300-%02d',1:10)))
  expect_identical(phase_keys('A2'),sprintf('range%d-rep%02d-binary-k100',rep(c(4L,6L,8L),each=3L),1:3))
  expect_error(phase_keys('C'),'phase')
})

test_that('schedule settings come from the archives and must equal the protocol', {
  root <- fake_archive()
  expect_identical(schedule_mcmc('A1','initial',root),list(nchain=2,nburn=3000,niter=5000,nthin=1))
  expect_identical(schedule_mcmc('B','long',root),list(nchain=4,nburn=6000,niter=12000,nthin=1))
  expect_identical(schedule_mcmc('B','pilot',root),list(nchain=2,nburn=200,niter=200,nthin=1))
  expect_identical(schedule_mcmc('A2','initial',root),list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L))
  expect_identical(schedule_mcmc('A2','long',root),list(nchain=4L,nburn=6000L,niter=12000L,nthin=1L))
  expect_identical(schedule_mcmc('A2','pilot',root),list(nchain=2L,nburn=200L,niter=200L,nthin=1L))
  expect_identical(short_mcmc(list(nchain=2L,nburn=3000L,niter=5000L,nthin=1L),50,50),
    list(nchain=2L,nburn=50L,niter=50L,nthin=1L))
  expect_error(schedule_mcmc('A1','initial',fake_archive(pr11_initial=list(nchain=2,nburn=1000,niter=5000,nthin=1))),
    'differ from the protocol')
  expect_error(schedule_mcmc('A2','initial',fake_archive(targeted_initial=list(nchain=2L,nburn=2999L,niter=5000L,nthin=1L))),
    'differ from the protocol')
  expect_error(schedule_mcmc('A1','medium',root),'schedule')
})

test_that('the real archived schedules equal the protocol numbers', {
  need_archives()
  for(phase in c('A1','A2','B')) {
    expect_equal(unlist(schedule_mcmc(phase,'initial',archives)),c(nchain=2,nburn=3000,niter=5000,nthin=1))
    expect_equal(unlist(schedule_mcmc(phase,'long',archives)),c(nchain=4,nburn=6000,niter=12000,nthin=1))
    expect_equal(unlist(schedule_mcmc(phase,'pilot',archives)),c(nchain=2,nburn=200,niter=200,nthin=1))
  }
})

test_that('the control schedule helper reproduces the recorded selections', {
  need_archives()
  a1 <- control_schedules('A1',archives)
  expect_identical(a1$key,phase_keys('A1'));expect_true(all(a1$schedule=='initial'))
  b <- control_schedules('B',archives)
  expect_identical(b$key,phase_keys('B'))
  expected_long <- c(sprintf('design-qnear_K6-sites300-%02d',c(2,4,7,8)),
    sprintf('design-qfar_K6-sites300-%02d',c(1,2,3,5,6,7,8,9)))
  expect_setequal(b$key[b$schedule=='long'],expected_long)
  expect_true(all(b$schedule[!b$key %in% expected_long]=='initial'))
  # The committed copies of the pr11 selection outputs agree.
  committed <- file.path(repo,'dev/simstudy/current-main-recheck/results')
  ls <- read.csv(file.path(committed,'long-selection.csv'));sm <- read.csv(file.path(committed,'selected-manifest.csv'))
  ab <- rbind(a1,b)
  expect_identical(ab$schedule=='long',ls$run_long[match(ab$key,ls$key)])
  expect_identical(ab$schedule,sm$schedule[match(ab$key,sm$key)])
  expect_identical(ab$control_fit,file.path(archives,'pr11-current-20260927',ab$schedule,paste0(ab$key,'-fit.rds')))
  a2 <- control_schedules('A2',archives)
  expect_identical(a2$key,phase_keys('A2'));expect_true(all(a2$schedule=='long'))
  expect_identical(a2$control_fit[a2$key=='range6-rep03-binary-k100'],
    file.path(archives,'spatial-targeted-20260927/long/range6-rep03-binary-k100-fit.rds'))
  others <- a2$key!='range6-rep03-binary-k100'
  expect_identical(a2$control_fit[others],
    file.path(archives,'spatial-amplitude-20260928/inverse_gamma/long',paste0(a2$key[others],'-fit.rds')))
  sa <- read.csv(file.path(repo,'dev/simstudy/spatial-amplitude-prior/results/fits.csv'))
  sa <- sa[sa$prior=='inverse_gamma',]
  expect_identical(a2$schedule,sa$phase[match(a2$key,sa$key)])
  expect_identical(basename(a2$control_fit),basename(sa$fit_file[match(a2$key,sa$key)]))
  expect_true(all(file.exists(c(ab$control_fit,a2$control_fit))))
  expect_true(all(grepl('^[0-9a-f]{32}$',a2$control_fit_md5)))
})

test_that('the control schedule helper refuses inconsistent selection records', {
  root <- tempfile('archives');pr11 <- file.path(root,'pr11-current-20260927');dir.create(pr11,recursive=TRUE)
  keys <- phase_keys('A1')
  sel <- data.frame(key=keys,schedule='initial',result_file=file.path(pr11,'initial',paste0(keys,'-result.rds')))
  saveRDS(sel,file.path(pr11,'selection.rds'))
  write.csv(data.frame(key=keys,run_long=FALSE),file.path(pr11,'long-selection.csv'),row.names=FALSE)
  expect_identical(control_schedules('A1',root)$schedule,rep('initial',20L))
  write.csv(data.frame(key=keys,run_long=c(TRUE,rep(FALSE,19L))),file.path(pr11,'long-selection.csv'),row.names=FALSE)
  expect_error(control_schedules('A1',root),'disagree')
  write.csv(data.frame(key=keys,run_long=FALSE),file.path(pr11,'long-selection.csv'),row.names=FALSE)
  saveRDS(sel[-1,],file.path(pr11,'selection.rds'))
  expect_error(control_schedules('A1',root),'missing')
})

test_that('fit paths follow the study layout', {
  expect_identical(fit_path('/s','A1',3,'pilot','jsdm-n0100-01'),'/s/fits/A1/sd3/pilot/jsdm-n0100-01-fit.rds')
  expect_identical(fit_path('/s','A2',1,'long','range6-rep01-binary-k100'),
    '/s/fits/A2/sd1/long/range6-rep01-binary-k100-fit.rds')
})

test_that('existing fits are resumed only when identical and are never overwritten', {
  f <- tempfile(fileext='-fit.rds')
  meta <- list(phase='A1',key='k',sd=3,schedule='pilot',mcmc=list(nchain=2,nburn=200,niter=200,nthin=1),
    input_md5='abc',fit_hashes=c(a='1'),script_hashes=c(b='2'))
  expect_identical(existing_fit_status(f,meta),'new')
  atomic_save(c(meta,list(fit=list(x=1),warnings=character())),f)
  expect_false(file.exists(paste0(f,'.tmp')))
  expect_identical(existing_fit_status(f,meta),'resume')
  before <- tools::md5sum(f)
  changed <- meta;changed$sd <- 5
  expect_error(existing_fit_status(f,changed),'Refusing to overwrite.*sd')
  changed <- meta;changed$fit_hashes <- c(a='9')
  expect_error(existing_fit_status(f,changed),'Refusing to overwrite.*fit_hashes')
  changed <- meta;changed$mcmc$niter <- 201
  expect_error(existing_fit_status(f,changed),'Refusing to overwrite.*mcmc')
  changed <- c(meta,list(extra='field'))
  expect_error(existing_fit_status(f,changed),'Refusing to overwrite.*extra')
  expect_identical(tools::md5sum(f),before)
})

test_that('archived fitting expressions are hash-checked and changed only by the intercept prior', {
  runner <- file.path(repo,'dev/simstudy/current-main-recheck/run.R')
  md5 <- unname(tools::md5sum(runner))
  x <- archived_fit_expression(runner,md5)
  expect_true(is.call(x));expect_identical(x[[1]],as.name('<-'));expect_identical(x[[2]],as.name('fit'))
  expect_error(archived_fit_expression(runner,'0123456789abcdef0123456789abcdef'),'md5')
  s <- with_intercept_prior(x)
  expect_false(any(grepl('runOccJSDM',deparse(s))))
  # Restoring the fitter symbol gives back the archived expression exactly.
  expect_identical(replace_in_call(s,as.name('.fitter'),quote(occJSDM::runOccJSDM)),x)
  capture <- function(...) list(...)
  set.seed(1);rng <- .Random.seed
  mcmc <- list(nchain=2,nburn=50,niter=50,nthin=1)
  jsdm_input <- list(data='DATA',n_factors=2L,n_lattrait=1L,covariates=c('c1','c2'),fit_rng=rng)
  jsdm_job <- list(family='jsdm',arm='sites100',priors=list())
  r <- evaluate_fit(s,jsdm_input,jsdm_job,mcmc,bindings=list(.fitter=make_fitter(3,target=capture)))
  expect_identical(r$fit,list('DATA',listParams=list(n_factors=2L,n_lattrait=1L),occCovariates=c('c1','c2'),
    collCovariates=NULL,spatCovariates=NULL,MCMCparams=mcmc,listPriors=list(sigma_b0=3)))
  expect_identical(r$rng_initial,rng)
  r0 <- evaluate_fit(s,jsdm_input,jsdm_job,mcmc,bindings=list(.fitter=make_fitter(NULL,target=capture)))
  expect_identical(r0$fit$listPriors,list())
  design_input <- list(sim=list(data_list='SIM'),scenario=list(d=2L,gt=1L,ncov_psi=2L),fit_rng=rng)
  design_job <- list(family='design',arm='sites300',priors=list(a_q=1,b_q=20))
  r <- evaluate_fit(s,design_input,design_job,mcmc,bindings=list(.fitter=make_fitter(5,target=capture)))
  expect_identical(r$fit,list('SIM',listParams=list(n_factors=2L,n_lattrait=1L),threshold=1,
    occCovariates=c('X_psi.EnvCov.1','X_psi.EnvCov.2'),collCovariates='X_theta',spatCovariates=NULL,
    MCMCparams=mcmc,listPriors=list(a_q=1,b_q=20,sigma_b0=5)))
  spatial <- file.path(repo,'dev/simstudy/spatial-amplitude-prior/run.R')
  y0 <- archived_fit_expression(spatial,unname(tools::md5sum(spatial)))
  y <- with_intercept_prior(y0)
  expect_false(any(grepl('\\bfun\\(',deparse(y),perl=TRUE)))
  expect_identical(replace_in_call(y,as.name('.fitter'),as.name('fun')),y0)
  spatial_input <- list(data=list(binary='BIN',low='LOW'),fit_rng=rng)
  spatial_job <- list(key='k',arm='binary',knots=100L)
  r <- evaluate_fit(y,spatial_input,spatial_job,mcmc,
    bindings=list(.fitter=make_fitter(2,target=capture),priors=list()))
  expect_identical(r$fit,list('BIN',listParams=list(n_factors=0L,n_lattrait=0L,n_supportpoints=100L),
    threshold=1,occCovariates='environment',collCovariates=NULL,
    spatCovariates=c('longitude','latitude'),MCMCparams=mcmc,listPriors=list(sigma_b0=2)))
  noisy <- function(...) {warning('first');warning('second');'done'}
  r <- evaluate_fit(y,spatial_input,spatial_job,mcmc,bindings=list(.fitter=make_fitter(2,target=noisy),priors=list()))
  expect_identical(r$warnings,c('first','second'));expect_identical(r$fit,'done')
  expect_error(with_intercept_prior(quote(fit <- other(x))),'no fitting call')
})

test_that('saved fits must record the requested intercept prior', {
  expect_true(check_fit_prior(list(infos=list(intercept_prior=list(mean=0,sd=3))),3))
  expect_error(check_fit_prior(list(infos=list(intercept_prior=list(mean=0,sd=1))),3),'intercept prior')
  expect_error(check_fit_prior(list(infos=list()),1),'intercept prior')
  expect_true(check_fit_prior(list(infos=list(intercept_prior=list(mean=0,sd=1))),NULL))
})
