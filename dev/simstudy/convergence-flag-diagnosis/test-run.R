# Standalone research tests for run.R, fixed-theta0.R and launch.R (Task 3 of
# PLAN.md in this directory). Run with `Rscript test-run.R` from this
# directory. Tests that need the study library, the pr11 archive or the
# community-5 input are skipped when they are absent; they read those files
# and write only to temporary directories.
library(testthat)
source('run.R')
source('launch.R')

archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
study <- file.path(archives,'convergence-diagnosis-20261001')
inputs_root <- file.path(archives,'intercept-prior-inputs')
have <- dir.exists(file.path(study,'library/occJSDM')) && dir.exists(file.path(archives,PR11_ARCHIVE)) && dir.exists(inputs_root)
need <- function() if(!have) skip('study library or archives not available')
if(have) {
  load_study_library(study)
  src <- pr11_job(archives,inputs_root)
  input <- readRDS(checked_input(src$input_file,src$input_md5))
  runner <- pr11_runner(SCRIPT_REPO,archives)
  expr <- with_intercept_prior(runner$expr)
}
tiny <- list(nchain=1,nburn=3,niter=5,nthin=4)

test_that('runs, chains, seeds and schedules are the frozen ones', {
  expect_identical(RUN_CHAINS[c('extended','a','b','c')],list(extended=1:16,a=1:8,b=1:8,c=1:8))
  expect_identical(chain_seed(1:16),20261001L+1:16)
  for(v in c('a','b','c')) expect_identical(vapply(1:8,function(k) run_spec(v,NULL,k)$seed,integer(1)),
    vapply(1:8,function(k) run_spec('extended',NULL,k)$seed,integer(1)))
  d <- run_spec('extended',NULL,3L)
  expect_identical(d$mcmc,list(nchain=1,nburn=10000,niter=10000,nthin=4))
  expect_identical(iterations_run(d$mcmc),50000)
  expect_identical(d$mcmc$niter*d$mcmc$nthin,40000)
  p <- run_spec('pilot','c',1L)
  expect_identical(p$variant,'c');expect_identical(p$schedule,'pilot')
  expect_identical(fit_file('/S',d),'/S/fits/extended/chain-03-fit.rds')
  expect_identical(fit_file('/S',p),'/S/fits/pilot/c/chain-01-fit.rds')
  expect_identical(fit_lock_path(fit_file('/S',d)),'/S/fits/extended/chain-03.lock')
})

test_that('argument parsing refuses unknown runs, chains outside the run and misplaced variants', {
  base <- c('--repo=/r','--study=/s/x','--inputs-root=/i')
  ok <- parse_run_args(c(base,'--run=b','--chain=8'))
  expect_identical(ok$spec$variant,'b');expect_identical(ok$archives,'/s')
  expect_error(parse_run_args(c(base,'--run=d','--chain=1')),'--run must be')
  expect_error(parse_run_args(c(base,'--run=a','--chain=9')),'--chain for run a')
  expect_error(parse_run_args(c(base,'--run=extended','--chain=17')),'--chain for run extended')
  expect_error(parse_run_args(c(base,'--run=a','--chain=1.5')),'whole number')
  expect_error(parse_run_args(c(base,'--run=pilot','--chain=1')),'needs --variant')
  expect_error(parse_run_args(c(base,'--run=a','--variant=b','--chain=1')),'only for --run=pilot')
  expect_error(parse_run_args(c(base,'--run=a')),'Missing --chain')
})

test_that('the installed runOccJSDM runs nburn + niter * nthin iterations and keeps niter draws', {
  need()
  body_text <- deparse(body(occJSDM::runOccJSDM))
  expect_true(any(grepl('for (iter in 1:(nburn + niter * nthin))',body_text,fixed=TRUE)))
  expect_true(any(grepl('if (iter > nburn & (iter - nburn)%%nthin == 0)',body_text,fixed=TRUE)))
  rec <- new.env()
  r <- evaluate_seeded(expr,input,src$job,tiny,bindings=list(.fitter=make_run_fitter('extended',input,rec)),seed=chain_seed(1))
  expect_identical(dim(r$fit$results_output$theta0_output),c(10L,5L,1L))
  expect_identical(dim(r$fit$results_output$jsdm_output$U_output),c(300L,2L,5L,1L))
  set.seed(chain_seed(1),kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
  expect_identical(r$rng_initial,.Random.seed)
  expect_identical(r$rng_initial[1],input$fit_rng[1])
})

test_that('drop_species removes species 6 everywhere and leaves every other species byte-identical', {
  need()
  derived <- drop_species(input)
  known <- vapply(SPECIES_ELEMENTS,function(e) element_path(e$path),character(1))
  expect_setequal(species_indexed_paths(input,10L),known)
  expect_identical(species_indexed_paths(derived,10L),character())
  expect_setequal(species_indexed_paths(derived,9L),known)
  expect_identical(derived$scenario$S,9L);expect_identical(derived$truth$datasettings$S,9L)
  expect_false('OTU_6' %in% colnames(derived$sim$data_list$OTU))
  slice <- function(x,margin,k) switch(as.character(margin),'0'=x[k],'1'=x[k,,drop=FALSE],'2'=x[,k,drop=FALSE])
  kept <- setdiff(1:10,6L)
  for(e in SPECIES_ELEMENTS) for(i in seq_along(kept)) {
    a <- slice(input[[e$path]],e$margin,kept[i]);b <- slice(derived[[e$path]],e$margin,i)
    expect_identical(serialize(b,NULL),serialize(a,NULL),label=paste(element_path(e$path),'species',kept[i]))
  }
  # The observation array column by column, with names.
  for(i in seq_along(kept)) {
    expect_identical(derived$sim$data_list$OTU[,i],input$sim$data_list$OTU[,kept[i]])
    expect_identical(colnames(derived$sim$data_list$OTU)[i],colnames(input$sim$data_list$OTU)[kept[i]])
  }
  # Nothing else changed: restoring the species elements and counts gives the input.
  restored <- derived
  for(e in SPECIES_ELEMENTS) restored[[e$path]] <- input[[e$path]]
  for(p in SPECIES_COUNTS) restored[[p]] <- input[[p]]
  expect_identical(restored,input)
  # The kept species' generating linear predictor is still B0 + X B + U L.
  info <- derived$sim$data_list$info;site <- info[!duplicated(info$Site),];site <- site[order(site$Site),]
  X <- scale(as.matrix(site[c('X_psi.EnvCov.1','X_psi.EnvCov.2')]))
  jp <- derived$sim$true_params$jsdmParams_true
  eta <- sweep(unname(X)%*%jp$B+jp$U%*%jp$L,2,jp$B0,'+')
  expect_lt(max(abs(eta-unname(jp$eta))),1e-10)
  expect_identical(drop_species(input),derived)
  bad <- input;bad$sim$extra <- 1:10
  expect_error(drop_species(bad),'differ from the known list')
})

test_that('the derived input is saved once and a different object is refused', {
  d <- tempfile();f <- file.path(d,'derived.rds')
  x <- list(a=1:3,b=letters)
  m1 <- save_or_check_derived(x,f);m2 <- save_or_check_derived(x,f)
  expect_identical(m1,m2);expect_identical(m1,unname(tools::md5sum(f)))
  expect_error(save_or_check_derived(list(a=1:4,b=letters),f),'differs')
  expect_false(dir.exists(paste0(f,'.lock')))
  unlink(d,recursive=TRUE)
})

test_that('the fixed-theta0 clone is the package sample_theta0 plus one statement', {
  need()
  expect_true(clone_matches_package())
  pkg <- package_function('sample_theta0')
  other <- pkg;body(other)[[3]] <- quote(z_all <- z[idx_z,,drop=FALSE])
  expect_false(clone_matches_package(other))
  found <- check_fixed_theta0_archive(file.path(SCRIPT_REPO,CLONE_ARCHIVE_DIR))
  h <- utils::read.csv(file.path(SCRIPT_REPO,CLONE_ARCHIVE_DIR,'hashes.csv'),colClasses='character')
  expect_identical(unname(found[h$item]),h$md5)
})

test_that('given the same RNG state the clone draws the package value for every other species', {
  need()
  pkg <- package_function('sample_theta0');fixed <- make_fixed_sample_theta0(6L,.038)
  for(seed in 1:5) for(b in c(20,100)) {
    set.seed(700+seed);n <- 40L;S <- 10L
    idx <- rep(seq_len(n),each=2L)
    z <- matrix(rbinom(n*S,1,.4),n,S);w <- matrix(rbinom(2*n*S,1,.3),2*n,S)
    set.seed(seed);a <- pkg(z,w,idx,1,b);ra <- .Random.seed
    set.seed(seed);f <- fixed(z,w,idx,1,b);rf <- .Random.seed
    expect_identical(f[-6],a[-6]);expect_identical(f[6],.038);expect_false(a[6]==.038)
    expect_identical(rf,ra)
  }
  expect_error(make_fixed_sample_theta0(6.5,.038),'whole number')
  expect_error(make_fixed_sample_theta0(6L,1),'in \\(0, 1\\)')
})

test_that('the fixed-theta0 fitter changes only the theta0 update', {
  need()
  run <- occJSDM::runOccJSDM;v <- fixed_theta0_value(input)
  expect_equal(v,.03800492,tolerance=1e-7)
  fitter <- make_fixed_theta0_fitter(6L,v)
  expect_identical(body(fitter),body(run));expect_identical(formals(fitter),formals(run))
  env <- environment(fitter)
  expect_identical(parent.env(env),asNamespace('occJSDM'))
  expect_identical(ls(env,all.names=TRUE),'sample_theta0')
  expect_identical(sort(ls(environment(env$sample_theta0),all.names=TRUE)),c('fixed_species','fixed_value'))
  expect_false(any(c('fixed_species','fixed_value') %in% all.names(body(run))))
  expect_true('sample_theta0' %in% all.names(body(run)))
  # The environment swap alone changes nothing: with the package's own update
  # the fit, warnings and final RNG state are identical to the extended path.
  base <- evaluate_seeded(expr,input,src$job,tiny,bindings=list(.fitter=make_run_fitter('extended',input,new.env())),seed=chain_seed(2))
  same <- evaluate_seeded(expr,input,src$job,tiny,bindings=list(.fitter=make_run_fitter('extended',input,new.env(),
    target=make_theta0_override_fitter(package_function('sample_theta0')))),seed=chain_seed(2))
  expect_identical(same$fit,base$fit);expect_identical(same$warnings,base$warnings)
  expect_identical(same$rng_final,base$rng_final)
  # With the clone, species 6 is held at the generating value and the others move.
  held <- evaluate_seeded(expr,input,src$job,tiny,bindings=list(.fitter=make_run_fitter('b',input,new.env())),seed=chain_seed(2))
  th <- held$fit$results_output$theta0_output
  expect_true(all(th[6,,1]==v))
  expect_true(all(apply(th[-6,,1],1,function(x) length(unique(x))==5L)))
})

test_that('the variant fitters pass listPriors exactly as each variant states', {
  capture <- function(...,listPriors) list(args=list(...),listPriors=listPriors)
  archived <- list(a_q=1,b_q=20)
  for(v in c('extended','c')) {
    rec <- new.env();out <- make_run_fitter(v,NULL,rec,target=capture)('data',threshold=1,listPriors=archived)
    expect_identical(out$listPriors,archived);expect_identical(rec$listPriors,archived)
    expect_identical(out$args,list('data',threshold=1))
  }
  rec <- new.env();out <- make_run_fitter('a',NULL,rec,target=capture)('data',listPriors=archived)
  expect_identical(out$listPriors,list(a_q=1,b_q=20,a_theta0=1,b_theta0=100))
  expect_identical(rec$listPriors,out$listPriors)
  expect_error(make_run_fitter('a',NULL,new.env(),target=capture)('data',listPriors=list(b_theta0=5)),'already set b_theta0')
})

test_that('the fitting statement is the md5-checked pr11 statement with only the fitting function replaced', {
  need()
  expect_identical(runner$md5,unname(tools::md5sum(file.path(SCRIPT_REPO,PR11_RUNNER))))
  back <- replace_in_call(expr,as.name('.fitter'),quote(occJSDM::runOccJSDM))
  expect_identical(back,runner$expr)
  expect_false(any(grepl('runOccJSDM',deparse(expr),fixed=TRUE)))
  expect_identical(src$job$priors,list(a_q=1,b_q=20))
})

test_that('a saved fit with different metadata is refused and an identical one resumes', {
  d <- tempfile();dir.create(d);dest <- file.path(d,'chain-01-fit.rds')
  meta <- list(run='a',chain=1L,seed=20261002L,mcmc=SCHEDULES$diagnostic)
  expect_identical(existing_fit_status(dest,meta),'new')
  saveRDS(c(meta,list(fit='x')),dest)
  expect_identical(existing_fit_status(dest,meta),'resume')
  changed <- meta;changed$seed <- 20261003L
  expect_error(existing_fit_status(dest,changed),'differing fields: seed')
  lock <- acquire_fit_lock(dest)
  expect_error(acquire_fit_lock(dest),'already exists')
  expect_true(release_fit_lock(lock))
  unlink(d,recursive=TRUE)
})

test_that('the library fingerprint check refuses a changed build', {
  need()
  csv <- file.path(SCRIPT_REPO,FINGERPRINT_FILE)
  expect_identical(check_library_fingerprint(study,csv),readLines(file.path(study,'source-revision.txt')))
  fp <- utils::read.csv(csv,colClasses='character');fp$md5[fp$file=='library/occJSDM/libs/occJSDM.so'] <- strrep('0',32)
  bad <- tempfile(fileext='.csv');utils::write.csv(fp,bad,row.names=FALSE)
  expect_error(check_library_fingerprint(study,bad),'changed library/occJSDM/libs/occJSDM.so')
  unlink(bad)
})

test_that('launcher run lists expand to single-chain jobs in order', {
  q <- parse_run_list('extended:1-16,a:1-8,b:1-8,c:1-8')
  expect_identical(nrow(q),40L)
  expect_identical(q$run[1:16],rep('extended',16));expect_identical(q$chain[1:16],1:16)
  expect_identical(unique(q$run),c('extended','a','b','c'))
  p <- parse_run_list('pilot-extended:1,pilot-a:1,pilot-b:1,pilot-c:1')
  expect_identical(p$run,rep('pilot',4));expect_identical(p$variant,c('extended','a','b','c'))
  expect_identical(parse_run_list('a:3,a:5-6')$chain,c(3L,5L,6L))
  expect_error(parse_run_list('a:1-9'),'chain')
  expect_error(parse_run_list('a:1,a:1'),'twice')
  expect_error(parse_run_list('z:1'),'Unknown run')
  expect_error(parse_run_list('a'),'RUN:CHAINS')
  args <- job_args(p[3,],'/r','/s','/i','/a')
  expect_identical(args,c('--repo=/r','--study=/s','--inputs-root=/i','--archives=/a','--run=pilot','--variant=b','--chain=1'))
  expect_identical(job_args(q[17,],'/r','/s','/i','/a')[5:6],c('--run=a','--chain=1'))
})

test_that('the launcher refuses to start while any fit lock exists', {
  d <- tempfile();dir.create(file.path(d,'fits','extended','chain-03.lock'),recursive=TRUE)
  expect_identical(fit_locks(d),file.path(d,'fits','extended','chain-03.lock'))
  q <- parse_run_list('a:1')
  out <- capture.output(res <- execute_jobs(q,d,max_procs=1L,command_fn=function(job) stop('must not start'),poll=.1))
  expect_identical(res$status,1L);expect_true(any(grepl('Refusing to start',out)))
  unlink(d,recursive=TRUE)
})

test_that('the launcher runs at most max_procs processes and records exit codes and a DONE file', {
  d <- tempfile();dir.create(d)
  q <- parse_run_list('a:1-5')
  # Chain 4 fails without a fit; the others save an (empty) fit file and succeed.
  fake <- function(job) list(command=file.path(R.home('bin'),'Rscript'),
    args=c('-e',sprintf('dir.create(dirname("%s"),recursive=TRUE,showWarnings=FALSE); if(%d!=4L) file.create("%s"); Sys.sleep(.6); cat("chain %d\\n"); quit(status=%d)',
      job$fit,job$chain,job$fit,job$chain,as.integer(job$chain==4L))))
  out <- capture.output(res <- execute_jobs(q,d,max_procs=2L,command_fn=fake,poll=.1,stamp='test'))
  expect_identical(res$status,1L);expect_identical(res$max_concurrent,2L)
  s <- utils::read.csv(file.path(res$dir,'summary.csv'),stringsAsFactors=FALSE)
  expect_identical(s$exit_code,c(0L,0L,0L,1L,0L))
  expect_identical(s$status,c(rep('completed',3),'failed','completed'))
  expect_identical(readLines(file.path(res$dir,'a-chain-04.exit')),'1')
  done <- readLines(file.path(res$dir,'DONE'))
  expect_true(all(c('exit_status: 1','completed: 4','failed: 1') %in% done))
  expect_true(any(grepl('chain 2',readLines(file.path(res$dir,'a-chain-02.log')))))
  expect_false(dir.exists(file.path(d,'logs','launcher.lock')))
  unlink(d,recursive=TRUE)
})
