library(testthat)
source('jobs.R');source('verify-helpers.R');source('launch.R');source('flags.R');source('select.R')

# Research tests for the phase launcher (launch.R), the frozen flag rules
# (flags.R) and the selection script (select.R). Run from this directory with
# `Rscript test-launch.R`. Tests that need the saved archives or the study
# archive are skipped when absent; set OCCJSDM_ARCHIVES to relocate them.
# The suite never writes into the live study archive and does not depend on
# its fits: command-line tests run against a temporary copy of its installed
# library and exported source, and the live study is only read (pilot fits).
here <- normalizePath('.')
repo <- normalizePath('../../..')
archives <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
inputs_root <- file.path(archives,'intercept-prior-inputs')
study <- file.path(archives,'intercept-prior-20260929')
have_archives <- all(dir.exists(file.path(archives,c(unname(ARCHIVES),'intercept-prior-inputs'))))
have_study <- have_archives && dir.exists(file.path(study,'library/occJSDM'))
need_archives <- function() if(!have_archives) skip('saved study archives not available')
need_study <- function() if(!have_study) skip('study archive not available')
rscript <- file.path(R.home('bin'),'Rscript')
run_cli <- function(script,args) {
  out <- suppressWarnings(system2(rscript,c(shQuote(script),shQuote(args)),stdout=TRUE,stderr=TRUE))
  list(status=attr(out,'status') %||% 0L,output=paste(out,collapse='\n'))
}
# A study holding only the installed library, exported source and revision of
# the live study (enough for the fingerprint check), with no fits or logs.
study_copy <- local({cache <- NULL;function() {
  if(is.null(cache)) {
    d <- tempfile('studycopy');dir.create(d)
    for(x in c('source','library')) stopifnot(file.copy(file.path(study,x),d,recursive=TRUE))
    stopifnot(file.copy(file.path(study,'source-revision.txt'),d));cache <<- normalizePath(d)
  }
  cache
}})
launch_args <- function(repo_dir,study_dir,...) c(paste0('--repo=',repo_dir),paste0('--study=',study_dir),
  paste0('--inputs-root=',inputs_root),paste0('--archives=',archives),'--phases=A1,A2','--sds=2,3,5','--max-procs=8',...)
# A copy of the files a launcher checkout needs, at another path.
copy_checkout <- function() {
  root <- tempfile('checkout');d <- file.path(root,'dev/simstudy/occupancy-intercept-prior')
  dir.create(file.path(d,'results'),recursive=TRUE)
  for(f in c('launch.R','jobs.R','verify-helpers.R','run.R','results/library-fingerprint.csv'))
    stopifnot(file.copy(file.path(here,f),file.path(d,f)))
  normalizePath(root)
}

base <- c('--repo=/r','--study=/s/results/study','--inputs-root=/in','--phases=A1,A2','--sds=2,3,5','--max-procs=8')
without <- function(name) base[!startsWith(base,paste0('--',name,'='))]

test_that('launcher arguments parse, default and validate', {
  a <- parse_launch_args(base)
  expect_identical(a$repo,'/r');expect_identical(a$study,'/s/results/study');expect_identical(a$inputs_root,'/in')
  expect_identical(a$phases,c('A1','A2'));expect_identical(a$sds,c(2,3,5));expect_identical(a$max_procs,8L)
  expect_false(a$dry_run);expect_null(a$only);expect_null(a$schedule_override);expect_null(a$queue_out)
  expect_identical(a$archives,'/s/results');expect_identical(a$poll,5)
  expect_true(parse_launch_args(c(base,'--dry-run'))$dry_run)
  b <- parse_launch_args(c(base,'--only=/k.txt','--schedule-override=long','--archives=/a','--queue-out=/q.csv','--poll-seconds=0.5'))
  expect_identical(b$only,'/k.txt');expect_identical(b$schedule_override,'long');expect_identical(b$archives,'/a')
  expect_identical(b$queue_out,'/q.csv');expect_identical(b$poll,0.5)
  expect_identical(parse_launch_args(c(without('phases'),'--phases=B'))$phases,'B')
  for(name in c('repo','study','inputs-root','phases','sds','max-procs'))
    expect_error(parse_launch_args(without(name)),paste0('Missing --',name))
  for(bad in c('A3','A1,A1','','A1,')) expect_error(parse_launch_args(c(without('phases'),paste0('--phases=',bad))),'--phases')
  # The control arm (SD 1) is reused, never refitted by the launcher.
  for(bad in c('1','4','2,2','2.5','')) expect_error(parse_launch_args(c(without('sds'),paste0('--sds=',bad))),'--sds')
  for(bad in c('0','9','x','2.5')) expect_error(parse_launch_args(c(without('max-procs'),paste0('--max-procs=',bad))),'--max-procs')
  expect_error(parse_launch_args(c(base,'--schedule-override=long')),'requires --only')
  expect_error(parse_launch_args(c(base,'--only=/k','--schedule-override=initial')),'--schedule-override')
  expect_error(parse_launch_args(c(base,'--poll-seconds=0')),'--poll-seconds')
  expect_error(parse_launch_args(c(base,'--dry-run=yes')),'Unknown option')
  expect_error(parse_launch_args(c(base,'--dry-run','--dry-run')),'Repeated')
  expect_error(parse_launch_args(c(base,'--workers=2')),'Unknown option')
  expect_error(parse_launch_args(c(base,'stray')),'Malformed')
})

test_that('the launcher invokes only the run.R of its own checkout, with the frozen fit scripts', {
  expect_identical(unname(FROZEN_FIT_MD5[c(RUN_SCRIPT,JOBS_SCRIPT,FINGERPRINT_FILE)]),
    c('e111cff6f079dcc783c7a2b32fa28eaf','d6dc59ff47a25de08fad60b8bde827d1','22d6191da2347a26009bdf15b6c31a4d'))
  expect_identical(unname(tools::md5sum(file.path(repo,names(FROZEN_FIT_MD5)))),unname(FROZEN_FIT_MD5))
  expect_identical(check_launch_repo(repo,file.path(here,'launch.R')),normalizePath(file.path(repo,RUN_SCRIPT)))
  expect_identical(check_launch_repo(repo),normalizePath(file.path(repo,RUN_SCRIPT)))
  expect_error(check_launch_repo(repo,file.path(tempdir(),'launch.R')),'is not the launcher of --repo')
  other <- copy_checkout()
  expect_identical(check_launch_repo(other,file.path(other,LAUNCH_SCRIPT)),file.path(other,RUN_SCRIPT))
  expect_error(check_launch_repo(other,file.path(here,'launch.R')),'is not the launcher of --repo')
  run <- file.path(other,RUN_SCRIPT);writeLines(c(readLines(run),'# edited'),run)
  expect_error(check_launch_repo(other),'Frozen fit script changed: dev/simstudy/occupancy-intercept-prior/run.R')
  unlink(run);expect_error(check_launch_repo(other),'run.R not found')
  # None of the new tooling enters a fit's script hashes.
  for(phase in PHASES) expect_false(any(grepl('launch|select|flags|check-flags',runner_script_files(phase))))
})

test_that('planning estimates order jobs longest first and simulate the concurrency cap', {
  expect_gt(expected_seconds('A2','long','range4-rep01-binary-k100'),expected_seconds('B','long','design-qfar_K6-sites300-01'))
  expect_gt(expected_seconds('B','initial','design-qfar_K6-sites300-01'),expected_seconds('A1','initial','jsdm-n0300-01'))
  expect_gt(expected_seconds('A1','initial','jsdm-n0300-01'),expected_seconds('A1','initial','jsdm-n0100-01'))
  expect_gt(expected_seconds('A1','long','jsdm-n0100-01'),expected_seconds('A1','initial','jsdm-n0100-01'))
  expect_identical(estimated_wall_seconds(c(10,10,10,10),2L),20)
  expect_identical(estimated_wall_seconds(c(9,5,4,3),2L),12)
  expect_identical(estimated_wall_seconds(c(3,3),8L),3)
  q <- order_queue(queue_frame(c('A1','A2','A1'),c(2,2,3),c('initial','long','initial'),
    c('jsdm-n0100-01','range4-rep01-binary-k100','jsdm-n0300-01'),c('initial','long','initial'),'/s'))
  expect_identical(q$key,c('range4-rep01-binary-k100','jsdm-n0300-01','jsdm-n0100-01'))
  expect_identical(q$order,1:3)
  expect_identical(q$fit[1],fit_path('/s','A2',2,'long','range4-rep01-binary-k100'))
  expect_identical(q$lock,fit_lock_path(q$fit))
})

test_that('each job runs the checkout\'s run.R with one worker and one key, in run.R\'s own interface', {
  run <- check_launch_repo(repo)
  q <- order_queue(queue_frame(c('A2','A1'),c(5,2),c('long','initial'),c('range8-rep03-binary-k100','jsdm-n0300-07'),
    c('long','initial'),'/study'))
  cmd <- run_r_command(run,repo,'/study','/inputs','/archives')
  for(i in 1:2) {
    x <- cmd(as.list(q[i,]))
    expect_identical(x$command,file.path(R.home('bin'),'Rscript'));expect_identical(x$args[1],run)
    a <- parse_run_args(x$args[-1])
    expect_identical(a$repo,repo);expect_identical(a$study,'/study');expect_identical(a$phase,q$phase[i])
    expect_identical(a$sd,q$sd[i]);expect_identical(a$schedule,q$schedule[i]);expect_identical(a$workers,1L)
    expect_identical(a$keys,q$key[i]);expect_identical(a$inputs_root,'/inputs');expect_identical(a$archives,'/archives')
    expect_identical(fit_path(a$study,a$phase,a$sd,a$schedule,a$keys),q$fit[i])
  }
})

test_that('the phase A queue is 60 A1 initial and 27 A2 long jobs at the control-selected schedules', {
  need_archives()
  q <- build_queue(c('A1','A2'),c(2,3,5),'/study',archives)
  expect_identical(nrow(q),87L);expect_identical(q$order,seq_len(87L))
  counts <- queue_counts(q)
  expect_identical(counts$phase,rep(c('A1','A2'),each=3L));expect_identical(counts$sd,rep(c(2,3,5),2L))
  expect_identical(counts$schedule,rep(c('initial','long'),each=3L));expect_identical(counts$jobs,rep(c(20L,9L),each=3L))
  for(phase in c('A1','A2')) {
    cs <- control_schedules(phase,archives)
    for(sd in c(2,3,5)) {
      r <- q[q$phase==phase & q$sd==sd,]
      expect_setequal(r$key,phase_keys(phase))
      expect_identical(r$schedule,cs$schedule[match(r$key,cs$key)]);expect_identical(r$control_schedule,r$schedule)
      expect_identical(r$fit,fit_path('/study',phase,sd,r$schedule,r$key))
    }
  }
  # Longest first: every A2 long fit before any A1 fit, 300-site before 100-site fits.
  expect_true(all(diff(q$expected_seconds)<=0))
  expect_lt(max(which(q$phase=='A2')),min(which(q$phase=='A1')))
  expect_lt(max(grep('-n0300-',q$key)),min(grep('-n0100-',q$key)))
  b <- build_queue('B',2,'/study',archives)
  expect_identical(as.vector(table(b$schedule)[c('initial','long')]),c(8L,12L))
  # A job list (select.R's long-keys.txt) restricts the queue; the override sets the schedule.
  f <- tempfile(fileext='.txt')
  write.csv(data.frame(phase=c('A1','A1'),sd=c(2,5),key=c('jsdm-n0100-03','jsdm-n0300-07')),f,row.names=FALSE)
  only <- read_only_jobs(f)
  r <- build_queue(c('A1','A2'),c(2,3,5),'/study',archives,only,'long')
  expect_identical(r$key,c('jsdm-n0300-07','jsdm-n0100-03'));expect_identical(r$sd,c(5,2))
  expect_identical(r$schedule,c('long','long'));expect_identical(r$control_schedule,c('initial','initial'))
  expect_identical(build_queue('A1',c(2,5),'/study',archives,only)$schedule,c('initial','initial'))
  bad <- function(phase,sd,key) {g <- tempfile();write.csv(data.frame(phase=phase,sd=sd,key=key),g,row.names=FALSE);read_only_jobs(g)}
  expect_error(build_queue(c('A1','A2'),c(2,3,5),'/study',archives,bad('A2',2,'range4-rep01-binary-k100'),'long'),
    'already at the long schedule')
  expect_error(build_queue('A1',2,'/study',archives,bad('A1',2,'jsdm-n9999-01')),'not a phase A1 key')
  expect_error(build_queue('A1',2,'/study',archives,bad('A1',3,'jsdm-n0100-01')),'not in --sds')
  expect_error(build_queue('A1',2,'/study',archives,bad('A2',2,'range4-rep01-binary-k100')),'not in --phases')
  expect_error(build_queue('A1',2,'/study',archives,bad(c('A1','A1'),2,c('jsdm-n0100-01','jsdm-n0100-01'))),'Duplicated')
  g <- tempfile();writeLines('phase,key\nA1,jsdm-n0100-01',g);expect_error(read_only_jobs(g),'phase, sd, key')
  # Amendment 1 (R23): A1 and A2 are launched separately; their queues are the two halves of the phase A queue.
  qa1 <- build_queue('A1',c(2,3,5),'/study',archives);qa2 <- build_queue('A2',c(2,3,5),'/study',archives)
  expect_identical(c(nrow(qa1),nrow(qa2)),c(60L,27L));expect_true(all(qa1$schedule=='initial'));expect_true(all(qa2$schedule=='long'))
  id <- function(x) paste(x$phase,x$sd,x$schedule,x$key)
  expect_setequal(c(id(qa1),id(qa2)),id(q))
  committed <- lapply(c('phase-a-queue.csv','phase-a1-queue.csv','phase-a2-queue.csv'),function(f)
    utils::read.csv(file.path(here,'results',f),stringsAsFactors=FALSE))
  expect_identical(c(nrow(committed[[2]]),nrow(committed[[3]])),c(60L,27L))
  expect_setequal(c(id(committed[[2]]),id(committed[[3]])),id(committed[[1]]))
  expect_identical(committed[[2]]$key,qa1$key);expect_identical(committed[[3]]$key,qa2$key)
  none <- bad(character(),numeric(),character());expect_identical(nrow(none),0L)
  for(override in list(NULL,'long')) {
    e <- build_queue(c('A1','A2'),c(2,3,5),'/study',archives,none,override)
    expect_identical(nrow(e),0L);expect_true(all(c(QUEUE_COLUMNS[QUEUE_COLUMNS!='fit_exists'],'lock') %in% names(e)))
    expect_identical(nrow(queue_counts(e)),0L)
  }
})

test_that('a dry run prints the queue and counts, writes the queue CSV and changes nothing in the study', {
  need_study();sc <- study_copy()
  before <- list.files(sc,recursive=TRUE,all.files=TRUE,include.dirs=TRUE)
  out <- tempfile(fileext='.csv')
  r <- run_cli(file.path(here,'launch.R'),launch_args(repo,sc,'--dry-run',paste0('--queue-out=',out)))
  expect_identical(r$status,0L)
  expect_match(r$output,'87 jobs');expect_match(r$output,'A1 +sd2 +initial +20');expect_match(r$output,'A2 +sd5 +long +9')
  expect_match(r$output,'Dry run: nothing started')
  q <- read.csv(out,stringsAsFactors=FALSE)
  expect_identical(names(q),QUEUE_COLUMNS);expect_identical(nrow(q),87L)
  expect_identical(q$fit[1],'fits/A2/sd2/long/range4-rep01-binary-k100-fit.rds')
  expect_false(any(q$fit_exists))
  expect_identical(list.files(sc,recursive=TRUE,all.files=TRUE,include.dirs=TRUE),before)
  # An empty job list (select.R found nothing to repeat) is an empty queue, not an error.
  empty <- tempfile(fileext='.txt');writeLines('"phase","sd","key"',empty)
  r <- run_cli(file.path(here,'launch.R'),launch_args(repo,sc,paste0('--only=',empty),'--schedule-override=long'))
  expect_identical(r$status,0L);expect_match(r$output,'Queue is empty; nothing to run')
  expect_false(dir.exists(file.path(sc,'logs')))
})

test_that('the command line refuses another checkout, changed frozen scripts and existing locks', {
  need_study();sc <- study_copy()
  other <- copy_checkout()
  r <- run_cli(file.path(here,'launch.R'),launch_args(other,sc,'--dry-run'))
  expect_identical(r$status,1L);expect_match(r$output,'is not the launcher of --repo')
  r <- run_cli(file.path(other,LAUNCH_SCRIPT),launch_args(other,sc,'--dry-run'))
  expect_identical(r$status,0L);expect_match(r$output,'Dry run: nothing started')
  run <- file.path(other,RUN_SCRIPT);writeLines(c(readLines(run),'# edited'),run)
  r <- run_cli(file.path(other,LAUNCH_SCRIPT),launch_args(other,sc,'--dry-run'))
  expect_identical(r$status,1L);expect_match(r$output,'Frozen fit script changed')
  unlink(run)
  r <- run_cli(file.path(other,LAUNCH_SCRIPT),launch_args(other,sc,'--dry-run'))
  expect_identical(r$status,1L);expect_match(r$output,'run.R not found')
  fake <- tempfile('study');dir.create(fake);fake <- normalizePath(fake)
  lock <- fit_lock_path(fit_path(fake,'A1',3,'initial','jsdm-n0300-07'));dir.create(lock,recursive=TRUE)
  for(extra in list('--dry-run',character())) {
    r <- run_cli(file.path(here,'launch.R'),launch_args(repo,fake,extra))
    expect_identical(r$status,1L);expect_match(r$output,'Refusing to start: 1 fit lock');expect_match(r$output,lock,fixed=TRUE)
  }
  expect_false(dir.exists(file.path(fake,'logs')));expect_true(dir.exists(lock))
})

# Stub fits: each touches a marker while it runs and records how many markers exist.
stub_command <- function(ts,behaviour=list(),sleep=0.5) function(job) {
  mode <- behaviour[[job$key]] %||% 'ok'
  run <- shQuote(file.path(ts,'running',job$key));dir <- shQuote(dirname(job$fit));fit <- shQuote(job$fit)
  q <- shQuote(file.path(dirname(job$fit),paste0(job$key,'-fit.QUARANTINE-20260929000000.rds')))
  pre <- sprintf('touch %s; ls %s | wc -l >> %s; sleep %s; rm %s;',run,shQuote(file.path(ts,'running')),
    shQuote(file.path(ts,'counts.txt')),sleep,run)
  body <- switch(mode,ok=sprintf('mkdir -p %s; touch %s; echo complete',dir,fit),
    resumed=sprintf('mkdir -p %s; touch %s; echo "t %s A1 sd2 initial already complete with identical settings; resumed"',dir,fit,job$key),
    fail='echo "ERROR: boom"; exit 1',quarantine=sprintf('mkdir -p %s; touch %s; exit 1',dir,q),nofit='echo done')
  list(command='/bin/sh',args=c('-c',paste(pre,body)))
}
stub_study <- function() {ts <- tempfile('study');dir.create(file.path(ts,'running'),recursive=TRUE);normalizePath(ts)}
stub_queue <- function(ts,keys) order_queue(queue_frame('A1',2,'initial',keys,'initial',ts))

test_that('no more than max-procs jobs run at once, each with its own log, and a clean queue exits 0', {
  ts <- stub_study();keys <- sprintf('jsdm-n0100-%02d',1:7)
  res <- execute_queue(stub_queue(ts,keys),ts,3L,stub_command(ts),poll=0.05,stamp='test')
  expect_identical(res$status,0L);expect_identical(res$max_concurrent,3L)
  counts <- as.integer(trimws(readLines(file.path(ts,'counts.txt'))))
  expect_length(counts,7L);expect_identical(max(counts),3L)
  expect_identical(res$dir,file.path(ts,'logs','launch-test'))
  s <- read.csv(file.path(res$dir,'summary.csv'),stringsAsFactors=FALSE)
  expect_identical(names(s),SUMMARY_COLUMNS);expect_identical(s$key,keys);expect_identical(s$status,rep('completed',7L))
  expect_identical(s$exit_code,rep(0L,7L));expect_true(all(s$elapsed_seconds>=0.4))
  expect_true(all(file.exists(file.path(res$dir,s$log))))
  done <- readLines(file.path(res$dir,'DONE'))
  expect_true('exit_status: 0' %in% done);expect_true('completed: 7' %in% done);expect_true('max_concurrent: 3' %in% done)
  expect_true(file.exists(file.path(res$dir,'queue.csv')));expect_true(file.exists(file.path(res$dir,'launcher.txt')))
  expect_false(dir.exists(file.path(ts,'logs','launcher.lock')))
  one <- execute_queue(stub_queue(ts,sprintf('jsdm-n0300-%02d',1:2)),ts,1L,stub_command(ts,sleep=0.1),poll=0.05,stamp='test')
  expect_identical(one$dir,paste0(file.path(ts,'logs','launch-test'),'-',Sys.getpid()))
  expect_identical(one$max_concurrent,1L)
})

test_that('statuses are completed, resumed, failed or quarantined, and any failure exits non-zero', {
  ts <- stub_study();keys <- sprintf('jsdm-n0100-%02d',1:5)
  modes <- list(`jsdm-n0100-02`='resumed',`jsdm-n0100-03`='fail',`jsdm-n0100-04`='quarantine',`jsdm-n0100-05`='nofit')
  res <- execute_queue(stub_queue(ts,keys),ts,8L,stub_command(ts,modes,sleep=0.1),poll=0.05,stamp='x')
  expect_identical(res$status,1L)
  s <- read.csv(file.path(res$dir,'summary.csv'),stringsAsFactors=FALSE)
  expect_identical(s$status,c('completed','resumed','failed','quarantined','failed'))
  expect_identical(s$exit_code,c(0L,0L,1L,1L,0L))
  done <- readLines(file.path(res$dir,'DONE'))
  for(line in c('exit_status: 1','jobs: 5','completed: 1','resumed: 1','failed: 2','quarantined: 1','not_run: 0'))
    expect_true(line %in% done,info=line)
  expect_identical(job_status(0L,'x resumed',tempfile(),character(),character()),'failed')
  expect_identical(job_status(NA_integer_,'',tempfile(),character(),character()),'failed')
})

test_that('an existing fit lock or a running launcher refuses the whole queue before anything starts', {
  ts <- stub_study();q <- stub_queue(ts,sprintf('jsdm-n0100-%02d',1:3))
  dir.create(q$lock[2],recursive=TRUE)
  expect_identical(existing_locks(q),q$lock[2])
  res <- execute_queue(q,ts,2L,stub_command(ts),poll=0.05,stamp='locked')
  expect_identical(res$status,1L);expect_null(res$dir)
  expect_false(file.exists(file.path(ts,'counts.txt')));expect_false(dir.exists(file.path(ts,'logs','launch-locked')))
  unlink(q$lock[2],recursive=TRUE)
  dir.create(file.path(ts,'logs','launcher.lock'),recursive=TRUE)
  res <- execute_queue(q,ts,2L,stub_command(ts),poll=0.05,stamp='second')
  expect_identical(res$status,1L);expect_false(file.exists(file.path(ts,'counts.txt')))
  expect_true(dir.exists(file.path(ts,'logs','launcher.lock')))
})

test_that('the launcher process exits non-zero when a job fails and zero when all succeed', {
  ts <- stub_study()
  driver <- function(mode) {
    f <- tempfile(fileext='.R')
    writeLines(c(sprintf('setwd(%s)',deparse(here)),"source('jobs.R');source('verify-helpers.R');source('launch.R')",
      sprintf("ts <- %s",deparse(ts)),
      sprintf("cmd <- function(job) list(command='/bin/sh',args=c('-c',%s))",
        deparse(if(mode=='fail') 'exit 3' else sprintf('mkdir -p "$(dirname %s)"; touch %s','"$FIT"','"$FIT"'))),
      "q <- order_queue(queue_frame('A1',2,'initial','jsdm-n0100-01','initial',ts))",
      "Sys.setenv(FIT=q$fit)",
      sprintf("res <- execute_queue(q,ts,1L,cmd,poll=0.05,stamp=%s)",deparse(mode)),"quit(save='no',status=res$status)"),f)
    system2(rscript,shQuote(f),stdout=FALSE,stderr=FALSE)
  }
  expect_identical(driver('fail'),1L)
  expect_identical(read.csv(file.path(ts,'logs/launch-fail/summary.csv'))$exit_code,3L)
  expect_identical(driver('ok'),0L)
})

test_that('the pr11 flag rule is the archived select.R expression, with its reasons', {
  expr <- pr11_flag_expression(repo)
  expect_true(is.call(expr))
  grid <- expand.grid(w=c(0L,2L),g=c(1.01,1.05,1.0500001,NA,Inf),e=c(1.01,1.06,NaN),u=c(0L,1L))
  for(i in seq_len(nrow(grid))) {
    x <- grid[i,];d <- list(max_group_rhat=x$g,max_element_rhat=x$e,unresolved_rhat=x$u)
    expected <- x$w>0 || x$u>0 || !is.finite(x$g) || !is.finite(x$e) || x$g>1.05 || x$e>1.05
    r <- pr11_rule(rep('w',x$w),d,expr)
    expect_identical(r$flagged,expected,info=paste(unlist(x),collapse=' '));expect_identical(length(r$reasons)>0L,expected)
  }
  r <- pr11_rule(c('a','b'),list(max_group_rhat=1.07,max_element_rhat=NA,unresolved_rhat=3L),expr)
  expect_identical(r$reasons,c('2 fitting warnings','max group Rhat 1.0700 > 1.05','max element Rhat unavailable',
    '3 non-finite or nonpositive Rhat'))
  # Diagnostics are combined exactly as current-main-recheck/run.R:79-81 did.
  sm <- function(x) if(!length(x) || all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)
  d <- pr11_diagnostics(list(groups=c(1.01,NA),elements=c(1.02,-1),additional=c(1.2,Inf)),sm)
  expect_identical(d,list(max_group_rhat=1.01,max_element_rhat=Inf,unresolved_rhat=3L))
  expect_identical(pr11_diagnostics(list(groups=numeric(),elements=1,additional=numeric()),sm)$max_group_rhat,NA_real_)
})

test_that('the spatial flag rule is the archived robust rule, which ignores native warnings', {
  need_archives()
  sc <- load_flag_scorers(repo,archives)
  ok <- list(groups=data.frame(metric=c('occupancy','occupancy','spatial_sd','range'),group=c('all','low','all','all'),
      rhat=c(1.01,1.01,1.01,NA),ess_mean=c(500,500,NA,NA)),
    species=data.frame(rhat=c(1.01,1.02),ess_mean=c(500,400)),
    elements=data.frame(metric=c('occupancy','intercept','range','spatial_sd'),rhat=c(1.01,1.01,NA,1.01)),
    spatial=list(diagnostics=data.frame(quantity=c('amplitude','field_mean_1'),rhat=1.01,ess_bulk=500,ess_median=500,
      ess_q025=500,ess_q975=500),field_diagnostics=data.frame(quantity=c('field_1','field_2'),rhat=1.01,ess_bulk=500,
      ess_median=500,ess_q025=500,ess_q975=500)))
  expect_identical(spatial_rule(ok,'Convergence issue in B_output',sc)$reasons,character())
  change <- function(table,column,row,value,sub=NULL) {x <- ok;if(is.null(sub)) x[[table]][[column]][row] <- value else
    x[[table]][[sub]][[column]][row] <- value;spatial_rule(x,character(),sc)}
  expect_identical(change('groups','rhat',2,1.06)$reasons,'group Rhat > 1.05')
  expect_identical(change('groups','ess_mean',1,99)$reasons,'occupancy group ESS < 100')
  expect_identical(change('species','ess_mean',2,99)$reasons,'occupancy species ESS < 100')
  expect_identical(change('species','rhat',1,1.2)$reasons,'species Rhat > 1.05')
  expect_identical(change('elements','rhat',2,NA)$reasons,'non-range element Rhat unavailable')
  expect_identical(change('elements','rhat',4,1.3)$reasons,'element Rhat > 1.05')
  expect_identical(change('spatial','ess_q025',1,50,'diagnostics')$reasons,'1 spatial summary traces fail rank/quantile checks')
  expect_identical(change('spatial','rhat',2,NA,'field_diagnostics')$reasons,'1 pointwise spatial fields fail rank/quantile checks')
  expect_true(change('spatial','rhat',2,NA,'field_diagnostics')$flagged)
})

test_that('the flag rules run on the three pilot fits from saved draws and report no occupancy error', {
  need_study()
  sc <- load_flag_scorers(repo,archives)
  for(phase in PHASES) {
    key <- if(phase=='A2') 'range6-rep01-binary-k100' else phase_keys(phase)[1]
    saved <- readRDS(fit_path(study,phase,3,'pilot',key))
    spec <- phase_jobs(phase,archives,inputs_root,key)[[1]]
    input <- readRDS(checked_input(spec$input_file,spec$input_md5))
    row <- fit_flags(phase,saved$fit,saved$warnings,input,sc)
    expect_identical(names(row),FLAG_COLUMNS);expect_identical(nrow(row),1L)
    expect_identical(row$rule,if(phase=='A2') 'spatial-robust-v1' else 'pr11')
    expect_identical(row$warnings,length(saved$warnings))
    expect_true(is.finite(row$max_group_rhat) && is.finite(row$max_element_rhat))
    # 200 retained draws per chain: every pilot is flagged. Mechanics only; never scored.
    expect_true(row$flagged);expect_true(nzchar(row$reasons))
    if(phase=='A2') expect_true(row$spatial_trace_flags>0L) else expect_true(is.na(row$spatial_trace_flags))
  }
  # The flag code never calls a scorer that computes estimation error.
  code <- c(readLines(file.path(here,'flags.R')),readLines(file.path(here,'select.R')))
  for(f in c('score_fit','score_spatial_fit','score_draw_block','score_current_fit','error_metrics',
    'make_design_scorer','load_scoring','add_group','score_amplitude_fit','summarise_field_draws'))
    expect_false(any(grepl(paste0('\\b',f,'\\('),code,perl=TRUE)),info=f)
})

test_that('control flags reproduce the archived diagnostics under each rule', {
  need_archives()
  sc <- load_flag_scorers(repo,archives)
  pr11 <- read.csv(file.path(archive_dir(archives,'pr11'),'long-selection.csv'),stringsAsFactors=FALSE)
  for(x in list(c('A1','jsdm-n0100-01'),c('B','design-qnear_K6-sites300-01'))) {
    saved <- readRDS(file.path(archive_dir(archives,'pr11'),'initial',paste0(x[2],'-fit.rds')))
    spec <- phase_jobs(x[1],archives,inputs_root,x[2])[[1]]
    row <- fit_flags(x[1],saved$fit,saved$warnings,readRDS(checked_input(spec$input_file,spec$input_md5)),sc)
    a <- pr11[pr11$key==x[2],]
    expect_identical(row$warnings,as.integer(a$warnings));expect_identical(row$unresolved_rhat,as.integer(a$unresolved_rhat))
    expect_lt(abs(row$max_group_rhat-a$max_group_rhat),1e-12);expect_lt(abs(row$max_element_rhat-a$max_element_rhat),1e-12)
    expect_identical(row$flagged,a$current_flag)
  }
  # An A2 inverse-gamma initial control whose archived robust reasons are not empty.
  key <- 'range6-rep03-binary-k100'
  saved <- readRDS(file.path(archive_dir(archives,'targeted'),'initial',paste0(key,'-fit.rds')))
  spec <- phase_jobs('A2',archives,inputs_root,key)[[1]]
  row <- fit_flags('A2',saved$fit,saved$warnings,readRDS(checked_input(spec$input_file,spec$input_md5)),sc)
  sel <- read.csv(file.path(archive_dir(archives,'amplitude'),'robust-v1/summary-binary-final/selection.csv'),stringsAsFactors=FALSE)
  expect_identical(row$reasons,sel$initial_ig_reasons[sel$key==key]);expect_true(row$flagged)
})

test_that('frozen tables are written once and never silently replaced', {
  d <- tempfile('sel');x <- data.frame(a=1:2,b=c('x','y'))
  f <- file.path(d,'t.csv');expect_identical(write_frozen_table(x,f),'written')
  expect_identical(write_frozen_table(x,f),'unchanged')
  expect_error(write_frozen_table(transform(x,a=3:4),f),'Refusing to replace')
  expect_identical(read.csv(f)$a,1:2)
  expect_false(any(grepl('\\.tmp$',list.files(d))))
})

test_that('the selection needs every first fit, lists the long repeats, and assembles each arm\'s selected fit', {
  need_archives()
  items <- selection_items(c('A1','A2'),c(2,3,5),'/nowhere',archives,repo)
  expect_identical(sum(items$role=='new'),87L);expect_identical(sum(items$role=='control'),29L)
  expect_identical(sum(items$role=='new' & items$schedule=='initial'),60L)
  expect_identical(sum(items$role=='control' & items$schedule=='long'),9L)
  ctl <- items[items$role=='control',]
  expect_true(all(ctl$sd==1));expect_true(all(nchar(ctl$expected_md5)==32L))
  expect_identical(ctl$fit_label[ctl$key=='range6-rep03-binary-k100'],'spatial-targeted-20260927/long/range6-rep03-binary-k100-fit.rds')
  expect_error(check_first_fits(items),'87 first fits are missing')
  # Amendment 1 (R23): A1 and A2 are selected in separate runs, each in its own directory.
  a1 <- selection_items('A1',c(2,3,5),'/nowhere',archives,repo);a2 <- selection_items('A2',c(2,3,5),'/nowhere',archives,repo)
  expect_identical(c(sum(a1$role=='new'),sum(a1$role=='control'),sum(a2$role=='new'),sum(a2$role=='control')),c(60L,20L,27L,9L))
  expect_identical(rbind(a1,a2),items)
  expect_identical(selection_dir('/s','A1'),'/s/selection/A1');expect_identical(selection_dir('/s','A2'),'/s/selection/A2')
  expect_identical(selection_dir('/s',c('A1','A2')),'/s/selection/A1-A2')
  # long-keys.txt is a launcher job list.
  sel <- data.frame(phase=c('A1','A1','A1'),sd=c(2,3,5),key=c('jsdm-n0100-01','jsdm-n0100-02','jsdm-n0300-03'),
    flagged=c(TRUE,FALSE,TRUE))
  f <- tempfile();write_frozen_table(long_keys(sel),f)
  expect_identical(read_only_jobs(f),data.frame(phase=c('A1','A1'),sd=c(2,5),key=c('jsdm-n0100-01','jsdm-n0300-03'),stringsAsFactors=FALSE))
  flag <- function(phase,sd,key,schedule,flagged,role='new') data.frame(role=role,phase=phase,sd=sd,key=key,
    schedule=schedule,fit=paste0(key,'-',schedule),fit_md5='m',rule='pr11',warnings=0L,max_group_rhat=1,
    max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA_integer_,spatial_field_flags=NA_integer_,
    flagged=flagged,reasons=if(flagged) 'x' else '',stringsAsFactors=FALSE)
  k1 <- 'jsdm-n0100-01';k2 <- 'jsdm-n0300-01'
  initial <- rbind(flag('A1',2,k1,'initial',TRUE),flag('A1',2,k2,'initial',FALSE))
  repeats <- flag('A1',2,k1,'long',TRUE)
  longs <- flag('A2',2,'r1','long',FALSE)
  controls <- rbind(flag('A1',1,k1,'initial',FALSE,'control'),flag('A1',1,k2,'initial',FALSE,'control'),
    flag('A2',1,'r1','long',TRUE,'control'))
  s <- final_selected(initial,repeats,longs,controls)
  expect_identical(s$selected_schedule[s$sd==2],c('long','initial','long'))
  expect_identical(s$first_schedule[s$sd==2],c('initial','initial','long'))
  expect_identical(s$long_repeat[s$sd==2],c(TRUE,FALSE,FALSE))
  expect_identical(s$flagged[s$sd==2],c(TRUE,FALSE,FALSE));expect_identical(s$fit[s$sd==2][1],'jsdm-n0100-01-long')
  # Amendment 1 (R21): A1 convergence counts are kept separately at 100 and 300 sites.
  c <- convergence_counts(s)
  expect_identical(names(c),c('phase','stratum','sd','role','selected_fits','selected_flagged','first_flagged','long_repeats'))
  expect_identical(c$stratum,c('n100','n300','n100','n300','all','all'));expect_identical(c$sd,c(1,1,2,2,1,2))
  a1 <- function(stratum,sd,col) c[[col]][c$phase=='A1' & c$stratum==stratum & c$sd==sd]
  expect_identical(a1('n100',2,'selected_flagged'),1L);expect_identical(a1('n300',2,'selected_flagged'),0L)
  expect_identical(a1('n100',1,'selected_flagged'),0L);expect_identical(a1('n100',2,'long_repeats'),1L)
  expect_identical(a1('n300',2,'long_repeats'),0L);expect_identical(a1('n100',2,'selected_fits'),1L)
  expect_identical(c$selected_flagged[c$phase=='A2' & c$sd==1],1L)
  expect_identical(size_stratum(c('jsdm-n0100-07','jsdm-n0300-07','range4-rep01-binary-k100','design-qfar_K6-sites300-01')),
    c('n100','n300','all','all'))
  expect_error(final_selected(initial,repeats[0,],longs,controls),'longer repeat is missing')
})

test_that('select.R before the first fits exist exits non-zero and writes nothing', {
  need_archives()
  empty <- tempfile('study');dir.create(empty);empty <- normalizePath(empty)
  sel <- function(...) run_cli(file.path(here,'select.R'),c(paste0('--repo=',repo),paste0('--study=',empty),
    paste0('--inputs-root=',inputs_root),paste0('--archives=',archives),'--sds=2,3,5',...))
  r <- sel('--phases=A1,A2','--mode=plan')
  expect_identical(r$status,1L);expect_match(r$output,'87 first fits are missing')
  r <- sel('--phases=A1','--mode=plan','--workers=2')
  expect_identical(r$status,1L);expect_match(r$output,'60 first fits are missing')
  r <- sel('--phases=A1','--mode=final')
  expect_identical(r$status,1L);expect_match(r$output,'60 first fits are missing')
  expect_identical(list.files(empty,recursive=TRUE,all.files=TRUE,include.dirs=TRUE),character())
})

# Stub flags for the end-to-end selection test: every metadata, md5, input and
# table round-trip check is real; only the Rhat computation is replaced.
stub_flags <- function(phase,fit,warnings,input,sc) data.frame(rule='stub',warnings=length(warnings),max_group_rhat=1,
  max_element_rhat=1,unresolved_rhat=0L,spatial_trace_flags=NA_integer_,spatial_field_flags=NA_integer_,
  flagged=isTRUE(fit$stub_flag),reasons=if(isTRUE(fit$stub_flag)) 'stub flag' else '',stringsAsFactors=FALSE)

test_that('select.R plan and final record a flagged initial fit, its single longer repeat and each arm\'s selected fit', {
  need_archives()
  ts <- tempfile('study');dir.create(ts);ts <- normalizePath(ts)
  specs <- phase_jobs('A1',archives,inputs_root)
  fake_fit <- function(key,schedule,flag) {
    rec <- list(phase='A1',key=key,sd=2,schedule=schedule,mcmc=schedule_mcmc('A1',schedule,archives),
      input_md5=specs[[key]]$input_md5,fit=list(infos=list(intercept_prior=list(mean=0,sd=2)),stub_flag=flag),
      warnings=character())
    f <- fit_path(ts,'A1',2,schedule,key);dir.create(dirname(f),recursive=TRUE,showWarnings=FALSE);saveRDS(rec,f);f
  }
  for(k in phase_keys('A1')) fake_fit(k,'initial',k=='jsdm-n0100-03')
  args <- function(mode) c(paste0('--repo=',repo),paste0('--study=',ts),paste0('--inputs-root=',inputs_root),
    paste0('--archives=',archives),'--phases=A1','--sds=2',paste0('--mode=',mode),'--workers=2')
  out <- file.path(ts,'selection','A1')
  expect_identical(select_main(args('plan'),flag_fn=stub_flags),0L)
  sel <- read_flag_table(file.path(out,'long-selection.csv'))
  expect_identical(nrow(sel),20L);expect_identical(sel$key[sel$flagged],'jsdm-n0100-03');expect_true(is.double(sel$sd))
  expect_identical(read_only_jobs(file.path(out,'long-keys.txt')),data.frame(phase='A1',sd=2,key='jsdm-n0100-03',stringsAsFactors=FALSE))
  expect_identical(nrow(read_flag_table(file.path(out,'long-fit-flags.csv'))),0L)
  ctl <- read_flag_table(file.path(out,'control-flags.csv'))
  expect_identical(nrow(ctl),20L);expect_true(all(ctl$sd==1));expect_false(any(ctl$flagged))
  tables <- file.path(out,c('long-selection.csv','long-keys.txt','long-fit-flags.csv','control-flags.csv'))
  before <- unname(tools::md5sum(tables))
  expect_identical(select_main(args('plan'),flag_fn=stub_flags),0L)
  expect_identical(unname(tools::md5sum(tables)),before)
  # Final mode needs the single longer repeat of every flagged initial fit.
  expect_error(select_main(args('final'),flag_fn=stub_flags),'Longer repeats missing: 1.*longer-repeat launcher')
  expect_false(file.exists(file.path(out,'selected-fits.csv')))
  q <- build_queue('A1',c(2,3,5),ts,archives,read_only_jobs(file.path(out,'long-keys.txt')),'long')
  expect_identical(q$fit,fit_path(ts,'A1',2,'long','jsdm-n0100-03'))
  fake_fit('jsdm-n0100-03','long',FALSE)
  expect_identical(select_main(args('final'),flag_fn=stub_flags),0L)
  expect_identical(nrow(read_flag_table(file.path(out,'repeat-flags.csv'))),1L)
  s <- utils::read.csv(file.path(out,'selected-fits.csv'),stringsAsFactors=FALSE)
  expect_identical(c(sum(s$role=='new'),sum(s$role=='control')),c(20L,20L))
  r <- s[s$role=='new' & s$key=='jsdm-n0100-03',]
  expect_identical(r$selected_schedule,'long');expect_true(r$long_repeat);expect_true(r$first_flagged);expect_false(r$flagged)
  expect_identical(r$fit,'fits/A1/sd2/long/jsdm-n0100-03-fit.rds')
  expect_true(all(s$selected_schedule[s$role=='new' & s$key!='jsdm-n0100-03']=='initial'))
  cc <- utils::read.csv(file.path(out,'convergence.csv'),stringsAsFactors=FALSE)
  expect_identical(cc$stratum,c('n100','n300','n100','n300'));expect_identical(cc$selected_flagged,rep(0L,4L))
  expect_identical(cc$long_repeats[cc$sd==2],c(1L,0L));expect_identical(cc$first_flagged[cc$sd==2],c(1L,0L))
  # A first fit changed after the selection was recorded is refused.
  fake_fit('jsdm-n0300-10','initial',TRUE)
  expect_error(select_main(args('final'),flag_fn=stub_flags),'changed after the selection was recorded')
})

# The documented detached form: nohup, a new session via perl's setsid, a pid
# file and an exit-status file, around a stub launcher run.
test_that('the documented detached wrapper starts a new session and records its pid and exit status', {
  if(!nzchar(Sys.which('perl'))) skip('perl not available')
  ts <- stub_study();logs <- file.path(ts,'logs');dir.create(logs)
  driver <- tempfile(fileext='.R')
  writeLines(c(sprintf('setwd(%s)',deparse(here)),"source('jobs.R');source('verify-helpers.R');source('launch.R')",
    sprintf('ts <- %s',deparse(ts)),
    "cmd <- function(job) list(command='/bin/sh',args=c('-c',sprintf('sleep 1; mkdir -p %s; touch %s',shQuote(dirname(job$fit)),shQuote(job$fit))))",
    "q <- order_queue(queue_frame('A1',2,'initial',c('jsdm-n0100-01','jsdm-n0100-02'),'initial',ts))",
    "res <- execute_queue(q,ts,2L,cmd,poll=0.05,stamp='detached');quit(save='no',status=res$status)"),driver)
  inner <- sprintf('echo $$ $(ps -o pgid= -p $$) > %s; %s %s; echo $? > %s',shQuote(file.path(logs,'ids')),
    shQuote(rscript),shQuote(driver),shQuote(file.path(logs,'stub.exit')))
  system(sprintf("nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c %s < /dev/null > %s 2>&1 & echo $! > %s",
    shQuote(inner),shQuote(file.path(logs,'stub.out')),shQuote(file.path(logs,'stub.pid'))))
  for(i in 1:600) {if(file.exists(file.path(logs,'stub.exit'))) break;Sys.sleep(0.1)}
  expect_identical(readLines(file.path(logs,'stub.exit')),'0')
  ids <- scan(file.path(logs,'ids'),quiet=TRUE);pid <- as.numeric(readLines(file.path(logs,'stub.pid')))
  expect_identical(ids[1],ids[2]);expect_identical(ids[1],pid)
  expect_false(identical(ids[2],as.numeric(system('ps -o pgid= -p $$',intern=TRUE))))
  expect_true('exit_status: 0' %in% readLines(file.path(logs,'launch-detached','DONE')))
})
