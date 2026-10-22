#!/usr/bin/env Rscript
# Launcher for the Hmsc refits (README.md in this directory).
#
#   Rscript launch.R --repo=REPO --root=EXTENSION_ARCHIVE --study=ARCHIVE
#     --max-procs=N [--jobs=all|N1,N2,...] [--poll-seconds=5] [--dry-run]
#
# Each selected job (results/jobs.csv, in job-number order) is one process,
#   Rscript REPO/dev/simstudy/hmsc-longer-runs/refit.R --repo=REPO --root=ROOT
#     --study=ARCHIVE --job-number=N
# with at most --max-procs (1 to 8) running at once. Each writes its own log
# and exit file under ARCHIVE/logs/launch-<timestamp>/. It refuses to start
# while any job lock exists under ARCHIVE/jobs or another launcher holds
# ARCHIVE/logs/launcher.lock, and records the md5s of common.R, refit.R,
# compare.R and README.md in launcher.txt. When every job has ended it writes
# summary.csv and DONE and exits 1 unless every job completed or resumed.
here <- local({
  a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  f <- NULL
  for(i in rev(seq_len(sys.nframe()))) {x <- sys.frame(i)$ofile;if(!is.null(x)) {f <- x;break}}
  dirname(normalizePath(f %||% sub('^--file=','',a[1])))
})
`%||%` <- function(x,y) if(is.null(x)) y else x
source(file.path(here,'common.R'))

RECORDED_FILES <- file.path(STUDY_REL,c('common.R','refit.R','compare.R','README.md'))
SUMMARY_COLUMNS <- c('order','label','job_number','job','seed','status','exit_code','elapsed_seconds','started','finished','pid','log','result')

parse_jobs <- function(text,jobs) {
  if(identical(text,'all')) q <- jobs else {
    if(!grepl('^[0-9]+(,[0-9]+)*$',text)) stop('--jobs must be all or a comma-separated list of job numbers')
    n <- as.integer(strsplit(text,',',fixed=TRUE)[[1]])
    absent <- setdiff(n,jobs$job_number)
    if(length(absent)) stop('Job number(s) not in the selection: ',paste(absent,collapse=', '))
    if(anyDuplicated(n)) stop('A job number is listed twice')
    q <- jobs[match(n,jobs$job_number),,drop=FALSE]
  }
  q <- q[order(q$job_number),,drop=FALSE]
  q$seed <- job_seed(q$job_number)
  q$label <- sprintf('job-%03d',q$job_number)
  q$order <- seq_len(nrow(q));rownames(q) <- NULL
  q
}

run_r_command <- function(repo,root,study) {
  script <- file.path(repo,STUDY_REL,'refit.R')
  function(job) list(command=file.path(R.home('bin'),'Rscript'),
    args=c(script,paste0('--repo=',repo),paste0('--root=',root),paste0('--study=',study),paste0('--job-number=',job$job_number)))
}

job_locks <- function(study) {
  d <- file.path(study,'jobs')
  if(!dir.exists(d)) return(character())
  x <- list.dirs(d,recursive=TRUE,full.names=TRUE)
  x[grepl('\\.lock$',x)]
}

launch_now <- function() format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')
atomic_write_csv <- function(x,path) {tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE);stopifnot(file.rename(tmp,path))}
atomic_write_lines <- function(x,path) {tmp <- paste0(path,'.tmp');writeLines(x,tmp);stopifnot(file.rename(tmp,path))}

job_status <- function(code,log,result) {
  if(is.na(code) || code!=0L || !file.exists(result)) return('failed')
  text <- if(file.exists(log)) readLines(log,warn=FALSE) else character()
  if(any(grepl(RESUMED_TEXT,text,fixed=TRUE))) 'resumed' else 'completed'
}

execute_jobs <- function(q,study,max_procs,command_fn,poll=5,stamp=format(Sys.time(),'%Y%m%d-%H%M%S'),recorded=character()) {
  stopifnot(nrow(q)>0L,max_procs %in% seq_len(8L))
  refused <- list(status=1L,dir=NULL,max_concurrent=0L)
  locks <- job_locks(study)
  if(length(locks)) {
    cat('Refusing to start: ',length(locks),' job lock(s) exist. Check that no refit.R process holds each one, then remove it by hand.\n',
      paste0('  ',locks,'\n'),sep='')
    return(refused)
  }
  logs <- file.path(study,'logs');dir.create(logs,recursive=TRUE,showWarnings=FALSE)
  guard <- file.path(logs,'launcher.lock')
  if(!dir.create(guard,showWarnings=FALSE)) {
    cat('Refusing to start: another launcher holds ',guard,'. Check that it is not running, then remove it by hand.\n',sep='')
    return(refused)
  }
  on.exit(unlink(guard,recursive=TRUE),add=TRUE)
  dir <- file.path(logs,paste0('launch-',stamp));if(file.exists(dir)) dir <- paste0(dir,'-',Sys.getpid())
  dir.create(dir)
  writeLines(c(paste('pid',Sys.getpid()),paste('since',launch_now()),paste('log',dir)),file.path(guard,'owner'))
  q$result <- file.path(study,'jobs',q$job,'result.rds')
  writeLines(c(paste('launcher_pid:',Sys.getpid()),paste('started:',launch_now()),paste('max_procs:',max_procs),
    paste('jobs:',nrow(q)),paste('R:',R.version.string),paste0('md5 ',names(recorded),': ',recorded)),file.path(dir,'launcher.txt'))
  utils::write.csv(q,file.path(dir,'queue.csv'),row.names=FALSE)
  res <- data.frame(order=q$order,label=q$label,job_number=q$job_number,job=q$job,seed=q$seed,status='queued',
    exit_code=NA_integer_,elapsed_seconds=NA_real_,started=NA_character_,finished=NA_character_,pid=NA_integer_,
    log=paste0(q$label,'.log'),result=q$result,stringsAsFactors=FALSE)
  running <- list();next_job <- 1L;max_seen <- 0L;error <- NULL
  finish <- function(i,code,started) {
    res$exit_code[i] <<- as.integer(code)
    res$elapsed_seconds[i] <<- round(as.numeric(difftime(Sys.time(),started,units='secs')),2)
    res$finished[i] <<- format(Sys.time(),'%Y-%m-%d %H:%M:%S')
    writeLines(as.character(code),file.path(dir,paste0(q$label[i],'.exit')))
    res$status[i] <<- job_status(as.integer(code),file.path(dir,res$log[i]),q$result[i])
    cat(launch_now(),res$status[i],q$label[i],q$job[i],'exit',code,'\n');flush.console()
  }
  start <- function(i) {
    job <- as.list(q[i,,drop=FALSE]);cmd <- command_fn(job);log <- file.path(dir,res$log[i])
    started <- Sys.time();res$started[i] <<- format(started,'%Y-%m-%d %H:%M:%S');res$status[i] <<- 'running'
    p <- tryCatch(processx::process$new(cmd$command,cmd$args,stdout=log,stderr='2>&1',cleanup=TRUE),
      error=function(e) {writeLines(paste('Launch error:',conditionMessage(e)),log);NULL})
    if(is.null(p)) {finish(i,NA_integer_,started);return(invisible())}
    res$pid[i] <<- p$get_pid()
    cat(launch_now(),'started',q$label[i],q$job[i],'pid',p$get_pid(),'\n');flush.console()
    running[[as.character(i)]] <<- list(process=p,started=started)
  }
  reap <- function() for(id in names(running)) {
    r <- running[[id]]
    if(!r$process$is_alive()) {r$process$wait();finish(as.integer(id),r$process$get_exit_status(),r$started);running[[id]] <<- NULL}
  }
  tryCatch({
    while(next_job<=nrow(q) || length(running)) {
      reap()
      while(length(running)<max_procs && next_job<=nrow(q)) {start(next_job);next_job <- next_job+1L}
      max_seen <- max(max_seen,length(running))
      atomic_write_csv(res[SUMMARY_COLUMNS],file.path(dir,'progress.csv'))
      if(length(running)) Sys.sleep(poll)
    }
  },error=function(e) {
    error <<- conditionMessage(e)
    cat(launch_now(),'LAUNCHER ERROR:',error,'; starting nothing further and waiting for running jobs\n');flush.console()
    while(length(running)) {tryCatch(reap(),error=function(e) running <<- list());if(length(running)) Sys.sleep(poll)}
  })
  res$status[res$status %in% c('queued','running')] <- 'not run'
  counts <- table(factor(res$status,levels=c('completed','resumed','failed','not run')))
  status <- if(is.null(error) && all(res$status %in% c('completed','resumed'))) 0L else 1L
  atomic_write_csv(res[SUMMARY_COLUMNS],file.path(dir,'summary.csv'))
  atomic_write_csv(res[SUMMARY_COLUMNS],file.path(dir,'progress.csv'))
  atomic_write_lines(c(paste('finished:',launch_now()),paste('exit_status:',status),paste('jobs:',nrow(res)),
    paste0(sub(' ','_',names(counts)),': ',as.integer(counts)),paste('max_concurrent:',max_seen),
    paste('launcher_error:',error %||% '')),file.path(dir,'DONE'))
  cat(launch_now(),'DONE:',paste(names(counts),as.integer(counts),collapse=', '),'; exit status',status,'\n')
  list(status=status,dir=dir,max_concurrent=max_seen)
}

launch_main <- function(args) {
  dry <- args=='--dry-run'
  o <- parse_options(args[!dry],known=c('repo','root','study','max-procs','jobs','poll-seconds'),required=c('repo','root','study','max-procs'))
  repo <- normalizePath(o$repo,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this launcher')
  root <- normalizePath(o$root,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  mp <- o$`max-procs`
  if(!grepl('^[0-9]+$',mp) || !as.integer(mp) %in% 1:8) stop('--max-procs must be an integer from 1 to 8')
  poll <- as.numeric(o$`poll-seconds` %||% '5');if(!is.finite(poll) || poll<=0) stop('--poll-seconds must be positive')
  check_fit_job(root);load_code(repo);load_hmsc(study)
  q <- parse_jobs(o$jobs %||% 'all',read_jobs(repo))
  recorded <- hash_files(file.path(repo,RECORDED_FILES),repo)
  cat('Queue:',nrow(q),'Hmsc refits, at most',mp,'at once; Hmsc',HMSC_VERSION,'from',file.path(study,'library'),'\n')
  cat(paste0('md5 ',names(recorded),': ',recorded,'\n'),sep='')
  if(any(dry)) {print(q[c('order','job_number','job','seed')],row.names=FALSE);cat('Dry run: nothing started.\n');return(0L)}
  execute_jobs(q,study,as.integer(mp),run_r_command(repo,root,study),poll=poll,recorded=recorded)$status
}

if(sys.nframe()==0L) {
  status <- tryCatch(launch_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
