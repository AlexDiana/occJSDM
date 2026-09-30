#!/usr/bin/env Rscript
# Launcher for the Task 3 diagnostic fits (README.md in this directory).
#
#   Rscript launch.R --repo=REPO --study=ARCHIVE --inputs-root=DIR --runs=LIST
#     --max-procs=N [--archives=DIR] [--poll-seconds=5] [--dry-run]
#
# LIST is comma-separated RUN:CHAINS items, RUN one of extended, a, b, c,
# pilot-extended, pilot-a, pilot-b, pilot-c and CHAINS one chain K or a range
# K1-K2, for example extended:1-16,a:1-8,b:1-8,c:1-8. Each chain is one
# process,
#   Rscript REPO/dev/simstudy/convergence-flag-diagnosis/run.R --repo=REPO
#     --study=ARCHIVE --inputs-root=DIR --archives=DIR --run=RUN --chain=K
# (pilots: --run=pilot --variant=V), started in list order with at most
# --max-procs (1 to 8) running at once. Each writes its own log and exit file,
# ARCHIVE/logs/launch-<timestamp>/<RUN>-chain-KK.log and .exit.
# It refuses to start while any fit lock (a *.lock directory) exists under
# ARCHIVE/fits, while another launcher holds ARCHIVE/logs/launcher.lock, when
# run.R, fixed-theta0.R, the library fingerprint or the clone archive differ
# from the frozen md5s below, or when the installed library does not match the
# fingerprint. When every job has ended it writes summary.csv and DONE (counts
# and exit status) and exits 1 unless every job completed or resumed.
# Start it detached; README.md gives the exact command.

`%||%` <- function(x,y) if(is.null(x)) y else x
here <- local({
  a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  f <- NULL
  for(i in rev(seq_len(sys.nframe()))) {x <- sys.frame(i)$ofile;if(!is.null(x)) {f <- x;break}}
  dirname(normalizePath(f %||% sub('^--file=','',a[1])))
})
source(file.path(here,'run.R'))

LAUNCH_SCRIPT <- file.path(STUDY_REL,'launch.R')
# Files hashed into every fit whose frozen versions the launcher insists on
# (set when the protocol was frozen; see README.md).
FROZEN_MD5 <- c(
  'dev/simstudy/convergence-flag-diagnosis/run.R'='TO-BE-SET',
  'dev/simstudy/convergence-flag-diagnosis/fixed-theta0.R'='TO-BE-SET',
  'dev/simstudy/convergence-flag-diagnosis/results/library-fingerprint.csv'='TO-BE-SET',
  'dev/simstudy/convergence-flag-diagnosis/results/fixed-theta0/hashes.csv'='TO-BE-SET')
LAUNCH_RUNS <- c(VARIANTS,paste0('pilot-',VARIANTS))
SUMMARY_COLUMNS <- c('order','label','run','variant','chain','seed','status','exit_code','elapsed_seconds',
  'started','finished','pid','log','fit')

# One row per single-chain job, in list order.
parse_run_list <- function(text) {
  items <- strsplit(text,',',fixed=TRUE)[[1]]
  rows <- lapply(items,function(it) {
    if(!grepl('^[a-z-]+:[0-9]+(-[0-9]+)?$',it)) stop('Run list items must be RUN:CHAINS, e.g. extended:1-16; got ',it)
    run <- sub(':.*$','',it);ch <- sub('^[^:]*:','',it)
    if(!run %in% LAUNCH_RUNS) stop('Unknown run ',run,'; use one of ',paste(LAUNCH_RUNS,collapse=', '))
    bounds <- as.integer(strsplit(ch,'-',fixed=TRUE)[[1]])
    chains <- if(length(bounds)==2L) seq(bounds[1],bounds[2]) else bounds
    pilot <- startsWith(run,'pilot-')
    base <- if(pilot) 'pilot' else run;variant <- if(pilot) sub('^pilot-','',run) else run
    for(k in chains) run_spec(base,if(pilot) variant else NULL,k)
    data.frame(run=base,variant=variant,chain=as.integer(chains),stringsAsFactors=FALSE)
  })
  q <- do.call(rbind,rows)
  dup <- duplicated(q[c('run','variant','chain')])
  if(any(dup)) stop('A chain is listed twice: ',q$run[dup][1],' ',q$variant[dup][1],' chain ',q$chain[dup][1])
  q$seed <- chain_seed(q$chain)
  q$label <- sprintf('%s-chain-%02d',ifelse(q$run=='pilot',paste0('pilot-',q$variant),q$run),q$chain)
  q$order <- seq_len(nrow(q));rownames(q) <- NULL
  q
}

job_args <- function(job,repo,study,inputs_root,archives)
  c(paste0('--repo=',repo),paste0('--study=',study),paste0('--inputs-root=',inputs_root),paste0('--archives=',archives),
    paste0('--run=',job$run),if(job$run=='pilot') paste0('--variant=',job$variant),paste0('--chain=',job$chain))

run_r_command <- function(repo,study,inputs_root,archives) {
  script <- file.path(repo,RUN_SCRIPT)
  function(job) list(command=file.path(R.home('bin'),'Rscript'),args=c(script,job_args(job,repo,study,inputs_root,archives)))
}

fit_locks <- function(study) {
  fits <- file.path(study,'fits')
  if(!dir.exists(fits)) return(character())
  d <- list.dirs(fits,recursive=TRUE,full.names=TRUE)
  d[grepl('\\.lock$',d)]
}

check_frozen <- function(repo) {
  files <- file.path(repo,names(FROZEN_MD5))
  found <- ifelse(file.exists(files),unname(tools::md5sum(files)),'missing')
  changed <- names(FROZEN_MD5)[found!=FROZEN_MD5]
  if(length(changed)) stop('Frozen fit file changed: ',paste(changed,collapse=', '),
    '; diagnostic fits must use the frozen run.R, fixed-theta0.R, fingerprint and clone archive')
  invisible(TRUE)
}

launch_now <- function() format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')
atomic_write_csv <- function(x,path) {tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE);stopifnot(file.rename(tmp,path))}
atomic_write_lines <- function(x,path) {tmp <- paste0(path,'.tmp');writeLines(x,tmp);stopifnot(file.rename(tmp,path))}

job_status <- function(code,log,fit) {
  if(is.na(code) || code!=0L || !file.exists(fit)) return('failed')
  text <- if(file.exists(log)) readLines(log,warn=FALSE) else character()
  if(any(grepl(RESUMED_TEXT,text,fixed=TRUE))) 'resumed' else 'completed'
}

# Run the jobs, at most max_procs at a time. Returns the exit status, the log
# directory (NULL if nothing started) and the largest concurrency reached.
execute_jobs <- function(q,study,max_procs,command_fn,poll=5,stamp=format(Sys.time(),'%Y%m%d-%H%M%S')) {
  stopifnot(nrow(q)>0L,max_procs %in% seq_len(8L))
  refused <- list(status=1L,dir=NULL,max_concurrent=0L)
  locks <- fit_locks(study)
  if(length(locks)) {
    cat('Refusing to start: ',length(locks),' fit lock(s) exist. A killed run may have left them: check that no run.R ',
      'process holds each one, then remove it by hand.\n',paste0('  ',locks,'\n'),sep='')
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
  q$fit <- vapply(seq_len(nrow(q)),function(i) fit_file(study,run_spec(q$run[i],if(q$run[i]=='pilot') q$variant[i] else NULL,
    q$chain[i])),character(1))
  writeLines(c(paste('launcher_pid:',Sys.getpid()),paste('started:',launch_now()),paste('max_procs:',max_procs),
    paste('jobs:',nrow(q)),paste('R:',R.version.string)),file.path(dir,'launcher.txt'))
  utils::write.csv(q,file.path(dir,'queue.csv'),row.names=FALSE)
  res <- data.frame(order=q$order,label=q$label,run=q$run,variant=q$variant,chain=q$chain,seed=q$seed,status='queued',
    exit_code=NA_integer_,elapsed_seconds=NA_real_,started=NA_character_,finished=NA_character_,pid=NA_integer_,
    log=paste0(q$label,'.log'),fit=q$fit,stringsAsFactors=FALSE)
  running <- list();next_job <- 1L;max_seen <- 0L;error <- NULL
  finish <- function(i,code,started) {
    res$exit_code[i] <<- as.integer(code)
    res$elapsed_seconds[i] <<- round(as.numeric(difftime(Sys.time(),started,units='secs')),2)
    res$finished[i] <<- format(Sys.time(),'%Y-%m-%d %H:%M:%S')
    writeLines(as.character(code),file.path(dir,paste0(q$label[i],'.exit')))
    res$status[i] <<- job_status(as.integer(code),file.path(dir,res$log[i]),q$fit[i])
    cat(launch_now(),res$status[i],q$label[i],'exit',code,'\n');flush.console()
  }
  start <- function(i) {
    job <- as.list(q[i,,drop=FALSE]);cmd <- command_fn(job);log <- file.path(dir,res$log[i])
    started <- Sys.time();res$started[i] <<- format(started,'%Y-%m-%d %H:%M:%S');res$status[i] <<- 'running'
    p <- tryCatch(processx::process$new(cmd$command,cmd$args,stdout=log,stderr='2>&1',cleanup=TRUE),
      error=function(e) {writeLines(paste('Launch error:',conditionMessage(e)),log);NULL})
    if(is.null(p)) {finish(i,NA_integer_,started);return(invisible())}
    res$pid[i] <<- p$get_pid()
    cat(launch_now(),'started',q$label[i],'pid',p$get_pid(),'\n');flush.console()
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
  o <- parse_options(args[!dry],known=c('repo','study','inputs-root','archives','runs','max-procs','poll-seconds'),
    required=c('repo','study','inputs-root','runs','max-procs'))
  repo <- normalizePath(o$repo,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this launcher')
  study <- normalizePath(o$study,mustWork=TRUE);inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  mp <- o$`max-procs`
  if(!grepl('^[0-9]+$',mp) || !as.integer(mp) %in% 1:8) stop('--max-procs must be an integer from 1 to 8')
  poll <- as.numeric(o$`poll-seconds` %||% '5');if(!is.finite(poll) || poll<=0) stop('--poll-seconds must be positive')
  q <- parse_run_list(o$runs)
  check_frozen(repo)
  revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
  cat('Queue:',nrow(q),'single-chain jobs, at most',mp,'at once; library revision',revision,'\n')
  cat('README.md md5',unname(tools::md5sum(file.path(repo,STUDY_REL,'README.md'))),'\n')
  if(any(dry)) {print(q,row.names=FALSE);cat('Dry run: nothing started.\n');return(0L)}
  execute_jobs(q,study,as.integer(mp),run_r_command(repo,study,inputs_root,archives),poll=poll)$status
}

if(sys.nframe()==0L) {
  status <- tryCatch(launch_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
