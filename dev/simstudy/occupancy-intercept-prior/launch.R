#!/usr/bin/env Rscript
# Detached phase launcher for the occupancy-intercept prior study.
#
#   Rscript launch.R --repo=REPO --study=STUDY --inputs-root=DIR --phases=A1,A2
#     --sds=2,3,5 --max-procs=8 [--dry-run] [--only=FILE] [--schedule-override=long]
#     [--archives=DIR] [--queue-out=FILE] [--poll-seconds=5]
#
# Builds one job per phase, community and prior SD, at the schedule of the
# control arm's selected fit (control_schedules in verify-helpers.R), orders
# the jobs longest first and runs at most --max-procs single-key processes
#   Rscript <repo>/dev/simstudy/occupancy-intercept-prior/run.R --repo=<repo> ...
#     --workers=1 --keys=<key>
# each with its own log under STUDY/logs/launch-<timestamp>/. It refuses to
# start while any fit lock of a queued key exists, or while another launcher
# holds STUDY/logs/launcher.lock. On completion it writes summary.csv and a DONE
# marker in the log directory and exits non-zero if any job did not complete.
# --only restricts the queue to the (phase, sd, key) rows of a CSV such as
# select.R's long-keys.txt; --schedule-override=long (which requires --only)
# fits those rows at the longer schedule for the single longer repeat.
# This file is not sourced or hashed by run.R. Start it detached, for example
#   nohup Rscript launch.R ... < /dev/null > STUDY/logs/launch.out 2>&1 &

LAUNCH_SCRIPT <- 'dev/simstudy/occupancy-intercept-prior/launch.R'
# run.R, jobs.R and the library fingerprint are hashed into every fit; the
# launcher refuses to start fits with any other version of them.
FROZEN_FIT_MD5 <- c(
  'dev/simstudy/occupancy-intercept-prior/run.R'='e111cff6f079dcc783c7a2b32fa28eaf',
  'dev/simstudy/occupancy-intercept-prior/jobs.R'='d6dc59ff47a25de08fad60b8bde827d1',
  'dev/simstudy/occupancy-intercept-prior/results/library-fingerprint.csv'='22d6191da2347a26009bdf15b6c31a4d')
# The control arm (SD 1) reuses archived fits and is never refitted here.
LAUNCH_SDS <- c(2,3,5)
OVERRIDE_SCHEDULES <- 'long'
QUEUE_COLUMNS <- c('order','phase','sd','schedule','key','control_schedule','expected_seconds','fit','fit_exists')
SUMMARY_COLUMNS <- c('order','key','phase','sd','schedule','status','exit_code','elapsed_seconds',
  'started','finished','pid','log','fit')
RESUMED_TEXT <- 'already complete with identical settings; resumed'

split_list <- function(text,name) {
  if(is.null(text) || !grepl('^[^,]+(,[^,]+)*$',text)) stop('--',name,' must be a comma-separated list')
  x <- strsplit(text,',',fixed=TRUE)[[1]]
  if(anyDuplicated(x)) stop('--',name,' repeats ',x[duplicated(x)][1])
  x
}

parse_launch_args <- function(args) {
  flag <- args=='--dry-run'
  if(sum(flag)>1L) stop('Repeated option: --dry-run')
  o <- parse_options(args[!flag],
    known=c('repo','study','inputs-root','archives','phases','sds','max-procs','only','schedule-override','queue-out','poll-seconds'),
    required=c('repo','study','inputs-root','phases','sds','max-procs'))
  phases <- split_list(o$phases,'phases')
  if(!all(phases %in% PHASES)) stop('--phases must name phases among ',paste(PHASES,collapse=', '))
  sds <- split_list(o$sds,'sds')
  if(!all(sds %in% as.character(LAUNCH_SDS))) stop('--sds must list prior SDs among ',paste(LAUNCH_SDS,collapse=', '),
    ' (the SD 1 control arm reuses archived fits)')
  mp <- o$`max-procs`
  if(!grepl('^[0-9]+$',mp) || !as.integer(mp) %in% seq_len(MAX_WORKERS))
    stop('--max-procs must be an integer from 1 to ',MAX_WORKERS)
  override <- o$`schedule-override`
  if(!is.null(override)) {
    if(!override %in% OVERRIDE_SCHEDULES) stop('--schedule-override must be one of ',paste(OVERRIDE_SCHEDULES,collapse=', '))
    if(is.null(o$only)) stop('--schedule-override requires --only (the single longer repeat is for listed jobs only)')
  }
  poll <- if(is.null(o$`poll-seconds`)) 5 else suppressWarnings(as.numeric(o$`poll-seconds`))
  if(length(poll)!=1L || !is.finite(poll) || poll<=0) stop('--poll-seconds must be a positive number')
  list(repo=o$repo,study=o$study,inputs_root=o$`inputs-root`,
    archives=o$archives %||% dirname(normalizePath(o$study,mustWork=FALSE)),
    phases=phases,sds=as.numeric(sds),max_procs=as.integer(mp),dry_run=any(flag),only=o$only,
    schedule_override=override,queue_out=o$`queue-out`,poll=poll)
}

# The run.R this launcher invokes is <repo>/RUN_SCRIPT, so run.R always sources
# and hashes the jobs.R of the same checkout. `script`, when known, is this
# launcher's own path, which must belong to the same checkout.
check_launch_repo <- function(repo,script=NULL) {
  repo <- normalizePath(repo,mustWork=TRUE)
  run <- file.path(repo,RUN_SCRIPT)
  if(!file.exists(run)) stop('run.R not found under --repo: ',run)
  if(!is.null(script)) {
    expected <- file.path(repo,LAUNCH_SCRIPT)
    if(!file.exists(script) || !identical(normalizePath(script),normalizePath(expected,mustWork=FALSE)))
      stop('This launcher (',script,') is not the launcher of --repo (',expected,'); run the launcher of the checkout whose run.R it starts')
  }
  files <- file.path(repo,names(FROZEN_FIT_MD5))
  found <- ifelse(file.exists(files),unname(tools::md5sum(files)),'missing')
  changed <- names(FROZEN_FIT_MD5)[found!=FROZEN_FIT_MD5]
  if(length(changed)) stop('Frozen fit script changed: ',paste(changed,collapse=', '),
    '; fits must use the frozen run.R, jobs.R and library fingerprint')
  normalizePath(run)
}

# Planning estimates (seconds per fit) from the archived runs recorded in the
# Task 2 report, used only to order the queue and to estimate wall time.
expected_seconds <- function(phase,schedule,key) {
  initial <- switch(phase,A1=if(grepl('-n0300-',key)) 33 else 16,B=140,A2=1450,stop('Unknown phase: ',phase))
  # The longer schedule runs 4 x 18,000 iterations against 2 x 8,000.
  switch(schedule,initial=initial,long=if(phase=='A2') 9400 else initial*4.5,stop('Unknown schedule: ',schedule))
}

# Longest-processing-time list scheduling: makespan at a concurrency cap.
estimated_wall_seconds <- function(seconds,max_procs) {
  free <- numeric(max_procs)
  for(s in sort(seconds,decreasing=TRUE)) {i <- which.min(free);free[i] <- free[i]+s}
  max(free)
}

queue_frame <- function(phase,sd,schedule,key,control_schedule,study) {
  n <- length(key)
  q <- data.frame(order=rep(NA_integer_,n),phase=rep(phase,length.out=n),sd=rep(as.numeric(sd),length.out=n),
    schedule=rep(schedule,length.out=n),key=key,control_schedule=rep(control_schedule,length.out=n),stringsAsFactors=FALSE)
  q$expected_seconds <- as.numeric(unlist(mapply(expected_seconds,q$phase,q$schedule,q$key,USE.NAMES=FALSE)))
  q$fit <- as.character(unlist(mapply(fit_path,study,q$phase,q$sd,q$schedule,q$key,USE.NAMES=FALSE)))
  q$lock <- fit_lock_path(q$fit)
  q
}

# Longest first; ties keep phase, SD and community order.
order_queue <- function(q) {
  if(!nrow(q)) {q$order <- integer();return(q)}
  rank <- function(p,k) match(k,phase_keys(p))
  within <- mapply(rank,q$phase,q$key,USE.NAMES=FALSE)
  q <- q[order(-q$expected_seconds,match(q$phase,PHASES),q$sd,within),,drop=FALSE]
  q$order <- seq_len(nrow(q));rownames(q) <- NULL
  q
}

read_only_jobs <- function(file) {
  x <- utils::read.csv(file,stringsAsFactors=FALSE,colClasses='character')
  if(!identical(names(x),c('phase','sd','key'))) stop('Job list ',file,' must have columns phase, sd, key')
  x$sd <- as.numeric(x$sd)
  if(anyNA(x$sd)) stop('Job list ',file,' has a non-numeric sd')
  x
}

build_queue <- function(phases,sds,study,archives,only=NULL,schedule_override=NULL) {
  control <- do.call(rbind,lapply(phases,control_schedules,archives=archives))
  jobs <- do.call(rbind,lapply(phases,function(p) data.frame(phase=p,sd=rep(sds,each=length(phase_keys(p))),
    key=rep(phase_keys(p),length(sds)),stringsAsFactors=FALSE)))
  if(!is.null(only)) {
    if(anyDuplicated(only)) stop('Duplicated job in the job list')
    bad <- !only$phase %in% phases
    if(any(bad)) stop('Job list phase ',only$phase[bad][1],' is not in --phases')
    bad <- !only$sd %in% sds
    if(any(bad)) stop('Job list SD ',only$sd[bad][1],' is not in --sds')
    for(i in seq_len(nrow(only))) if(!only$key[i] %in% phase_keys(only$phase[i]))
      stop('Job list key ',only$key[i],' is not a phase ',only$phase[i],' key')
    jobs <- only
  }
  cs <- control$schedule[match(paste(jobs$phase,jobs$key),paste(control$phase,control$key))]
  stopifnot(!anyNA(cs))
  schedule <- cs
  if(!is.null(schedule_override)) {
    same <- cs==schedule_override
    if(any(same)) stop('Job ',jobs$phase[same][1],' ',jobs$key[same][1],' is already at the ',schedule_override,
      ' schedule; the protocol allows no further escalation')
    schedule <- rep(schedule_override,nrow(jobs))
  }
  if(!nrow(jobs)) return(order_queue(queue_frame(character(),numeric(),character(),character(),character(),study)))
  order_queue(queue_frame(jobs$phase,jobs$sd,schedule,jobs$key,cs,study))
}

queue_counts <- function(q) {
  if(!nrow(q)) return(data.frame(phase=character(),sd=numeric(),schedule=character(),jobs=integer()))
  a <- stats::aggregate(list(jobs=q$key),q[c('phase','sd','schedule')],length)
  a <- a[order(match(a$phase,PHASES),a$sd,a$schedule),];rownames(a) <- NULL
  a$jobs <- as.integer(a$jobs);a
}

existing_locks <- function(q) q$lock[dir.exists(q$lock)]

relative_to <- function(path,root) {
  prefix <- paste0(normalizePath(root,mustWork=FALSE),'/')
  ifelse(startsWith(path,prefix),substring(path,nchar(prefix)+1L),path)
}
queue_table <- function(q,study) {
  x <- q;x$fit_exists <- file.exists(q$fit);x$fit <- relative_to(q$fit,study)
  x[QUEUE_COLUMNS]
}

# The command each job runs: one run.R process, one worker, one key.
run_r_command <- function(run_script,repo,study,inputs_root,archives) {
  force(run_script)
  function(job) list(command=file.path(R.home('bin'),'Rscript'),args=c(run_script,paste0('--repo=',repo),
    paste0('--study=',study),paste0('--phase=',job$phase),paste0('--sd=',format(job$sd)),
    paste0('--schedule=',job$schedule),'--workers=1',paste0('--keys=',job$key),
    paste0('--inputs-root=',inputs_root),paste0('--archives=',archives)))
}

quarantine_files <- function(fit) {
  key <- sub('-fit\\.rds$','',basename(fit))
  list.files(dirname(fit),pattern=paste0('^',gsub('([.+])','\\\\\\1',key),'-fit\\.QUARANTINE-.*\\.rds$'))
}

# completed: exit 0 with the fit saved; resumed: run.R found an identical fit;
# quarantined: a new quarantine file for the key; anything else failed.
job_status <- function(exit_code,log_text,fit,quarantine_before,quarantine_after) {
  if(is.na(exit_code)) return('failed')
  if(exit_code==0L) {
    if(!file.exists(fit)) return('failed')
    return(if(any(grepl(RESUMED_TEXT,log_text,fixed=TRUE))) 'resumed' else 'completed')
  }
  if(length(setdiff(quarantine_after,quarantine_before))) 'quarantined' else 'failed'
}

atomic_write_lines <- function(lines,path) {tmp <- paste0(path,'.tmp');writeLines(lines,tmp);stopifnot(file.rename(tmp,path))}
atomic_write_csv <- function(x,path) {tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE);stopifnot(file.rename(tmp,path))}
launch_now <- function() format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')

refuse_locks <- function(q) {
  locks <- existing_locks(q)
  if(!length(locks)) return(FALSE)
  cat('Refusing to start: ',length(locks),if(length(locks)>1L) ' fit locks exist' else ' fit lock exists',' for queued keys. ',
    'A killed run may have left them: check that no run.R process holds each one, then remove it by hand.\n',
    paste0('  ',locks,'\n'),sep='')
  TRUE
}

# Run the queue with at most max_procs processes. Returns the exit status, the
# log directory (NULL if nothing started) and the largest concurrency reached.
execute_queue <- function(q,study,max_procs,command_fn,poll=5,stamp=format(Sys.time(),'%Y%m%d-%H%M%S')) {
  stopifnot(nrow(q)>0L,max_procs %in% seq_len(MAX_WORKERS),!anyNA(q$order))
  refused <- list(status=1L,dir=NULL,max_concurrent=0L)
  if(refuse_locks(q)) return(refused)
  logs <- file.path(study,'logs');dir.create(logs,recursive=TRUE,showWarnings=FALSE)
  guard <- file.path(logs,'launcher.lock')
  if(!dir.create(guard,showWarnings=FALSE)) {
    owner <- file.path(guard,'owner')
    cat('Refusing to start: another launcher holds ',guard,' (',
      if(file.exists(owner)) paste(readLines(owner,warn=FALSE),collapse=', ') else 'owner unknown',
      '). Check that it is not running, then remove the directory by hand.\n',sep='')
    return(refused)
  }
  on.exit(unlink(guard,recursive=TRUE),add=TRUE)
  dir <- file.path(logs,paste0('launch-',stamp))
  if(file.exists(dir)) dir <- paste0(dir,'-',Sys.getpid())
  dir.create(dir)
  writeLines(c(paste('pid',Sys.getpid()),paste('host',Sys.info()[['nodename']]),paste('since',launch_now()),paste('log',dir)),
    file.path(guard,'owner'))
  writeLines(c(paste('launcher_pid:',Sys.getpid()),paste('started:',launch_now()),paste('max_procs:',max_procs),
    paste('jobs:',nrow(q)),paste('R:',R.version.string),paste('processx:',utils::packageVersion('processx'))),
    file.path(dir,'launcher.txt'))
  utils::write.csv(queue_table(q,study),file.path(dir,'queue.csv'),row.names=FALSE)
  res <- data.frame(order=q$order,key=q$key,phase=q$phase,sd=q$sd,schedule=q$schedule,status='queued',
    exit_code=NA_integer_,elapsed_seconds=NA_real_,started=NA_character_,finished=NA_character_,pid=NA_integer_,
    log=sprintf('%03d-%s-sd%s-%s-%s.log',q$order,q$phase,vapply(q$sd,format,''),q$schedule,q$key),
    fit=relative_to(q$fit,study),stringsAsFactors=FALSE)
  running <- list();next_job <- 1L;max_seen <- 0L;error <- NULL
  start_job <- function(i) {
    job <- as.list(q[i,,drop=FALSE]);log <- file.path(dir,res$log[i])
    before <- quarantine_files(q$fit[i])
    cmd <- command_fn(job)
    started <- Sys.time()
    p <- tryCatch(processx::process$new(cmd$command,cmd$args,stdout=log,stderr='2>&1',cleanup=TRUE),
      error=function(e) {writeLines(paste('Launch error:',conditionMessage(e)),log);NULL})
    res$started[i] <<- format(started,'%Y-%m-%d %H:%M:%S');res$status[i] <<- 'running'
    if(is.null(p)) {finish_job(i,NA_integer_,started,before);return(invisible())}
    res$pid[i] <<- p$get_pid()
    cat(launch_now(),'started',res$log[i],'pid',p$get_pid(),'\n');flush.console()
    running[[as.character(i)]] <<- list(process=p,started=started,before=before)
  }
  finish_job <- function(i,code,started,before) {
    log <- file.path(dir,res$log[i])
    text <- if(file.exists(log)) readLines(log,warn=FALSE) else character()
    res$exit_code[i] <<- as.integer(code)
    res$elapsed_seconds[i] <<- round(as.numeric(difftime(Sys.time(),started,units='secs')),2)
    res$finished[i] <<- format(Sys.time(),'%Y-%m-%d %H:%M:%S')
    res$status[i] <<- job_status(as.integer(code),text,q$fit[i],before,quarantine_files(q$fit[i]))
    cat(launch_now(),res$status[i],res$log[i],'exit',code,'\n');flush.console()
  }
  reap <- function() for(id in names(running)) {
    r <- running[[id]]
    if(!r$process$is_alive()) {
      r$process$wait();finish_job(as.integer(id),r$process$get_exit_status(),r$started,r$before)
      running[[id]] <<- NULL
    }
  }
  progress <- function() atomic_write_csv(res,file.path(dir,'progress.csv'))
  tryCatch({
    while(next_job<=nrow(q) || length(running)) {
      reap()
      while(length(running)<max_procs && next_job<=nrow(q)) {start_job(next_job);next_job <- next_job+1L}
      max_seen <- max(max_seen,length(running))
      progress()
      if(length(running)) Sys.sleep(poll)
    }
  },error=function(e) {
    error <<- conditionMessage(e)
    cat(launch_now(),'LAUNCHER ERROR:',error,'; starting nothing further and waiting for running jobs\n');flush.console()
    while(length(running)) {tryCatch(reap(),error=function(e) running <<- list());if(length(running)) Sys.sleep(poll)}
  })
  res$status[res$status %in% c('queued','running')] <- 'not run'
  counts <- table(factor(res$status,levels=c('completed','resumed','failed','quarantined','not run')))
  status <- if(is.null(error) && all(res$status %in% c('completed','resumed'))) 0L else 1L
  atomic_write_csv(res[SUMMARY_COLUMNS],file.path(dir,'summary.csv'))
  progress()
  atomic_write_lines(c(paste('finished:',launch_now()),paste('exit_status:',status),paste('jobs:',nrow(res)),
    paste0(sub(' ','_',names(counts)),': ',as.integer(counts)),paste('max_concurrent:',max_seen),
    paste('launcher_error:',error %||% '')),file.path(dir,'DONE'))
  cat(launch_now(),'DONE:',paste(names(counts),as.integer(counts),collapse=', '),'; exit status',status,'\n')
  list(status=status,dir=dir,max_concurrent=max_seen)
}

launch_main <- function(args,script=NULL) {
  tryCatch({
    opt <- parse_launch_args(args)
    repo <- normalizePath(opt$repo,mustWork=TRUE)
    run_script <- check_launch_repo(repo,script)
    study <- normalizePath(opt$study,mustWork=TRUE);archives <- normalizePath(opt$archives,mustWork=TRUE)
    inputs_root <- normalizePath(opt$inputs_root,mustWork=TRUE)
    only <- if(is.null(opt$only)) NULL else read_only_jobs(opt$only)
    q <- build_queue(opt$phases,opt$sds,study,archives,only,opt$schedule_override)
    if(!nrow(q)) {cat('Queue is empty; nothing to run.\n');return(0L)}
    counts <- queue_counts(q)
    cat('Queue:',nrow(q),'jobs, at most',opt$max_procs,'concurrent; planning estimate',
      round(estimated_wall_seconds(q$expected_seconds,opt$max_procs)/3600,1),'h wall time\n')
    cat(paste0('  ',counts$phase,' sd',format(counts$sd),' ',counts$schedule,' ',counts$jobs,'\n'),sep='')
    cat('Existing fits (run.R resumes an identical fit and refuses any other):',sum(file.exists(q$fit)),'\n')
    if(refuse_locks(q)) return(1L)
    revision <- check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
    cat('Library fingerprint matches revision',revision,'\n')
    if(!is.null(opt$queue_out)) utils::write.csv(queue_table(q,study),opt$queue_out,row.names=FALSE)
    if(opt$dry_run) {
      print(queue_table(q,study)[c('order','phase','sd','schedule','key','expected_seconds','fit_exists')],row.names=FALSE)
      cat('Dry run: nothing started.\n')
      return(0L)
    }
    res <- execute_queue(q,study,opt$max_procs,run_r_command(run_script,repo,study,inputs_root,archives),poll=opt$poll)
    res$status
  },error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]))
  repo_arg <- sub('^--repo=','',grep('^--repo=',args,value=TRUE))
  if(length(repo_arg)!=1L) {cat('ERROR: missing or repeated --repo\n');quit(save='no',status=1)}
  expected <- file.path(normalizePath(repo_arg,mustWork=FALSE),LAUNCH_SCRIPT)
  # Check this before sourcing anything from --repo.
  if(!identical(script,normalizePath(expected,mustWork=FALSE))) {
    cat('ERROR: This launcher (',script,') is not the launcher of --repo (',expected,
      '); run the launcher of the checkout whose run.R it starts\n',sep='')
    quit(save='no',status=1)
  }
  here <- dirname(script)
  source(file.path(here,'jobs.R'));source(file.path(here,'verify-helpers.R'))
  quit(save='no',status=launch_main(args,script))
}
