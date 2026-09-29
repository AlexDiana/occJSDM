#!/usr/bin/env Rscript
# Convergence selection for the occupancy-intercept prior study, recorded
# before any error is compared. Applies the frozen flag rules of flags.R to
# saved fits, identically for the new arms and the reused controls.
#
#   Rscript select.R --repo=REPO --study=STUDY --inputs-root=DIR --phases=A1,A2
#     --sds=2,3,5 --mode=plan|final [--workers=N] [--archives=DIR] [--out=DIR]
#
# plan, once the launcher's DONE marker reports exit status 0, writes to OUT
# (default STUDY/selection/<phases>, for example STUDY/selection/A1; see
# AMENDMENT-1.md), in this order:
#   long-selection.csv  every new-arm first fit at the initial schedule, flagged or not
#   long-keys.txt       the flagged ones, as a launcher job list for --only
#   long-fit-flags.csv  new-arm first fits already at the longer schedule; no
#                       further escalation, so their flags are retained
#   control-flags.csv   each control's selected fit under the same rule
# final, once the single longer repeats are complete, writes
#   repeat-flags.csv    each longer repeat
#   selected-fits.csv   each arm's selected fit per community, with its flags
#   convergence.csv     selected fits still flagged, per phase, A1 site count and arm
# A table that exists is never replaced: a recomputation must be byte-identical.
# Nothing here computes or reads an occupancy error. Not hashed into fits.

read_flag_table <- function(path) {
  if(!file.exists(path)) stop('Required selection table missing: ',path)
  utils::read.csv(path,stringsAsFactors=FALSE,colClasses=c(role='character',phase='character',key='character',
    schedule='character',fit='character',fit_md5='character',rule='character',reasons='character'))
}

empty_flag_table <- function() {
  x <- data.frame(role=character(),phase=character(),sd=numeric(),key=character(),schedule=character(),fit=character(),
    fit_md5=character(),rule=character(),warnings=integer(),max_group_rhat=numeric(),max_element_rhat=numeric(),
    unresolved_rhat=integer(),spatial_trace_flags=integer(),spatial_field_flags=integer(),flagged=logical(),
    reasons=character(),stringsAsFactors=FALSE)
  stopifnot(identical(names(x),c(ITEM_COLUMNS,FLAG_COLUMNS)));x
}

select_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','phases','sds','mode','workers','out'),
    required=c('repo','study','inputs-root','phases','sds','mode'))
  phases <- split_list(o$phases,'phases')
  if(!all(phases %in% PHASES)) stop('--phases must name phases among ',paste(PHASES,collapse=', '))
  sds <- split_list(o$sds,'sds')
  if(!all(sds %in% as.character(LAUNCH_SDS))) stop('--sds must list prior SDs among ',paste(LAUNCH_SDS,collapse=', '))
  sds <- as.numeric(sds)
  if(!o$mode %in% c('plan','final')) stop('--mode must be plan or final')
  workers <- parse_workers(o$workers %||% '1')
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  out <- o$out %||% selection_dir(study,phases)
  sc <- load_flag_scorers(repo,archives)
  items <- selection_items(phases,sds,study,archives,repo)
  check_first_fits(items)
  flags <- function(rows) {x <- evaluate_items(rows,sc,archives,inputs_root,workers);if(is.null(x)) empty_flag_table() else x}
  report <- function(name,status,x) cat(name,status,':',nrow(x),'fits,',sum(x$flagged),'flagged\n')
  new <- items$role=='new'
  if(o$mode=='plan') {
    # The long-repeat selection is fixed first, from convergence flags alone.
    initial <- flags(items[new & items$schedule=='initial',,drop=FALSE])
    report('long-selection.csv',write_frozen_table(initial,file.path(out,'long-selection.csv')),initial)
    write_frozen_table(long_keys(initial),file.path(out,'long-keys.txt'))
    if(nrow(initial)) print(initial[initial$flagged,c('phase','sd','key','reasons')],row.names=FALSE)
    longs <- flags(items[new & items$schedule=='long',,drop=FALSE])
    report('long-fit-flags.csv',write_frozen_table(longs,file.path(out,'long-fit-flags.csv')),longs)
    controls <- flags(items[!new,,drop=FALSE])
    report('control-flags.csv',write_frozen_table(controls,file.path(out,'control-flags.csv')),controls)
  } else {
    initial <- read_flag_table(file.path(out,'long-selection.csv'))
    longs <- read_flag_table(file.path(out,'long-fit-flags.csv'))
    controls <- read_flag_table(file.path(out,'control-flags.csv'))
    id <- function(x) paste(x$role,x$phase,x$sd,x$key,x$schedule)
    recorded <- rbind(initial,longs,controls)
    if(!setequal(id(recorded),id(items)) || anyDuplicated(id(recorded)))
      stop('The recorded plan tables do not cover exactly the requested phases and SDs; run --mode=plan first')
    now_md5 <- unname(tools::md5sum(items$fit[match(id(recorded),id(items))]))
    if(!identical(now_md5,recorded$fit_md5)) stop('A fit changed after the selection was recorded: ',
      recorded$fit[now_md5!=recorded$fit_md5][1])
    flagged <- initial[initial$flagged,,drop=FALSE]
    repeats <- items[0,]
    if(nrow(flagged)) {
      fits <- as.character(mapply(fit_path,study,flagged$phase,flagged$sd,'long',flagged$key,USE.NAMES=FALSE))
      repeats <- data.frame(role='new',phase=flagged$phase,sd=flagged$sd,key=flagged$key,schedule='long',fit=fits,
        fit_label=path_within(fits,study),expected_md5=NA_character_,stringsAsFactors=FALSE)
      check_first_fits(repeats)
    }
    repeat_flags <- flags(repeats)
    report('repeat-flags.csv',write_frozen_table(repeat_flags,file.path(out,'repeat-flags.csv')),repeat_flags)
    selected <- final_selected(initial,repeat_flags,longs,controls)
    report('selected-fits.csv',write_frozen_table(selected,file.path(out,'selected-fits.csv')),selected)
    counts <- convergence_counts(selected)
    write_frozen_table(counts,file.path(out,'convergence.csv'))
    print(counts,row.names=FALSE)
  }
  files <- c(file.path(repo,'dev/simstudy/occupancy-intercept-prior',c('select.R','flags.R','jobs.R','verify-helpers.R','launch.R')))
  provenance <- data.frame(kind=c(rep('script',length(files)),rep('scorer',length(sc$hashes)),'args','time'),
    name=c(path_within(files,repo),names(sc$hashes),'args','finished'),
    value=c(unname(tools::md5sum(files)),unname(sc$hashes),paste(args,collapse=' '),format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')))
  utils::write.csv(provenance,file.path(out,paste0('provenance-',o$mode,'-',format(Sys.time(),'%Y%m%d%H%M%S'),'.csv')),row.names=FALSE)
  cat('Selection',o$mode,'recorded in',out,'\n')
  0L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  repo_arg <- sub('^--repo=','',grep('^--repo=',commandArgs(trailingOnly=TRUE),value=TRUE))
  # The scripts sourced here and the scorer files read from --repo must be one checkout.
  if(length(repo_arg)!=1L || !identical(normalizePath(file.path(repo_arg,'dev/simstudy/occupancy-intercept-prior'),mustWork=FALSE),here)) {
    cat('ERROR: run the select.R of the --repo checkout (',here,' is not under --repo)\n',sep='');quit(save='no',status=1)
  }
  for(f in c('jobs.R','verify-helpers.R','launch.R','flags.R')) source(file.path(here,f))
  status <- tryCatch(select_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
