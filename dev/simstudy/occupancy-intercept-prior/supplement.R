#!/usr/bin/env Rscript
# The documented manual record of R27 (AMENDMENT-2.md): the valid selected fits
# of an SD that select.R could not record, because its check_first_fits refuses
# while any first fit of a requested SD is missing. Frozen select.R and flags.R
# are not edited; this script applies their frozen functions instead.
#
#   Rscript supplement.R --repo=REPO --study=STUDY --inputs-root=DIR --phase=A1|A2
#     [--archives=DIR] [--workers=N] [--out=DIR]
#
# Reads two files the controller writes in STUDY/selection/<phase>:
#   manual-missing.csv   phase, sd, key, reason: the communities without a
#                        valid selected fit, with the log or quarantine file
#   manual-selected.csv  phase, sd, key, schedule, fit, fit_md5: every other
#                        community of each such SD, with the fit the selection
#                        rule selects (fit relative to STUDY)
# An SD in manual-selected.csv must be absent from selected-fits.csv, have at
# least one manual-missing row, and cover every community of the phase exactly
# once with its manual-missing rows. For every listed fit, and for the first fit
# of a longer repeat, it applies flags.R's evaluate_item with the frozen flag
# rule fit_flags, exactly as select.R does, and checks the README's selection
# rule (R12): an initial-schedule fit is selected only if unflagged; a longer
# repeat of an initial-schedule community only if its first fit is flagged; a
# first fit at the longer schedule as it is. It writes, once (a recomputation
# must be identical), OUT/<phase>/manual-flags.csv (the flag rows, in the
# format of select.R's flag tables) and OUT/<phase>/manual-selection.csv (the
# format of selected-fits.csv, plus first_fit and first_fit_md5), and then
# scores every listed fit with analysis.R's score_item. OUT defaults to
# STUDY/summary. Like analysis.R it refuses to run before the other
# sub-phase's selection exists (AMENDMENT-1.md, R29). Such an SD fails criterion 4 of its sub-phase (R27); its
# criteria 1 to 3 are reported descriptively. Not sourced or hashed by run.R.

read_manual_selected <- function(sel_dir,phase) {
  path <- file.path(sel_dir,'manual-selected.csv')
  if(!file.exists(path)) return(NULL)
  x <- utils::read.csv(path,stringsAsFactors=FALSE,colClasses=c(phase='character',sd='numeric',key='character',
    schedule='character',fit='character',fit_md5='character'))
  if(!identical(names(x),c('phase','sd','key','schedule','fit','fit_md5')) || any(x$phase!=phase) ||
     !all(x$key %in% phase_keys(phase)) || !all(x$sd %in% setdiff(ARMS,1)) || !all(x$schedule %in% c('initial','long')) ||
     anyDuplicated(x[c('sd','key')]) || any(!grepl('^[0-9a-f]{32}$',x$fit_md5)) ||
     !identical(x$fit,sprintf('fits/%s/sd%s/%s/%s-fit.rds',phase,format(x$sd),x$schedule,x$key)))
    stop('Malformed manual record ',path,': columns phase, sd, key, schedule, fit, fit_md5; new-arm SDs, phase ',phase,
      ' keys, the protocol fit path of each schedule and an md5 per fit')
  x
}

# The rules of the manual record, given the selected-fits.csv SDs and the
# control schedule of every community.
check_supplement <- function(phase,ms,missing,selected_sds,control_schedule) {
  keys <- phase_keys(phase)
  for(sd in unique(ms$sd)) {
    if(sd %in% selected_sds) stop('SD ',sd,' is in selected-fits.csv; the manual record is only for an SD select.R could not record')
    k <- ms$key[ms$sd==sd];gone <- missing$key[missing$sd==sd]
    if(!length(gone)) stop('SD ',sd,' has no manual-missing.csv row; without one the SD belongs in select.R')
    if(length(intersect(k,gone))) stop('Community both selected and recorded as missing for SD ',sd,': ',intersect(k,gone)[1])
    if(!setequal(c(k,gone),keys) || length(c(k,gone))!=length(keys)) stop('SD ',sd,' does not cover every community of phase ',phase,' exactly once')
    s <- ms[ms$sd==sd,,drop=FALSE];cs <- control_schedule[s$key]
    if(any(cs=='long' & s$schedule!='long')) stop('SD ',sd,': a community with a longer-schedule control is fitted at the longer schedule only')
  }
  invisible(TRUE)
}

# Flag items of select.R's form: each listed fit, and the first fit (at the
# initial schedule) of every longer repeat of an initial-schedule community.
supplement_flag_items <- function(phase,ms,study,control_schedule) {
  item <- function(sd,key,schedule,fit,md5,kind) data.frame(role='new',phase=phase,sd=sd,key=key,schedule=schedule,
    fit=file.path(study,fit),fit_label=fit,expected_md5=md5,kind=kind,stringsAsFactors=FALSE)
  rows <- lapply(seq_len(nrow(ms)),function(i) item(ms$sd[i],ms$key[i],ms$schedule[i],ms$fit[i],ms$fit_md5[i],'selected'))
  rep <- ms[ms$schedule=='long' & control_schedule[ms$key]=='initial',,drop=FALSE]
  rows <- c(rows,lapply(seq_len(nrow(rep)),function(i) item(rep$sd[i],rep$key[i],'initial',
    sprintf('fits/%s/sd%s/initial/%s-fit.rds',phase,format(rep$sd[i]),rep$key[i]),NA_character_,'first')))
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# The selection table, after applying the selection rule to the flags.
supplement_table <- function(phase,ms,flags,control_schedule) {
  one <- function(i) {
    r <- ms[i,];f <- flags[flags$sd==r$sd & flags$key==r$key & flags$schedule==r$schedule,,drop=FALSE]
    if(nrow(f)!=1L || !identical(f$fit_md5,r$fit_md5)) stop('No flags of the listed fit ',r$fit)
    cs <- control_schedule[[r$key]];first <- NULL
    if(r$schedule=='initial' && isTRUE(f$flagged)) stop('Listed initial fit ',r$fit,' is flagged; the selection rule selects its single longer repeat')
    if(r$schedule=='long' && cs=='initial') {
      first <- flags[flags$sd==r$sd & flags$key==r$key & flags$schedule=='initial',,drop=FALSE]
      if(nrow(first)!=1L || !isTRUE(first$flagged)) stop('Listed longer repeat ',r$fit,' has no flagged first fit; the selection rule keeps an unflagged first fit')
    }
    data.frame(phase=phase,sd=r$sd,key=r$key,role='new',first_schedule=if(is.null(first)) r$schedule else 'initial',
      selected_schedule=r$schedule,long_repeat=!is.null(first),fit=r$fit,fit_md5=r$fit_md5,rule=f$rule,
      first_flagged=if(is.null(first)) f$flagged else first$flagged,flagged=f$flagged,reasons=f$reasons,
      first_fit=if(is.null(first)) NA_character_ else first$fit,first_fit_md5=if(is.null(first)) NA_character_ else first$fit_md5,
      stringsAsFactors=FALSE)
  }
  x <- do.call(rbind,lapply(seq_len(nrow(ms)),one))
  x <- x[order(x$sd,match(x$key,phase_keys(phase))),,drop=FALSE];rownames(x) <- NULL;x
}

control_schedule_of <- function(phase,archives) {cs <- control_schedules(phase,archives);stats::setNames(cs$schedule,cs$key)}

# flag_fn and score_fn are fit_flags and score_item except in tests.
supplement_main <- function(args,flag_fn=fit_flags,score_fn=score_item) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','phase','workers','out'),
    required=c('repo','study','inputs-root','phase'))
  phase <- o$phase;if(!phase %in% SCORED_PHASES) stop('--phase must be A1 or A2')
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE);inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  workers <- parse_workers(o$workers %||% '1');out <- o$out %||% file.path(study,'summary')
  sel_dir <- selection_dir(study,phase)
  # It scores fits, so AMENDMENT-1's rule holds here too (R29): not before the other sub-phase's selection exists.
  other <- file.path(selection_dir(study,setdiff(SCORED_PHASES,phase)),'selected-fits.csv')
  if(!file.exists(other)) stop('AMENDMENT-1.md: no new-arm occupancy error is computed before both selections exist, and ',other,' is missing')
  ms <- read_manual_selected(sel_dir,phase);if(is.null(ms)) stop('No manual record ',file.path(sel_dir,'manual-selected.csv'))
  missing <- read_manual_missing(sel_dir,phase)
  sf <- file.path(sel_dir,'selected-fits.csv');selected_sds <- if(file.exists(sf)) unique(read_selected(sf)$sd) else numeric()
  cs <- control_schedule_of(phase,archives)
  check_supplement(phase,ms,missing,selected_sds,cs)
  items <- supplement_flag_items(phase,ms,study,cs)
  flags <- evaluate_items(items,load_flag_scorers(repo,archives),archives,inputs_root,workers,flag_fn)
  table <- supplement_table(phase,ms,flags,cs)
  d <- file.path(out,phase);dir.create(d,recursive=TRUE,showWarnings=FALSE)
  cat('manual-flags.csv',write_frozen_table(flags,file.path(d,'manual-flags.csv')),'\n')
  cat('manual-selection.csv',write_frozen_table(table,file.path(d,'manual-selection.csv')),'\n')
  sc <- load_metric_scorers(repo,archives)
  scored <- items;scored$expected_md5 <- flags$fit_md5[match(paste(scored$sd,scored$key,scored$schedule),paste(flags$sd,flags$key,flags$schedule))]
  for(i in seq_len(nrow(scored))) score_fn(scored[i,,drop=FALSE],sc,archives,inputs_root,out)
  cat('Recorded and scored',nrow(ms),'manually selected fits of phase',phase,'(SD',paste(unique(ms$sd),collapse=', '),
    ') in',d,'; each such SD fails criterion 4 (R27)\n')
  0L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  repo_arg <- sub('^--repo=','',grep('^--repo=',commandArgs(trailingOnly=TRUE),value=TRUE))
  if(length(repo_arg)!=1L || !identical(normalizePath(file.path(repo_arg,'dev/simstudy/occupancy-intercept-prior'),mustWork=FALSE),here)) {
    cat('ERROR: run the supplement.R of the --repo checkout (',here,' is not under --repo)\n',sep='');quit(save='no',status=1)
  }
  for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R')) source(file.path(here,f))
  status <- tryCatch(supplement_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
