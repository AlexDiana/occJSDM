#!/usr/bin/env Rscript
# Verification that flags.R reproduces the archived convergence diagnostics.
#
#   Rscript check-flags.R --repo=REPO --archives=DIR --inputs-root=DIR [--workers=N] [--out=FILE]
#
# Recomputes the flags of all 70 reused control fits in results/control-provenance.csv
# from their saved draws and compares them with what each archive recorded when
# it selected them: pr11 long-selection.csv (A1 and B initial fits: warnings,
# maximum group and element Rhat, unresolved Rhat and the flag), pr11
# selected-manifest.csv (B longer fits), and spatial-amplitude robust-v1
# selection.csv (A2 initial inverse-gamma reasons) and fits.csv (A2 selected
# longer inverse-gamma reasons). Rhat must agree within 1e-12, counts, flags and
# reason strings exactly. Writes results/flag-crosscheck.csv by default and
# exits non-zero on any disagreement. Uses no occupancy error: of the files it
# reads, only the spatial-amplitude fits.csv holds any (archived control
# occupancy columns), and only its reason strings are used.

check_main <- function(args) {
  o <- parse_options(args,known=c('repo','archives','inputs-root','workers','out'),required=c('repo','archives','inputs-root'))
  repo <- normalizePath(o$repo,mustWork=TRUE);archives <- normalizePath(o$archives,mustWork=TRUE)
  inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE);workers <- parse_workers(o$workers %||% '1')
  out <- o$out %||% file.path(repo,'dev/simstudy/occupancy-intercept-prior/results/flag-crosscheck.csv')
  sc <- load_flag_scorers(repo,archives)
  prov <- utils::read.csv(file.path(repo,CONTROL_PROVENANCE),stringsAsFactors=FALSE)
  items <- data.frame(role='control',phase=prov$phase,sd=1,key=prov$key,schedule=prov$schedule,
    fit=file.path(archives,prov$control_fit),fit_label=prov$control_fit,expected_md5=prov$file_md5,stringsAsFactors=FALSE)
  flags <- evaluate_items(items,sc,archives,inputs_root,workers)
  pr11 <- archive_dir(archives,'pr11')
  initial <- utils::read.csv(file.path(pr11,'long-selection.csv'),stringsAsFactors=FALSE)
  selected <- utils::read.csv(file.path(pr11,'selected-manifest.csv'),stringsAsFactors=FALSE)
  final <- file.path(archive_dir(archives,'amplitude'),'robust-v1/summary-binary-final')
  a2_initial <- utils::read.csv(file.path(final,'selection.csv'),stringsAsFactors=FALSE,colClasses=c(initial_ig_reasons='character'))
  a2_long <- utils::read.csv(file.path(final,'fits.csv'),stringsAsFactors=FALSE,colClasses=c(reasons='character'))
  a2_long <- a2_long[a2_long$prior=='inverse_gamma',]
  rows <- lapply(seq_len(nrow(flags)),function(i) {
    f <- flags[i,];base <- data.frame(phase=f$phase,key=f$key,schedule=f$schedule,selected=prov$selected[i],rule=f$rule,
      flagged=f$flagged,reasons=f$reasons,stringsAsFactors=FALSE)
    if(f$rule=='pr11') {
      a <- if(f$schedule=='initial') initial[initial$key==f$key,] else selected[selected$key==f$key & selected$schedule=='long',]
      stopifnot(nrow(a)==1L)
      d <- list(max_group_rhat=a$max_group_rhat,max_element_rhat=a$max_element_rhat,unresolved_rhat=as.integer(a$unresolved_rhat))
      af <- pr11_rule(rep('w',a$warnings),d,sc$flag_expr)
      if(f$schedule=='initial') stopifnot(identical(af$flagged,a$current_flag))
      gd <- abs(f$max_group_rhat-a$max_group_rhat);ed <- abs(f$max_element_rhat-a$max_element_rhat)
      cbind(base,source=if(f$schedule=='initial') 'pr11 long-selection.csv' else 'pr11 selected-manifest.csv',
        archived_flagged=af$flagged,archived_reasons=paste(af$reasons,collapse='; '),group_rhat_difference=gd,element_rhat_difference=ed,
        match=identical(f$flagged,af$flagged) && f$warnings==a$warnings && f$unresolved_rhat==a$unresolved_rhat && gd<=1e-12 && ed<=1e-12)
    } else {
      reasons <- if(f$schedule=='initial') a2_initial$initial_ig_reasons[a2_initial$key==f$key] else a2_long$reasons[a2_long$key==f$key]
      stopifnot(length(reasons)==1L)
      cbind(base,source=if(f$schedule=='initial') 'robust-v1 selection.csv' else 'robust-v1 fits.csv',archived_flagged=nzchar(reasons),
        archived_reasons=reasons,group_rhat_difference=NA_real_,element_rhat_difference=NA_real_,
        match=identical(f$reasons,reasons) && identical(f$flagged,nzchar(reasons)))
    }
  })
  x <- do.call(rbind,rows)
  utils::write.csv(x,out,row.names=FALSE)
  cat('Recomputed flags of',nrow(x),'control fits;',sum(x$match),'match the archived record;',
    sum(x$flagged),'flagged. Largest Rhat difference',format(max(c(x$group_rhat_difference,x$element_rhat_difference),na.rm=TRUE)),'\n')
  if(!all(x$match)) {print(x[!x$match,],row.names=FALSE);return(1L)}
  0L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  for(f in c('jobs.R','verify-helpers.R','launch.R','flags.R')) source(file.path(here,f))
  status <- tryCatch(check_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
