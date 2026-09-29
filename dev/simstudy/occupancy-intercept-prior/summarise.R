#!/usr/bin/env Rscript
# Tables and gate of the phase A scorer (README "Outcomes and gate";
# AMENDMENT-1.md R21; AMENDMENT-2.md R27).
#
#   Rscript summarise.R --repo=REPO --study=STUDY --phase=A1|A2 --arms=1[,2,3,5]
#     [--archives=DIR] [--inputs-root=DIR] [--out=DIR]
#   Rscript summarise.R --repo=REPO --study=STUDY --phase=A [--out=DIR]
#
# Reads only analysis.R's score records (OUT/<phase>/scores, OUT defaults to
# STUDY/summary) and the selection tables; loads no fit. Every record must be
# of the fit the selection names (md5 and path) and made by the current scorer
# code. Writes to OUT/<phase>:
#   band-error.csv      arm, stratum, band or group and community: signed error E,
#                       |E| and mean absolute cell error, in percentage points
#   coverage.csv        the same keys: 95% interval coverage
#   mae.csv             overall MAE M per community (scope primary is the gate's)
#   b0.csv              B0 bias, absolute bias and coverage per community (logit scale)
#   convergence.csv     selected fits and flags per stratum and arm, from the
#                       selection (controls: control-flags.csv); before a
#                       selection exists, controls only, from results/flag-crosscheck.csv
#   summary-means.csv   across-community means per arm, stratum, scope and group
#   provenance.csv      scorer hashes, arguments, records read
# and, when the control and a new arm are both requested:
#   gate.csv            one row per SD: the numbers of every criterion, PASS or FAIL
#   gate-detail.csv     the same, one row per criterion and component
#   schedule-matched.csv  criteria 1 to 3 with each new arm's first fits (R12), descriptive
# When the gates of both A1 and A2 exist, or with --phase=A, OUT/A/gate.csv
# holds the A1 and A2 rows and one phase A row per SD (passes A1 and A2).
# Not sourced or hashed by run.R.

GATE_BANDS <- c('low','middle','high')
# "At most" thresholds are compared with this allowance so that a difference
# equal to the threshold in exact arithmetic (for example coverage .53 - .50)
# is not failed by binary rounding. It is far below any attainable nonzero
# difference in these outcomes; strict comparisons (criterion 1) use none.
GATE_TOL <- 1e-9

gate_config <- function(phase) switch(phase,
  A1=list(strata=c('n100','n300'),c1_stratum='pooled',gated='low'),
  A2=list(strata='all',c1_stratum='all',gated=c('low','rare_1_5pct')),
  stop('No gate is defined here for phase ',phase))

ceil_two_thirds <- function(n) as.integer((2L*n+2L)%/%3L)

detail_row <- function(phase,sd,criterion,component,stratum,communities=NA,needed=NA,improved=NA,value_new=NA,
  value_control=NA,statistic=NA,threshold=NA,missing_new=NA,pass)
  data.frame(phase=phase,sd=sd,criterion=criterion,component=component,stratum=stratum,communities=as.integer(communities),
    needed=as.integer(needed),improved=as.integer(improved),value_new=as.numeric(value_new),value_control=as.numeric(value_control),
    statistic=as.numeric(statistic),threshold=as.numeric(threshold),missing_new=as.integer(missing_new),pass=pass,stringsAsFactors=FALSE)

# One quantity of one group and stratum for the fits valid in both arms (the
# same communities for both arms; a group empty in a community is left out).
paired <- function(x,sd,group,stratum,column) {
  d <- x[x$group==group & x$stratum==stratum & x$cells>0,,drop=FALSE]
  a <- d[d$sd==sd,,drop=FALSE];b <- d[d$sd==1,,drop=FALSE];keys <- intersect(a$key,b$key)
  list(new=a[[column]][match(keys,a$key)],control=b[[column]][match(keys,b$key)],n=length(keys))
}

# Criterion 1 for one gated group. A community's value is the mean |E| of its
# fits in the pooled strata: in A1 its n100 and n300 fits (R21), in A2 its one
# fit. Improved means strictly smaller than the control's; at least
# ceil(2n/3) of the n communities with the group in both arms, and a strictly
# smaller across-community mean.
criterion1 <- function(x,phase,sd,group) {
  cfg <- gate_config(phase)
  d <- x[x$group==group & x$cells>0 & x$stratum %in% cfg$strata,,drop=FALSE]
  value <- function(a) {
    y <- d[d$sd==a,,drop=FALSE]
    v <- vapply(split(y,y$community),function(z)
      if(setequal(z$stratum,cfg$strata) && !anyDuplicated(z$stratum)) mean(z$abs_signed_error) else NA_real_,numeric(1))
    v[!is.na(v)]
  }
  vn <- value(sd);vc <- value(1);ids <- intersect(names(vn),names(vc))
  n <- length(ids);needed <- ceil_two_thirds(n);improved <- sum(vn[ids]<vc[ids])
  mn <- if(n) mean(vn[ids]) else NA_real_;mc <- if(n) mean(vc[ids]) else NA_real_
  detail_row(phase,sd,'c1',group,cfg$c1_stratum,n,needed,improved,mn,mc,mn-mc,
    pass=n>0L && improved>=needed && mn<mc)
}

# Criterion 2 at one stratum: each band's mean |E| worse by at most 1 point and
# the mean overall MAE worse by at most 0.5 point. A band with no cells in any
# community cannot be worse and passes.
criterion2 <- function(x,phase,sd,stratum) {
  one <- function(group,column,component,limit) {
    p <- paired(x,sd,group,stratum,column)
    mn <- if(p$n) mean(p$new) else NA_real_;mc <- if(p$n) mean(p$control) else NA_real_
    detail_row(phase,sd,'c2',component,stratum,p$n,value_new=mn,value_control=mc,statistic=mn-mc,threshold=limit,
      pass=p$n==0L || mn-mc<=limit+GATE_TOL)
  }
  rbind(do.call(rbind,lapply(GATE_BANDS,function(g) one(g,'abs_signed_error',g,1))),one('all','mean_abs_cell_error','mae',.5))
}

# Criterion 3 at one stratum: each band's mean coverage lower by at most 0.03.
criterion3 <- function(x,phase,sd,stratum) do.call(rbind,lapply(GATE_BANDS,function(g) {
  p <- paired(x,sd,g,stratum,'coverage')
  mn <- if(p$n) mean(p$new) else NA_real_;mc <- if(p$n) mean(p$control) else NA_real_
  detail_row(phase,sd,'c3',g,stratum,p$n,value_new=mn,value_control=mc,statistic=mc-mn,threshold=.03,
    pass=p$n==0L || mc-mn<=.03+GATE_TOL)
}))

# Criterion 4 at one stratum: flagged selected fits of the new arm at most the
# control's (R20), and no community without a valid selected fit (R25, R27).
criterion4 <- function(conv,phase,sd,stratum,valid_new) {
  keys <- phase_keys(phase);keys <- keys[stratum_of(phase,keys)==stratum]
  missing <- length(setdiff(keys,valid_new))
  fc <- conv$selected_flagged[conv$sd==1 & conv$stratum==stratum]
  if(length(fc)!=1L || is.na(fc)) stop('No control convergence count for phase ',phase,' stratum ',stratum)
  fn <- conv$selected_flagged[conv$sd==sd & conv$stratum==stratum];fn <- if(length(fn)==1L) fn else NA_real_
  detail_row(phase,sd,'c4','flags',stratum,value_new=fn,value_control=fc,statistic=fn-fc,threshold=0,missing_new=missing,
    pass=missing==0L && isTRUE(fn<=fc))
}

gate_row <- function(d,phase,sd,descriptive=FALSE) {
  row <- list(sd=sd,subphase=phase)
  put <- function(prefix,names,values) row[paste0(prefix,names)] <<- values
  for(i in seq_len(nrow(d))) {
    r <- d[i,,drop=FALSE]
    switch(r$criterion,
      c1=put(paste0('c1_',r$component,'_'),c('communities','improved','needed','mean_new','mean_control','pass'),
        list(r$communities,r$improved,r$needed,r$value_new,r$value_control,r$pass)),
      c2=put(paste0('c2_',r$component,'_',r$stratum,'_'),c('new','control','diff','pass'),list(r$value_new,r$value_control,r$statistic,r$pass)),
      c3=put(paste0('c3_',r$component,'_',r$stratum,'_'),c('new','control','drop','pass'),list(r$value_new,r$value_control,r$statistic,r$pass)),
      c4=put(paste0('c4_',r$stratum,'_'),c('flagged_new','flagged_control','missing_new','pass'),
        list(r$value_new,r$value_control,r$missing_new,r$pass)))
  }
  for(k in c('c1','c2','c3','c4')) row[[paste0(k,'_pass')]] <- if(any(d$criterion==k)) all(d$pass[d$criterion==k]) else NA
  row$result <- if(descriptive) 'DESCRIPTIVE' else if(all(d$pass)) 'PASS' else 'FAIL'
  as.data.frame(row,stringsAsFactors=FALSE)
}

# One new arm against the control within one sub-phase. scores: primary rows
# of both arms (sd, key, stratum, community, group, cells, abs_signed_error,
# mean_abs_cell_error, coverage). conv: selected fits and flags per stratum and
# arm; NULL for a descriptive comparison without criterion 4. missing: the
# manual record (sd, key) of communities without a valid selected fit.
gate_subphase <- function(phase,scores,conv,sd,missing=NULL) {
  cfg <- gate_config(phase);keys <- phase_keys(phase)
  x <- scores[scores$scope=='primary' & scores$sd %in% c(1,sd),,drop=FALSE]
  if(!setequal(unique(x$key[x$sd==1]),keys)) stop('The control arm lacks a valid selected fit for some community of phase ',phase)
  gone <- if(is.null(missing)) character() else missing$key[missing$sd==sd]
  valid <- setdiff(unique(x$key[x$sd==sd]),gone)
  if(length(setdiff(valid,keys))) stop('Unknown community for SD ',sd,': ',setdiff(valid,keys)[1])
  x <- x[x$sd==1 | x$key %in% valid,,drop=FALSE]
  d <- do.call(rbind,c(lapply(cfg$gated,function(g) criterion1(x,phase,sd,g)),
    lapply(cfg$strata,function(h) criterion2(x,phase,sd,h)),lapply(cfg$strata,function(h) criterion3(x,phase,sd,h)),
    if(!is.null(conv)) lapply(cfg$strata,function(h) criterion4(conv,phase,sd,h,valid))))
  rownames(d) <- NULL
  list(detail=d,row=gate_row(d,phase,sd,descriptive=is.null(conv)))
}

rbind_fill <- function(xs) {
  xs <- Filter(Negate(is.null),xs);cols <- unique(unlist(lapply(xs,names)))
  x <- do.call(rbind,lapply(xs,function(x) {for(k in setdiff(cols,names(x))) x[[k]] <- NA;x[cols]}));rownames(x) <- NULL;x
}

# Phase A rows: an SD passes phase A only if it passes A1 and A2 (R9); the
# smallest passing SD is marked.
combine_phase_a <- function(a1,a2) {
  if(!setequal(a1$sd,a2$sd) || anyDuplicated(a1$sd) || anyDuplicated(a2$sd) || any(a1$subphase!='A1') || any(a2$subphase!='A2'))
    stop('Phase A needs both an A1 and an A2 gate row for each SD')
  sds <- sort(a1$sd);r1 <- a1$result[match(sds,a1$sd)];r2 <- a2$result[match(sds,a2$sd)];pass <- r1=='PASS' & r2=='PASS'
  a <- data.frame(sd=sds,subphase='A',a1_result=r1,a2_result=r2,result=ifelse(pass,'PASS','FAIL'),
    smallest_passing=pass & sds==(if(any(pass)) min(sds[pass]) else -Inf),stringsAsFactors=FALSE)
  rbind_fill(list(a1,a2,a))
}

# ---------------------------------------------------------------------------
# Records and selection tables
# ---------------------------------------------------------------------------

read_record <- function(path,md5,label,hashes,role,sd,kind) {
  if(!file.exists(path)) stop('Not scored: ',path,'; run analysis.R first')
  r <- readRDS(path)
  if(!identical(r$version,SCORE_VERSION)) stop('Score record ',path,' has version ',r$version %||% 'none')
  if(!identical(r$fit_md5,md5) || !identical(r$fit_label,label) || !identical(r$role,role) ||
     !identical(r$sd,as.numeric(sd)) || !identical(r$kind,kind))
    stop('Score record ',path,' is not of the selected fit ',label,' (md5 ',md5,')')
  if(!identical(r$hashes,hashes)) stop('Score record ',path,' was made by other scorer code; rescore with analysis.R')
  r
}

same_table <- function(a,b,cols) all(vapply(cols,function(k)
  if(is.numeric(a[[k]]) || is.numeric(b[[k]]) || is.logical(a[[k]])) identical(as.numeric(a[[k]]),as.numeric(b[[k]])) else
    identical(as.character(a[[k]]),as.character(b[[k]])),logical(1)))

convergence_table <- function(phase,arms,sel_dir,repo,missing) {
  keys <- phase_keys(phase);strata <- unique(stratum_of(phase,keys))
  cols <- c('phase','stratum','sd','role','selected_fits','selected_flagged','first_flagged','long_repeats')
  f <- file.path(sel_dir,'convergence.csv')
  if(file.exists(f)) {
    cv <- utils::read.csv(f,stringsAsFactors=FALSE)
    again <- convergence_counts(read_selected(file.path(sel_dir,'selected-fits.csv')))
    if(nrow(cv)!=nrow(again) || !same_table(cv,again,cols)) stop('convergence.csv disagrees with selected-fits.csv in ',sel_dir)
    cf <- utils::read.csv(file.path(sel_dir,'control-flags.csv'),stringsAsFactors=FALSE)
    for(h in strata) if(!identical(as.numeric(sum(cf$flagged[stratum_of(phase,cf$key)==h])),
      as.numeric(cv$selected_flagged[cv$sd==1 & cv$stratum==h]))) stop('control-flags.csv disagrees with convergence.csv in ',sel_dir)
    source <- 'selection'
  } else {
    if(any(arms!=1)) stop('No selection recorded in ',sel_dir,': only the control arm can be summarised')
    fc <- utils::read.csv(file.path(repo,ANALYSIS_REL,'results/flag-crosscheck.csv'),stringsAsFactors=FALSE)
    fc <- fc[fc$phase==phase & fc$selected,,drop=FALSE]
    if(!setequal(fc$key,keys) || anyDuplicated(fc$key) || !all(fc$match)) stop('results/flag-crosscheck.csv does not cover the phase ',phase,' controls')
    cv <- do.call(rbind,lapply(strata,function(h) {z <- fc[stratum_of(phase,fc$key)==h,]
      data.frame(phase=phase,stratum=h,sd=1,role='control',selected_fits=nrow(z),selected_flagged=sum(z$flagged),
        first_flagged=NA_integer_,long_repeats=NA_integer_,stringsAsFactors=FALSE)}))
    source <- 'results/flag-crosscheck.csv (flags.R on the control fits; control selection not yet recorded for this phase)'
  }
  cv <- cv[cv$sd %in% arms,cols,drop=FALSE]
  for(sd in setdiff(arms,unique(cv$sd))) cv <- rbind(cv,data.frame(phase=phase,stratum=strata,sd=sd,role='new',selected_fits=0L,
    selected_flagged=NA_integer_,first_flagged=NA_integer_,long_repeats=NA_integer_,stringsAsFactors=FALSE))
  cv$expected_fits <- vapply(cv$stratum,function(h) sum(stratum_of(phase,keys)==h),integer(1))
  cv$manual_missing <- vapply(seq_len(nrow(cv)),function(i) sum(missing$sd==cv$sd[i] & stratum_of(phase,missing$key)==cv$stratum[i]),integer(1))
  cv$missing_fits <- cv$expected_fits-cv$selected_fits
  cv$source <- source;rownames(cv) <- NULL
  cv[order(cv$sd,cv$stratum),,drop=FALSE]
}

# Records of the requested arms, each checked against the fit it must be of.
collect_records <- function(phase,arms,study,archives,repo,out,hashes,sel_dir,missing) {
  recs <- list()
  if(1 %in% arms) {
    ctl <- control_items(phase,study,archives,repo)
    for(i in seq_len(nrow(ctl))) recs[[length(recs)+1L]] <- read_record(score_path(out,phase,1,ctl$key[i],ctl$schedule[i]),
      ctl$expected_md5[i],ctl$fit_label[i],hashes,'control',1,'selected')
  }
  for(sd in setdiff(arms,1)) {
    x <- new_items(phase,sd,study,sel_dir)
    for(i in seq_len(nrow(x))) recs[[length(recs)+1L]] <- read_record(score_path(out,phase,sd,x$key[i],x$schedule[i]),
      x$expected_md5[i],x$fit_label[i],hashes,'new',sd,x$kind[i])
    gone <- setdiff(phase_keys(phase),x$key[x$kind=='selected'])
    if(length(gone)) cat('R27: SD',sd,'has no valid selected fit for',length(gone),'communities:',paste(gone,collapse=', '),'\n')
  }
  recs
}

record_table <- function(recs,part) {
  meta <- function(r) data.frame(phase=r$phase,sd=r$sd,arm=r$role,kind=r$kind,key=r$key,stratum=r$stratum,community=r$community,
    schedule=r$schedule,fit=r$fit_label,fit_md5=r$fit_md5,stringsAsFactors=FALSE)
  x <- do.call(rbind,lapply(recs,function(r) cbind(meta(r)[rep(1L,nrow(r[[part]])),,drop=FALSE],r[[part]],row.names=NULL)))
  rownames(x) <- NULL;x
}

summary_means <- function(g) {
  by <- interaction(g$sd,g$stratum,g$scope,g$group,drop=TRUE,lex.order=TRUE)
  x <- do.call(rbind,lapply(split(g,by),function(z) data.frame(phase=z$phase[1],sd=z$sd[1],arm=z$arm[1],stratum=z$stratum[1],
    scope=z$scope[1],group=z$group[1],communities=nrow(z),mean_signed_error=mean(z$signed_error),
    mean_abs_signed_error=mean(z$abs_signed_error),mean_abs_cell_error=mean(z$mean_abs_cell_error),mean_coverage=mean(z$coverage),
    stringsAsFactors=FALSE)))
  rownames(x) <- NULL;x
}

write_table <- function(x,dir,name) {
  path <- file.path(dir,name);tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE);stopifnot(file.rename(tmp,path));path
}

summarise_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','archives','inputs-root','phase','arms','out'),required=c('repo','study','phase'))
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  out <- o$out %||% file.path(study,'summary')
  if(o$phase=='A') {
    g <- lapply(SCORED_PHASES,function(p) {f <- file.path(out,p,'gate.csv');if(!file.exists(f)) stop('Missing ',f)
      utils::read.csv(f,stringsAsFactors=FALSE)})
    x <- combine_phase_a(g[[1]],g[[2]]);dir.create(file.path(out,'A'),recursive=TRUE,showWarnings=FALSE)
    cat('Wrote',write_table(x,file.path(out,'A'),'gate.csv'),'\n');print(x[x$subphase=='A',],row.names=FALSE)
    return(0L)
  }
  if(!o$phase %in% SCORED_PHASES) stop('--phase must be A1, A2 or A')
  if(is.null(o$arms)) stop('Missing --arms')
  phase <- o$phase;arms <- parse_arms(o$arms);new_sds <- setdiff(arms,1)
  sel_dir <- selection_dir(study,phase)
  missing <- if(dir.exists(sel_dir)) read_manual_missing(sel_dir,phase) else read_manual_missing(tempfile(),phase)
  hashes <- current_score_hashes(repo,archives)
  recs <- collect_records(phase,arms,study,archives,repo,out,hashes,sel_dir,missing)
  selected <- Filter(function(r) r$kind=='selected',recs)
  # Per-community tables hold every scored fit, labelled by kind: 'selected',
  # and 'first' for a first fit replaced by its longer repeat (sensitivity only).
  every <- record_table(recs,'groups');groups <- every[every$kind=='selected',,drop=FALSE]
  meta <- c('phase','sd','arm','kind','key','stratum','community','schedule','fit','fit_md5')
  d <- file.path(out,phase);dir.create(d,recursive=TRUE,showWarnings=FALSE)
  write_table(every[c(meta,'scope','group','cells','truth_mean','estimate_mean','signed_error','abs_signed_error','mean_abs_cell_error')],d,'band-error.csv')
  write_table(every[c(meta,'scope','group','cells','coverage')],d,'coverage.csv')
  m <- every[every$group=='all',c(meta,'scope','cells','mean_abs_cell_error')];names(m)[names(m)=='mean_abs_cell_error'] <- 'mae'
  write_table(m,d,'mae.csv')
  write_table(record_table(recs,'b0'),d,'b0.csv')
  conv <- convergence_table(phase,arms,sel_dir,repo,missing)
  write_table(conv,d,'convergence.csv')
  write_table(summary_means(groups),d,'summary-means.csv')
  if(1 %in% arms && length(new_sds)) {
    rows <- list();details <- list();matched <- list()
    for(sd in new_sds) {
      g <- gate_subphase(phase,groups,conv,sd,missing)
      rows[[length(rows)+1L]] <- g$row;details[[length(details)+1L]] <- g$detail
      # Schedule-matched sensitivity (R12): the first fit wherever a longer repeat replaced it.
      first <- Filter(function(r) r$kind=='first' && r$sd==sd,recs)
      sm <- groups
      if(length(first)) {
        f <- record_table(first,'groups');sm <- rbind(sm[!(sm$sd==sd & sm$key %in% f$key),,drop=FALSE],f[names(sm)])
      }
      s <- gate_subphase(phase,sm,NULL,sd,missing)$detail
      s$first_fits_used <- length(first);matched[[length(matched)+1L]] <- s
    }
    gate <- rbind_fill(rows)
    write_table(gate,d,'gate.csv');write_table(do.call(rbind,details),d,'gate-detail.csv')
    write_table(do.call(rbind,matched),d,'schedule-matched.csv')
    print(gate[c('sd','subphase',grep('_pass$',names(gate),value=TRUE),'result')],row.names=FALSE)
    other <- file.path(out,setdiff(SCORED_PHASES,phase),'gate.csv')
    if(file.exists(other)) {
      o2 <- utils::read.csv(other,stringsAsFactors=FALSE)
      if(setequal(o2$sd,gate$sd)) {
        both <- if(phase=='A1') combine_phase_a(gate,o2) else combine_phase_a(o2,gate)
        dir.create(file.path(out,'A'),recursive=TRUE,showWarnings=FALSE);write_table(both,file.path(out,'A'),'gate.csv')
        print(both[both$subphase=='A',],row.names=FALSE)
      }
    }
  }
  prov <- data.frame(kind=c(rep('hash',length(hashes)),'args','finished',rep('record',length(recs))),
    name=c(names(hashes),'args','time',vapply(recs,function(r) paste(r$sd,r$kind,r$key,r$schedule),'')),
    value=c(unname(hashes),paste(args,collapse=' '),format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),vapply(recs,function(r) r$fit_md5,'')),
    stringsAsFactors=FALSE)
  write_table(prov,d,'provenance.csv')
  cat('Summarised',length(selected),'selected fits of phase',phase,'(arms',paste(arms,collapse=', '),') in',d,'\n')
  0L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  repo_arg <- sub('^--repo=','',grep('^--repo=',commandArgs(trailingOnly=TRUE),value=TRUE))
  if(length(repo_arg)!=1L || !identical(normalizePath(file.path(repo_arg,'dev/simstudy/occupancy-intercept-prior'),mustWork=FALSE),here)) {
    cat('ERROR: run the summarise.R of the --repo checkout (',here,' is not under --repo)\n',sep='');quit(save='no',status=1)
  }
  for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R')) source(file.path(here,f))
  status <- tryCatch(summarise_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
