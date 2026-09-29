#!/usr/bin/env Rscript
# Tables and gate of the scorer (README "Outcomes and gate"; AMENDMENT-1.md
# R21; AMENDMENT-2.md R24 and R27).
#
#   Rscript summarise.R --repo=REPO --study=STUDY --phase=A1|A2|B --arms=1[,2,3,5]
#     [--archives=DIR] [--inputs-root=DIR] [--out=DIR]
#   Rscript summarise.R --repo=REPO --study=STUDY --phase=A|final [--out=DIR]
#
# Reads only the score records of analysis.R (A1, A2) and analysis-b.R (B)
# (OUT/<phase>/scores, OUT defaults to STUDY/summary) and the selection tables;
# loads no fit. Every record must be of the fit the selection names (md5 and
# path) and made by the current scorer code of its phase. Writes to OUT/<phase>:
#   band-error.csv      arm, stratum, band or group and community: signed error E,
#                       |E| and mean absolute cell error, in percentage points
#   coverage.csv        the same keys: 95% interval coverage
#   mae.csv             overall MAE M per community (scope primary is the gate's)
#   b0.csv              B0 bias, absolute bias and coverage per community (logit scale)
#   convergence.csv     selected fits and flags per stratum and arm, from the
#                       selection (controls: control-flags.csv); for an SD that
#                       select.R could not record, from supplement.R's manual
#                       record (R27), or none; before a selection exists,
#                       controls only, from results/flag-crosscheck.csv. In
#                       phase B the strata are the contamination levels qnear
#                       and qfar (R24), counted from selected-fits.csv and
#                       checked against select.R's single-stratum totals
#   summary-means.csv   across-community means per arm, stratum, scope and group
#   b0-means.csv        across-community means of the B0 outcomes per arm and stratum
#   beta-theta.csv      phase B only: collection intercept (beta_theta) bias,
#                       absolute bias and coverage, and the B0 and collection
#                       intercept posterior correlation, per community
#   beta-theta-means.csv  phase B only: their across-community means per arm and stratum
#   provenance.csv      scorer hashes, arguments, records read
# and, when the control and a new arm are both requested:
#   gate.csv            one row per SD: the numbers of every criterion and the
#                       result: PASS, FAIL, or UNDETERMINED when no criterion
#                       fails but one is evaluated over no community (NA)
#   gate-detail.csv     the same, one row per criterion and component, with
#                       tolerance_decisive: whether GATE_TOL decided an "at most" comparison
#   schedule-matched.csv  criteria 1 to 3 with each new arm's first fits (R12), descriptive
# When the gates of both A1 and A2 exist, or with --phase=A, OUT/A/gate.csv
# holds the A1 and A2 rows and one phase A row per SD (passes A1 and A2).
# Phase B's gate is separate: OUT/B/gate.csv, one phase B row per SD. With
# --phase=final, OUT/final/decision.csv combines OUT/A/gate.csv and
# OUT/B/gate.csv: per SD the phase A and phase B results and the protocol's
# consequence (README "Gate"). No verdict is printed; it is read from gate.csv
# or decision.csv after verify.R has passed for the same phase.
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
  # R24: criterion 1 pooled over the two contamination levels of each replicate,
  # the no-harm criteria at each level; the low band only (R22).
  B=list(strata=c('qnear','qfar'),c1_stratum='pooled',gated='low'),
  stop('No gate is defined here for phase ',phase))

ceil_two_thirds <- function(n) as.integer((2L*n+2L)%/%3L)

detail_row <- function(phase,sd,criterion,component,stratum,communities=NA,needed=NA,improved=NA,value_new=NA,
  value_control=NA,statistic=NA,threshold=NA,missing_new=NA,pass,tolerance_decisive=NA)
  data.frame(phase=phase,sd=sd,criterion=criterion,component=component,stratum=stratum,communities=as.integer(communities),
    needed=as.integer(needed),improved=as.integer(improved),value_new=as.numeric(value_new),value_control=as.numeric(value_control),
    statistic=as.numeric(statistic),threshold=as.numeric(threshold),missing_new=as.integer(missing_new),pass=as.logical(pass),
    tolerance_decisive=as.logical(tolerance_decisive),stringsAsFactors=FALSE)

# An "at most" comparison over n communities: NA when there is nothing to
# compare, and whether GATE_TOL decided it (the difference lies in (limit, limit + GATE_TOL]).
at_most <- function(n,difference,limit) {
  if(n==0L) return(list(pass=NA,decisive=NA))
  list(pass=difference<=limit+GATE_TOL,decisive=difference>limit && difference<=limit+GATE_TOL)
}

# One quantity of one group and stratum for the fits valid in both arms (the
# same communities for both arms; a group empty in a community is left out).
paired <- function(x,sd,group,stratum,column) {
  d <- x[x$group==group & x$stratum==stratum & x$cells>0,,drop=FALSE]
  a <- d[d$sd==sd,,drop=FALSE];b <- d[d$sd==1,,drop=FALSE];keys <- intersect(a$key,b$key)
  list(new=a[[column]][match(keys,a$key)],control=b[[column]][match(keys,b$key)],n=length(keys))
}

# Criterion 1 for one gated group. A community's value is the mean |E| of its
# fits in the pooled strata: in A1 its n100 and n300 fits (R21), in B its qnear
# and qfar fits (R24), in A2 its one fit. Improved means strictly smaller than the control's; at least
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
    pass=if(n==0L) NA else improved>=needed && mn<mc)
}

# Criterion 2 at one stratum: each band's mean |E| worse by at most 1 point and
# the mean overall MAE worse by at most 0.5 point. Over no community (a band
# with no cells anywhere, or an arm without valid fits, R27) it is NA, never PASS.
criterion2 <- function(x,phase,sd,stratum) {
  one <- function(group,column,component,limit) {
    p <- paired(x,sd,group,stratum,column)
    mn <- if(p$n) mean(p$new) else NA_real_;mc <- if(p$n) mean(p$control) else NA_real_;a <- at_most(p$n,mn-mc,limit)
    detail_row(phase,sd,'c2',component,stratum,p$n,value_new=mn,value_control=mc,statistic=mn-mc,threshold=limit,
      pass=a$pass,tolerance_decisive=a$decisive)
  }
  rbind(do.call(rbind,lapply(GATE_BANDS,function(g) one(g,'abs_signed_error',g,1))),one('all','mean_abs_cell_error','mae',.5))
}

# Criterion 3 at one stratum: each band's mean coverage lower by at most 0.03
# (NA over no community).
criterion3 <- function(x,phase,sd,stratum) do.call(rbind,lapply(GATE_BANDS,function(g) {
  p <- paired(x,sd,g,stratum,'coverage')
  mn <- if(p$n) mean(p$new) else NA_real_;mc <- if(p$n) mean(p$control) else NA_real_;a <- at_most(p$n,mc-mn,.03)
  detail_row(phase,sd,'c3',g,stratum,p$n,value_new=mn,value_control=mc,statistic=mc-mn,threshold=.03,
    pass=a$pass,tolerance_decisive=a$decisive)
}))

# Criterion 4 at one stratum: flagged selected fits of the new arm at most the
# control's (R20), and no community without a valid selected fit (R25, R27).
criterion4 <- function(conv,phase,sd,stratum,valid_new) {
  keys <- phase_keys(phase);keys <- keys[phase_stratum(phase,keys)==stratum]
  missing <- length(setdiff(keys,valid_new))
  fc <- conv$selected_flagged[conv$sd==1 & conv$stratum==stratum]
  if(length(fc)!=1L || is.na(fc)) stop('No control convergence count for phase ',phase,' stratum ',stratum)
  fn <- conv$selected_flagged[conv$sd==sd & conv$stratum==stratum];fn <- if(length(fn)==1L) fn else NA_real_
  detail_row(phase,sd,'c4','flags',stratum,value_new=fn,value_control=fc,statistic=fn-fc,threshold=0,missing_new=missing,
    pass=missing==0L && isTRUE(fn<=fc))
}

# FAIL if any criterion fails; otherwise UNDETERMINED if any is NA (evaluated
# over no community); otherwise PASS. Only PASS passes.
gate_result <- function(pass,descriptive=FALSE) {
  if(descriptive) return('DESCRIPTIVE')
  if(any(!pass,na.rm=TRUE)) 'FAIL' else if(anyNA(pass)) 'UNDETERMINED' else 'PASS'
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
  row$result <- gate_result(d$pass,descriptive)
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
  both <- intersect(gone,unique(x$key[x$sd==sd]))
  if(length(both)) stop('Community both scored and recorded as missing for SD ',sd,': ',both[1])
  valid <- unique(x$key[x$sd==sd])
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

read_record <- function(path,md5,label,hashes,role,sd,kind,version=SCORE_VERSION) {
  if(!file.exists(path)) stop('Not scored: ',path,'; run analysis.R (phase B: analysis-b.R) first')
  r <- readRDS(path)
  if(!identical(r$version,version)) stop('Score record ',path,' has version ',r$version %||% 'none')
  if(!identical(r$fit_md5,md5) || !identical(r$fit_label,label) || !identical(r$role,role) ||
     !identical(r$sd,as.numeric(sd)) || !identical(r$kind,kind))
    stop('Score record ',path,' is not of the selected fit ',label,' (md5 ',md5,')')
  if(!identical(r$hashes,hashes)) stop('Score record ',path,' was made by other scorer code; rescore with analysis.R (phase B: analysis-b.R)')
  r
}

# The record version and scorer hashes of a phase: analysis.R's for A1 and A2,
# analysis-b.R's for B.
phase_scorer <- function(phase,repo,archives) if(phase=='B') list(version=B_SCORE_VERSION,hashes=b_score_hashes(repo,archives)) else
  list(version=SCORE_VERSION,hashes=current_score_hashes(repo,archives))

# convergence_counts (flags.R) by the phase's own strata: in phase B the
# contamination levels (R24), which select.R's convergence.csv does not
# separate, since its size_stratum gives every phase B key the stratum 'all'.
stratum_counts <- function(selected,phase) {
  selected$stratum <- phase_stratum(phase,selected$key)
  groups <- split(selected,interaction(selected$phase,selected$stratum,selected$sd,drop=TRUE,lex.order=TRUE))
  x <- do.call(rbind,lapply(groups,function(g) data.frame(phase=g$phase[1],stratum=g$stratum[1],sd=g$sd[1],role=g$role[1],
    selected_fits=nrow(g),selected_flagged=sum(g$flagged),first_flagged=if(g$role[1]=='control') NA_integer_ else sum(g$first_flagged),
    long_repeats=if(g$role[1]=='control') NA_integer_ else sum(g$long_repeat),stringsAsFactors=FALSE)))
  x <- x[order(x$sd,x$stratum),,drop=FALSE];rownames(x) <- NULL
  x
}

# Per-level counts must add up to select.R's totals for every arm.
check_level_totals <- function(levels,totals,sel_dir) {
  cols <- c('selected_fits','selected_flagged','first_flagged','long_repeats')
  for(sd in unique(c(levels$sd,totals$sd))) {
    a <- levels[levels$sd==sd,,drop=FALSE];b <- totals[totals$sd==sd,,drop=FALSE]
    same <- nrow(b)==1L && all(vapply(cols,function(k) {x <- a[[k]];y <- b[[k]]
      if(all(is.na(x)) && is.na(y)) TRUE else !anyNA(x) && !is.na(y) && sum(x)==y},logical(1)))
    if(!same) stop('convergence.csv in ',sel_dir,' does not hold the totals of the per-level counts for SD ',sd)
  }
  invisible(TRUE)
}

same_table <- function(a,b,cols) all(vapply(cols,function(k)
  if(is.numeric(a[[k]]) || is.numeric(b[[k]]) || is.logical(a[[k]])) identical(as.numeric(a[[k]]),as.numeric(b[[k]])) else
    identical(as.character(a[[k]]),as.character(b[[k]])),logical(1)))

convergence_table <- function(phase,arms,sel_dir,repo,missing,supplement=NULL) {
  keys <- phase_keys(phase);strata <- unique(phase_stratum(phase,keys))
  cols <- c('phase','stratum','sd','role','selected_fits','selected_flagged','first_flagged','long_repeats')
  f <- file.path(sel_dir,'convergence.csv')
  if(file.exists(f)) {
    cv <- utils::read.csv(f,stringsAsFactors=FALSE)
    selected <- read_selected(file.path(sel_dir,'selected-fits.csv'))
    again <- convergence_counts(selected)
    if(nrow(cv)!=nrow(again) || !same_table(cv,again,cols)) stop('convergence.csv disagrees with selected-fits.csv in ',sel_dir)
    source <- 'selection'
    if(phase=='B') {
      # R24: criterion 4 at each contamination level, counted from selected-fits.csv.
      levels <- stratum_counts(selected,phase);check_level_totals(levels,cv,sel_dir)
      cv <- levels;source <- 'selection (per contamination level, R24)'
    }
    cf <- utils::read.csv(file.path(sel_dir,'control-flags.csv'),stringsAsFactors=FALSE)
    for(h in strata) if(!identical(as.numeric(sum(cf$flagged[phase_stratum(phase,cf$key)==h])),
      as.numeric(cv$selected_flagged[cv$sd==1 & cv$stratum==h]))) stop('control-flags.csv disagrees with convergence.csv in ',sel_dir)
  } else {
    if(any(arms!=1)) stop('No selection recorded in ',sel_dir,': only the control arm can be summarised')
    fc <- utils::read.csv(file.path(repo,ANALYSIS_REL,'results/flag-crosscheck.csv'),stringsAsFactors=FALSE)
    fc <- fc[fc$phase==phase & fc$selected,,drop=FALSE]
    if(!setequal(fc$key,keys) || anyDuplicated(fc$key) || !all(fc$match)) stop('results/flag-crosscheck.csv does not cover the phase ',phase,' controls')
    cv <- do.call(rbind,lapply(strata,function(h) {z <- fc[phase_stratum(phase,fc$key)==h,]
      data.frame(phase=phase,stratum=h,sd=1,role='control',selected_fits=nrow(z),selected_flagged=sum(z$flagged),
        first_flagged=NA_integer_,long_repeats=NA_integer_,stringsAsFactors=FALSE)}))
    source <- 'results/flag-crosscheck.csv (flags.R on the control fits; control selection not yet recorded for this phase)'
  }
  cv <- cv[cv$sd %in% arms,cols,drop=FALSE];cv$source <- if(nrow(cv)) source else character()
  # An SD that select.R could not record (R27): its supplement's selected fits, or none.
  if(!is.null(supplement) && nrow(supplement)) {
    sup <- (if(phase=='B') stratum_counts(supplement,phase) else convergence_counts(supplement))[cols]
    sup$source <- 'manual supplement (R27)';cv <- rbind(cv,sup)
  }
  for(sd in setdiff(arms,unique(cv$sd))) cv <- rbind(cv,data.frame(phase=phase,stratum=strata,sd=sd,role='new',selected_fits=0L,
    selected_flagged=NA_integer_,first_flagged=NA_integer_,long_repeats=NA_integer_,source='absent (R27)',stringsAsFactors=FALSE))
  for(sd in setdiff(arms,1)) for(h in setdiff(strata,cv$stratum[cv$sd==sd])) cv <- rbind(cv,data.frame(phase=phase,stratum=h,sd=sd,
    role='new',selected_fits=0L,selected_flagged=NA_integer_,first_flagged=NA_integer_,long_repeats=NA_integer_,source='absent (R27)',
    stringsAsFactors=FALSE))
  cv$expected_fits <- vapply(cv$stratum,function(h) sum(phase_stratum(phase,keys)==h),integer(1))
  cv$manual_missing <- vapply(seq_len(nrow(cv)),function(i) sum(missing$sd==cv$sd[i] & phase_stratum(phase,missing$key)==cv$stratum[i]),integer(1))
  cv$missing_fits <- cv$expected_fits-cv$selected_fits
  cv <- cv[order(cv$sd,cv$stratum),c(cols,'expected_fits','manual_missing','missing_fits','source'),drop=FALSE];rownames(cv) <- NULL
  cv
}

# Items of one new arm: select.R's selected-fits.csv, or, for an SD that
# select.R could not record (R27), supplement.R's manual-selection.csv, which
# must list exactly the fits of the controller's manual-selected.csv.
arm_items <- function(phase,sd,study,sel_dir,out) {
  sf <- file.path(sel_dir,'selected-fits.csv');mf <- file.path(out,phase,'manual-selection.csv')
  in_sel <- file.exists(sf) && {z <- read_selected(sf);any(z$role=='new' & z$sd==sd)}
  ms <- if(file.exists(mf)) read_selected(mf) else NULL;in_sup <- !is.null(ms) && any(ms$sd==sd)
  if(in_sel && in_sup) stop('SD ',sd,' is both in selected-fits.csv and in the manual supplement ',mf)
  if(in_sel) return(list(items=new_items(phase,sd,study,sel_dir),table=NULL))
  empty <- new_items_empty()
  if(!in_sup) return(list(items=empty,table=NULL))
  x <- ms[ms$sd==sd,,drop=FALSE];m <- read_manual_selected(sel_dir,phase);m <- m[m$sd==sd,,drop=FALSE]
  if(!setequal(paste(x$key,x$fit,x$fit_md5),paste(m$key,m$fit,m$fit_md5)) || nrow(x)!=nrow(m))
    stop('manual-selection.csv does not list the fits of manual-selected.csv for SD ',sd,'; rerun supplement.R')
  item <- function(key,schedule,fit,md5,kind) data.frame(role='new',phase=phase,sd=sd,key=key,schedule=schedule,
    fit=file.path(study,fit),fit_label=fit,expected_md5=md5,kind=kind,stringsAsFactors=FALSE)
  rows <- c(lapply(seq_len(nrow(x)),function(i) item(x$key[i],x$selected_schedule[i],x$fit[i],x$fit_md5[i],'selected')),
    lapply(which(x$long_repeat %in% TRUE),function(i) item(x$key[i],'initial',as.character(x$first_fit[i]),as.character(x$first_fit_md5[i]),'first')))
  list(items=do.call(rbind,rows),table=x)
}
new_items_empty <- function() data.frame(role=character(),phase=character(),sd=numeric(),key=character(),schedule=character(),
  fit=character(),fit_label=character(),expected_md5=character(),kind=character(),stringsAsFactors=FALSE)

# Records of the requested arms, each checked against the fit it must be of,
# and the supplement tables of the SDs recorded by hand (R27).
collect_records <- function(phase,arms,study,archives,repo,out,hashes,sel_dir,missing,version=SCORE_VERSION) {
  recs <- list();sup <- list()
  if(1 %in% arms) {
    ctl <- control_items(phase,study,archives,repo)
    for(i in seq_len(nrow(ctl))) recs[[length(recs)+1L]] <- read_record(score_path(out,phase,1,ctl$key[i],ctl$schedule[i]),
      ctl$expected_md5[i],ctl$fit_label[i],hashes,'control',1,'selected',version)
  }
  for(sd in setdiff(arms,1)) {
    a <- arm_items(phase,sd,study,sel_dir,out);x <- a$items;sup[[length(sup)+1L]] <- a$table
    for(i in seq_len(nrow(x))) recs[[length(recs)+1L]] <- read_record(score_path(out,phase,sd,x$key[i],x$schedule[i]),
      x$expected_md5[i],x$fit_label[i],hashes,'new',sd,x$kind[i],version)
    gone <- setdiff(phase_keys(phase),x$key[x$kind=='selected'])
    if(length(gone)) cat('R27: SD',sd,'has no valid selected fit for',length(gone),'communities:',paste(gone,collapse=', '),'\n')
  }
  sup <- Filter(Negate(is.null),sup)
  list(records=recs,supplement=if(length(sup)) do.call(rbind,sup) else NULL)
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

# B0 bias and coverage averaged over species within a community (b0.csv), then
# over communities (README "Secondary outcomes").
b0_means <- function(b) {
  by <- interaction(b$sd,b$stratum,drop=TRUE,lex.order=TRUE)
  x <- do.call(rbind,lapply(split(b,by),function(z) data.frame(phase=z$phase[1],sd=z$sd[1],arm=z$arm[1],stratum=z$stratum[1],
    communities=nrow(z),mean_b0_bias=mean(z$b0_bias),mean_b0_abs_bias=mean(z$b0_abs_bias),mean_b0_coverage=mean(z$b0_coverage),
    stringsAsFactors=FALSE)))
  rownames(x) <- NULL;x
}

# The collection intercept outcomes averaged over communities (phase B).
beta_theta_means <- function(t) {
  by <- interaction(t$sd,t$stratum,drop=TRUE,lex.order=TRUE)
  x <- do.call(rbind,lapply(split(t,by),function(z) data.frame(phase=z$phase[1],sd=z$sd[1],arm=z$arm[1],stratum=z$stratum[1],
    communities=nrow(z),mean_bt_bias=mean(z$bt_bias),mean_bt_abs_bias=mean(z$bt_abs_bias),mean_bt_coverage=mean(z$bt_coverage),
    mean_b0_bt_correlation=mean(z$b0_bt_correlation),stringsAsFactors=FALSE)))
  rownames(x) <- NULL;x
}

# The final decision per SD (README "Gate"): an SD is recommended only if it
# passes phase A and phase B, and then only the smallest such SD; one that
# passes phase A but fails phase B is a binary-only improvement and does not
# change the default; only PASS passes, so an UNDETERMINED phase B does not
# change the default either; phase B runs only for SDs that pass phase A. Any
# default change is made in a separate reviewed pull request.
# a: the phase A rows (sd, result); b: the phase B rows (sd, result).
final_decision <- function(a,b) {
  if(anyDuplicated(a$sd) || anyDuplicated(b$sd)) stop('One phase A and at most one phase B result per SD')
  sds <- sort(a$sd);ra <- a$result[match(sds,a$sd)];rb <- b$result[match(sds,b$sd)]
  if(length(setdiff(b$sd,a$sd))) stop('A phase B result for an SD without a phase A result')
  if(any(ra!='PASS' & !is.na(rb))) stop('A phase B result exists for an SD that did not pass phase A; the protocol runs phase B only for SDs that pass it')
  if(any(ra=='PASS' & is.na(rb))) stop('No phase B result for SD ',paste(sds[ra=='PASS' & is.na(rb)],collapse=', '),', which passed phase A')
  rb[is.na(rb)] <- 'NOT RUN'
  both <- ra=='PASS' & rb=='PASS'
  consequence <- ifelse(ra!='PASS','fails phase A: not adopted (phase B not run)',
    ifelse(rb=='PASS','passes phases A and B',
    ifelse(rb=='FAIL','binary-only improvement (passes phase A, fails phase B): does not change the default',
      'phase B undetermined (passes phase A; phase B neither passes nor fails): not adopted, does not change the default')))
  decision <- if(any(both)) sprintf('sigma_b0 = %s, the smallest SD passing phases A and B, is recommended as the new default; the change is made in a separate reviewed pull request',
    format(min(sds[both]))) else 'no SD passes phases A and B: the default stays at sigma_b0 = 1'
  data.frame(sd=sds,phase_a_result=ra,phase_b_result=rb,passes_both=both,recommended=both & sds==(if(any(both)) min(sds[both]) else -Inf),
    consequence=consequence,study_decision=decision,stringsAsFactors=FALSE)
}

write_table <- function(x,dir,name) {
  path <- file.path(dir,name);tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE);stopifnot(file.rename(tmp,path));path
}

summarise_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','archives','inputs-root','phase','arms','out'),required=c('repo','study','phase'))
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  out <- o$out %||% file.path(study,'summary')
  if(o$phase=='final') {
    need <- function(p) {f <- file.path(out,p,'gate.csv');if(!file.exists(f)) stop('Missing ',f);utils::read.csv(f,stringsAsFactors=FALSE)}
    a <- need('A');a <- a[a$subphase=='A',,drop=FALSE];b <- need('B');b <- b[b$subphase=='B',,drop=FALSE]
    x <- final_decision(a[c('sd','result')],b[c('sd','result')]);dir.create(file.path(out,'final'),recursive=TRUE,showWarnings=FALSE)
    cat('Wrote',write_table(x,file.path(out,'final'),'decision.csv'),'\n')
    return(0L)
  }
  if(o$phase=='A') {
    g <- lapply(SCORED_PHASES,function(p) {f <- file.path(out,p,'gate.csv');if(!file.exists(f)) stop('Missing ',f)
      utils::read.csv(f,stringsAsFactors=FALSE)})
    x <- combine_phase_a(g[[1]],g[[2]]);dir.create(file.path(out,'A'),recursive=TRUE,showWarnings=FALSE)
    cat('Wrote',write_table(x,file.path(out,'A'),'gate.csv'),'\n')
    return(0L)
  }
  if(!o$phase %in% c(SCORED_PHASES,'B')) stop('--phase must be A1, A2, A, B or final')
  if(is.null(o$arms)) stop('Missing --arms')
  phase <- o$phase;arms <- parse_arms(o$arms);new_sds <- setdiff(arms,1)
  sel_dir <- selection_dir(study,phase)
  missing <- if(dir.exists(sel_dir)) read_manual_missing(sel_dir,phase) else read_manual_missing(tempfile(),phase)
  scorer <- phase_scorer(phase,repo,archives);hashes <- scorer$hashes
  collected <- collect_records(phase,arms,study,archives,repo,out,hashes,sel_dir,missing,scorer$version);recs <- collected$records
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
  b0 <- record_table(recs,'b0');write_table(b0,d,'b0.csv')
  conv <- convergence_table(phase,arms,sel_dir,repo,missing,collected$supplement)
  write_table(conv,d,'convergence.csv')
  write_table(summary_means(groups),d,'summary-means.csv')
  write_table(b0_means(b0[b0$kind=='selected',,drop=FALSE]),d,'b0-means.csv')
  if(phase=='B') {
    bt <- record_table(recs,'theta');write_table(bt,d,'beta-theta.csv')
    write_table(beta_theta_means(bt[bt$kind=='selected',,drop=FALSE]),d,'beta-theta-means.csv')
  }
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
    # No verdict is printed: it is read from gate.csv once verify.R has passed (SCORING.md).
    cat('Wrote',file.path(d,c('gate.csv','gate-detail.csv','schedule-matched.csv')),sep='\n  ');cat('\n')
    other <- file.path(out,setdiff(SCORED_PHASES,phase),'gate.csv')
    if(phase %in% SCORED_PHASES && file.exists(other)) {
      o2 <- utils::read.csv(other,stringsAsFactors=FALSE)
      if(setequal(o2$sd,gate$sd)) {
        both <- if(phase=='A1') combine_phase_a(gate,o2) else combine_phase_a(o2,gate)
        dir.create(file.path(out,'A'),recursive=TRUE,showWarnings=FALSE)
        cat('Wrote',write_table(both,file.path(out,'A'),'gate.csv'),'\n')
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
  for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R','analysis-b.R','supplement.R')) source(file.path(here,f))
  status <- tryCatch(summarise_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
