#!/usr/bin/env Rscript
# Phase A scorer of the occupancy-intercept prior study: per-community outcomes
# from the saved draws of each arm's selected fit (README "Outcomes and gate",
# AMENDMENT-1.md, AMENDMENT-2.md). The same code scores every arm, including the
# reused controls (R10).
#
#   Rscript analysis.R --repo=REPO --study=STUDY --inputs-root=DIR --phase=A1|A2
#     --arms=1[,2,3,5] [--workers=N] [--archives=DIR] [--out=DIR]
#     [--allow-before-both-selections=TRUE]
#
# Arm 1 is the control: its selected fit per community is the archived one of
# results/control-provenance.csv, checked against control_schedules(). A new
# arm's selected fits are those of STUDY/selection/<phase>/selected-fits.csv,
# with every fit's md5 checked; where a flagged first fit was replaced by its
# longer repeat, the first fit is scored too, for the schedule-matched
# sensitivity (R12). New arms are refused before the phase's selection exists
# and, as AMENDMENT-1.md requires, before both A1 and A2 selections exist,
# unless a dated amendment permits otherwise and the option says so.
#
# One record per scored fit is written to OUT/<phase>/scores/sd<sd>/<key>-<schedule>.rds
# (OUT defaults to STUDY/summary): the cell table (truth, posterior-mean
# estimate, 95% interval, containment), the group table and the B0 table. A
# record of the same fit md5 and scorer hashes is resumed, not recomputed.
# summarise.R builds the study tables and the gate from these records.
#
# Definitions are taken from the archived scorers after md5 checks
# (flags.R's checked_scorer_file), not copied:
#   A1  jsdm-sample-size-recheck/helpers.R: score_fit (probability per draw,
#       lines 77-78; posterior mean, 82; truth, 59-61; band errors, 62-92) and
#       error_metrics (45-47). The per-draw statements are parsed out of
#       score_fit and evaluated here to keep every draw for the intervals.
#   A2  spatial-targeted-recheck/score.R: reconstruct_spatial_draws (31-47),
#       spatial_groups (50-56).
#   Both score.R's score_draw_block (58-77): posterior mean, central 95%
#       interval from type-7 quantiles of pooled draws, inclusive containment
#       (62-66, 73); also used for B0 (131, 138-139).
# Only the cell sets (README "Bands and groups") and the aggregation into
# percentage points are written here. Not sourced or hashed by run.R.

ANALYSIS_REL <- 'dev/simstudy/occupancy-intercept-prior'
SCORE_VERSION <- 'intercept-prior-phase-a-v1'
METRIC_FILES <- c('dev/simstudy/jsdm-sample-size-recheck/helpers.R','dev/simstudy/spatial-targeted-recheck/score.R')
SCORED_PHASES <- c('A1','A2')
ARMS <- c(1,2,3,5)
A1_ORIGINAL_SITES <- 100L

# ---------------------------------------------------------------------------
# Archived definitions
# ---------------------------------------------------------------------------

current_score_hashes <- function(repo,archives) {
  scripts <- file.path(repo,ANALYSIS_REL,c('analysis.R','jobs.R','flags.R','verify-helpers.R'))
  archived <- vapply(METRIC_FILES,checked_scorer_file,character(1),repo=repo,archives=archives)
  h <- c(hash_files(scripts,repo),hash_files(archived,repo))
  h[order(names(h))]
}

load_metric_scorers <- function(repo,archives) {
  paths <- vapply(METRIC_FILES,checked_scorer_file,character(1),repo=repo,archives=archives)
  sc <- new.env(parent=globalenv())
  extract_definitions(paths[[1]],c('error_metrics','trace_stats','score_fit'),sc)
  extract_definitions(paths[[2]],c('trace_diagnostics','independent_bases','reconstruct_spatial_draws',
    'spatial_groups','score_draw_block'),sc)
  # The per-draw statements of score_fit (helpers.R:77-78), verbatim.
  eta <- nested_assignments(list(body(sc$score_fit)),as.name('eta'))
  p <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('plogis')),
    nested_assignments(list(body(sc$score_fit)),as.name('p')))
  if(length(eta)!=1L || length(p)!=1L) stop('Expected one per-draw eta and probability statement in score_fit')
  sc$a1_eta <- eta[[1]];sc$a1_probability <- p[[1]]
  sc$score_hashes <- current_score_hashes(repo,archives)
  sc
}

# ---------------------------------------------------------------------------
# Cells: truth, posterior-mean estimate and interval of every scored cell
# ---------------------------------------------------------------------------

make_cells <- function(phase,truth,estimate,covered,n_original=nrow(truth),target_prevalence=NULL,
  lower=NULL,upper=NULL,b0=NULL,checks=list()) {
  truth <- unname(as.matrix(truth))
  reshape <- function(x) {x <- unname(x);if(is.null(dim(x))) {stopifnot(length(x)==length(truth));matrix(x,nrow(truth))} else as.matrix(x)}
  estimate <- reshape(estimate);covered <- reshape(covered)
  stopifnot(phase %in% SCORED_PHASES,identical(dim(truth),dim(estimate)),identical(dim(truth),dim(covered)),
    is.logical(covered),!anyNA(covered),all(is.finite(truth)),all(is.finite(estimate)),n_original<=nrow(truth))
  if(phase=='A2') stopifnot(length(target_prevalence)==ncol(truth))
  shape <- function(x) if(is.null(x)) matrix(NA_real_,nrow(truth),ncol(truth)) else {stopifnot(length(x)==length(truth));matrix(x,nrow(truth))}
  list(phase=phase,n=nrow(truth),S=ncol(truth),n_original=as.integer(n_original),truth=truth,estimate=estimate,
    lower=shape(lower),upper=shape(upper),covered=covered,target_prevalence=target_prevalence,b0=b0,checks=checks)
}

# Posterior mean, central 95% interval and containment per cell (score.R:58-77).
interval_cells <- function(draws,truth,sc) {
  e <- sc$score_draw_block(draws,truth,'cell')$elements
  data.frame(truth=e$truth,estimate=e$estimate,lower=e$lower,upper=e$upper,covered=e$covered)
}

b0_cells <- function(B0_output,truth,sc) {
  e <- sc$score_draw_block(B0_output,truth,'intercept')$elements
  data.frame(species=seq_along(truth),truth=e$truth,estimate=e$estimate,bias=e$bias,lower=e$lower,upper=e$upper,covered=e$covered)
}

# Every retained draw of every cell, by the per-draw statements of score_fit.
# The evaluated expressions are only ever the two statements parsed out of the
# tracked archived scorer after its md5 matched the record made when it scored
# the controls (load_metric_scorers), never text from fits, inputs or other data.
a1_probability_draws <- function(fit,sc) {
  j <- fit$results_output$jsdm_output;X <- fit$X_psi
  S <- dim(j$B0_output)[1];ni <- dim(j$B0_output)[2];nc <- dim(j$B0_output)[3]
  env <- new.env(parent=globalenv());env$X <- X;env$j <- j
  out <- array(NA_real_,c(nrow(X)*S,ni,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    env$it <- it;env$ch <- ch
    env$eta <- eval(sc$a1_eta,env)
    out[,it,ch] <- as.vector(eval(sc$a1_probability,env))
  }
  out
}

# Phase A1: the archived scorer gives truth, estimate and band errors; the
# retained draws give the intervals. Both must agree exactly on the estimate.
a1_cells <- function(fit,input,sc) {
  r <- sc$score_fit(fit,input)
  n <- input$n;S <- ncol(input$data$OTU)
  stopifnot(identical(dim(r$truth),c(as.integer(n),S)),identical(dim(r$estimate),dim(r$truth)))
  iv <- interval_cells(a1_probability_draws(fit,sc),as.vector(r$truth),sc)
  difference <- max(abs(iv$estimate-as.vector(r$estimate)))
  if(!is.finite(difference) || difference>=1e-12) stop('Draw reconstruction disagrees with score_fit: ',difference)
  cells <- make_cells('A1',truth=r$truth,estimate=r$estimate,covered=iv$covered,n_original=A1_ORIGINAL_SITES,
    lower=iv$lower,upper=iv$upper,b0=b0_cells(fit$results_output$jsdm_output$B0_output,input$truth$B0,sc))
  # The group errors must be the archived scorer's own (helpers.R:89-92).
  g <- group_table(cells,sc);gd <- 0
  for(k in seq_len(nrow(r$groups))) {
    a <- r$groups[k,];if(a$cells==0) next
    scope <- if(a$scope=='original100') 'primary' else 'allsites'
    if(scope=='allsites' && n<=A1_ORIGINAL_SITES) next
    m <- g[g$scope==scope & g$group==a$group,,drop=FALSE]
    if(nrow(m)!=1L || m$cells!=a$cells) stop('Cell sets disagree with score_fit for ',a$scope,' ',a$group)
    gd <- max(gd,abs(m$signed_error-100*a$bias),abs(m$mean_abs_cell_error-100*a$mae))
  }
  if(gd>1e-10) stop('Group errors disagree with score_fit: ',gd)
  cells$checks <- list(estimate_vs_score_fit=difference,groups_vs_score_fit=gd)
  cells
}

# Phase A2: the archived spatial reconstruction, binary arm, with the checks
# of score_spatial_fit (score.R:83-86).
a2_cells <- function(fit,input,sc,knots=100L) {
  t <- input$truth;js <- fit$results_output$jsdm_output
  n <- input$settings$n;S <- input$settings$S
  stopifnot(fit$infos$ps==knots,fit$infos$n_factors==0L,identical(fit$infos$speciesNames,colnames(input$data$binary$OTU)),
    fit$infos$model=='binary',max(abs(fit$Xs-t$Xs))<1e-12,max(abs(fit$X_psi-t$X))<1e-12,identical(dim(t$psi),c(as.integer(n),as.integer(S))))
  iv <- interval_cells(sc$reconstruct_spatial_draws(fit)$probability,t$psi,sc)
  make_cells('A2',truth=t$psi,estimate=iv$estimate,covered=iv$covered,target_prevalence=t$target_prevalence,
    lower=iv$lower,upper=iv$upper,b0=b0_cells(js$B0_output,t$B0,sc))
}

community_cells <- function(phase,fit,input,sc) switch(phase,A1=a1_cells(fit,input,sc),A2=a2_cells(fit,input,sc),
  stop('No scorer for phase ',phase))

# ---------------------------------------------------------------------------
# Groups and community outcomes
# ---------------------------------------------------------------------------

# Cell sets (README "Bands and groups"): truth bands below 0.2, 0.2 to 0.8
# inclusive and above 0.8. A1: the original sites are primary (helpers.R:62);
# fits with more sites are also scored on all sites, descriptively; the rare
# group is every fitted site of each species whose mean truth over all fitted
# sites is below 0.2 (run_nonspatial_recheck_balanced.R:129-134). A2: all sites,
# with the prevalence groups of score.R's spatial_groups and the gated rare
# group of the 1% and 5% species together (R18).
cell_masks <- function(cells,sc) {
  n <- cells$n;S <- cells$S;tv <- as.vector(cells$truth)
  bands <- function(ix) list(all=ix,low=ix[tv[ix]<.2],middle=ix[tv[ix]>=.2 & tv[ix]<=.8],high=ix[tv[ix]>.8])
  add <- function(masks,scope,sets) c(masks,lapply(names(sets),function(g) list(scope=scope,group=g,idx=sets[[g]])))
  if(cells$phase=='A1') {
    sites <- function(s) unlist(lapply(s,function(k) (k-1L)*n+seq_len(cells$n_original)))
    masks <- add(list(),'primary',bands(sites(seq_len(S))))
    rare <- which(colMeans(cells$truth)<.2)
    masks <- add(masks,'primary',list(rare_below_20pct=unlist(lapply(rare,function(s) (s-1L)*n+seq_len(n)))))
    if(n>cells$n_original) masks <- add(masks,'allsites',bands(seq_len(n*S)))
    return(masks)
  }
  g <- sc$spatial_groups(cells$truth,cells$target_prevalence)
  stopifnot(setequal(names(g),c('all','low','medium','high',paste0('prevalence_',c(1,5,25,75),'pct'))))
  masks <- add(list(),'primary',bands(seq_len(n*S)))
  # The bands written above are exactly score.R's (medium there is middle here).
  names_of <- vapply(masks,function(m) m$group,'')
  for(k in c('low','medium','high')) stopifnot(identical(which(g[[k]]),masks[[match(if(k=='medium') 'middle' else k,names_of)]]$idx))
  prevalence <- paste0('prevalence_',c(1,5,25,75),'pct')
  add(masks,'primary',c(stats::setNames(lapply(prevalence,function(k) which(g[[k]])),prevalence),
    list(rare_1_5pct=which(g$prevalence_1pct | g$prevalence_5pct))))
}

# One row per non-empty cell set: signed error E and mean absolute cell error
# in percentage points (error_metrics, helpers.R:45-47), |E|, and coverage.
group_table <- function(cells,sc) {
  tv <- as.vector(cells$truth);ev <- as.vector(cells$estimate);cv <- as.vector(cells$covered)
  rows <- lapply(cell_masks(cells,sc),function(m) {
    if(!length(m$idx)) return(NULL)
    e <- sc$error_metrics(tv[m$idx],ev[m$idx])
    data.frame(scope=m$scope,group=m$group,cells=length(m$idx),truth_mean=e[['truth']],estimate_mean=e[['estimate']],
      signed_error=100*e[['bias']],abs_signed_error=abs(100*e[['bias']]),mean_abs_cell_error=100*e[['mae']],
      coverage=mean(cv[m$idx]),stringsAsFactors=FALSE)
  })
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# B0 on the fitted (logit) scale, averaged over species (README "Secondary outcomes").
b0_summary <- function(b0) {
  if(is.null(b0)) return(data.frame(species=0L,b0_bias=NA_real_,b0_abs_bias=NA_real_,b0_coverage=NA_real_))
  data.frame(species=nrow(b0),b0_bias=mean(b0$bias),b0_abs_bias=mean(abs(b0$bias)),b0_coverage=mean(b0$covered))
}

community_of <- function(phase,key) if(phase=='A1') sub('^jsdm-n0[13]00-','',key) else key
stratum_of <- function(phase,key) if(phase=='A1') size_stratum(key) else rep('all',length(key))

score_path <- function(root,phase,sd,key,schedule)
  file.path(root,phase,'scores',paste0('sd',format(sd)),paste0(key,'-',schedule,'.rds'))

score_record <- function(item,cells,sc,hashes,input_md5,fit_md5=item$expected_md5,warnings=NA_integer_) {
  n <- cells$n;S <- cells$S
  list(version=SCORE_VERSION,phase=item$phase,sd=as.numeric(item$sd),role=item$role,kind=item$kind %||% 'selected',
    key=item$key,schedule=item$schedule,stratum=stratum_of(item$phase,item$key),community=community_of(item$phase,item$key),
    fit_label=item$fit_label,fit_md5=fit_md5,input_md5=input_md5,warnings=warnings,hashes=hashes,
    groups=group_table(cells,sc),b0=b0_summary(cells$b0),b0_species=cells$b0,
    cells=data.frame(site=rep(seq_len(n),S),species=rep(seq_len(S),each=n),truth=as.vector(cells$truth),
      estimate=as.vector(cells$estimate),lower=as.vector(cells$lower),upper=as.vector(cells$upper),covered=as.vector(cells$covered)),
    checks=cells$checks,scored_at=format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'))
}

# ---------------------------------------------------------------------------
# Which fits: selected controls and new-arm selections
# ---------------------------------------------------------------------------

control_items <- function(phase,study,archives,repo) {
  x <- selection_items(phase,numeric(0),study,archives,repo)
  x <- x[x$role=='control',,drop=FALSE];x$kind <- 'selected';rownames(x) <- NULL
  if(!identical(sort(x$key),sort(phase_keys(phase)))) stop('Control items do not cover phase ',phase)
  x
}

read_selected <- function(path) {
  if(!file.exists(path)) stop('Selection table missing: ',path)
  utils::read.csv(path,stringsAsFactors=FALSE,colClasses=c(phase='character',sd='numeric',key='character',role='character',
    first_schedule='character',selected_schedule='character',fit='character',fit_md5='character',rule='character',reasons='character'))
}

# The manual record of communities without a valid selected fit (AMENDMENT-2, R27).
read_manual_missing <- function(sel_dir,phase) {
  path <- file.path(sel_dir,'manual-missing.csv')
  if(!file.exists(path)) return(data.frame(phase=character(),sd=numeric(),key=character(),reason=character(),stringsAsFactors=FALSE))
  x <- utils::read.csv(path,stringsAsFactors=FALSE,colClasses=c(phase='character',sd='numeric',key='character',reason='character'))
  if(!identical(names(x),c('phase','sd','key','reason')) || any(x$phase!=phase) || !all(x$key %in% phase_keys(phase)) ||
     !all(x$sd %in% setdiff(ARMS,1)) || any(!nzchar(x$reason)) || anyDuplicated(x[c('sd','key')]))
    stop('Malformed manual record ',path,': columns phase, sd, key, reason; new-arm SDs and phase ',phase,' keys only')
  x
}

new_items <- function(phase,sd,study,sel_dir) {
  sf <- read_selected(file.path(sel_dir,'selected-fits.csv'))
  missing <- read_manual_missing(sel_dir,phase)
  x <- sf[sf$role=='new' & sf$phase==phase & sf$sd==sd,,drop=FALSE]
  if(anyDuplicated(x$key) || !all(x$key %in% phase_keys(phase))) stop('Malformed selection for SD ',sd,' in ',sel_dir)
  both <- intersect(x$key,missing$key[missing$sd==sd])
  if(length(both)) stop('Community both selected and recorded as missing for SD ',sd,': ',both[1])
  item <- function(r,schedule,fit,md5,kind) data.frame(role='new',phase=phase,sd=sd,key=r$key,schedule=schedule,
    fit=file.path(study,fit),fit_label=fit,expected_md5=md5,kind=kind,stringsAsFactors=FALSE)
  rows <- lapply(seq_len(nrow(x)),function(i) item(x[i,],x$selected_schedule[i],x$fit[i],x$fit_md5[i],'selected'))
  repeated <- x[x$long_repeat %in% TRUE,,drop=FALSE]
  if(nrow(repeated)) {
    first <- utils::read.csv(file.path(sel_dir,'long-selection.csv'),stringsAsFactors=FALSE,
      colClasses=c(sd='numeric',key='character',schedule='character',fit='character',fit_md5='character'))
    for(i in seq_len(nrow(repeated))) {
      f <- first[first$sd==sd & first$key==repeated$key[i] & first$role=='new',,drop=FALSE]
      if(nrow(f)!=1L || !isTRUE(f$flagged) || f$schedule!=repeated$first_schedule[i]) stop('No recorded flagged first fit for ',repeated$key[i])
      rows[[length(rows)+1L]] <- item(repeated[i,],f$schedule,f$fit,f$fit_md5,'first')
    }
  }
  if(!length(rows)) return(data.frame(role=character(),phase=character(),sd=numeric(),key=character(),schedule=character(),
    fit=character(),fit_label=character(),expected_md5=character(),kind=character(),stringsAsFactors=FALSE))
  out <- do.call(rbind,rows);rownames(out) <- NULL;out
}

# Score one fit after checking it is the fit the item names (flags.R's
# evaluate_item: md5, job record, schedule, prior, input md5).
score_item <- function(item,sc,archives,inputs_root,out_root) {
  path <- score_path(out_root,item$phase,item$sd,item$key,item$schedule)
  if(file.exists(path)) {
    old <- readRDS(path);md5 <- unname(tools::md5sum(item$fit))
    if(identical(old$fit_md5,md5) && (is.na(item$expected_md5) || identical(md5,item$expected_md5)) &&
       identical(old$hashes,sc$score_hashes) && identical(old$role,item$role) && identical(old$sd,as.numeric(item$sd)) &&
       identical(old$kind,item$kind %||% 'selected')) {
      cat(format(Sys.time()),item$phase,paste0('sd',item$sd),item$key,item$schedule,'already scored; resumed\n');flush.console()
      return(path)
    }
  }
  captured <- NULL
  scorer <- function(phase,fit,warnings,input,sc) {captured <<- community_cells(phase,fit,input,sc);data.frame(n_warnings=length(warnings))}
  row <- evaluate_item(item,sc,archives,inputs_root,flag_fn=scorer)
  spec <- phase_jobs(item$phase,archives,inputs_root,item$key)[[1]]
  rec <- score_record(item,captured,sc,sc$score_hashes,input_md5=spec$input_md5,fit_md5=row$fit_md5,warnings=row$n_warnings)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);atomic_save(rec,path)
  cat(format(Sys.time()),item$phase,paste0('sd',item$sd),item$key,item$schedule,'scored\n');flush.console()
  path
}

parse_arms <- function(text) {
  a <- strsplit(text %||% '',',',fixed=TRUE)[[1]]
  if(!length(a) || any(!a %in% as.character(ARMS)) || anyDuplicated(a))
    stop('--arms must list arms among ',paste(ARMS,collapse=', '),' (1 is the control), without repeats')
  sort(as.numeric(a))
}

analysis_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','phase','arms','workers','out','allow-before-both-selections'),
    required=c('repo','study','inputs-root','phase','arms'))
  if(!o$phase %in% SCORED_PHASES) stop('--phase must be A1 or A2; phase B is not scored by this script')
  arms <- parse_arms(o$arms);workers <- parse_workers(o$workers %||% '1')
  allow <- o$`allow-before-both-selections` %||% 'FALSE'
  if(!allow %in% c('TRUE','FALSE')) stop('--allow-before-both-selections must be TRUE or FALSE')
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  out <- o$out %||% file.path(study,'summary');phase <- o$phase
  sel_dir <- selection_dir(study,phase);new_sds <- setdiff(arms,1)
  if(length(new_sds)) {
    f <- file.path(sel_dir,'selected-fits.csv')
    if(!file.exists(f)) stop('No new-arm occupancy error is computed before select.R --mode=final has written ',f)
    other <- file.path(selection_dir(study,setdiff(SCORED_PHASES,phase)),'selected-fits.csv')
    if(!file.exists(other) && allow!='TRUE') stop('AMENDMENT-1.md: no new-arm occupancy error is computed before both the A1 and ',
      'the A2 selections exist, and ',other,' is missing. Use --allow-before-both-selections=TRUE only under a dated amendment that permits it.')
  }
  sc <- load_metric_scorers(repo,archives)
  items <- list()
  if(1 %in% arms) {
    ctl <- control_items(phase,study,archives,repo)
    cf <- file.path(sel_dir,'control-flags.csv')
    if(file.exists(cf)) {
      recorded <- utils::read.csv(cf,stringsAsFactors=FALSE,colClasses=c(key='character',fit='character',fit_md5='character'))
      i <- match(ctl$key,recorded$key)
      if(anyNA(i) || !identical(recorded$fit[i],ctl$fit_label) || !identical(recorded$fit_md5[i],ctl$expected_md5))
        stop('The recorded control selection ',cf,' disagrees with results/control-provenance.csv')
    }
    items[[1]] <- ctl
  }
  for(sd in new_sds) items[[length(items)+1L]] <- new_items(phase,sd,study,sel_dir)
  items <- do.call(rbind,items)
  cat('Scoring',nrow(items),'fits of phase',phase,'(arms',paste(arms,collapse=', '),') into',file.path(out,phase),'\n');flush.console()
  work <- function(i) tryCatch(score_item(items[i,,drop=FALSE],sc,archives,inputs_root,out),
    error=function(e) structure(conditionMessage(e),class='score-error'))
  res <- if(workers==1L) lapply(seq_len(nrow(items)),work) else
    parallel::mclapply(seq_len(nrow(items)),work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  bad <- !vapply(res,function(x) is.character(x) && !inherits(x,'score-error') && length(x)==1L,logical(1))
  prov <- data.frame(kind=c(rep('hash',length(sc$score_hashes)),'args','finished','fits','failed'),
    name=c(names(sc$score_hashes),'args','time','scored','scored'),
    value=c(unname(sc$score_hashes),paste(args,collapse=' '),format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),nrow(items),sum(bad)))
  dir.create(file.path(out,phase,'scores'),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(prov,file.path(out,phase,'scores',paste0('provenance-analysis-',format(Sys.time(),'%Y%m%d%H%M%S'),'.csv')),row.names=FALSE)
  if(any(bad)) {
    cat('FAILED',sum(bad),'of',nrow(items),'fits:\n')
    for(i in which(bad)) cat(' ',items$key[i],paste0('sd',items$sd[i]),items$schedule[i],':',paste(as.character(res[[i]]),collapse=' '),'\n')
    return(1L)
  }
  cat('Scored',nrow(items),'fits.\n')
  0L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  repo_arg <- sub('^--repo=','',grep('^--repo=',commandArgs(trailingOnly=TRUE),value=TRUE))
  # The scripts sourced here and the archived scorers read from --repo must be one checkout.
  if(length(repo_arg)!=1L || !identical(normalizePath(file.path(repo_arg,ANALYSIS_REL),mustWork=FALSE),here)) {
    cat('ERROR: run the analysis.R of the --repo checkout (',here,' is not under --repo)\n',sep='');quit(save='no',status=1)
  }
  for(f in c('jobs.R','verify-helpers.R','flags.R')) source(file.path(here,f))
  status <- tryCatch(analysis_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
