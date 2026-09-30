#!/usr/bin/env Rscript
# Phase B scorer of the occupancy-intercept prior study: per-community outcomes
# of the two-stage fits from the saved draws of each arm's selected fit (README
# "Outcomes and gate", "Secondary outcomes"; AMENDMENT-2.md R24 and R27). The
# same code scores every arm, including the reused pr11 controls (R10).
#
#   Rscript analysis-b.R --repo=REPO --study=STUDY --inputs-root=DIR --phase=B
#     --arms=1[,2,3,5] [--workers=N] [--archives=DIR] [--out=DIR]
#
# This file extends analysis.R without editing it: analysis.R's md5 is one of
# the scorer hashes of every phase A score record, so any edit to it would
# invalidate all phase A records. The phase A scorer is therefore untouched;
# this script sources analysis.R for the phase-independent pieces (the
# selection readers, the control items, the interval and B0 blocks and the
# record path) and adds the phase B definitions below. Phase B records carry
# their own version and hash set (b_score_hashes), which includes both files.
#
# Arm 1 is the control: its selected fit per community is the archived pr11 one
# of results/control-provenance.csv. A new arm's selected fits are those of
# STUDY/selection/B/selected-fits.csv (select.R --phases=B --mode=final), and,
# where a flagged first fit was replaced by its single longer repeat, the first
# fit too, for the schedule-matched sensitivity (R12). New arms are refused
# before that file exists (README: no occupancy error before select-final).
#
# One record per scored fit is written to OUT/B/scores/sd<sd>/<key>-<schedule>.rds
# (OUT defaults to STUDY/summary), resumed when the fit md5 and scorer hashes
# are unchanged.
#
# Definitions are taken from the archived scorers after md5 checks
# (flags.R's checked_scorer_file), not copied:
#   current-main-recheck/helpers.R: load_scoring and score_current_fit, the
#     exact scorer the pr11 archive applied to these controls. load_scoring
#     builds the design scorer from nonspatial-bias-recheck/
#     run_nonspatial_recheck_balanced.R (score_fit, lines 59-149) and
#     nonspatial-design-recheck/design_helpers.R (make_design_scorer, 117-137,
#     which adds the original-site groups), with the all-chain trace summary.
#     It supplies the truth (lines 71-73), the posterior-mean estimate (118),
#     the band and rare-group errors (122-134, and the original-site groups),
#     the B0 truth and the collection intercept truth (78-85), and checks the
#     design, the truth and the package's posterior means on the way.
#   run_nonspatial_recheck_balanced.R: the per-draw statements of score_fit
#     (lines 114-115), parsed out and evaluated here to keep the draws of the
#     scored cells for their intervals.
#   spatial-targeted-recheck/score.R: score_draw_block (58-77), the central
#     95% interval from type-7 quantiles of pooled draws with inclusive
#     containment, used for cells, B0 and the collection intercept, as in phase A.
#   jsdm-sample-size-recheck/helpers.R: error_metrics (45-47).
# Only the scored cell sets (README "Bands and groups", phase B) and the
# aggregation into percentage points are written here. Not sourced or hashed
# by run.R.

B_SCORE_VERSION <- 'intercept-prior-phase-b-v1'
B_ANALYSIS_FILES <- c('analysis.R','analysis-b.R','jobs.R','flags.R','verify-helpers.R')
B_METRIC_FILES <- c(helpers='dev/simstudy/current-main-recheck/helpers.R',
  balanced='dev/simstudy/nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R',
  design='dev/simstudy/nonspatial-design-recheck/design_helpers.R',
  jsdm='dev/simstudy/jsdm-sample-size-recheck/helpers.R',
  score='dev/simstudy/spatial-targeted-recheck/score.R')
B_KEY_PATTERN <- '^design-(qnear|qfar)_K6-sites300-([0-9]{2})$'

# ---------------------------------------------------------------------------
# Strata and communities (AMENDMENT-2.md, R24)
# ---------------------------------------------------------------------------

# qnear and qfar of one replicate share their generating occupancy: one
# generating community (the replicate) at two contamination levels (strata).
b_stratum <- function(key) {
  if(!all(grepl(B_KEY_PATTERN,key))) stop('Not a phase B key: ',key[!grepl(B_KEY_PATTERN,key)][1])
  sub(B_KEY_PATTERN,'\\1',key)
}
b_community <- function(key) {b_stratum(key);sub(B_KEY_PATTERN,'\\2',key)}
phase_stratum <- function(phase,key) if(phase=='B') b_stratum(key) else stratum_of(phase,key)
phase_community <- function(phase,key) if(phase=='B') b_community(key) else community_of(phase,key)

# ---------------------------------------------------------------------------
# Archived definitions
# ---------------------------------------------------------------------------

b_score_hashes <- function(repo,archives) {
  scripts <- file.path(repo,ANALYSIS_REL,B_ANALYSIS_FILES)
  archived <- vapply(B_METRIC_FILES,checked_scorer_file,character(1),repo=repo,archives=archives)
  h <- c(hash_files(scripts,repo),hash_files(archived,repo))
  h[order(names(h))]
}

load_b_scorers <- function(repo,archives) {
  paths <- vapply(B_METRIC_FILES,checked_scorer_file,character(1),repo=repo,archives=archives)
  sc <- new.env(parent=globalenv())
  extract_definitions(paths[['helpers']],c('extract_functions','load_scoring','score_current_fit'),sc)
  # load_scoring reads the three scorer files below --repo, whose md5s were just checked.
  sc$scoring <- sc$load_scoring(repo)
  extract_definitions(paths[['score']],c('trace_diagnostics','score_draw_block'),sc)
  extract_definitions(paths[['jsdm']],'error_metrics',sc)
  # The per-draw statements of score_fit (run_nonspatial_recheck_balanced.R:114-115), verbatim.
  body_of <- list(body(sc$scoring$design$score_fit))
  eta <- nested_assignments(body_of,as.name('eta'));p <- nested_assignments(body_of,quote(psi_draws[,it,ch]))
  if(length(eta)!=1L || length(p)!=1L) stop('Expected one per-draw eta and probability statement in the phase B score_fit')
  sc$b_eta <- eta[[1]];sc$b_probability <- p[[1]]
  sc$score_hashes <- b_score_hashes(repo,archives)
  sc
}

# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------

# Scored cells (README, phase B): the original sites of every species, and
# every fitted site of each rare species (mean truth over all fitted sites
# below 0.2), in the column-major order of the sites-by-species matrix.
b_scored_cells <- function(truth,n_original) {
  n <- nrow(truth);S <- ncol(truth)
  original <- unlist(lapply(seq_len(S),function(s) (s-1L)*n+seq_len(n_original)))
  rare <- unlist(lapply(which(colMeans(truth)<.2),function(s) (s-1L)*n+seq_len(n)))
  sort(unique(as.integer(c(original,rare))))
}

# truth and estimate for every fitted cell; lower, upper and covered for the
# scored cells (NA elsewhere). b0 and theta: per-species tables.
b_make_cells <- function(truth,estimate,covered,n_original,lower=NULL,upper=NULL,b0=NULL,theta=NULL,checks=list()) {
  truth <- unname(as.matrix(truth));n <- nrow(truth);S <- ncol(truth)
  shape <- function(x,na) if(is.null(x)) matrix(na,n,S) else {x <- unname(x);stopifnot(length(x)==n*S);matrix(x,n,S)}
  estimate <- shape(estimate,NA_real_);covered <- shape(covered,NA)
  stopifnot(is.logical(covered),all(is.finite(truth)),all(is.finite(estimate)),n_original<=n)
  scored <- b_scored_cells(truth,n_original)
  if(anyNA(covered[scored])) stop('A scored cell has no interval')
  list(phase='B',n=n,S=S,n_original=as.integer(n_original),truth=truth,estimate=estimate,lower=shape(lower,NA_real_),
    upper=shape(upper,NA_real_),covered=covered,scored=scored,b0=b0,theta=theta,checks=checks)
}

# The draws of the given cells, by the per-draw statements of score_fit. The
# evaluated expressions are only ever the two statements parsed out of the
# tracked archived scorer after its md5 matched its record (load_b_scorers),
# never text from fits, inputs or other data.
b_probability_draws <- function(fit,sc,cells) {
  jo <- fit$results_output$jsdm_output;ni <- dim(jo$B0_output)[2];nc <- dim(jo$B0_output)[3]
  env <- new.env(parent=globalenv());env$fit <- fit;env$jo <- jo
  out <- array(NA_real_,c(length(cells),ni,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    env$it <- it;env$ch <- ch
    env$eta <- eval(sc$b_eta,env)
    out[,it,ch] <- eval(sc$b_probability,env)[cells]
  }
  out
}

# Phase B: the archived design scorer gives the truth, the estimate, the group
# errors and the B0 and collection intercept truths; the retained draws give
# the intervals, B0 and collection intercept summaries and their correlation.
b_cells <- function(fit,input,job,sc) {
  r <- sc$score_current_fit(fit,input,job,sc$scoring)
  n <- input$scenario$n;S <- input$scenario$S;n0 <- input$design$original_sites
  stopifnot(identical(dim(r$truth),c(as.integer(n),as.integer(S))),identical(dim(r$estimate),dim(r$truth)),
    is.numeric(n0),length(n0)==1L,n0>=1,n0<=n)
  idx <- b_scored_cells(r$truth,n0)
  iv <- interval_cells(b_probability_draws(fit,sc,idx),as.vector(r$truth)[idx],sc)
  difference <- max(abs(iv$estimate-as.vector(r$estimate)[idx]))
  if(!is.finite(difference) || difference>=1e-12) stop('Draw reconstruction disagrees with the design scorer: ',difference)
  full <- function(v,na) {x <- rep(na,n*S);x[idx] <- v;x}
  ro <- fit$results_output;jo <- ro$jsdm_output
  el <- function(metric) {e <- r$elements[r$elements$metric==metric,,drop=FALSE];stopifnot(identical(e$element,seq_len(S)));e}
  b0 <- b0_cells(jo$B0_output,el('B0')$truth,sc)
  bt_draws <- array(ro$beta_theta_output[1,,,,drop=FALSE],dim(jo$B0_output))
  theta <- b0_cells(bt_draws,el('collection_intercept')$truth,sc)
  # Pearson correlation of B0 and the collection intercept over the pooled draws of each species.
  theta$correlation <- vapply(seq_len(S),function(s) stats::cor(as.vector(jo$B0_output[s,,]),as.vector(bt_draws[s,,])),numeric(1))
  estimate_check <- max(abs(c(b0$estimate-el('B0')$estimate,theta$estimate-el('collection_intercept')$estimate)))
  if(!is.finite(estimate_check) || estimate_check>=1e-12) stop('B0 or collection intercept estimates disagree with the design scorer: ',estimate_check)
  cells <- b_make_cells(r$truth,r$estimate,full(iv$covered,NA),n0,lower=full(iv$lower,NA_real_),upper=full(iv$upper,NA_real_),
    b0=b0,theta=theta)
  # The group errors must be the archived scorer's own: its original-site bands
  # (design_helpers.R:123-127) and its rare group over all fitted sites
  # (run_nonspatial_recheck_balanced.R:129-134).
  g <- b_group_table(cells,sc);gd <- 0
  pairs <- rbind(c('occupancy_original_sites','all','all'),c('occupancy_original_sites','low','low'),
    c('occupancy_original_sites','medium','middle'),c('occupancy_original_sites','high','high'),c('occupancy','rare_below_20pct','rare_below_20pct'))
  for(k in seq_len(nrow(pairs))) {
    a <- r$groups[r$groups$metric==pairs[k,1] & r$groups$group==pairs[k,2],,drop=FALSE];m <- g[g$group==pairs[k,3],,drop=FALSE]
    if(nrow(a)>1L || nrow(m)!=nrow(a) || (nrow(a) && m$cells!=a$n_elements)) stop('Cell sets disagree with the design scorer for ',pairs[k,3])
    if(nrow(a)) gd <- max(gd,abs(m$signed_error-100*a$bias),abs(m$mean_abs_cell_error-100*a$mae))
  }
  if(!setequal(g$group,pairs[pairs[,3] %in% g$group,3]) || gd>1e-10) stop('Group errors disagree with the design scorer: ',gd)
  cells$checks <- list(estimate_vs_score_fit=difference,groups_vs_score_fit=gd,parameters_vs_score_fit=estimate_check)
  cells
}

# ---------------------------------------------------------------------------
# Groups and community outcomes
# ---------------------------------------------------------------------------

# Cell sets (README "Bands and groups"): the truth bands below 0.2, 0.2 to 0.8
# inclusive and above 0.8 of the original sites; the rare group, every fitted
# site of each species whose mean truth over all fitted sites is below 0.2
# (descriptive, R22). There is no all-site scope in phase B.
b_cell_masks <- function(cells) {
  n <- cells$n;S <- cells$S;tv <- as.vector(cells$truth)
  original <- unlist(lapply(seq_len(S),function(s) (s-1L)*n+seq_len(cells$n_original)))
  rare <- unlist(lapply(which(colMeans(cells$truth)<.2),function(s) (s-1L)*n+seq_len(n)))
  sets <- list(all=original,low=original[tv[original]<.2],middle=original[tv[original]>=.2 & tv[original]<=.8],
    high=original[tv[original]>.8],rare_below_20pct=rare)
  lapply(names(sets),function(g) list(scope='primary',group=g,idx=as.integer(sets[[g]])))
}

# One row per non-empty cell set, as analysis.R's group_table.
b_group_table <- function(cells,sc) {
  tv <- as.vector(cells$truth);ev <- as.vector(cells$estimate);cv <- as.vector(cells$covered)
  rows <- lapply(b_cell_masks(cells),function(m) {
    if(!length(m$idx)) return(NULL)
    e <- sc$error_metrics(tv[m$idx],ev[m$idx])
    data.frame(scope=m$scope,group=m$group,cells=length(m$idx),truth_mean=e[['truth']],estimate_mean=e[['estimate']],
      signed_error=100*e[['bias']],abs_signed_error=abs(100*e[['bias']]),mean_abs_cell_error=100*e[['mae']],
      coverage=mean(cv[m$idx]),stringsAsFactors=FALSE)
  })
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# The collection intercept (beta_theta) bias and coverage and the B0 and
# collection intercept posterior correlation, averaged over species
# (README "Secondary outcomes", phase B only).
b_theta_summary <- function(theta) {
  if(is.null(theta)) return(data.frame(species=0L,bt_bias=NA_real_,bt_abs_bias=NA_real_,bt_coverage=NA_real_,b0_bt_correlation=NA_real_))
  data.frame(species=nrow(theta),bt_bias=mean(theta$bias),bt_abs_bias=mean(abs(theta$bias)),bt_coverage=mean(theta$covered),
    b0_bt_correlation=mean(theta$correlation))
}

b_score_record <- function(item,cells,sc,hashes,input_md5,fit_md5=item$expected_md5,warnings=NA_integer_) {
  n <- cells$n;S <- cells$S
  list(version=B_SCORE_VERSION,phase=item$phase,sd=as.numeric(item$sd),role=item$role,kind=item$kind %||% 'selected',
    key=item$key,schedule=item$schedule,stratum=b_stratum(item$key),community=b_community(item$key),
    fit_label=item$fit_label,fit_md5=fit_md5,input_md5=input_md5,warnings=warnings,hashes=hashes,
    groups=b_group_table(cells,sc),b0=b0_summary(cells$b0),b0_species=cells$b0,theta=b_theta_summary(cells$theta),theta_species=cells$theta,
    cells=data.frame(site=rep(seq_len(n),S),species=rep(seq_len(S),each=n),truth=as.vector(cells$truth),
      estimate=as.vector(cells$estimate),lower=as.vector(cells$lower),upper=as.vector(cells$upper),covered=as.vector(cells$covered)),
    checks=cells$checks,scored_at=format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'))
}

# Score one fit after checking it is the fit the item names (flags.R's
# evaluate_item: md5, job record, schedule, prior, input md5).
b_score_item <- function(item,sc,archives,inputs_root,out_root) {
  if(!identical(item$phase,'B')) stop('b_score_item scores phase B only')
  path <- score_path(out_root,item$phase,item$sd,item$key,item$schedule)
  if(file.exists(path)) {
    old <- readRDS(path);md5 <- unname(tools::md5sum(item$fit))
    if(identical(old$version,B_SCORE_VERSION) && identical(old$fit_md5,md5) && (is.na(item$expected_md5) || identical(md5,item$expected_md5)) &&
       identical(old$hashes,sc$score_hashes) && identical(old$role,item$role) && identical(old$sd,as.numeric(item$sd)) &&
       identical(old$kind,item$kind %||% 'selected')) {
      cat(format(Sys.time()),item$phase,paste0('sd',item$sd),item$key,item$schedule,'already scored; resumed\n');flush.console()
      return(path)
    }
  }
  spec <- phase_jobs(item$phase,archives,inputs_root,item$key)[[1]]
  captured <- NULL
  scorer <- function(phase,fit,warnings,input,sc) {captured <<- b_cells(fit,input,spec$job,sc);data.frame(n_warnings=length(warnings))}
  row <- evaluate_item(item,sc,archives,inputs_root,flag_fn=scorer)
  rec <- b_score_record(item,captured,sc,sc$score_hashes,input_md5=spec$input_md5,fit_md5=row$fit_md5,warnings=row$n_warnings)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);atomic_save(rec,path)
  cat(format(Sys.time()),item$phase,paste0('sd',item$sd),item$key,item$schedule,'scored\n');flush.console()
  path
}

analysis_b_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','phase','arms','workers','out'),
    required=c('repo','study','inputs-root','phase','arms'))
  if(!identical(o$phase,'B')) stop('--phase must be B: analysis-b.R scores phase B only (analysis.R scores A1 and A2)')
  arms <- parse_arms(o$arms);workers <- parse_workers(o$workers %||% '1')
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  out <- o$out %||% file.path(study,'summary');phase <- 'B'
  sel_dir <- selection_dir(study,phase);new_sds <- setdiff(arms,1)
  if(length(new_sds)) {
    f <- file.path(sel_dir,'selected-fits.csv')
    if(!file.exists(f)) stop('No new-arm occupancy error is computed before select.R --phases=B --mode=final has written ',f)
  }
  sc <- load_b_scorers(repo,archives)
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
  cat('Scoring',nrow(items),'fits of phase B (arms',paste(arms,collapse=', '),') into',file.path(out,phase),'\n');flush.console()
  work <- function(i) tryCatch(b_score_item(items[i,,drop=FALSE],sc,archives,inputs_root,out),
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
  if(length(repo_arg)!=1L || !identical(normalizePath(file.path(repo_arg,'dev/simstudy/occupancy-intercept-prior'),mustWork=FALSE),here)) {
    cat('ERROR: run the analysis-b.R of the --repo checkout (',here,' is not under --repo)\n',sep='');quit(save='no',status=1)
  }
  for(f in c('jobs.R','verify-helpers.R','flags.R','analysis.R')) source(file.path(here,f))
  status <- tryCatch(analysis_b_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
