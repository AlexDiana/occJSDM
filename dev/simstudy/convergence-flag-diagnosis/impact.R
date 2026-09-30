#!/usr/bin/env Rscript
# Task 4 of PLAN.md in this directory: the other four flagged fits, and the
# impact of the community 5 diagnosis on the conclusions that used its pr11
# fit. It reads saved fits, inputs and archived tables; it fits nothing and
# writes only results/impact/ (and, with --cache, one scratch file). Archives
# and the current-main-recheck and occupancy-intercept-prior directories are
# read, never written.
#
#   Rscript impact.R [--out=DIR] [--cache=FILE]
#
# --cache=FILE keeps the scored fits between runs while the tables are being
# developed; the committed results come from a run without it.
#
# Step 1 (brief): the Task 1 labels of the other four flagged fits and the
#   sources of their pr11 flags (results/anatomy/ and current-main-recheck
#   results/). No species of theirs is separated, so none gets an extended run.
# Step 2: community 5 of the high-contamination 300-site design scored four
#   ways: the pr11 pooled estimate, reproduced from the saved fit; the pooled
#   means of the pr11 chain groups {1,3} (near-truth) and {2,4} (mirror); and
#   the pooled means of the 16-chain extended run of Task 3. Each is substituted
#   for the published community 5 scores and the affected comparisons of
#   current-main-recheck/REPORT.md are recomputed.
# Step 3: the intercept-prior phase B gate recomputed with the extended run in
#   place of the flagged control fit of community 5 (qfar), and that control
#   counted as unflagged if the pr11 flag rule, applied to the 16 extended
#   chains bound into one fit (validated first on the pr11 fit), leaves it
#   unflagged; a descriptive sensitivity, not a re-run of the gate.
#
# Definitions are the archived ones, loaded after md5 checks, not copied:
#   current-main-recheck/helpers.R (load_scoring, score_current_fit) and the
#     scorer files it builds from, through the intercept-prior phase B loader
#     load_b_scorers (analysis-b.R), which checks their md5s against the pr11
#     settings record (flags.R checked_scorer_file); the intercept-prior
#     scripts are checked against that study's results/provenance-B.csv.
#   current-main-recheck/summarise.R: normalise_groups and mean_ci (md5 checked
#     against current-main-recheck/results/source-hashes.csv). Its inline
#     arm, design-contrast and code-version aggregation is restated in
#     arm_summary(), design_contrast() and version_change(); test-impact.R
#     proves them by reproducing chain-sensitivity-summary.csv,
#     design-paired-summary.csv and paired-summary.csv.
#   The design scorer's original-site groups (design_helpers.R:117-127) are
#     restated in recheck_groups() with the archived error_metrics
#     (jsdm-sample-size-recheck/helpers.R:45-47); the proof is the reproduction
#     of the published community 5 rows from the saved fit.
#   Phase B cells and groups: analysis-b.R (b_cells, b_probability_draws,
#     b_scored_cells, b_make_cells, b_group_table), analysis.R (interval_cells)
#     and spatial-targeted-recheck/score.R (score_draw_block); the gate:
#     occupancy-intercept-prior/summarise.R (gate_subphase and its criteria).
# Inputs are read through the path remap with md5 refusal (jobs.R); nothing
# under /Users/douglasyu/Documents is read.

`%||%` <- function(x,y) if(is.null(x)) y else x

IMPACT_DIR <- local({
  here <- NULL
  for(i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if(!is.null(f)) {here <- dirname(normalizePath(f));break}
  }
  if(is.null(here)) {
    a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
    if(length(a)) here <- dirname(normalizePath(sub('^--file=','',a[1])))
  }
  here %||% getwd()
})
# REPO, ARCHIVES, PR11_ARCHIVE, INPUTS_ROOT, RECHECK_RESULTS, FLAGGED_KEYS,
# REFERENCE_KEYS, selected_fit() and the jobs.R remap come from anatomy.R.
source(file.path(IMPACT_DIR,'anatomy.R'))

TARGET_KEY <- 'design-qfar_K6-sites300-05'
PR11_FIT_MD5 <- '07b2d9eb4356b7973f87485a0f8bcd09'
SOURCE_INPUT_MD5 <- '574ed4df46b9c1bd227c79d2185a7fcb'
ANATOMY_MD5 <- 'bd98ff3f40cad972f96bdeca5a07e206'  # frozen in AMENDMENT-1.md
EXTENDED_CHAINS <- 16L
RECHECK_DIR <- file.path(REPO,'dev/simstudy/current-main-recheck')
IP_REL <- 'dev/simstudy/occupancy-intercept-prior'
IP_DIR <- file.path(REPO,IP_REL)
IP_RESULTS <- file.path(IP_DIR,'results')
IP_ARCHIVE <- file.path(ARCHIVES,'intercept-prior-20260929')
DIAG_ARCHIVE <- file.path(ARCHIVES,'convergence-diagnosis-20261001')
IP_SCRIPTS <- c('jobs.R','verify-helpers.R','flags.R','analysis.R','analysis-b.R','summarise.R')
ANATOMY_RESULTS <- file.path(IMPACT_DIR,'results/anatomy')
IMPACT_OUT <- file.path(IMPACT_DIR,'results/impact')
PP <- 100  # probability scale to percentage points

# ---- Archived definitions ----------------------------------------------------

# The intercept-prior scripts, each checked against that study's scorer
# provenance record where it has one (summarise.R has none; the reproduction of
# gate-detail-B.csv in test-impact.R is its proof), sourced for their
# definitions only (each file's command-line block is guarded).
load_ip <- function() {
  prov <- utils::read.csv(file.path(IP_RESULTS,'provenance-B.csv'),stringsAsFactors=FALSE)
  recorded <- stats::setNames(prov$value[prov$kind=='hash'],prov$name[prov$kind=='hash'])
  ip <- new.env(parent=globalenv())
  for(f in IP_SCRIPTS) {
    path <- file.path(IP_DIR,f);rel <- file.path(IP_REL,f)
    if(rel %in% names(recorded) && !identical(unname(tools::md5sum(path)),recorded[[rel]]))
      stop('Intercept-prior script differs from its provenance record: ',rel)
    sys.source(path,envir=ip)
  }
  ip
}

load_recheck_summaries <- function(ip) {
  path <- file.path(RECHECK_DIR,'summarise.R')
  x <- utils::read.csv(file.path(RECHECK_RESULTS,'source-hashes.csv'),stringsAsFactors=FALSE)
  recorded <- ip$recorded_hash(stats::setNames(x$md5,x$path),'dev/simstudy/current-main-recheck/summarise.R')
  if(!identical(unname(tools::md5sum(path)),recorded)) stop('current-main-recheck/summarise.R differs from its record')
  env <- new.env(parent=globalenv())
  ip$extract_definitions(path,c('normalise_groups','mean_ci'),env)
  env
}

impact_setup <- function(archives=ARCHIVES,inputs_root=INPUTS_ROOT) {
  if(!identical(unname(tools::md5sum(file.path(IMPACT_DIR,'anatomy.R'))),ANATOMY_MD5))
    stop('anatomy.R differs from the md5 frozen in AMENDMENT-1.md')
  ip <- load_ip()
  sc <- ip$load_b_scorers(REPO,archives)
  list(ip=ip,sc=sc,recheck=load_recheck_summaries(ip),flag=load_flag_rule(ip,archives),archives=archives,inputs_root=inputs_root)
}

# ---- The pr11 flag rule ------------------------------------------------------------

# The pr11 flag rule for two-stage fits (current-main-recheck/select.R:17-19 on
# the diagnostics of its run.R:79-81) as the intercept-prior study applies it
# to saved fits: flags.R's fit_flags for phase B (b_rhats, pr11_diagnostics and
# the archived flag expression), with the archived trace summary, loaded by
# load_flag_scorers after its md5 checks and validated there on all 70 control
# fits. The fitting warnings the rule counts are those of the package's own
# computeDiagnostics (R/diagnostics.R), which runOccJSDM calls on its final
# output; its definition is parsed from the source the extended run's library
# was installed from, after checking that file against the pr11 source record
# (the two revisions' files are identical).
load_flag_rule <- function(ip,archives=ARCHIVES) {
  fr <- ip$load_flag_scorers(REPO,archives)
  path <- file.path(DIAG_ARCHIVE,'source/R/diagnostics.R')
  x <- utils::read.csv(file.path(RECHECK_RESULTS,'source-hashes.csv'),stringsAsFactors=FALSE)
  recorded <- ip$recorded_hash(stats::setNames(x$md5,x$path),'source-main/R/diagnostics.R')
  if(!identical(unname(tools::md5sum(path)),recorded)) stop('R/diagnostics.R of the diagnosis library differs from the pr11 source record')
  env <- new.env(parent=globalenv());ip$extract_definitions(path,'computeDiagnostics',env)
  fr$computeDiagnostics <- env$computeDiagnostics;fr$diagnostics_file <- path;fr$diagnostics_md5 <- recorded
  fr
}

# The fitting warnings computeDiagnostics issues on a results_output, its
# printed summary silenced.
package_warnings <- function(ro,computeDiagnostics) {
  w <- character()
  withCallingHandlers(suppressMessages(computeDiagnostics(ro)),
    warning=function(x) {w <<- c(w,conditionMessage(x));invokeRestart('muffleWarning')})
  w
}

# Single-chain results_output lists bound along the chain dimension into the
# multi-chain results_output runOccJSDM would return for these chains: every
# array whose last two dimensions are (iterations, 1), at both levels. The
# per-fit posterior summaries (z, psi, w and theta outputs) and WAIC have no
# chain dimension and are left out; computeDiagnostics does not read them.
bind_chains <- function(parts) {
  K <- length(parts);ni <- dim(parts[[1]]$jsdm_output$B0_output)[2]
  level <- function(get) {
    first <- get(parts[[1]])
    keep <- names(first)[vapply(first,function(a) {d <- dim(a);length(d)>=2L && d[length(d)]==1L && d[length(d)-1L]==ni},logical(1))]
    out <- list()
    for(nm in keep) {
      d <- dim(first[[nm]])
      xs <- lapply(parts,function(p) {a <- get(p)[[nm]];if(!identical(dim(a),d)) stop('Chains differ in the dimensions of ',nm);a})
      d[length(d)] <- K;out[[nm]] <- array(unlist(xs,use.names=FALSE),d)
    }
    out
  }
  ro <- level(function(p) p[names(p)!='jsdm_output'])
  ro$jsdm_output <- level(function(p) p$jsdm_output)
  ro
}

read_recheck <- function(name) utils::read.csv(file.path(RECHECK_RESULTS,name),stringsAsFactors=FALSE)

# ---- Current-main-recheck definitions ----------------------------------------

# The original-site groups of the design scorer (design_helpers.R:117-127):
# the first n0 sites of every species, and within them the truth bands below
# 0.2, 0.2 to 0.8 inclusive and above 0.8, each only when non-empty; errors of
# the posterior-mean estimate by the archived error_metrics.
recheck_groups <- function(truth,estimate,n0,error_metrics) {
  truth <- as.matrix(truth);estimate <- as.matrix(estimate)
  n <- nrow(truth);S <- ncol(truth)
  stopifnot(identical(dim(truth),dim(estimate)),length(n0)==1L,n0>=1,n0<=n)
  base <- unlist(lapply(seq_len(S),function(j) (j-1L)*n+seq_len(n0)))
  bt <- as.vector(truth)[base];be <- as.vector(estimate)[base]
  sets <- list(all=seq_along(base),low=which(bt<.2),medium=which(bt>=.2 & bt<=.8),high=which(bt>.8))
  rows <- lapply(names(sets),function(g) {
    ix <- sets[[g]];if(!length(ix)) return(NULL)
    e <- error_metrics(bt[ix],be[ix])
    data.frame(metric='occupancy_original_sites',group=g,n_elements=length(ix),truth=e[['truth']],
      estimate=e[['estimate']],bias=e[['bias']],mae=e[['mae']],rmse=e[['rmse']],stringsAsFactors=FALSE)
  })
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# The published community scores with one community's rows of one version
# replaced by a substitute estimate's groups. The cell sets and truth must be
# those of the rows replaced; chain diagnostics of a substitute are undefined.
substitute_scores <- function(scores,key,version,groups,metric='occupancy_original_sites') {
  hit <- which(scores$key==key & scores$version==version & scores$metric==metric)
  if(!length(hit) || nrow(groups)!=length(hit) || anyDuplicated(groups$group) || !setequal(groups$group,scores$group[hit]))
    stop('Substitute groups differ from the rows they replace for ',key)
  g <- groups[match(scores$group[hit],groups$group),,drop=FALSE]
  if(any(g$n_elements!=scores$n_elements[hit])) stop('Substitute cells differ from the rows they replace for ',key)
  if(max(abs(g$truth-scores$truth[hit]))>1e-12) stop('Substitute truth differs from the rows it replaces for ',key)
  for(k in c('estimate','bias','mae','rmse')) scores[[k]][hit] <- g[[k]]
  for(k in intersect(c('rhat','ess_mean','mcse','chain_gap'),names(scores))) scores[[k]][hit] <- NA_real_
  scores
}

primary_rows <- function(scores,version,family,scenario,arm,group) {
  p <- scores[scores$version==version & scores$family==family & scores$scenario==scenario & scores$arm==arm &
    scores$metric=='occupancy_original_sites' & scores$group==group,,drop=FALSE]
  if(nrow(p)!=10L || !setequal(p$replicate,1:10)) stop('Expected ten communities for ',paste(version,scenario,arm,group))
  p
}

# Ten-community means of one arm (chain-sensitivity.R's selected_mae and selected_bias).
arm_summary <- function(scores,version,scenario,arm,group,family='design') {
  p <- primary_rows(scores,version,family,scenario,arm,group)
  c(mae=mean(p$mae),bias=mean(p$bias))
}

# Reference-design MAE minus alternative-design MAE per community, its mean and
# 95% t interval over the ten communities, and how many communities improve
# (summarise.R's design_pairs and design_summary).
design_contrast <- function(scores,version,scenario,reference,alternative,group,mean_ci,family='design') {
  d <- merge(primary_rows(scores,version,family,scenario,reference,group),primary_rows(scores,version,family,scenario,alternative,group),
    by=c('family','scenario','replicate','group'),suffixes=c('_reference','_alternative'))
  stopifnot(nrow(d)==10L,max(abs(d$truth_reference-d$truth_alternative))<1e-12,all(d$n_elements_reference==d$n_elements_alternative))
  reduction <- d$mae_reference-d$mae_alternative;ci <- mean_ci(reduction)
  data.frame(reduction_mae=ci[['mean']],lower=ci[['lower']],upper=ci[['upper']],improved=sum(reduction>0),communities=10L)
}

# Current minus archived code, paired by community (summarise.R's paired and summary_rows).
version_change <- function(scores,scenario,arm,group,mean_ci,family='design') {
  keys <- c('key','family','arm','scenario','replicate','metric','group')
  p <- merge(primary_rows(scores,'historical',family,scenario,arm,group),primary_rows(scores,'current',family,scenario,arm,group),
    by=keys,suffixes=c('_historical','_current'))
  stopifnot(nrow(p)==10L,max(abs(p$truth_historical-p$truth_current))<1e-12)
  m <- mean_ci(p$mae_current-p$mae_historical);b <- mean_ci(p$bias_current-p$bias_historical)
  data.frame(historical_mae=mean(p$mae_historical),current_mae=mean(p$mae_current),delta_mae_mean=m[['mean']],
    delta_mae_lower=m[['lower']],delta_mae_upper=m[['upper']],historical_bias=mean(p$bias_historical),
    current_bias=mean(p$bias_current),delta_bias_mean=b[['mean']],delta_bias_lower=b[['lower']],delta_bias_upper=b[['upper']])
}

# The published observed-chain sensitivity range (chain-sensitivity.R): for an
# arm, of the ten-community mean MAE or signed error; for a design contrast
# ('reference:alternative'), of the MAE reduction.
sensitivity_range <- function(kind,scenario,what,group,column) {
  if(kind=='arm') {
    s <- read_recheck('chain-sensitivity-summary.csv')
    r <- s[s$family=='design' & s$scenario==scenario & s$arm==what & s$group==group,,drop=FALSE]
    stopifnot(nrow(r)==1L,column %in% c('mae','bias'))
    return(c(min=r[[paste0('minimum_',column)]],max=r[[paste0('maximum_',column)]]))
  }
  stopifnot(kind=='contrast',column=='reduction_mae')
  arms <- strsplit(what,':',fixed=TRUE)[[1]]
  s <- read_recheck('design-chain-sensitivity.csv')
  r <- s[s$version=='current' & s$family=='design' & s$scenario==scenario & s$reference==arms[1] & s$alternative==arms[2] & s$group==group,,drop=FALSE]
  stopifnot(nrow(r)==1L)
  c(min=r$minimum_reduction,max=r$maximum_reduction)
}

# ---- Phase B definitions -------------------------------------------------------

# The interval columns interval_cells() (analysis.R) takes from score_draw_block().
iv_from_block <- function(e) data.frame(truth=e$truth,estimate=e$estimate,lower=e$lower,upper=e$upper,covered=e$covered)

# Per-community group table of the phase B selected fits, as summarise.R passes
# it to the gate (record_table of the selected records' groups), rebuilt from
# the archive's band-error.csv and coverage.csv after md5 checks against the
# study's audit record (results/audit-source-B.csv).
phase_b_groups <- function(archive=IP_ARCHIVE) {
  audit <- utils::read.csv(file.path(IP_RESULTS,'audit-source-B.csv'),stringsAsFactors=FALSE)
  read <- function(name) {
    f <- file.path(archive,'summary/B',name);rec <- audit$md5[audit$file==file.path('summary/B',name)]
    if(length(rec)!=1L || !identical(unname(tools::md5sum(f)),rec)) stop('Archived table differs from its audit record: ',f)
    utils::read.csv(f,stringsAsFactors=FALSE,colClasses=c(community='character',fit_md5='character'))
  }
  be <- read('band-error.csv');cv <- read('coverage.csv')
  meta <- c('phase','sd','arm','kind','key','stratum','community','schedule','fit','fit_md5','scope','group','cells')
  g <- merge(be,cv,by=meta,sort=FALSE)
  if(nrow(g)!=nrow(be) || nrow(g)!=nrow(cv)) stop('band-error.csv and coverage.csv do not pair up')
  g <- g[g$kind=='selected',,drop=FALSE]
  g <- g[order(g$sd,g$key,match(g$group,c('all','low','middle','high','rare_below_20pct'))),,drop=FALSE]
  rownames(g) <- NULL;g
}

phase_b_convergence <- function() {
  f <- file.path(IP_RESULTS,'convergence-B.csv')
  m <- utils::read.csv(file.path(IP_RESULTS,'evidence-manifest.csv'),stringsAsFactors=FALSE)
  if(!identical(unname(tools::md5sum(f)),m$md5[m$file=='convergence-B.csv'])) stop('convergence-B.csv differs from its evidence record')
  utils::read.csv(f,stringsAsFactors=FALSE)
}

# gate_subphase (intercept-prior summarise.R) for each new SD against the control.
phase_b_gate <- function(groups,conv,ip,sds=c(2,3,5)) {
  res <- lapply(sds,function(sd) ip$gate_subphase('B',groups,conv,sd,NULL))
  d <- do.call(rbind,lapply(res,`[[`,'detail'));rownames(d) <- NULL
  list(detail=d,rows=ip$rbind_fill(lapply(res,`[[`,'row')))
}

# The group rows of one community of one arm replaced by a substitute's.
substitute_phase_b <- function(groups,key,new,sd=1,label=NULL) {
  hit <- which(groups$sd==sd & groups$key==key)
  if(!length(hit) || nrow(new)!=length(hit) || anyDuplicated(new$group) || !setequal(new$group,groups$group[hit]))
    stop('Substitute groups differ from the rows they replace for ',key)
  n <- new[match(groups$group[hit],new$group),,drop=FALSE]
  if(any(n$cells!=groups$cells[hit]) || any(n$scope!=groups$scope[hit])) stop('Substitute cells differ from the rows they replace for ',key)
  if(max(abs(n$truth_mean-groups$truth_mean[hit]))>1e-12) stop('Substitute truth differs from the rows it replaces for ',key)
  for(k in c('estimate_mean','signed_error','abs_signed_error','mean_abs_cell_error','coverage')) groups[[k]][hit] <- n[[k]]
  if(!is.null(label)) {groups$fit[hit] <- label;groups$fit_md5[hit] <- NA_character_}
  groups
}

# The control counted with one fewer flagged selected fit at one contamination level.
converged_convergence <- function(conv,stratum,fewer=1L) {
  at <- conv$sd==1 & conv$stratum==stratum
  if(sum(at)!=1L || conv$selected_flagged[at]<fewer) stop('No control count to lower at ',stratum)
  conv$selected_flagged[at] <- conv$selected_flagged[at]-as.integer(fewer)
  conv
}

# Each new arm counted with one fewer flagged selected fit at the key's
# contamination level, for each arm whose selected fit of the key is flagged
# (results/selection-selected-fits-B.csv, md5 checked against the study's
# evidence record): the bound on criterion 4 if every arm's fit of the
# community were as converged as the extended run.
arms_unflagged_convergence <- function(conv,key) {
  f <- file.path(IP_RESULTS,'selection-selected-fits-B.csv')
  m <- utils::read.csv(file.path(IP_RESULTS,'evidence-manifest.csv'),stringsAsFactors=FALSE)
  if(!identical(unname(tools::md5sum(f)),m$md5[m$file=='selection-selected-fits-B.csv'])) stop('selection-selected-fits-B.csv differs from its evidence record')
  sel <- utils::read.csv(f,stringsAsFactors=FALSE);sel <- sel[sel$key==key & sel$role=='new',,drop=FALSE]
  stratum <- sub('^design-(qnear|qfar)_K6-sites300-[0-9]{2}$','\\1',key)
  for(sd in sel$sd[sel$flagged %in% TRUE]) {
    at <- conv$sd==sd & conv$stratum==stratum
    if(sum(at)!=1L || conv$selected_flagged[at]<1L) stop('No count to lower for SD ',sd,' at ',stratum)
    conv$selected_flagged[at] <- conv$selected_flagged[at]-1L
  }
  conv
}

# ---- The other four flagged fits ----------------------------------------------

# Species of one row of current-main-recheck flagged-diagnostics.csv, from the
# element layouts of the design scorer: per-species blocks (score_fit blocks
# B0, theta0, collection intercept and slope), primer-by-species blocks p, q and
# q_nominal (column-major, two primers), environment slopes (B_output,
# covariate-by-species, two covariates; helpers.R diagnostic_block) and
# original-site probabilities (base_idx, the original sites of each species in
# turn; design_helpers.R:121). Groups and sigma_h belong to no species.
element_species <- function(block,metric,element,n0=100L,n_cov=2L,n_primer=2L) {
  element <- as.integer(element)
  if(identical(block,'groups') || is.na(element)) return(NA_integer_)
  if(identical(block,'elements')) {
    if(metric %in% c('B0','theta0','collection_intercept','collection_slope')) return(element)
    if(metric %in% c('p','q','q_nominal')) return(as.integer((element-1L)%/%n_primer+1L))
    stop('Unknown element metric ',metric)
  }
  if(identical(block,'additional_diagnostics')) return(switch(metric,
    environment_slope=as.integer((element-1L)%/%n_cov+1L),
    original_probability=as.integer((element-1L)%/%as.integer(n0)+1L),
    sigma_h=NA_integer_,stop('Unknown additional diagnostic ',metric)))
  stop('Unknown diagnostic block ',block)
}

count_label <- function(x) if(!length(x)) '' else {t <- table(x);paste(sprintf('%s %d',names(t),as.integer(t)),collapse='; ')}
split_list <- function(x) if(is.na(x) || !nzchar(x)) character() else strsplit(x,';',fixed=TRUE)[[1]]

# One row per fit and species: the Task 1 label (from its non-loading
# quantities, ruling R6) with the quantities behind it, and the sources of the
# fit's pr11 flag (current-main-recheck diagnostics and warnings).
other_flagged_fits <- function(keys=setdiff(FLAGGED_KEYS,TARGET_KEY)) {
  text_cols <- c('separated_quantities','drifting_quantities','separated_loadings','drifting_loadings')
  labels <- utils::read.csv(file.path(ANATOMY_RESULTS,'species-labels.csv'),stringsAsFactors=FALSE,
    colClasses=stats::setNames(rep('character',length(text_cols)),text_cols))
  for(k in text_cols) labels[[k]][is.na(labels[[k]])] <- ''
  diag <- read_recheck('diagnostic-summary.csv');flags <- read_recheck('flagged-diagnostics.csv');warn <- read_recheck('warnings.csv')
  rows <- list()
  for(key in keys) {
    l <- labels[labels$key==key,,drop=FALSE];d <- diag[diag$key==key,,drop=FALSE]
    if(nrow(l)!=10L || nrow(d)!=1L) stop('Expected ten species labels and one diagnostic row for ',key)
    f <- flags[flags$key==key,,drop=FALSE]
    f$species <- if(nrow(f)) mapply(element_species,f$block,f$metric,f$element,USE.NAMES=FALSE) else integer()
    w <- warn$warning[warn$key==key]
    community <- f[is.na(f$species),,drop=FALSE]
    by_species <- vapply(sort(unique(f$species[!is.na(f$species)])),function(s)
      sprintf('species %d: %s',s,count_label(f$metric[f$species %in% s])),character(1))
    source <- c(if(nrow(f)) sprintf('%d scored diagnostics with Rhat above 1.05 (%s%s)',nrow(f),paste(by_species,collapse='; '),
      if(nrow(community)) paste0('; community groups: ',count_label(community$metric)) else '') else 'no scored diagnostic above Rhat 1.05',
      if(length(w)) sprintf('%d package warning%s (%s)',length(w),if(length(w)>1L) 's' else '',paste(w,collapse=' | ')) else 'no package warning')
    separated <- any(l$label=='separated')
    for(s in seq_len(10L)) {
      ls <- l[l$species==s,,drop=FALSE];fs <- f[f$species %in% s,,drop=FALSE]
      rows[[length(rows)+1L]] <- data.frame(key=key,arm=d$arm,scenario=d$scenario,replicate=d$replicate,schedule=d$schedule,species=s,
        task1_label=ls$label,separated_quantities=ls$separated_quantities,drifting_quantities=ls$drifting_quantities,
        n_drifting_quantities=length(split_list(ls$drifting_quantities)),n_separated_quantities=length(split_list(ls$separated_quantities)),
        loading_labels=paste0(if(nzchar(ls$separated_loadings)) paste0('separated: ',ls$separated_loadings) else '',
          if(nzchar(ls$drifting_loadings)) paste0('drifting: ',ls$drifting_loadings) else ''),
        max_separation=ls$max_separation,max_rhat_all_chains=ls$max_rhat,
        flagged_diagnostics=nrow(fs),flagged_metrics=count_label(fs$metric),
        max_flagged_rhat=if(nrow(fs)) max(fs$rhat) else NA_real_,
        fit_package_warnings=length(w),package_warnings=paste(w,collapse=' | '),
        fit_max_group_rhat=d$max_group_rhat,fit_max_element_rhat=d$max_element_rhat,
        fit_flag_sources=paste(source,collapse='; '),fit_any_separated=separated,extended_run_needed=separated,
        stringsAsFactors=FALSE)
    }
  }
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# ---- Chain groups ---------------------------------------------------------------

# The pr11 chain groups of Task 3: the chains the anchored classifier labels
# mirror in-sample (results/modes/anchor-validation.csv), and the rest; checked
# against Task 1's grouping (results/anatomy/chain-groups.csv, species 6).
pr11_chain_groups <- function() {
  a <- utils::read.csv(file.path(IMPACT_DIR,'results/modes/anchor-validation.csv'),stringsAsFactors=FALSE)
  a <- a[a$check=='in-sample share of draws labelled mirror',,drop=FALSE]
  stopifnot(identical(sort(as.integer(a$chain)),1:4))
  mirror <- sort(as.integer(a$chain[a$value>.5]));near <- sort(as.integer(a$chain[a$value<.5]))
  g <- utils::read.csv(file.path(ANATOMY_RESULTS,'chain-groups.csv'),stringsAsFactors=FALSE)
  g <- g[g$key==TARGET_KEY & g$species==6L,,drop=FALSE]
  if(!identical(sort(as.integer(g$chain[g$chain_group==2L])),near) || !identical(sort(as.integer(g$chain[g$chain_group==1L])),mirror))
    stop('Task 1 and Task 3 chain groups disagree')
  list(near_truth=near,mirror=mirror)
}

# ---- Scoring the saved fits (slow) -----------------------------------------------

# The saved pr11 fit of community 5, scored by the archived phase B cell scorer
# (b_cells), which runs the current-main-recheck scorer inside it; the latter's
# result is kept as well. Per-chain posterior means come from the archived
# per-draw statements (b_probability_draws), the chain-group intervals from
# score_draw_block on the pooled draws of the group's chains.
score_pr11 <- function(ctx,key=TARGET_KEY) {
  ip <- ctx$ip;sc <- ctx$sc
  spec <- ip$phase_jobs('B',ctx$archives,ctx$inputs_root,key)[[1]]
  fit_file <- selected_fit(key);md5 <- unname(tools::md5sum(fit_file))
  ver <- read_recheck('verification-selected.csv')
  if(!identical(md5,PR11_FIT_MD5) || !identical(md5,ver$fit_md5[ver$key==key])) stop('pr11 fit md5 differs from its verification record: ',fit_file)
  saved <- readRDS(fit_file)
  if(!identical(saved$job[names(spec$job)],spec$job)) stop('The saved fit job record differs from the archived pr11 job: ',key)
  input <- readRDS(ip$checked_input(spec$input_file,spec$input_md5))
  log_step('pr11 fit and input read (md5s match); scoring with the archived scorers')
  captured <- NULL;original <- sc$score_current_fit
  sc$score_current_fit <- function(...) {captured <<- original(...);captured}
  cells <- tryCatch(ip$b_cells(saved$fit,input,spec$job,sc),finally=sc$score_current_fit <- original)
  r <- captured;if(is.null(r)) stop('The design scorer was not run')
  log_step('pr11 flag rule and package diagnostics on the saved fit')
  flags <- list(row=ip$fit_flags('B',saved$fit,saved$warnings,input,ctx$flag),warnings_saved=saved$warnings,
    warnings_recomputed=package_warnings(saved$fit$results_output,ctx$flag$computeDiagnostics))
  recheck_rows <- ctx$recheck$normalise_groups(r,spec$job)
  truth <- r$truth;n <- nrow(truth);S <- ncol(truth);n0 <- as.integer(input$design$original_sites)
  log_step('archived scorers done; reconstructing per-chain draws')
  draws <- ip$b_probability_draws(saved$fit,sc,seq_len(n*S))
  nc <- dim(draws)[3];ni <- dim(draws)[2]
  chain_means <- vapply(seq_len(nc),function(ch) rowMeans(draws[,,ch]),numeric(n*S))
  pooled_check <- max(abs(rowMeans(chain_means)-as.vector(r$estimate)))
  if(!(pooled_check<1e-12)) stop('Per-chain means do not pool to the scorer estimate: ',pooled_check)
  idx <- ip$b_scored_cells(truth,n0)
  kept <- draws[idx,,,drop=FALSE];rm(draws);invisible(gc())
  tv <- as.vector(truth)[idx]
  groups <- pr11_chain_groups()
  block <- function(chains) {
    b <- sc$score_draw_block(kept[,,chains,drop=FALSE],tv,'cell')$elements
    list(iv=iv_from_block(b),max_cell_rhat=max(b$rhat,na.rm=TRUE))
  }
  log_step('interval blocks: all four chains')
  all4 <- block(seq_len(nc))
  b_iv <- data.frame(truth=as.vector(cells$truth)[idx],estimate=as.vector(cells$estimate)[idx],lower=as.vector(cells$lower)[idx],
    upper=as.vector(cells$upper)[idx],covered=as.vector(cells$covered)[idx])
  interval_check <- max(abs(unlist(all4$iv[c('truth','estimate','lower','upper')])-unlist(b_iv[c('truth','estimate','lower','upper')])))
  if(!(interval_check<1e-12) || !identical(all4$iv$covered,b_iv$covered)) stop('Interval blocks disagree with b_cells: ',interval_check)
  sources <- list(pr11_pooled=list(chains=seq_len(nc),draws=ni*nc,estimate=r$estimate,iv=b_iv,max_cell_rhat=all4$max_cell_rhat))
  for(g in names(groups)) {
    log_step(paste('interval blocks: chain group',g))
    ch <- groups[[g]];b <- block(ch)
    est <- matrix(rowMeans(chain_means[,ch,drop=FALSE]),n,S)
    if(!(max(abs(b$iv$estimate-est[idx]))<1e-12)) stop('Chain-group estimate disagrees with its draws: ',g)
    sources[[g]] <- list(chains=ch,draws=ni*length(ch),estimate=est,iv=b$iv,max_cell_rhat=b$max_cell_rhat)
  }
  list(key=key,fit_file=fit_file,fit_md5=md5,input_md5=spec$input_md5,truth=truth,n0=n0,idx=idx,chain_means=chain_means,
    chain_groups=groups,sources=sources,recheck_rows=recheck_rows,phase_b_cells=cells,X_psi=saved$fit$X_psi,X_theta=saved$fit$X_theta,flags=flags,
    species=saved$fit$infos$speciesNames,n_warnings=length(saved$warnings),
    checks=c(pooled_vs_scorer=pooled_check,intervals_vs_b_cells=interval_check,unlist(cells$checks)))
}

# The 16 single-chain fits of the extended run (Task 3), each checked against
# its recorded md5, the community 5 input and the pr11 design; the pooled
# estimate is the mean of the chain means (equal draws per chain), the
# intervals and cell Rhats come from score_draw_block on all pooled draws. The
# chains are then bound into one 16-chain fit, on which the package
# diagnostics give the fitting warnings and the pr11 flag rule is applied, as
# for a fit that had run the 16 chains together.
score_extended <- function(ctx,base,chains=seq_len(EXTENDED_CHAINS)) {
  prov <- utils::read.csv(file.path(IMPACT_DIR,'results/extended/provenance.csv'),stringsAsFactors=FALSE)
  prov <- prov[prov$run=='extended',,drop=FALSE]
  n <- nrow(base$truth);S <- ncol(base$truth);idx <- base$idx
  chain_means <- matrix(NA_real_,n*S,length(chains));kept <- NULL;psi_check <- 0;files <- list();parts <- list();shell <- NULL
  for(j in seq_along(chains)) {
    row <- prov[prov$chain==chains[j],,drop=FALSE];if(nrow(row)!=1L) stop('No provenance row for extended chain ',chains[j])
    file <- file.path(ctx$archives,row$fit_file);md5 <- unname(tools::md5sum(file))
    if(!identical(md5,row$fit_md5)) stop('Extended fit md5 differs from its record: ',file)
    saved <- readRDS(file);fit <- saved$fit
    if(!identical(saved$input_md5,SOURCE_INPUT_MD5) || !identical(saved$run,'extended') || !identical(as.integer(saved$chain),chains[j]) ||
       length(saved$warnings)!=0L || !identical(fit$infos$speciesNames,base$species) || !identical(dim(fit$X_psi),dim(base$X_psi)) ||
       max(abs(unname(fit$X_psi)-unname(base$X_psi)))>1e-12 || !identical(dim(fit$X_theta),dim(base$X_theta)) ||
       max(abs(unname(fit$X_theta)-unname(base$X_theta)))>1e-12 || dim(fit$results_output$jsdm_output$B0_output)[3]!=1L)
      stop('Extended fit is not the expected chain of community 5: ',file)
    d <- ctx$ip$b_probability_draws(fit,ctx$sc,seq_len(n*S))
    m <- rowMeans(d[,,1]);psi_check <- max(psi_check,max(abs(m-as.vector(fit$results_output$psi_output))))
    if(!(psi_check<1e-10)) stop('Reconstructed draws disagree with the stored posterior mean: ',file)
    if(is.null(kept)) kept <- array(NA_real_,c(length(idx),dim(d)[2],length(chains)))
    kept[,,j] <- d[idx,,1];chain_means[,j] <- m
    ro <- fit$results_output;parts[[j]] <- ro[setdiff(names(ro),c('z_output','psi_output','w_output','theta_output'))]
    if(is.null(shell)) shell <- fit[setdiff(names(fit),'results_output')]
    rm(ro)
    files[[j]] <- data.frame(chain=chains[j],fit_file=row$fit_file,fit_md5=md5,draws=dim(d)[2])
    rm(saved,fit,d);invisible(gc())
    log_step(sprintf('extended chain %d read and reconstructed',chains[j]))
  }
  estimate <- matrix(rowMeans(chain_means),n,S)
  b <- ctx$sc$score_draw_block(kept,as.vector(base$truth)[idx],'cell')$elements
  est_check <- max(abs(b$estimate-estimate[idx]))
  if(!(est_check<1e-12)) stop('Pooled extended estimate disagrees with its draws: ',est_check)
  rm(kept);invisible(gc())
  log_step('pr11 flag rule and package diagnostics on the',length(chains),'bound chains')
  combined <- c(shell,list(results_output=bind_chains(parts)));rm(parts);invisible(gc())
  if(!identical(dim(combined$results_output$jsdm_output$B0_output)[3],length(chains))) stop('The bound fit does not hold every chain')
  spec <- ctx$ip$phase_jobs('B',ctx$archives,ctx$inputs_root,TARGET_KEY)[[1]]
  if(!identical(spec$input_md5,SOURCE_INPUT_MD5)) stop('Unexpected input record for ',TARGET_KEY)
  input <- readRDS(ctx$ip$checked_input(spec$input_file,spec$input_md5))
  w <- package_warnings(combined$results_output,ctx$flag$computeDiagnostics)
  flags <- list(row=ctx$ip$fit_flags('B',combined,w,input,ctx$flag),warnings=w)
  rm(combined);invisible(gc())
  list(chains=chains,draws=sum(vapply(files,function(f) f$draws,numeric(1))),estimate=estimate,iv=iv_from_block(b),
    max_cell_rhat=max(b$rhat,na.rm=TRUE),files=do.call(rbind,files),flags=flags,
    checks=c(psi_vs_stored=psi_check,pooled_vs_draws=est_check))
}

# Both scorers' groups for one source estimate of community 5: the
# current-main-recheck original-site groups, and the phase B groups with the
# source's own 95% intervals.
source_scores <- function(src,base,ctx) {
  truth <- base$truth;n <- nrow(truth);S <- ncol(truth);idx <- base$idx
  if(!identical(src$iv$truth,as.vector(truth)[idx])) stop('Source intervals are not of the scored cells')
  full <- function(v,na) {x <- rep(na,n*S);x[idx] <- v;x}
  cells <- ctx$ip$b_make_cells(truth,src$estimate,full(src$iv$covered,NA),base$n0,lower=full(src$iv$lower,NA_real_),
    upper=full(src$iv$upper,NA_real_))
  list(recheck=recheck_groups(truth,src$estimate,base$n0,ctx$sc$error_metrics),phase_b=ctx$ip$b_group_table(cells,ctx$sc))
}

log_step <- function(...) {cat(format(Sys.time(),'%H:%M:%S'),...,'\n');flush.console()}

# ---- Tables -----------------------------------------------------------------------

SOURCE_LABELS <- c(pr11_pooled='pr11 fit, all four chains pooled (as published)',
  near_truth='pr11 fit, chains 1 and 3 (near-truth mode)',mirror='pr11 fit, chains 2 and 4 (mirror mode)',
  extended='extended run of Task 3, 16 chains pooled')

community5_table <- function(scored,sources) {
  rows <- lapply(names(sources),function(s) {
    rg <- scored[[s]]$recheck;pb <- scored[[s]]$phase_b
    pb$group[pb$group=='middle'] <- 'medium';pb <- pb[match(rg$group,pb$group),,drop=FALSE]
    if(any(pb$cells!=rg$n_elements) || max(abs(PP*rg$bias-pb$signed_error),abs(PP*rg$mae-pb$mean_abs_cell_error))>1e-10)
      stop('The two scorers disagree on the cells or errors of ',s)
    data.frame(source=s,description=SOURCE_LABELS[[s]],chains=paste(sources[[s]]$chains,collapse=';'),n_chains=length(sources[[s]]$chains),
      draws=sources[[s]]$draws,group=rg$group,cells=rg$n_elements,truth_mean=rg$truth,estimate_mean=rg$estimate,
      bias=rg$bias,mae=rg$mae,rmse=rg$rmse,signed_error_pp=pb$signed_error,mean_abs_cell_error_pp=pb$mean_abs_cell_error,
      coverage=pb$coverage,max_cell_rhat_original_sites=sources[[s]]$max_cell_rhat,stringsAsFactors=FALSE)
  })
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# The comparisons of current-main-recheck/REPORT.md that use community 5 of
# the high-contamination 300-site design, and the four-sample one beside them.
# Each statement is quoted verbatim from REPORT.md with its bold markup
# removed (' ... ' joins non-adjacent sentences of one paragraph); quantity
# names the number the row recomputes.
COMPARISONS <- list(
  list(id='sites300_mae',section='Do the practical conclusions change?',kind='arm',arm='sites300',group='all',column='mae',
    quantity='overall MAE, 300 sites, high contamination (16.97)',
    statement='Using 300 sites gives 15.20 and 16.97 points.',
    rule='a reported value, not itself a conclusion'),
  list(id='baseline_to_sites300',section='Do the practical conclusions change?',kind='contrast',reference='baseline',alternative='sites300',group='all',
    quantity='MAE reduction from the baseline to 300 sites, high contamination, and communities improving',
    statement='Using 300 sites gives 15.20 and 16.97 points. Every community improves relative to the baseline in each of these comparisons.',
    rule='holds if the mean reduction is positive and all 10 communities improve; reverses if the mean reduction is not positive'),
  list(id='baseline_to_field4',section='Do the practical conclusions change?',kind='contrast',reference='baseline',alternative='field4',group='all',
    quantity='MAE reduction from the baseline to four field samples, high contamination, and communities improving',
    statement=paste('In the two-stage model, four field samples reduce error from 17.97 to 15.26 points under low contamination and from 20.52 to 17.23 under high contamination.',
      '... Every community improves relative to the baseline in each of these comparisons.'),
    rule='as baseline_to_sites300; does not use the 300-site fit of community 5, so unchanged by construction'),
  list(id='field4_to_sites300',section='Do the practical conclusions change?',kind='contrast',reference='field4',alternative='sites300',group='all',
    quantity='MAE reduction from four field samples to 300 sites, high contamination (0.26, -0.90 to 1.42)',
    statement=paste('The additional advantage of 300 sites over four field samples is small: 0.06 points under low contamination and 0.26 under high contamination.',
      'Their paired intervals span zero (-0.98 to 1.11 and -0.90 to 1.42 points). These data do not identify a clear winner between those designs.'),
    rule='holds while the 95% interval spans zero; reverses if it excludes zero (a clear winner); the sign of the mean is reported separately'),
  list(id='sites300_low_bias',section='Do the practical conclusions change?',kind='arm',arm='sites300',group='low',column='bias',
    quantity='mean signed error, true probability below 0.2, 300 sites, high contamination',
    statement='Low probabilities remain too high and high probabilities remain too low.',
    rule='holds if the mean signed error is positive'),
  list(id='sites300_high_bias',section='Do the practical conclusions change?',kind='arm',arm='sites300',group='high',column='bias',
    quantity='mean signed error, true probability above 0.8, 300 sites, high contamination',
    statement='Low probabilities remain too high and high probabilities remain too low.',
    rule='holds if the mean signed error is negative'),
  list(id='sites300_code_change',section='How much did probability accuracy change?',kind='version',arm='sites300',group='all',
    quantity='current minus archived overall MAE, 300 sites, high contamination (0.249, -0.197 to 0.695)',
    statement='Every paired 95% interval for these eleven changes includes zero.',
    rule='holds while the 95% interval includes zero'),
  list(id='two_stage_largest_code_change',section='How much did probability accuracy change?',kind='version_max',group='all',
    arms=c('baseline','field4','sites300','knownU'),scenarios=c('qnear_K6','qfar_K6'),
    quantity='largest absolute current minus archived overall MAE change over the eight two-stage arms (0.249)',
    statement='Across the ten-community summaries, overall mean absolute error changed by at most 0.025 percentage points in binary JSDM and 0.249 points in the two-stage comparisons.',
    rule='a reported value, not itself a conclusion; note names the arm attaining it'))
SCENARIO <- 'qfar_K6'

# The largest absolute change of a set of arm changes, and where it occurs.
largest_change <- function(changes) {
  i <- which.max(abs(changes$delta));list(value=PP*abs(changes$delta[i]),lower=NA_real_,upper=NA_real_,improved=NA_integer_,
    note=paste(changes$scenario[i],changes$arm[i]))
}

comparison_value <- function(cmp,scores,mean_ci) {
  if(cmp$kind=='arm') {
    a <- arm_summary(scores,'current',SCENARIO,cmp$arm,cmp$group)
    return(list(value=PP*a[[cmp$column]],lower=NA_real_,upper=NA_real_,improved=NA_integer_,note=NA_character_))
  }
  if(cmp$kind=='contrast') {
    d <- design_contrast(scores,'current',SCENARIO,cmp$reference,cmp$alternative,cmp$group,mean_ci)
    return(list(value=PP*d$reduction_mae,lower=PP*d$lower,upper=PP*d$upper,improved=as.integer(d$improved),note=NA_character_))
  }
  if(cmp$kind=='version_max') {
    g <- expand.grid(scenario=cmp$scenarios,arm=cmp$arms,stringsAsFactors=FALSE)
    g$delta <- vapply(seq_len(nrow(g)),function(i) version_change(scores,g$scenario[i],g$arm[i],cmp$group,mean_ci)$delta_mae_mean,numeric(1))
    return(largest_change(g))
  }
  v <- version_change(scores,SCENARIO,cmp$arm,cmp$group,mean_ci)
  list(value=PP*v$delta_mae_mean,lower=PP*v$delta_mae_lower,upper=PP*v$delta_mae_upper,improved=NA_integer_,note=NA_character_)
}

published_value <- function(cmp) {
  if(cmp$kind=='arm') {
    s <- read_recheck('chain-sensitivity-summary.csv')
    r <- s[s$family=='design' & s$scenario==SCENARIO & s$arm==cmp$arm & s$group==cmp$group,]
    return(list(value=PP*r[[paste0('selected_',cmp$column)]],lower=NA_real_,upper=NA_real_,improved=NA_integer_,note=NA_character_))
  }
  if(cmp$kind=='contrast') {
    s <- read_recheck('design-paired-summary.csv')
    r <- s[s$version=='current' & s$family=='design' & s$scenario==SCENARIO & s$reference==cmp$reference & s$alternative==cmp$alternative & s$group==cmp$group,]
    return(list(value=PP*r$reduction_mae,lower=PP*r$lower,upper=PP*r$upper,improved=as.integer(r$improved),note=NA_character_))
  }
  s <- read_recheck('paired-summary.csv')
  s <- s[s$comparison=='selected' & s$family=='design' & s$metric=='occupancy_original_sites' & s$group==cmp$group,]
  if(cmp$kind=='version_max') {
    s <- s[s$scenario %in% cmp$scenarios & s$arm %in% cmp$arms,];stopifnot(nrow(s)==length(cmp$scenarios)*length(cmp$arms))
    return(largest_change(data.frame(scenario=s$scenario,arm=s$arm,delta=s$delta_mae_mean)))
  }
  r <- s[s$scenario==SCENARIO & s$arm==cmp$arm,]
  list(value=PP*r$delta_mae_mean,lower=PP*r$delta_mae_lower,upper=PP*r$delta_mae_upper,improved=NA_integer_,note=NA_character_)
}

statement_holds <- function(cmp,v) switch(cmp$id,
  sites300_mae=,two_stage_largest_code_change=NA,
  baseline_to_sites300=,baseline_to_field4=v$value>0 && identical(v$improved,10L),
  field4_to_sites300=v$lower<0 && v$upper>0,
  sites300_low_bias=v$value>0,
  sites300_high_bias=v$value<0,
  sites300_code_change=v$lower<=0 && v$upper>=0,
  stop('No rule for ',cmp$id))

direction_reverses <- function(cmp,v) switch(cmp$id,
  sites300_mae=,two_stage_largest_code_change=NA,
  baseline_to_sites300=,baseline_to_field4=!(v$value>0),
  field4_to_sites300=!(v$lower<0 && v$upper>0),
  sites300_low_bias=!(v$value>0),
  sites300_high_bias=!(v$value<0),
  sites300_code_change=!(v$lower<=0 && v$upper>=0),
  stop('No rule for ',cmp$id))

affected_table <- function(scores_by_source,mean_ci) {
  rows <- list()
  for(cmp in COMPARISONS) {
    pub <- published_value(cmp)
    range <- switch(cmp$kind,arm=PP*sensitivity_range('arm',SCENARIO,cmp$arm,cmp$group,cmp$column),
      contrast=PP*sensitivity_range('contrast',SCENARIO,paste(cmp$reference,cmp$alternative,sep=':'),cmp$group,'reduction_mae'),
      version=,version_max=c(min=NA_real_,max=NA_real_))
    for(s in c('published',names(scores_by_source))) {
      v <- if(s=='published') pub else comparison_value(cmp,scores_by_source[[s]],mean_ci)
      holds <- statement_holds(cmp,v);rev <- direction_reverses(cmp,v)
      rows[[length(rows)+1L]] <- data.frame(comparison=cmp$id,report_section=cmp$section,quantity=cmp$quantity,statement=cmp$statement,
        source=s,unit='percentage points',value=v$value,lower=v$lower,upper=v$upper,improved=v$improved,note=v$note,
        change_from_published=v$value-pub$value,sensitivity_min=range[['min']],sensitivity_max=range[['max']],
        inside_sensitivity_range=if(is.na(range[['min']])) NA else v$value>=range[['min']]-1e-12 && v$value<=range[['max']]+1e-12,
        statement_holds=holds,reverses=rev,rule=cmp$rule,stringsAsFactors=FALSE)
    }
  }
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

# The published gate beside the gate with the converged control (converged)
# and, for criterion 4, with every arm's own flagged community 5 fit also
# counted as unflagged (arms_unflagged).
phase_b_table <- function(published,converged,arms_unflagged) {
  key <- c('sd','criterion','component','stratum')
  a <- published$detail;b <- converged$detail;u <- arms_unflagged$detail
  tag <- function(d) paste(d$sd,d$criterion,d$component,d$stratum)
  if(!identical(tag(a),tag(b)) || !identical(tag(a),tag(u))) stop('Gate details do not pair up')
  if(any(a$value_new!=b$value_new,na.rm=TRUE)) stop('A new arm value changed')
  if(any(u$criterion!='c4' & (u$pass!=b$pass | u$value_new!=b$value_new),na.rm=TRUE)) stop('Counting the arms unflagged changed more than criterion 4')
  x <- data.frame(a[key],communities=a$communities,needed=a$needed,improved_published=a$improved,improved_converged=b$improved,
    value_new=a$value_new,value_new_arms_unflagged=u$value_new,value_control_published=a$value_control,value_control_converged=b$value_control,
    statistic_published=a$statistic,statistic_converged=b$statistic,threshold=a$threshold,
    pass_published=a$pass,pass_converged=b$pass,pass_changed=a$pass!=b$pass,pass_arms_unflagged=u$pass,result_published=NA_character_,
    result_converged=NA_character_,result_arms_unflagged=NA_character_,stringsAsFactors=FALSE)
  gate <- lapply(seq_len(nrow(published$rows)),function(i) {
    p <- published$rows[i,];q <- converged$rows[i,];z <- arms_unflagged$rows[i,];stopifnot(p$sd==q$sd,p$sd==z$sd)
    do.call(rbind,lapply(c('c1','c2','c3','c4','result'),function(k) {
      pass <- function(r) if(k=='result') NA else r[[paste0(k,'_pass')]]
      res <- function(r) if(k=='result') r$result else NA_character_
      data.frame(sd=p$sd,criterion=if(k=='result') 'gate' else k,component=if(k=='result') 'result' else 'all components',stratum='both',
        communities=NA_integer_,needed=NA_integer_,improved_published=NA_integer_,improved_converged=NA_integer_,value_new=NA_real_,
        value_new_arms_unflagged=NA_real_,value_control_published=NA_real_,value_control_converged=NA_real_,statistic_published=NA_real_,
        statistic_converged=NA_real_,threshold=NA_real_,pass_published=pass(p),pass_converged=pass(q),
        pass_changed=if(k=='result') p$result!=q$result else pass(p)!=pass(q),pass_arms_unflagged=pass(z),
        result_published=res(p),result_converged=res(q),result_arms_unflagged=res(z),stringsAsFactors=FALSE)
    }))
  })
  x <- rbind(x,do.call(rbind,gate));rownames(x) <- NULL
  x[order(x$sd,match(x$criterion,c('c1','c2','c3','c4','gate'))),,drop=FALSE]
}

write_csv <- function(x,dir,name) {
  path <- file.path(dir,name);tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE);stopifnot(file.rename(tmp,path));path
}

# ---- Main -----------------------------------------------------------------------------

impact_main <- function(args) {
  opt <- function(name,default=NULL) {
    hit <- grep(paste0('^--',name,'='),args,value=TRUE)
    if(length(hit)) sub(paste0('^--',name,'='),'',hit[length(hit)]) else default
  }
  bad <- setdiff(sub('=.*$','',sub('^--','',args)),c('out','cache'))
  if(length(bad)) stop('Unknown option: --',bad[1])
  out <- opt('out',IMPACT_OUT);cache <- opt('cache')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  ctx <- impact_setup()
  log_step('archived definitions loaded (md5s match)')
  heavy <- if(!is.null(cache) && file.exists(cache)) {log_step('reading cache',cache);readRDS(cache)} else {
    p <- score_pr11(ctx);e <- score_extended(ctx,p)
    h <- list(pr11=p,extended=e);if(!is.null(cache)) saveRDS(h,cache);h
  }
  p <- heavy$pr11;e <- heavy$extended
  sources <- c(p$sources,list(extended=e))
  scored <- lapply(sources,source_scores,base=p,ctx=ctx)
  checks <- list()
  add_check <- function(name,value,tolerance,pass=value<tolerance) checks[[length(checks)+1L]] <<-
    data.frame(check=name,max_abs_difference=value,tolerance=tolerance,pass=pass,stringsAsFactors=FALSE)

  # Step 2: reproduction first, then substitution.
  scores <- read_recheck('community-scores.csv')
  published5 <- scores[scores$key==TARGET_KEY & scores$version=='current',,drop=FALSE]
  m <- merge(published5,p$recheck_rows,by=c('metric','group'),suffixes=c('_published','_fit'))
  if(nrow(m)!=nrow(published5)) stop('The scorer groups do not pair with the published rows')
  num <- c('n_elements','truth','estimate','bias','mae','rmse','rhat','ess_mean','mcse','chain_gap')
  add_check('community 5 current rows of community-scores.csv (every metric and group, every numeric column) from the pr11 fit',
    max(vapply(num,function(k) max(abs(m[[paste0(k,'_published')]]-m[[paste0(k,'_fit')]])),numeric(1))),1e-10)
  pub_band <- published5[published5$metric=='occupancy_original_sites',];rg <- scored$pr11_pooled$recheck
  pub_band <- pub_band[match(rg$group,pub_band$group),]
  add_check('community 5 original-site band errors (recheck_groups on the reproduced estimate) against community-scores.csv',
    max(abs(c(rg$truth-pub_band$truth,rg$estimate-pub_band$estimate,rg$bias-pub_band$bias,rg$mae-pub_band$mae,rg$rmse-pub_band$rmse))),1e-10)
  ch <- read_recheck('flagged-chain-scores.csv');ch <- ch[ch$key==TARGET_KEY,];gap <- 0
  for(k in 0:ncol(p$chain_means)) {
    est <- if(k==0L) p$sources$pr11_pooled$estimate else matrix(p$chain_means[,k],nrow(p$truth))
    g <- recheck_groups(p$truth,est,p$n0,ctx$sc$error_metrics);z <- ch[ch$chain==k,];z <- z[match(g$group,z$group),]
    gap <- max(gap,abs(g$bias-z$bias),abs(g$mae-z$mae))
  }
  add_check('community 5 per-chain and pooled band errors against flagged-chain-scores.csv',gap,1e-10)
  fd <- read_recheck('flagged-diagnostics.csv')
  fd <- fd[fd$key==TARGET_KEY & fd$metric=='original_probability',]
  add_check('largest original-site cell Rhat of the pr11 fit (score_draw_block) against flagged-diagnostics.csv',
    abs(p$sources$pr11_pooled$max_cell_rhat-max(fd$rhat)),1e-10)
  add_check('per-chain means pooled against the design scorer estimate',p$checks[['pooled_vs_scorer']],1e-12)
  add_check('pooled four-chain intervals against b_cells',p$checks[['intervals_vs_b_cells']],1e-12)
  add_check('extended chains: reconstructed means against each fit\'s stored psi_output',e$checks[['psi_vs_stored']],1e-10)
  add_check('extended run: pooled estimate against its pooled draws',e$checks[['pooled_vs_draws']],1e-12)

  # The pr11 flag rule: validated on the pr11 fit, then applied to the 16 bound chains.
  fp <- p$flags$row;cf <- utils::read.csv(file.path(IP_RESULTS,'selection-control-flags-B.csv'),stringsAsFactors=FALSE)
  cf <- cf[cf$key==TARGET_KEY,,drop=FALSE];ds <- read_recheck('diagnostic-summary.csv');ds <- ds[ds$key==TARGET_KEY,,drop=FALSE]
  gap <- max(abs(c(fp$max_group_rhat-cf$max_group_rhat,fp$max_element_rhat-cf$max_element_rhat)))
  add_check('pr11 flag rule on the pr11 fit against the intercept-prior control flags (selection-control-flags-B.csv): Rhats, warnings, flag and reasons',
    gap,1e-12,pass=gap<1e-12 && fp$warnings==cf$warnings && fp$unresolved_rhat==cf$unresolved_rhat && identical(fp$flagged,cf$flagged) &&
      identical(fp$reasons,cf$reasons))
  gap <- max(abs(c(fp$max_group_rhat-ds$max_group_rhat,fp$max_element_rhat-ds$max_element_rhat)))
  add_check('pr11 flag rule on the pr11 fit against current-main-recheck diagnostic-summary.csv: Rhats, warnings, unresolved',
    gap,1e-12,pass=gap<1e-12 && fp$warnings==ds$warnings && fp$unresolved_rhat==ds$unresolved_rhat)
  add_check('package diagnostics recomputed on the pr11 fit give its saved fitting warnings',NA_real_,NA_real_,
    pass=identical(sort(p$flags$warnings_recomputed),sort(p$flags$warnings_saved)) && length(p$flags$warnings_saved)==7L)
  fe <- e$flags$row
  flag_rows <- rbind(
    data.frame(fit='pr11 fit (published record)',source='results/selection-control-flags-B.csv (intercept-prior)',chains=4L,
      draws_per_chain=12000L,warnings=as.integer(cf$warnings),warning_messages=NA_character_,max_group_rhat=cf$max_group_rhat,
      max_element_rhat=cf$max_element_rhat,unresolved_rhat=as.integer(cf$unresolved_rhat),flagged=cf$flagged,reasons=cf$reasons,stringsAsFactors=FALSE),
    data.frame(fit='pr11 fit (recomputed)',source='fit_flags on the saved fit; warnings recomputed by computeDiagnostics',chains=4L,
      draws_per_chain=12000L,warnings=as.integer(fp$warnings),warning_messages=paste(sort(p$flags$warnings_recomputed),collapse=' | '),
      max_group_rhat=fp$max_group_rhat,max_element_rhat=fp$max_element_rhat,unresolved_rhat=as.integer(fp$unresolved_rhat),flagged=fp$flagged,
      reasons=fp$reasons,stringsAsFactors=FALSE),
    data.frame(fit='extended run, 16 chains bound into one fit',source='fit_flags on the bound chains; warnings by computeDiagnostics on the bound output',
      chains=length(e$chains),draws_per_chain=as.integer(e$files$draws[1]),warnings=as.integer(fe$warnings),
      warning_messages=paste(e$flags$warnings,collapse=' | '),max_group_rhat=fe$max_group_rhat,max_element_rhat=fe$max_element_rhat,
      unresolved_rhat=as.integer(fe$unresolved_rhat),flagged=fe$flagged,reasons=fe$reasons,stringsAsFactors=FALSE))
  mean_ci <- ctx$recheck$mean_ci
  subs <- lapply(scored,function(s) substitute_scores(scores,TARGET_KEY,'current',s$recheck))
  for(cmp in COMPARISONS) {
    a <- published_value(cmp);b <- comparison_value(cmp,scores,mean_ci);c <- comparison_value(cmp,subs$pr11_pooled,mean_ci)
    add_check(paste('published',cmp$id,'recomputed from community-scores.csv'),max(abs(unlist(a[1:3])-unlist(b[1:3])),na.rm=TRUE),1e-8,
      pass=max(abs(unlist(a[1:3])-unlist(b[1:3])),na.rm=TRUE)<1e-8 && identical(a$improved,b$improved))
    add_check(paste('published',cmp$id,'recomputed with community 5 reproduced from the fit'),max(abs(unlist(a[1:3])-unlist(c[1:3])),na.rm=TRUE),1e-8,
      pass=max(abs(unlist(a[1:3])-unlist(c[1:3])),na.rm=TRUE)<1e-8 && identical(a$improved,c$improved))
  }
  # Every other row of the substituted score tables is the published one.
  for(s in names(subs)) {
    hit <- subs[[s]]$key==TARGET_KEY & subs[[s]]$version=='current' & subs[[s]]$metric=='occupancy_original_sites'
    add_check(paste('substitution',s,'leaves every other community-score row unchanged'),NA_real_,NA_real_,
      pass=sum(hit)==4L && identical(subs[[s]][!hit,],scores[!hit,]))
  }

  # Step 3.
  groups <- phase_b_groups();conv <- phase_b_convergence()
  gate <- phase_b_gate(groups,conv,ctx$ip)
  pubd <- utils::read.csv(file.path(IP_RESULTS,'gate-detail-B.csv'),stringsAsFactors=FALSE)
  nd <- c('communities','needed','improved','value_new','value_control','statistic','threshold','missing_new')
  add_check('phase B gate-detail-B.csv recomputed from the archived per-community tables',
    max(abs(as.matrix(gate$detail[nd])-as.matrix(pubd[nd])),na.rm=TRUE),1e-12,
    pass=max(abs(as.matrix(gate$detail[nd])-as.matrix(pubd[nd])),na.rm=TRUE)<1e-12 && identical(gate$detail$pass,pubd$pass) &&
      identical(gate$rows$result,utils::read.csv(file.path(IP_RESULTS,'gate-B.csv'),stringsAsFactors=FALSE)$result))
  ctl <- groups[groups$sd==1 & groups$key==TARGET_KEY,];pb <- scored$pr11_pooled$phase_b;ctl <- ctl[match(pb$group,ctl$group),]
  cols <- c('cells','truth_mean','estimate_mean','signed_error','abs_signed_error','mean_abs_cell_error','coverage')
  add_check('phase B control groups of community 5 (band-error.csv, coverage.csv) from the pr11 fit',
    max(abs(as.matrix(pb[cols])-as.matrix(ctl[cols]))),1e-10)
  g2 <- substitute_phase_b(groups,TARGET_KEY,scored$extended$phase_b,label='extended run (Task 3), 16 chains')
  # The control counts as unflagged only if the pr11 rule leaves the extended run unflagged.
  conv2 <- converged_convergence(conv,'qfar',fewer=if(isTRUE(fe$flagged)) 0L else 1L)
  gate2 <- phase_b_gate(g2,conv2,ctx$ip)
  gate3 <- phase_b_gate(g2,arms_unflagged_convergence(conv2,TARGET_KEY),ctx$ip)
  pbt <- phase_b_table(gate,gate2,gate3)

  write_csv(other_flagged_fits(),out,'other-flagged-fits.csv')
  write_csv(community5_table(scored,sources),out,'community5-scores.csv')
  write_csv(affected_table(subs[c('pr11_pooled','near_truth','mirror','extended')],mean_ci),out,'affected-comparisons.csv')
  write_csv(pbt,out,'intercept-prior-phase-b.csv')
  write_csv(flag_rows,out,'flag-rule.csv')
  chk <- do.call(rbind,checks)
  write_csv(chk,out,'reproduction.csv')
  hashed <- c(file.path(IMPACT_DIR,c('impact.R','test-impact.R','anatomy.R')),file.path(IP_DIR,IP_SCRIPTS),
    file.path(RECHECK_DIR,c('helpers.R','summarise.R','chain-sensitivity.R')),
    file.path(REPO,c(ctx$ip$B_METRIC_FILES)))
  hashed <- c(unique(hashed),ctx$flag$diagnostics_file)
  label <- ifelse(startsWith(hashed,paste0(REPO,'/')),sub(paste0('^',REPO,'/'),'',hashed),
    paste0('archives/',sub(paste0('^',ARCHIVES,'/'),'',hashed)))
  write_csv(data.frame(file=label,md5=unname(tools::md5sum(hashed))),out,'source-hashes.csv')
  rel <- function(f) sub(paste0('^',ARCHIVES,'/'),'',f)
  prov <- rbind(data.frame(source='pr11_pooled',chain=NA_integer_,fit_file=rel(p$fit_file),fit_md5=p$fit_md5,draws=p$sources$pr11_pooled$draws,
      input_md5=p$input_md5,stringsAsFactors=FALSE),
    data.frame(source='extended',chain=e$files$chain,fit_file=e$files$fit_file,fit_md5=e$files$fit_md5,draws=e$files$draws,input_md5=SOURCE_INPUT_MD5,
      stringsAsFactors=FALSE))
  write_csv(prov,out,'provenance.csv')
  failed <- chk$check[!chk$pass]
  if(length(failed)) {cat('FAILED checks:\n',paste(' ',failed,collapse='\n'),'\n');return(1L)}
  log_step('wrote',out,'; all',nrow(chk),'checks pass')
  0L
}

if(sys.nframe()==0L) {
  status <- tryCatch(impact_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
