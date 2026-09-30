# Frozen convergence-flag rules for the occupancy-intercept prior study and the
# selection bookkeeping built on them. Every flag is computed from saved draws,
# identically for the control and new arms. No occupancy error is computed or
# read here, and none is used. Generating truth enters in two ways: in every
# phase it defines which cells form a probability band or prevalence group, as
# in the archived scorers; in A2 the archived robust rule also screens the
# field_projection traces, which project each field draw onto the centred true
# spatial field, and robust_field_draws also returns a field-recovery summary
# against the true field, which spatial_tables discards unread. Not sourced or
# hashed by run.R. Source jobs.R and verify-helpers.R first. See AMENDMENT-1.md.
#
# Phases A1 and B use the current-main recheck rule (select.R:17-19 there, with
# diagnostics combined as its run.R:79-81): flag a fit with any fitting warning,
# maximum group or element Rhat above 1.05, or any non-finite or nonpositive
# Rhat. Group and element Rhat are those of the archived scorer for the family
# (jsdm-sample-size-recheck/helpers.R for A1; the nonspatial-bias-recheck scorer
# as extended by nonspatial-design-recheck/design_helpers.R and
# current-main-recheck/helpers.R for B), reproduced here without their error
# columns. Phase A2 uses the spatial-amplitude robust-v1 rule (robust.R:100-110
# there): the spatial-targeted diagnostic rule without native warnings, with the
# amplitude diagnosed by rank/quantile methods, plus the rank/quantile screens of
# the spatial summary traces and pointwise fields.

FLAG_COLUMNS <- c('rule','warnings','max_group_rhat','max_element_rhat','unresolved_rhat',
  'spatial_trace_flags','spatial_field_flags','flagged','reasons')
ITEM_COLUMNS <- c('role','phase','sd','key','schedule','fit','fit_md5')
CONTROL_PROVENANCE <- 'dev/simstudy/occupancy-intercept-prior/results/control-provenance.csv'
PR11_SELECT <- 'dev/simstudy/current-main-recheck/select.R'

# Archived files the flag code reads definitions from or reproduces, each with
# the record of its md5 made when it was last used to select or score fits.
SCORER_FILES <- c(
  'dev/simstudy/current-main-recheck/select.R'='repo-source-hashes',
  'dev/simstudy/current-main-recheck/helpers.R'='pr11-settings',
  'dev/simstudy/jsdm-sample-size-recheck/helpers.R'='pr11-settings',
  'dev/simstudy/nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R'='pr11-settings',
  'dev/simstudy/nonspatial-design-recheck/design_helpers.R'='pr11-settings',
  'dev/simstudy/spatial-targeted-recheck/score.R'='amplitude-settings',
  'dev/simstudy/spatial-targeted-recheck/analysis.R'='amplitude-settings',
  'dev/simstudy/spatial-amplitude-prior/robust.R'='robust-provenance',
  'dev/simstudy/spatial-amplitude-prior/rescore.R'='robust-provenance')

recorded_scorer_md5 <- function(file,record,repo,archives) {
  hashes <- switch(record,
    'repo-source-hashes'={x <- utils::read.csv(file.path(repo,'dev/simstudy/current-main-recheck/results/source-hashes.csv'),
      stringsAsFactors=FALSE);stats::setNames(x$md5,x$path)},
    'pr11-settings'=readRDS(file.path(archive_dir(archives,'pr11'),'initial/settings.rds'))$source_hashes,
    'amplitude-settings'=readRDS(file.path(archive_dir(archives,'amplitude'),'inverse_gamma/long/settings.rds'))$script_hashes,
    'robust-provenance'=readRDS(file.path(archive_dir(archives,'amplitude'),'robust-v1/summary-binary-final/provenance.rds'))$script_hashes,
    stop('Unknown md5 record: ',record))
  recorded_hash(hashes,file)
}

checked_scorer_file <- function(file,repo,archives) {
  path <- file.path(repo,file);recorded <- recorded_scorer_md5(file,SCORER_FILES[[file]],repo,archives)
  found <- unname(tools::md5sum(path))
  if(!identical(found,recorded)) stop('Archived scorer changed: ',file,' has md5 ',found,', recorded ',recorded)
  path
}

is_assignment <- function(e,lhs) is.call(e) && length(e)==3L && identical(e[[1]],as.name('<-')) && identical(e[[2]],lhs)

# Top-level `name <- value` statements of a file, evaluated in env. Every eval
# in this file runs code parsed only from a tracked archived research script
# whose md5 was first checked against its record (checked_scorer_file,
# pr11_flag_expression), never text from fits, inputs or other data.
extract_definitions <- function(path,names,env) {
  exprs <- parse(path,keep.source=FALSE)
  for(name in names) {
    hits <- Filter(function(e) is_assignment(e,as.name(name)),exprs)
    if(length(hits)!=1L) stop('Expected one definition of ',name,' in ',path,', found ',length(hits))
    eval(hits[[1]],envir=env)
  }
  invisible(env)
}

# Right-hand sides of every assignment `lhs <- value` nested anywhere in exprs.
nested_assignments <- function(exprs,lhs) {
  hits <- list()
  visit <- function(e) {
    if(!is.call(e)) return(invisible())
    if(is_assignment(e,lhs)) hits[[length(hits)+1L]] <<- e[[3]]
    for(i in seq_along(e)) if(!is_missing_arg(e,i)) visit(e[[i]])
  }
  for(e in exprs) visit(e)
  hits
}

# The archived flag expression of current-main-recheck/select.R, verbatim.
pr11_flag_expression <- function(repo) {
  path <- file.path(repo,PR11_SELECT)
  x <- utils::read.csv(file.path(repo,'dev/simstudy/current-main-recheck/results/source-hashes.csv'),stringsAsFactors=FALSE)
  recorded <- recorded_hash(stats::setNames(x$md5,x$path),PR11_SELECT)
  if(!identical(unname(tools::md5sum(path)),recorded)) stop('Archived selection rule changed: ',path)
  hits <- nested_assignments(parse(path,keep.source=FALSE),as.name('flag'))
  if(length(hits)!=1L) stop('Expected one flag rule in ',path,', found ',length(hits))
  hits[[1]]
}

load_flag_scorers <- function(repo,archives) {
  paths <- vapply(names(SCORER_FILES),checked_scorer_file,character(1),repo=repo,archives=archives)
  pr11 <- new.env(parent=globalenv())
  extract_definitions(paths[['dev/simstudy/jsdm-sample-size-recheck/helpers.R']],'trace_stats',pr11)
  extract_definitions(paths[['dev/simstudy/current-main-recheck/helpers.R']],'safe_max',pr11)
  # current-main-recheck/helpers.R:19-24 replaces the design scorer's trace
  # summary by one that compares every chain; take that definition verbatim.
  loader <- Filter(function(e) is_assignment(e,as.name('load_scoring')),
    parse(paths[['dev/simstudy/current-main-recheck/helpers.R']],keep.source=FALSE))
  override <- nested_assignments(loader,quote(design$trace_summary))
  if(length(override)!=1L) stop('Expected one design trace summary override in current-main-recheck/helpers.R')
  trace_summary <- eval(override[[1]],pr11)
  spatial <- new.env(parent=globalenv())
  extract_definitions(paths[['dev/simstudy/spatial-targeted-recheck/score.R']],
    c('trace_diagnostics','independent_bases','reconstruct_spatial_draws','spatial_groups'),spatial)
  extract_definitions(paths[['dev/simstudy/spatial-targeted-recheck/analysis.R']],'diagnostic_reasons',spatial)
  extract_definitions(paths[['dev/simstudy/spatial-amplitude-prior/robust.R']],
    c('robust_trace_diagnostics','robust_flag_rows','robust_field_draws','robust_reasons'),spatial)
  list(flag_expr=pr11_flag_expression(repo),trace_stats=pr11$trace_stats,trace_summary=trace_summary,
    safe_max=pr11$safe_max,spatial=spatial,hashes=c(stats::setNames(unname(tools::md5sum(paths)),names(paths))))
}

# current-main-recheck/run.R:79-81.
pr11_diagnostics <- function(rhats,safe_max) {
  all_rhat <- c(rhats$groups,rhats$elements,rhats$additional)
  list(max_group_rhat=safe_max(rhats$groups),max_element_rhat=safe_max(c(rhats$elements,rhats$additional)),
    unresolved_rhat=sum(!is.finite(all_rhat) | all_rhat<=0))
}

# Evaluates the archived expression on x$warnings and d, and states why.
pr11_rule <- function(warnings,d,flag_expr) {
  env <- new.env(parent=baseenv());env$x <- list(warnings=warnings);env$d <- d
  flagged <- eval(flag_expr,env)
  reasons <- character();w <- length(warnings)
  if(w) reasons <- c(reasons,sprintf('%d fitting warning%s',w,if(w==1L) '' else 's'))
  for(nm in c('group','element')) {
    v <- d[[paste0('max_',nm,'_rhat')]]
    if(!is.finite(v)) reasons <- c(reasons,sprintf('max %s Rhat unavailable',nm)) else
      if(v>1.05) reasons <- c(reasons,sprintf('max %s Rhat %.4f > 1.05',nm,v))
  }
  if(d$unresolved_rhat>0L) reasons <- c(reasons,sprintf('%d non-finite or nonpositive Rhat',d$unresolved_rhat))
  if(!is.logical(flagged) || length(flagged)!=1L || is.na(flagged) || flagged!=(length(reasons)>0L))
    stop('Flag reasons disagree with the archived rule')
  list(flagged=flagged,reasons=reasons)
}

rhat_of <- function(f) function(x) unname(f(x)['rhat'])
block_rhats <- function(arr,ni,nc,rh) {
  dd <- dim(arr);np <- prod(head(dd,-2L));flat <- array(arr,c(np,ni,nc))
  vapply(seq_len(np),function(k) rh(flat[k,,]),numeric(1))
}

# Group, element and additional Rhat of the A1 (binary JSDM) scorer,
# jsdm-sample-size-recheck/helpers.R:53-99 with current-main-recheck/helpers.R:94-95.
a1_rhats <- function(fit,input,trace_stats) {
  stopifnot(fit$infos$model=='binary',fit$infos$ps==0,
    is.null(fit$results_output$p_output),is.null(fit$results_output$q_output))
  j <- fit$results_output$jsdm_output
  n <- input$n;S <- ncol(input$data$OTU);ni <- dim(j$B0_output)[2];nc <- dim(j$B0_output)[3]
  X <- fit$X_psi;tp <- input$truth
  true_eta <- sweep(X%*%tp$B+tp$U%*%tp$L,2,tp$B0,'+')
  stopifnot(max(abs(true_eta-tp$eta))<1e-10)
  truth <- plogis(true_eta)
  original <- unlist(lapply(seq_len(S),function(s)(s-1L)*n+seq_len(100L)))
  masks <- list()
  for(scope in c('original100','allsites')) {
    base <- if(scope=='original100') original else seq_len(n*S)
    for(band in c('all','low','middle','high')) {
      p <- truth[base]
      pick <- switch(band,all=rep(TRUE,length(p)),low=p<.2,middle=p>=.2&p<=.8,high=p>.8)
      masks[[paste(scope,band,sep=':')]] <- base[pick]
    }
  }
  traces <- array(0,c(length(masks),ni,nc))
  original_draws <- array(0,c(100L*S,ni,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    eta <- sweep(X%*%j$B_output[,,it,ch]+j$U_output[,,it,ch]%*%j$L_output[,,it,ch],2,j$B0_output[,it,ch],'+')
    p <- plogis(eta)
    original_draws[,it,ch] <- p[original]
    traces[,it,ch] <- vapply(masks,function(ix) mean(p[ix]),numeric(1))
  }
  rh <- rhat_of(trace_stats)
  list(groups=vapply(seq_along(masks),function(k) rh(traces[k,,]),numeric(1)),
    elements=c(block_rhats(j$B0_output,ni,nc,rh),block_rhats(j$B_output,ni,nc,rh),block_rhats(original_draws,ni,nc,rh)),
    additional=rh(j$sigmah_output))
}

# Group, element and additional Rhat of the B (two-stage design) scorer:
# nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R:59-134, with the
# original-site groups of nonspatial-design-recheck/design_helpers.R:119-129
# and the additional diagnostics and trace summary of
# current-main-recheck/helpers.R:19-41.
b_rhats <- function(fit,input,trace_summary) {
  sim <- input$sim;scenario <- input$scenario;tp <- sim$true_params;jp <- tp$jsdmParams_true
  ro <- fit$results_output;jo <- ro$jsdm_output
  S <- scenario$S;n <- scenario$n;ni <- dim(jo$B0_output)[2];nc <- dim(jo$B0_output)[3]
  stopifnot(fit$infos$ps==0,fit$infos$model==scenario$model,identical(fit$infos$speciesNames,colnames(sim$data_list$OTU)))
  info <- sim$data_list$info
  site_info <- info[!duplicated(info$Site),,drop=FALSE];site_info <- site_info[order(site_info$Site),,drop=FALSE]
  rawx <- as.matrix(site_info[paste0('X_psi.EnvCov.',1:scenario$ncov_psi)])
  stopifnot(max(abs(unname(scale(rawx))-unname(fit$X_psi)))<1e-10)
  true_eta <- sweep(fit$X_psi %*% jp$B + jp$U %*% jp$L,2,jp$B0,'+')
  stopifnot(max(abs(true_eta-jp$eta))<1e-10)
  true_psi <- plogis(true_eta)
  if(scenario$model=='two_stage') {
    sample_info <- info[!duplicated(info[c('Site','Sample')]),,drop=FALSE]
    sample_info <- sample_info[order(sample_info$Site,sample_info$Sample),,drop=FALSE]
  } else sample_info <- info[order(info$Site,seq_len(nrow(info))),,drop=FALSE]
  rawt <- sample_info$X_theta
  stopifnot(max(abs(unname(cbind(1,scale(rawt)))-unname(fit$X_theta)))<1e-10)
  true_bt <- tp$beta_theta_true
  true_bt[1,] <- true_bt[1,]+mean(rawt)*true_bt[2,];true_bt[2,] <- sd(rawt)*true_bt[2,]
  blocks <- list(B0=list(jo$B0_output,length(jp$B0)),collection_intercept=list(ro$beta_theta_output[1,,,,drop=FALSE],ncol(true_bt)),
    collection_slope=list(ro$beta_theta_output[2,,,,drop=FALSE],ncol(true_bt)),theta0=list(ro$theta0_output,length(input$truth$params$theta0)))
  if(scenario$model=='two_stage') blocks <- c(blocks,list(p=list(ro$p_output,length(tp$p_true)),
    q=list(ro$q_output,length(tp$q_true)),q_nominal=list(ro$q_output,length(tp$q_true))))
  rh <- rhat_of(trace_summary)
  groups <- numeric();elements <- numeric()
  for(nm in names(blocks)) {
    len <- blocks[[nm]][[2]];x <- array(blocks[[nm]][[1]],c(len,ni,nc))
    elements <- c(elements,vapply(seq_len(len),function(i) rh(x[i,,]),numeric(1)))
    groups <- c(groups,rh(apply(x,c(2,3),mean)))
    if(nm=='collection_slope') for(sg in c(-1,0,1)) {
      ix <- which(sign(true_bt[2,])==sg)
      if(length(ix)) groups <- c(groups,rh(apply(x[ix,,,drop=FALSE],c(2,3),mean)))
    }
  }
  # Occupancy groups: all fitted sites, the three truth bands, the original
  # sites and their bands, and the prevalence groups (empty groups skipped).
  tpv <- as.vector(true_psi)
  n0 <- input$design$original_sites
  base_idx <- unlist(lapply(seq_len(S),function(j)((j-1L)*n)+seq_len(n0)))
  base_truth <- tpv[base_idx]
  band <- function(p,b) switch(b,low=p<.2,medium=p>=.2 & p<=.8,high=p>.8)
  masks <- list(all=seq_along(tpv))
  for(b in c('low','medium','high')) masks[[b]] <- which(band(tpv,b))
  masks$original_all <- base_idx
  for(b in c('low','medium','high')) masks[[paste0('original_',b)]] <- base_idx[which(band(base_truth,b))]
  prevalence <- colMeans(true_psi)
  for(g in c('rare_below_20pct','common_20pct_or_more')) {
    species <- if(g=='rare_below_20pct') which(prevalence<.2) else which(prevalence>=.2)
    masks[[g]] <- unlist(lapply(species,function(s) ((s-1L)*n+1L):(s*n)))
  }
  masks <- masks[lengths(masks)>0L]
  psi_traces <- array(NA_real_,c(length(masks),ni,nc));theta_trace <- matrix(NA_real_,ni,nc)
  original_draws <- array(NA_real_,c(length(base_idx),ni,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    eta <- sweep(fit$X_psi%*%jo$B_output[,,it,ch]+jo$U_output[,,it,ch]%*%jo$L_output[,,it,ch],2,jo$B0_output[,it,ch],'+')
    psi <- as.vector(plogis(eta))
    psi_traces[,it,ch] <- vapply(masks,function(ix) mean(psi[ix]),numeric(1))
    theta_trace[it,ch] <- mean(as.vector(plogis(fit$X_theta%*%ro$beta_theta_output[,,it,ch])))
    original_draws[,it,ch] <- psi[base_idx]
  }
  groups <- c(groups,vapply(seq_along(masks),function(k) rh(psi_traces[k,,]),numeric(1)),rh(theta_trace))
  list(groups=groups,elements=elements,additional=c(block_rhats(jo$B_output,ni,nc,rh),
    block_rhats(original_draws,ni,nc,rh),block_rhats(jo$sigmah_output,ni,nc,rh)))
}

# Diagnostic tables of the A2 (spatial binary) scorer without error columns:
# spatial-targeted-recheck/score.R:79-151 (occupancy groups, cells, blocks and
# species), with the rank/quantile amplitude substitution of
# spatial-amplitude-prior/rescore.R:42-53 and the spatial screens of
# spatial-amplitude-prior/robust.R:38-86. Those screens use the true spatial
# field for the field_projection traces; the field-recovery summary that
# robust_field_draws also returns is discarded unread.
spatial_tables <- function(fit,input,sp) {
  t <- input$truth;js <- fit$results_output$jsdm_output
  n <- input$settings$n;S <- input$settings$S;ni <- dim(js$B0_output)[2];nc <- dim(js$B0_output)[3]
  stopifnot(fit$infos$ps==100L,fit$infos$n_factors==0L,identical(fit$infos$speciesNames,colnames(input$data$binary$OTU)),
    fit$infos$model=='binary',max(abs(fit$Xs-t$Xs))<1e-12,max(abs(fit$X_psi-t$X))<1e-12)
  x <- sp$reconstruct_spatial_draws(fit)$probability
  td <- sp$trace_diagnostics
  rows <- function(arr,len,bni,bnc) as.data.frame(t(vapply(seq_len(len),function(i) td(matrix(arr[i,,],bni,bnc)),numeric(6))))
  groups <- list(data.frame(metric='occupancy',group='all',as.list(td(apply(x,c(2,3),mean)))))
  masks <- sp$spatial_groups(t$psi,t$target_prevalence)
  for(g in setdiff(names(masks),'all')) {
    idx <- which(masks[[g]]);if(!length(idx)) next
    groups[[length(groups)+1L]] <- data.frame(metric='occupancy',group=g,as.list(td(apply(x[idx,,,drop=FALSE],c(2,3),mean))))
  }
  elements <- list(data.frame(metric='occupancy',element=seq_len(n*S),rows(x,n*S,ni,nc)))
  blocks <- list(intercept=list(js$B0_output,length(t$B0)),environment_slope=list(js$B_output,length(t$B)),
    range=list(matrix(fit$infos$l_s_grid[js$idx_ls_output],ni,nc),length(t$range)),
    spatial_sd=list(js$sigmabs_output,length(t$field_sd)))
  for(nm in names(blocks)) {
    draws <- blocks[[nm]][[1]];len <- blocks[[nm]][[2]]
    bni <- dim(draws)[length(dim(draws))-1L];bnc <- tail(dim(draws),1);xb <- array(draws,c(len,bni,bnc))
    elements[[length(elements)+1L]] <- data.frame(metric=nm,element=seq_len(len),rows(xb,len,bni,bnc))
    groups[[length(groups)+1L]] <- data.frame(metric=nm,group='all',as.list(td(apply(xb,c(2,3),mean))))
  }
  species <- do.call(rbind,lapply(seq_len(S),function(s) {
    idx <- ((s-1L)*n+1L):(s*n);data.frame(species=s,as.list(td(apply(x[idx,,,drop=FALSE],c(2,3),mean))))
  }))
  spatial <- sp$robust_field_draws(js$Bs_output,js$idx_ls_output,sp$independent_bases(fit),t$field,js$sigmabs_output)
  diag_cols <- c('rhat','ess_bulk','ess_median','ess_q025','ess_q975','mcse_median','mcse_q025','mcse_q975',
    'chain_median_gap','constant_chain')
  out <- list(groups=do.call(rbind,groups),species=species,elements=do.call(rbind,elements),
    spatial=list(diagnostics=spatial$diagnostics[c('quantity',diag_cols)],
      field_diagnostics=spatial$field_diagnostics[c('quantity','site','species',diag_cols)]))
  amplitude <- out$spatial$diagnostics[1,]
  stopifnot(identical(amplitude$quantity,'amplitude'))
  for(tb in c('groups','elements')) {
    id <- out[[tb]]$metric=='spatial_sd';stopifnot(sum(id)==1L)
    for(nm in c('ess_mean','mcse','chain_gap')) out[[tb]][id,nm] <- NA_real_
    for(nm in c('rhat','ess_bulk','constant_chain')) out[[tb]][id,nm] <- amplitude[[nm]]
  }
  out
}

spatial_rule <- function(tables,warnings,sc) {
  reasons <- sc$spatial$robust_reasons(c(tables,list(warnings=warnings)))
  list(flagged=length(reasons)>0L,reasons=reasons)
}

finite_max <- function(x) if(!length(x) || !any(is.finite(x))) NA_real_ else max(x[is.finite(x)])

# One row of FLAG_COLUMNS for a saved fit.
fit_flags <- function(phase,fit,warnings,input,sc) {
  if(phase=='A2') {
    tables <- spatial_tables(fit,input,sc$spatial)
    r <- spatial_rule(tables,warnings,sc)
    e <- tables$elements
    return(data.frame(rule='spatial-robust-v1',warnings=length(warnings),max_group_rhat=finite_max(tables$groups$rhat),
      max_element_rhat=finite_max(c(e$rhat,tables$species$rhat)),unresolved_rhat=sum(!is.finite(e$rhat[e$metric!='range'])),
      spatial_trace_flags=sum(sc$spatial$robust_flag_rows(tables$spatial$diagnostics)),
      spatial_field_flags=sum(sc$spatial$robust_flag_rows(tables$spatial$field_diagnostics)),
      flagged=r$flagged,reasons=paste(r$reasons,collapse='; '),stringsAsFactors=FALSE))
  }
  rhats <- switch(phase,A1=a1_rhats(fit,input,sc$trace_stats),B=b_rhats(fit,input,sc$trace_summary),stop('Unknown phase: ',phase))
  d <- pr11_diagnostics(rhats,sc$safe_max)
  r <- pr11_rule(warnings,d,sc$flag_expr)
  data.frame(rule='pr11',warnings=length(warnings),max_group_rhat=d$max_group_rhat,max_element_rhat=d$max_element_rhat,
    unresolved_rhat=d$unresolved_rhat,spatial_trace_flags=NA_integer_,spatial_field_flags=NA_integer_,
    flagged=r$flagged,reasons=paste(r$reasons,collapse='; '),stringsAsFactors=FALSE)
}

path_within <- function(path,root) {
  prefix <- paste0(normalizePath(root,mustWork=FALSE),'/')
  ifelse(startsWith(path,prefix),substring(path,nchar(prefix)+1L),path)
}

# Every fit the selection reads: each new arm's first fit (at the schedule of
# the control's selected fit) and each control's selected fit.
selection_items <- function(phases,sds,study,archives,repo) {
  prov <- utils::read.csv(file.path(repo,CONTROL_PROVENANCE),stringsAsFactors=FALSE)
  rows <- list()
  for(phase in phases) {
    cs <- control_schedules(phase,archives)
    if(phase=='A2' && any(cs$schedule!='long')) stop('A2 first fits at the initial schedule are not part of this protocol')
    for(sd in sds) {
      fits <- as.character(mapply(fit_path,study,phase,sd,cs$schedule,cs$key,USE.NAMES=FALSE))
      rows[[length(rows)+1L]] <- data.frame(role='new',phase=phase,sd=sd,key=cs$key,schedule=cs$schedule,fit=fits,
        fit_label=path_within(fits,study),expected_md5=NA_character_,stringsAsFactors=FALSE)
    }
    p <- prov[prov$phase==phase & prov$selected,];i <- match(cs$key,p$key)
    if(anyNA(i) || !identical(p$schedule[i],cs$schedule)) stop('Control provenance disagrees with the recorded selection for phase ',phase)
    fits <- file.path(archives,p$control_fit[i])
    if(!identical(normalizePath(fits,mustWork=FALSE),normalizePath(cs$control_fit,mustWork=FALSE)))
      stop('Control provenance names other fits than the recorded selection for phase ',phase)
    if(phase=='A2' && !identical(p$file_md5[i],cs$control_fit_md5)) stop('A2 control fit md5s disagree with fits.csv')
    rows[[length(rows)+1L]] <- data.frame(role='control',phase=phase,sd=1,key=cs$key,schedule=cs$schedule,fit=fits,
      fit_label=p$control_fit[i],expected_md5=p$file_md5[i],stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}

# `repeats` names the fits as the single longer repeats of final mode.
check_first_fits <- function(items,repeats=FALSE) {
  new <- items$role=='new'
  missing <- items$fit[new & !file.exists(items$fit)]
  if(length(missing)) {
    if(repeats) stop('Longer repeats missing: ',length(missing),' (for example ',missing[1],
      '); run the final selection only after the longer-repeat launcher has written DONE with exit status 0')
    stop(length(missing),' first fits are missing (for example ',missing[1],
      '); run the selection only after the launcher has written DONE with exit status 0')
  }
  locks <- fit_lock_path(items$fit[new]);locks <- locks[dir.exists(locks)]
  if(length(locks)) stop(length(locks),' fit locks exist (for example ',locks[1],'); a fit may still be running')
  absent <- items$fit[!new & !file.exists(items$fit)]
  if(length(absent)) stop('Control fit missing: ',absent[1])
  invisible(TRUE)
}

# Flags of one saved fit, after checking that it is the fit the item names.
# flag_fn is fit_flags except in tests of the surrounding bookkeeping.
evaluate_item <- function(item,sc,archives,inputs_root,flag_fn=fit_flags) {
  sd <- as.numeric(item$sd)
  spec <- phase_jobs(item$phase,archives,inputs_root,item$key)[[1]]
  md5 <- unname(tools::md5sum(item$fit))
  if(!is.na(item$expected_md5) && !identical(md5,item$expected_md5))
    stop('Control fit md5 mismatch for ',item$fit_label,': recorded ',item$expected_md5,', found ',md5)
  saved <- readRDS(item$fit)
  if(item$role=='control') {
    check_protocol_mcmc(saved$mcmc,item$schedule,item$fit_label)
    fields <- if(item$phase=='A2') c('key','community','grid_index','replicate','arm','knots','input_md5') else names(spec$job)
    if(!identical(saved$job[fields],spec$job[fields])) stop('Control fit job record differs from the archived job: ',item$fit_label)
  } else {
    same <- c(phase=identical(saved$phase,item$phase),key=identical(saved$key,item$key),sd=identical(saved$sd,sd),
      schedule=identical(saved$schedule,item$schedule),mcmc=identical(saved$mcmc,schedule_mcmc(item$phase,item$schedule,archives)),
      input_md5=identical(saved$input_md5,spec$input_md5))
    if(!all(same)) stop('Fit ',item$fit_label,' is not the expected fit; differing: ',paste(names(same)[!same],collapse=', '))
    check_fit_prior(saved$fit,sd)
  }
  input <- readRDS(checked_input(spec$input_file,spec$input_md5))
  cbind(data.frame(role=item$role,phase=item$phase,sd=sd,key=item$key,schedule=item$schedule,
    fit=item$fit_label,fit_md5=md5,stringsAsFactors=FALSE),flag_fn(item$phase,saved$fit,saved$warnings,input,sc))
}

evaluate_items <- function(items,sc,archives,inputs_root,workers=1L,flag_fn=fit_flags) {
  if(!nrow(items)) return(NULL)
  work <- function(i) tryCatch(evaluate_item(items[i,,drop=FALSE],sc,archives,inputs_root,flag_fn),
    error=function(e) structure(conditionMessage(e),class='flag-error'))
  out <- if(workers==1L) lapply(seq_len(nrow(items)),work) else
    parallel::mclapply(seq_len(nrow(items)),work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  bad <- !vapply(out,is.data.frame,logical(1))
  if(any(bad)) stop('Flag computation failed for ',sum(bad),' fits: ',
    paste(items$fit_label[bad],vapply(out[bad],function(x) paste(as.character(x),collapse=' '),''),sep=': ',collapse=' | '))
  x <- do.call(rbind,out);rownames(x) <- NULL;x
}

# Write a table once; an identical recomputation is accepted, any other refused.
write_frozen_table <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- paste0(path,'.tmp');utils::write.csv(x,tmp,row.names=FALSE)
  if(file.exists(path)) {
    same <- identical(unname(tools::md5sum(tmp)),unname(tools::md5sum(path)));unlink(tmp)
    if(!same) stop('Refusing to replace ',path,': the recomputed table differs from the recorded one')
    return('unchanged')
  }
  stopifnot(file.rename(tmp,path))
  'written'
}

# The launcher job list of initial fits that need the single longer repeat.
long_keys <- function(selection) {
  x <- selection[selection$flagged,c('phase','sd','key'),drop=FALSE];rownames(x) <- NULL;x
}

# Each arm's selected fit: a flagged initial fit is replaced by its single
# longer repeat (whatever the repeat's flags); a first fit already at the
# longer schedule is kept; a control keeps its archived selected fit.
final_selected <- function(initial,repeats,longs,controls) {
  id <- function(x) paste(x$phase,x$sd,x$key)
  flagged <- initial[initial$flagged,,drop=FALSE]
  extra <- setdiff(id(repeats),id(flagged))
  if(length(extra)) stop('Unexpected long repeat: ',extra[1])
  if(nrow(repeats) && any(repeats$schedule!='long')) stop('A long repeat is not at the long schedule')
  row <- function(first,chosen,role,long_repeat) data.frame(phase=first$phase,sd=first$sd,key=first$key,role=role,
    first_schedule=if(role=='control') NA_character_ else first$schedule,selected_schedule=chosen$schedule,
    long_repeat=long_repeat,fit=chosen$fit,fit_md5=chosen$fit_md5,rule=chosen$rule,
    first_flagged=if(role=='control') NA else first$flagged,flagged=chosen$flagged,reasons=chosen$reasons,stringsAsFactors=FALSE)
  out <- list()
  for(i in seq_len(nrow(initial))) {
    first <- initial[i,,drop=FALSE]
    if(first$flagged) {
      k <- repeats[id(repeats)==id(first),,drop=FALSE]
      if(nrow(k)!=1L) stop('The single longer repeat is missing for ',id(first))
      out[[length(out)+1L]] <- row(first,k,'new',TRUE)
    } else out[[length(out)+1L]] <- row(first,first,'new',FALSE)
  }
  for(i in seq_len(nrow(longs))) out[[length(out)+1L]] <- row(longs[i,,drop=FALSE],longs[i,,drop=FALSE],'new',FALSE)
  for(i in seq_len(nrow(controls))) out[[length(out)+1L]] <- row(controls[i,,drop=FALSE],controls[i,,drop=FALSE],'control',NA)
  x <- do.call(rbind,out)
  x <- x[order(match(x$phase,PHASES),x$sd),,drop=FALSE];rownames(x) <- NULL
  x
}

# Default output directory of select.R: one per set of phases selected
# together (Amendment 1, R23), for example STUDY/selection/A1.
selection_dir <- function(study,phases) file.path(study,'selection',paste(phases,collapse='-'))

# Amendment 1 (R21): A1 no-harm criteria apply at 100 and 300 sites separately.
size_stratum <- function(key) ifelse(grepl('^jsdm-n0100-',key),'n100',ifelse(grepl('^jsdm-n0300-',key),'n300','all'))

# Per phase, stratum and arm (sd 1 is the control): selected fits still flagged
# after the single longer repeat, the gate's convergence criterion.
convergence_counts <- function(selected) {
  selected$stratum <- size_stratum(selected$key)
  groups <- split(selected,interaction(selected$phase,selected$stratum,selected$sd,drop=TRUE,lex.order=TRUE))
  x <- do.call(rbind,lapply(groups,function(g) data.frame(phase=g$phase[1],stratum=g$stratum[1],sd=g$sd[1],role=g$role[1],
    selected_fits=nrow(g),selected_flagged=sum(g$flagged),first_flagged=if(g$role[1]=='control') NA_integer_ else sum(g$first_flagged),
    long_repeats=if(g$role[1]=='control') NA_integer_ else sum(g$long_repeat),stringsAsFactors=FALSE)))
  x <- x[order(match(x$phase,PHASES),x$sd,x$stratum),,drop=FALSE];rownames(x) <- NULL
  x
}
