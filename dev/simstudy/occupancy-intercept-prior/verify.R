#!/usr/bin/env Rscript
# Independent audit of the phase A scorer, from the saved posterior draws.
#
#   Rscript verify.R --repo=REPO --study=STUDY --inputs-root=DIR --phase=A1|A2
#     --arms=1[,2,3,5] [--archives=DIR] [--summary=DIR] [--out=DIR] [--workers=N]
#   Rscript verify.R --repo=REPO --study=STUDY --phase=A [--summary=DIR] [--out=DIR]
#
# The fits to audit are read here from their own records: the controls from
# results/control-provenance.csv; a new arm from select.R's selected-fits.csv,
# or, for an SD select.R could not record (R27), from the controller's
# manual-selected.csv, with manual-missing.csv for the rest (AMENDMENT-2). The
# summary must hold exactly these fits (SUMMARY defaults to STUDY/summary). For
# each fit this script re-reads the saved fit and its input and recomputes,
# with code written here and nothing else:
#   the per-draw occupancy probabilities (A1: inverse logit of X B + U L + B0;
#   A2: the same with the spatial field, from its own kernel basis, checked
#   against the package's native basis), the truth (A1 from the generating
#   parameters, A2 the input's psi, checked against its parameters), posterior
#   means, type-7 0.025 and 0.975 quantiles, containment, the bands and groups,
#   E, |E|, mean absolute cell error, coverage, overall MAE and B0 bias and
#   coverage; and the convergence diagnostics of the phase's flag rule.
# From these per-community values it then recomputes every across-community
# table: summary-means.csv, b0-means.csv, convergence.csv (flag counts from its
# own diagnostics) and, with the control and a new arm, gate criteria 1 to 4
# (gate-detail.csv, gate.csv) and the schedule-matched sensitivity. It requires
# summarise.R's tables to agree (VERIFY_TOL_* below) and the stored diagnostics
# to be reproduced: the selection's flag tables when a selection exists (and
# supplement.R's manual-flags.csv), and for the controls also the archived
# result files. With --phase=A it recomputes the phase A rows from its own A1
# and A2 gate results. It never loads the scorer's code (neither the study's
# analysis, summary or supplement scripts nor the archived scorers); jobs.R
# supplies only paths, option parsing, md5-checked input lookup and the library
# fingerprint check. Writes to OUT/<phase> (OUT defaults to STUDY/verify):
# verify-arms-<arms>.csv (per fit), tables-arms-<arms>.csv (every table
# check), gate-audit-arms-<arms>.csv (its own gate) and audit-source.csv
# (md5 of this script, jobs.R and every table it read), and exits 1 on any
# disagreement.

# Tolerances. The summaries are means of at most 48,000 draws of probabilities
# in [0, 1], aggregated here along other routes (matrix sums in place of mean,
# a hand-written type-7 quantile, cell sets built afresh). Their rounding
# differences are of order 1e-14 in probability (the spatial-amplitude audit
# found at most 1.7e-14 for group means), so 1e-12 in probability, 1e-10 in
# percentage points, leaves a hundredfold margin while being far below any
# difference that could matter. Coverage is a proportion of identical
# containment decisions, so it must agree to rounding (1e-12).
# Rank-based diagnostics are not continuous in the last bit of the draws: with
# an even number of draws the two middle draws are tied in the folded
# statistic, and rounding decides whether that tie survives, which moves Rhat
# by up to about 1e-6 (seen when group means were first taken with colMeans
# rather than mean). As in the spatial-amplitude audit (verify-robust.R), the
# strict diagnostic check therefore rebuilds the draws in the archive's exact
# numerical representation (plogis, per-draw mean, the Cholesky and solve
# basis, the robust rule's blocked field products) and requires Rhat to 1e-12
# (as check-flags.R) and mean-based ESS to relative 1e-9; the package's native
# basis is checked separately against that basis (1e-10, as score.R did).
VERIFY_TOL_POINTS <- 1e-10
VERIFY_TOL_COVERAGE <- 1e-12
VERIFY_TOL_RHAT <- 1e-12
VERIFY_TOL_ESS_RELATIVE <- 1e-9
VERIFY_TOL_BASIS <- 1e-10

audit_logistic <- function(x) stats::plogis(x)
# Per-draw mean over a set of cells (rows) of a cells-by-draws matrix.
draw_means <- function(d,rows) vapply(seq_len(ncol(d)),function(k) mean(d[rows,k]),numeric(1))

# Type-7 sample quantile of a sorted vector: index 1 + (N - 1) p, linear
# interpolation between the two neighbouring order statistics.
audit_q7 <- function(s,p) {
  index <- 1+(length(s)-1)*p;lo <- floor(index);h <- index-lo
  if(h==0) s[lo] else (1-h)*s[lo]+h*s[lo+1L]
}

# Cells from a draws matrix (cells by pooled draws) or array (cells, iterations, chains).
audit_cells_from_draws <- function(phase,draws,truth,target_prevalence=NULL,n_original=nrow(truth)) {
  d <- if(length(dim(draws))==3L) matrix(draws,dim(draws)[1]) else draws
  truth <- unname(as.matrix(truth));stopifnot(nrow(d)==length(truth))
  N <- ncol(d);estimate <- as.vector(d%*%rep(1,N))/N
  q <- t(apply(d,1,function(x) {s <- sort(x);c(audit_q7(s,.025),audit_q7(s,.975))}))
  tv <- as.vector(truth)
  list(phase=phase,n=nrow(truth),S=ncol(truth),n_original=as.integer(n_original),truth=truth,
    estimate=matrix(estimate,nrow(truth)),lower=matrix(q[,1],nrow(truth)),upper=matrix(q[,2],nrow(truth)),
    covered=matrix(q[,1]<=tv & tv<=q[,2],nrow(truth)),target_prevalence=target_prevalence)
}

# Bands (below 0.2, 0.2 to 0.8 inclusive, above 0.8), the A1 rare group (every
# fitted site of each species with mean truth below 0.2 over all fitted sites),
# the A1 all-site scope, and the A2 prevalence groups and 1% plus 5% group.
audit_groups <- function(cells) {
  tv <- as.vector(cells$truth);ev <- as.vector(cells$estimate);cv <- as.vector(cells$covered)
  site <- rep(seq_len(cells$n),cells$S);species <- rep(seq_len(cells$S),each=cells$n)
  sets <- list()
  add <- function(scope,group,keep) sets[[length(sets)+1L]] <<- list(scope=scope,group=group,idx=which(keep))
  bands <- function(scope,keep) {add(scope,'all',keep);add(scope,'low',keep & tv<.2);add(scope,'middle',keep & tv>=.2 & tv<=.8)
    add(scope,'high',keep & tv>.8)}
  if(cells$phase=='A1') {
    bands('primary',site<=cells$n_original)
    mean_truth <- vapply(seq_len(cells$S),function(s) mean(tv[species==s]),numeric(1))
    add('primary','rare_below_20pct',species %in% which(mean_truth<.2))
    if(cells$n>cells$n_original) bands('allsites',rep(TRUE,length(tv)))
  } else {
    bands('primary',rep(TRUE,length(tv)))
    tp <- cells$target_prevalence
    for(p in c(1,5,25,75)) add('primary',paste0('prevalence_',p,'pct'),abs(tp[species]-p/100)<1e-10)
    add('primary','rare_1_5pct',abs(tp[species]-.01)<1e-10 | abs(tp[species]-.05)<1e-10)
  }
  rows <- lapply(sets,function(s) {
    if(!length(s$idx)) return(NULL)
    e <- ev[s$idx]-tv[s$idx]
    data.frame(scope=s$scope,group=s$group,cells=length(s$idx),signed_error=100*mean(e),abs_signed_error=abs(100*mean(e)),
      mean_abs_cell_error=100*mean(abs(e)),coverage=mean(cv[s$idx]),stringsAsFactors=FALSE)
  })
  x <- do.call(rbind,rows);rownames(x) <- NULL;x
}

audit_b0 <- function(B0_output,truth) {
  d <- matrix(B0_output,dim(B0_output)[1]);stopifnot(nrow(d)==length(truth))
  q <- t(apply(d,1,function(x) {s <- sort(x);c(audit_q7(s,.025),audit_q7(s,.975))}))
  bias <- as.vector(d%*%rep(1,ncol(d)))/ncol(d)-truth
  c(species=length(truth),b0_bias=mean(bias),b0_abs_bias=mean(abs(bias)),b0_coverage=mean(q[,1]<=truth & truth<=q[,2]))
}

# Phase A1: X B + U L + B0 per draw and for the truth (the input's generating
# B0, B, U and L with the fitted design, checked against its eta).
audit_a1_cells <- function(fit,input) {
  js <- fit$results_output$jsdm_output;X <- fit$X_psi;n <- nrow(X)
  S <- dim(js$B0_output)[1];ni <- dim(js$B0_output)[2];nc <- dim(js$B0_output)[3];P <- ncol(X)
  stopifnot(n==input$n,S==ncol(input$data$OTU),fit$infos$ps==0L)
  tp <- input$truth;d <- ncol(tp$U)
  eta <- X%*%tp$B+tp$U%*%tp$L+matrix(tp$B0,n,S,byrow=TRUE)
  stopifnot(max(abs(eta-tp$eta))<1e-10)
  draws <- matrix(NA_real_,n*S,ni*nc)
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    e <- X%*%matrix(js$B_output[,,it,ch],P,S)+matrix(js$U_output[,,it,ch],n)%*%matrix(js$L_output[,,it,ch],ncol(js$U_output),S)+
      matrix(js$B0_output[,it,ch],n,S,byrow=TRUE)
    draws[,(ch-1L)*ni+it] <- audit_logistic(as.vector(e))
  }
  cells <- audit_cells_from_draws('A1',draws,audit_logistic(eta),n_original=100L)
  cells$draws <- draws;cells$ni <- ni;cells$nc <- nc
  cells
}

# The spatial basis of every grid range: the squared-exponential kernel from
# sites to knots times the inverse transposed Cholesky factor of the knot
# kernel with the package's 1e-5 jitter, K_sx L^-T with L lower triangular,
# solved as the archive did so that the draws are bit-identical to it.
audit_a2_bases <- function(fit) {
  xs <- fit$infos$list_Xs;k <- xs$X_tilde;s <- xs$X_s
  d2 <- function(a,b) outer(a[,1],b[,1],'-')^2+outer(a[,2],b[,2],'-')^2
  lapply(fit$infos$l_s_grid,function(l) {
    L <- t(chol(exp(-d2(k,k)/(2*l^2))+diag(1e-5,nrow(k))))
    t(solve(L,t(exp(-d2(s,k)/(2*l^2)))))[xs$Xs_index,,drop=FALSE]
  })
}

# The field draws of the robust rule's screens, species by species, with each
# chain's iterations at one grid range multiplied in blocks of 100 as that rule does.
audit_a2_field <- function(fit,bases) {
  js <- fit$results_output$jsdm_output;n <- nrow(fit$Xs);ps <- fit$infos$ps
  S <- dim(js$B0_output)[1];ni <- dim(js$B0_output)[2];nc <- dim(js$B0_output)[3]
  out <- matrix(NA_real_,n*S,ni*nc)
  for(s in seq_len(S)) for(ch in seq_len(nc)) for(g in sort(unique(js$idx_ls_output[,ch]))) {
    ii <- which(js$idx_ls_output[,ch]==g)
    for(block in split(ii,ceiling(seq_along(ii)/100)))
      out[(s-1L)*n+seq_len(n),(ch-1L)*ni+block] <- bases[[g]]%*%matrix(js$Bs_output[,s,block,ch,drop=FALSE],ps,length(block))
  }
  out
}

# Per-draw probability of every cell of a binary spatial fit.
audit_a2_draws <- function(fit,bases=audit_a2_bases(fit)) {
  js <- fit$results_output$jsdm_output;n <- nrow(fit$Xs);ps <- fit$infos$ps;P <- ncol(fit$X_psi)
  S <- dim(js$B0_output)[1];ni <- dim(js$B0_output)[2];nc <- dim(js$B0_output)[3]
  stopifnot(fit$infos$n_factors==0L)
  out <- matrix(NA_real_,n*S,ni*nc)
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    f <- bases[[js$idx_ls_output[it,ch]]]%*%matrix(js$Bs_output[,,it,ch],ps,S)
    e <- f+fit$X_psi%*%matrix(js$B_output[,,it,ch],P,S)+matrix(js$B0_output[,it,ch],n,S,byrow=TRUE)
    out[,(ch-1L)*ni+it] <- audit_logistic(as.vector(e))
  }
  out
}

audit_a2_cells <- function(fit,input,native=NULL) {
  t <- input$truth;n <- input$settings$n;S <- input$settings$S
  stopifnot(fit$infos$ps==100L,fit$infos$model=='binary',identical(fit$infos$speciesNames,colnames(input$data$binary$OTU)),
    max(abs(fit$X_psi-t$X))<1e-12,max(abs(fit$Xs-t$Xs))<1e-12)
  truth <- audit_logistic(matrix(t$B0,n,S,byrow=TRUE)+t$X%*%t(t$B)+t$field)
  truth_difference <- max(abs(truth-t$psi))
  if(truth_difference>1e-12) stop('Input psi disagrees with its generating parameters: ',truth_difference)
  if(max(abs(colMeans(t$psi)-t$target_prevalence))>1e-10) stop('Input psi does not have the designed prevalences')
  bases <- audit_a2_bases(fit)
  basis_difference <- if(is.null(native)) NA_real_ else native(fit,bases)
  x <- audit_a2_draws(fit,bases)
  cells <- audit_cells_from_draws('A2',x,t$psi,t$target_prevalence)
  cells$draws <- x;cells$field <- audit_a2_field(fit,bases)
  cells$ni <- dim(fit$results_output$jsdm_output$B0_output)[2];cells$nc <- dim(fit$results_output$jsdm_output$B0_output)[3]
  cells$truth_difference <- truth_difference;cells$basis_difference <- basis_difference
  cells
}

# ---------------------------------------------------------------------------
# Convergence diagnostics, recomputed
# ---------------------------------------------------------------------------

audit_rhat <- function(x) suppressWarnings(as.numeric(posterior::rhat(x)))
audit_safe_max <- function(x) if(!length(x) || all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)
audit_finite_max <- function(x) if(!length(x) || !any(is.finite(x))) NA_real_ else max(x[is.finite(x)])
as_chains <- function(v,ni,nc) matrix(v,ni,nc)

# A1 and B rule (current-main-recheck select.R:17-19 with run.R:79-81): groups
# are the mean probability over each band of the original and all fitted
# sites; elements B0, B and every original-site probability; plus sigma_h.
audit_diag_a1 <- function(fit,cells,warnings) {
  js <- fit$results_output$jsdm_output;ni <- cells$ni;nc <- cells$nc;n <- cells$n;S <- cells$S
  tv <- as.vector(cells$truth);site <- rep(seq_len(n),S);original <- which(site<=100L)
  group <- list()
  for(base in list(original,seq_len(n*S))) {
    p <- tv[base]
    for(m in list(base,base[p<.2],base[p>=.2 & p<=.8],base[p>.8])) group[[length(group)+1L]] <- draw_means(cells$draws,m)
  }
  groups <- vapply(group,function(v) audit_rhat(as_chains(v,ni,nc)),numeric(1))
  block <- function(a) {m <- matrix(a,prod(head(dim(a),-2L)));vapply(seq_len(nrow(m)),function(k) audit_rhat(as_chains(m[k,],ni,nc)),numeric(1))}
  elements <- c(block(js$B0_output),block(js$B_output),vapply(original,function(k) audit_rhat(as_chains(cells$draws[k,],ni,nc)),numeric(1)))
  additional <- audit_rhat(js$sigmah_output)
  all <- c(groups,elements,additional)
  mg <- audit_safe_max(groups);me <- audit_safe_max(c(elements,additional));un <- sum(!is.finite(all) | all<=0)
  flagged <- warnings>0L || un>0L || !is.finite(mg) || !is.finite(me) || mg>1.05 || me>1.05
  list(summary=list(warnings=warnings,max_group_rhat=mg,max_element_rhat=me,unresolved_rhat=un,spatial_trace_flags=NA,
    spatial_field_flags=NA,flagged=flagged),vectors=list(groups=groups,elements=elements,additional=additional))
}

# Mean-based trace diagnostics of the spatial archive: unavailable when a chain
# is constant or there are fewer than four iterations.
audit_td <- function(x) {
  if(any(apply(x,2,sd)==0) || nrow(x)<4L) return(c(rhat=NA_real_,ess_mean=NA_real_))
  c(rhat=audit_rhat(x),ess_mean=suppressWarnings(as.numeric(posterior::ess_mean(x))))
}
# Rank and quantile diagnostics of the robust-v1 rule: unavailable when a chain is constant.
audit_rd <- function(x) {
  stopifnot(all(is.finite(x)))
  if(any(apply(x,2,function(z) all(z==z[1])))) return(c(rhat=NA_real_,ess_bulk=NA_real_,ess_median=NA_real_,ess_q025=NA_real_,ess_q975=NA_real_))
  q <- function(p) suppressWarnings(as.numeric(posterior::ess_quantile(x,probs=p)))
  c(rhat=audit_rhat(x),ess_bulk=suppressWarnings(as.numeric(posterior::ess_bulk(x))),ess_median=q(.5),ess_q025=q(.025),ess_q975=q(.975))
}
audit_screen_fails <- function(m) apply(m,1,function(r) any(!is.finite(r)) || r[['rhat']]>1.05 || any(r[-1]<100))

# A2 rule (spatial-amplitude robust.R:100-110): occupancy groups, species and
# cells, the B0, slope, range and amplitude blocks (amplitude by rank
# diagnostics), and the rank/quantile screens of the amplitude, per-species
# field mean, RMS and projection on the centred true field, and every
# pointwise field entry. Native warnings do not flag.
audit_diag_a2 <- function(fit,input,cells,warnings) {
  js <- fit$results_output$jsdm_output;ni <- cells$ni;nc <- cells$nc;n <- cells$n;S <- cells$S
  tv <- as.vector(cells$truth);species <- rep(seq_len(S),each=n);tp <- cells$target_prevalence
  occ_sets <- list(all=seq_along(tv),low=which(tv<.2),medium=which(tv>=.2 & tv<=.8),high=which(tv>.8))
  for(p in c(1,5,25,75)) occ_sets[[paste0('prevalence_',p,'pct')]] <- which(abs(tp[species]-p/100)<1e-10)
  occ_sets <- occ_sets[lengths(occ_sets)>0L]
  occ <- t(vapply(occ_sets,function(ix) audit_td(as_chains(draw_means(cells$draws,ix),ni,nc)),numeric(2)))
  amp <- audit_rd(js$sigmabs_output)
  block_list <- list(intercept=js$B0_output,environment_slope=js$B_output,
    range=array(fit$infos$l_s_grid[js$idx_ls_output],c(1L,ni,nc)))
  block_group <- list();block_element <- list()
  for(nm in names(block_list)) {
    m <- matrix(block_list[[nm]],prod(head(dim(block_list[[nm]]),-2L)))
    block_element[[nm]] <- vapply(seq_len(nrow(m)),function(k) audit_td(as_chains(m[k,],ni,nc))[['rhat']],numeric(1))
    block_group[[nm]] <- audit_td(as_chains(draw_means(m,seq_len(nrow(m))),ni,nc))[['rhat']]
  }
  groups <- c(stats::setNames(occ[,'rhat'],paste0('occupancy:',rownames(occ))),
    stats::setNames(unlist(block_group),paste0(names(block_group),':all')),'spatial_sd:all'=amp[['rhat']])
  cell_rhat <- vapply(seq_along(tv),function(k) audit_td(as_chains(cells$draws[k,],ni,nc))[['rhat']],numeric(1))
  elements <- c(cell_rhat,unlist(block_element),amp[['rhat']])
  element_metric <- c(rep('occupancy',length(tv)),rep(names(block_element),lengths(block_element)),'spatial_sd')
  sp <- t(vapply(seq_len(S),function(s) audit_td(as_chains(draw_means(cells$draws,which(species==s)),ni,nc)),numeric(2)))
  # Spatial screens.
  tf <- input$truth$field;tc <- sweep(tf,2,colMeans(tf),'-');denom <- colSums(tc^2)
  traces <- list(amplitude=js$sigmabs_output)
  fm <- fr <- fp <- list()
  for(s in seq_len(S)) {
    f <- cells$field[species==s,,drop=FALSE]
    fm[[s]] <- as_chains(colMeans(f),ni,nc);fr[[s]] <- as_chains(sqrt(colMeans(f^2)),ni,nc)
    fp[[s]] <- as_chains(colSums(f*tc[,s])/denom[s],ni,nc)
  }
  traces <- c(traces,stats::setNames(fm,paste0('field_mean_',seq_len(S))),stats::setNames(fr,paste0('field_rms_',seq_len(S))),
    stats::setNames(fp,paste0('field_projection_',seq_len(S))))
  trace_diag <- t(vapply(traces,audit_rd,numeric(5)))
  field_diag <- t(vapply(seq_len(n*S),function(k) audit_rd(as_chains(cells$field[k,],ni,nc)),numeric(5)))
  trace_flags <- sum(audit_screen_fails(trace_diag));field_flags <- sum(audit_screen_fails(field_diag))
  unresolved <- sum(!is.finite(elements[element_metric!='range']))
  reasons <- c(any(groups>1.05,na.rm=TRUE),any(sp[,'rhat']>1.05,na.rm=TRUE),any(elements>1.05,na.rm=TRUE),
    any(occ[,'ess_mean']<100,na.rm=TRUE),any(sp[,'ess_mean']<100,na.rm=TRUE),unresolved>0L,trace_flags>0L,field_flags>0L)
  list(summary=list(warnings=warnings,max_group_rhat=audit_finite_max(groups),max_element_rhat=audit_finite_max(c(elements,sp[,'rhat'])),
      unresolved_rhat=unresolved,spatial_trace_flags=trace_flags,spatial_field_flags=field_flags,flagged=any(reasons)),
    vectors=list(groups=groups,occupancy_ess=occ[,'ess_mean'],elements=elements,element_metric=element_metric,species=sp,
      trace_diag=trace_diag,field_diag=field_diag))
}

# ---------------------------------------------------------------------------
# Stored diagnostics
# ---------------------------------------------------------------------------

# The selection's flag tables (every fit select.R evaluated), if recorded.
selection_flag_rows <- function(sel_dir) {
  files <- file.path(sel_dir,c('long-selection.csv','repeat-flags.csv','long-fit-flags.csv','control-flags.csv'))
  x <- lapply(files[file.exists(files)],function(f) utils::read.csv(f,stringsAsFactors=FALSE,
    colClasses=c(role='character',sd='numeric',key='character',fit='character',fit_md5='character',reasons='character')))
  if(!length(x)) return(NULL)
  do.call(rbind,x)
}

compare_summary <- function(mine,stored) {
  num <- function(a,b) if(is.na(a) && is.na(b)) 0 else if(is.na(a) || is.na(b)) Inf else abs(a-b)
  counts <- c('warnings','unresolved_rhat','spatial_trace_flags','spatial_field_flags')
  same_count <- vapply(counts,function(k) identical(as.numeric(mine[[k]]),as.numeric(stored[[k]])) ||
    (is.na(mine[[k]]) && is.na(stored[[k]])),logical(1))
  list(group=num(mine$max_group_rhat,stored$max_group_rhat),element=num(mine$max_element_rhat,stored$max_element_rhat),
    counts=all(same_count),flag=identical(as.logical(mine$flagged),as.logical(stored$flagged)))
}

vector_difference <- function(a,b,relative=FALSE) {
  a <- as.numeric(a);b <- as.numeric(b)
  if(length(a)!=length(b) || !identical(is.na(a),is.na(b))) return(Inf)
  ok <- !is.na(a);if(!any(ok)) return(0)
  d <- abs(a[ok]-b[ok]);if(relative) d <- d/pmax(abs(b[ok]),1)
  max(d)
}

# Stored per-element diagnostics of the archived control results.
archive_comparison <- function(phase,f,archives,diag) {
  if(phase=='A1') {
    rf <- file.path(archives,sub('-fit[.]rds$','-result.rds',f$fit));r <- readRDS(rf)
    stored <- list(warnings=length(r$warnings),max_group_rhat=r$diagnostics$max_group_rhat,max_element_rhat=r$diagnostics$max_element_rhat,
      unresolved_rhat=r$diagnostics$unresolved_rhat,spatial_trace_flags=NA,spatial_field_flags=NA,
      flagged=length(r$warnings)>0L || r$diagnostics$unresolved_rhat>0L || !is.finite(r$diagnostics$max_group_rhat) ||
        !is.finite(r$diagnostics$max_element_rhat) || r$diagnostics$max_group_rhat>1.05 || r$diagnostics$max_element_rhat>1.05)
    rh <- max(vector_difference(diag$vectors$groups,r$groups$rhat),vector_difference(diag$vectors$elements,r$elements$rhat),
      vector_difference(diag$vectors$additional,r$additional_diagnostics$rhat))
    return(list(source=paste('archive',basename(rf)),stored=stored,rhat=rh,ess=0,crossings=0L))
  }
  fits <- utils::read.csv(file.path(archive_dir(archives,'amplitude'),'robust-v1/summary-binary-final/fits.csv'),stringsAsFactors=FALSE,
    colClasses=c(reasons='character'))
  row <- fits[fits$prior=='inverse_gamma' & fits$key==f$key,,drop=FALSE]
  if(nrow(row)!=1L || !identical(row$fit_md5,f$fit_md5) || !identical(unname(tools::md5sum(row$result_file)),row$result_md5))
    stop('No archived robust-v1 result of this control fit: ',f$key)
  r <- readRDS(row$result_file)
  rows_screen <- function(d) sum(apply(as.matrix(d[c('rhat','ess_bulk','ess_median','ess_q025','ess_q975')]),1,function(z)
    any(!is.finite(z)) || z[1]>1.05 || any(z[-1]<100)))
  el <- r$elements
  stored <- list(warnings=length(r$warnings),max_group_rhat=audit_finite_max(r$groups$rhat),
    max_element_rhat=audit_finite_max(c(el$rhat,r$species$rhat)),unresolved_rhat=sum(!is.finite(el$rhat[el$metric!='range'])),
    spatial_trace_flags=rows_screen(r$spatial$diagnostics),spatial_field_flags=rows_screen(r$spatial$field_diagnostics),
    flagged=length(r$reasons)>0L)
  v <- diag$vectors
  gkey <- paste(r$groups$metric,r$groups$group,sep=':')
  if(!setequal(gkey,names(v$groups))) return(list(source='archive robust-v1',stored=stored,rhat=Inf,ess=Inf,crossings=NA_integer_))
  occ <- r$groups$metric=='occupancy'
  td <- r$spatial$diagnostics[match(rownames(v$trace_diag),r$spatial$diagnostics$quantity),]
  fd <- r$spatial$field_diagnostics
  cols <- c('rhat','ess_bulk','ess_median','ess_q025','ess_q975')
  rh <- max(vector_difference(v$groups,r$groups$rhat[match(names(v$groups),gkey)]),
    vector_difference(v$elements,el$rhat[order(match(el$metric,unique(v$element_metric)),el$element)]),
    vector_difference(v$species[,'rhat'],r$species$rhat),vector_difference(v$trace_diag[,'rhat'],td$rhat),
    vector_difference(v$field_diag[,'rhat'],fd$rhat))
  ess <- max(vector_difference(v$occupancy_ess,r$groups$ess_mean[occ][match(names(v$occupancy_ess),r$groups$group[occ])],TRUE),
    vector_difference(v$species[,'ess_mean'],r$species$ess_mean,TRUE),
    vector_difference(v$trace_diag[,cols[-1]],as.matrix(td[cols[-1]]),TRUE),vector_difference(v$field_diag[,cols[-1]],as.matrix(fd[cols[-1]]),TRUE))
  crossings <- sum(audit_screen_fails(v$trace_diag)!=audit_screen_fails(as.matrix(td[cols])))+
    sum(audit_screen_fails(v$field_diag)!=audit_screen_fails(as.matrix(fd[cols])))
  list(source=paste('archive',basename(row$result_file)),stored=stored,rhat=rh,ess=ess,crossings=crossings)
}

# ---------------------------------------------------------------------------
# One fit
# ---------------------------------------------------------------------------

# Stratum and community of a key, stated here afresh: the A1 fits at 100 and
# 300 sites of replicate r form generating community r (R21); A2 keys are communities.
audit_stratum <- function(phase,key) if(phase=='A1') ifelse(grepl('^jsdm-n0100-',key),'n100',ifelse(grepl('^jsdm-n0300-',key),'n300',NA)) else rep('all',length(key))
audit_community <- function(phase,key) if(phase=='A1') sub('^jsdm-n0[13]00-','',key) else key

audit_fit <- function(f,phase,study,archives,inputs_root,tabs,flag_rows,native) {
  path <- if(f$arm=='control') file.path(archives,f$fit) else file.path(study,f$fit)
  if(!identical(unname(tools::md5sum(path)),f$fit_md5)) stop('Fit md5 differs from the summary: ',f$fit)
  saved <- readRDS(path);fit <- saved$fit;warnings <- length(saved$warnings)
  spec <- phase_jobs(phase,archives,inputs_root,f$key)[[1]]
  input <- readRDS(checked_input(spec$input_file,spec$input_md5))
  if(phase=='A1') {cells <- audit_a1_cells(fit,input);diag <- audit_diag_a1(fit,cells,warnings)}
  else {cells <- audit_a2_cells(fit,input,native);diag <- audit_diag_a2(fit,input,cells,warnings)}
  g <- audit_groups(cells)
  pick <- function(x) x[x$sd==f$sd & x$kind==f$kind & x$key==f$key,,drop=FALSE]
  be <- pick(tabs$band);cv <- pick(tabs$coverage);ma <- pick(tabs$mae);b0 <- pick(tabs$b0)
  id <- function(x) paste(x$scope,x$group)
  same_sets <- setequal(id(g),id(be)) && nrow(g)==nrow(be) && setequal(id(g),id(cv)) && nrow(cv)==nrow(g)
  i <- match(id(g),id(be));j <- match(id(g),id(cv))
  diff <- function(a,b) if(anyNA(a) || anyNA(b) || length(a)!=length(b)) Inf else if(!length(a)) 0 else max(abs(a-b))
  all_g <- g[g$group=='all',];k <- match(all_g$scope,ma$scope)
  b <- audit_b0(fit$results_output$jsdm_output$B0_output,input$truth$B0)
  stored_sel <- if(!is.null(flag_rows)) flag_rows[flag_rows$role==f$arm & flag_rows$sd==f$sd & flag_rows$key==f$key &
    flag_rows$fit_md5==f$fit_md5,,drop=FALSE] else NULL
  if(!is.null(stored_sel) && nrow(stored_sel)>1L) stored_sel <- stored_sel[1,,drop=FALSE]
  sel <- if(!is.null(stored_sel) && nrow(stored_sel)==1L) compare_summary(diag$summary,stored_sel) else NULL
  arch <- if(f$arm=='control') archive_comparison(phase,f,archives,diag) else NULL
  ac <- if(!is.null(arch)) compare_summary(diag$summary,arch$stored) else NULL
  out <- data.frame(phase=phase,sd=f$sd,arm=f$arm,kind=f$kind,key=f$key,schedule=f$schedule,fit_md5=f$fit_md5,
    groups_compared=nrow(g),group_sets_equal=same_sets,cells_equal=same_sets && identical(as.integer(g$cells),as.integer(be$cells[i])),
    max_diff_signed_error=diff(g$signed_error,be$signed_error[i]),max_diff_abs_signed_error=diff(g$abs_signed_error,be$abs_signed_error[i]),
    max_diff_mean_abs_cell_error=diff(g$mean_abs_cell_error,be$mean_abs_cell_error[i]),max_diff_coverage=diff(g$coverage,cv$coverage[j]),
    max_diff_mae=diff(all_g$mean_abs_cell_error,ma$mae[k]),
    max_diff_b0=if(nrow(b0)==1L && b0$species==b[['species']]) max(abs(c(b0$b0_bias,b0$b0_abs_bias,b0$b0_coverage)-b[-1])) else Inf,
    truth_difference=cells$truth_difference %||% NA_real_,native_basis_difference=cells$basis_difference %||% NA_real_,
    selection_diagnostics=!is.null(sel),sel_group_rhat_diff=sel$group %||% NA_real_,sel_element_rhat_diff=sel$element %||% NA_real_,
    sel_counts_equal=sel$counts %||% NA,sel_flag_equal=sel$flag %||% NA,
    archive_source=arch$source %||% NA_character_,archive_group_rhat_diff=ac$group %||% NA_real_,archive_element_rhat_diff=ac$element %||% NA_real_,
    archive_counts_equal=ac$counts %||% NA,archive_flag_equal=ac$flag %||% NA,archive_rhat_vector_diff=arch$rhat %||% NA_real_,
    archive_ess_relative_diff=arch$ess %||% NA_real_,archive_threshold_crossings=arch$crossings %||% NA_integer_,
    flagged=diag$summary$flagged,max_group_rhat=diag$summary$max_group_rhat,max_element_rhat=diag$summary$max_element_rhat,stringsAsFactors=FALSE)
  ok <- function(x,tol) !is.na(x) && x<=tol
  out$pass <- out$group_sets_equal && out$cells_equal && ok(out$max_diff_signed_error,VERIFY_TOL_POINTS) &&
    ok(out$max_diff_abs_signed_error,VERIFY_TOL_POINTS) && ok(out$max_diff_mean_abs_cell_error,VERIFY_TOL_POINTS) &&
    ok(out$max_diff_mae,VERIFY_TOL_POINTS) && ok(out$max_diff_coverage,VERIFY_TOL_COVERAGE) && ok(out$max_diff_b0,VERIFY_TOL_POINTS) &&
    (phase!='A2' || ok(out$native_basis_difference,VERIFY_TOL_BASIS)) &&
    (!out$selection_diagnostics || (ok(out$sel_group_rhat_diff,VERIFY_TOL_RHAT) && ok(out$sel_element_rhat_diff,VERIFY_TOL_RHAT) &&
      isTRUE(out$sel_counts_equal) && isTRUE(out$sel_flag_equal))) &&
    (is.null(arch) || (ok(out$archive_group_rhat_diff,VERIFY_TOL_RHAT) && ok(out$archive_element_rhat_diff,VERIFY_TOL_RHAT) &&
      isTRUE(out$archive_counts_equal) && isTRUE(out$archive_flag_equal) && ok(out$archive_rhat_vector_diff,VERIFY_TOL_RHAT) &&
      ok(out$archive_ess_relative_diff,VERIFY_TOL_ESS_RELATIVE) && identical(out$archive_threshold_crossings,0L))) &&
    (out$selection_diagnostics || !is.null(arch))
  meta <- data.frame(sd=f$sd,arm=f$arm,kind=f$kind,key=f$key,stratum=audit_stratum(phase,f$key),community=audit_community(phase,f$key),
    stringsAsFactors=FALSE)
  list(row=out,groups=cbind(meta[rep(1L,nrow(g)),,drop=FALSE],g,row.names=NULL),
    b0=cbind(meta,as.data.frame(as.list(b))),flagged=diag$summary$flagged)
}

# ---------------------------------------------------------------------------
# Which fits: the selection, read here from its own records
# ---------------------------------------------------------------------------

# The selected fit of every community of every requested arm, the first fits of
# longer repeats, and the communities without a valid selected fit (R27).
# Controls: results/control-provenance.csv. New arms: select.R's
# selected-fits.csv, or, for an SD select.R could not record, the controller's
# manual-selected.csv; manual-missing.csv lists the rest.
audit_selection <- function(phase,arms,repo,study,summary) {
  keys <- phase_keys(phase)
  prov <- utils::read.csv(file.path(repo,'dev/simstudy/occupancy-intercept-prior/results/control-provenance.csv'),stringsAsFactors=FALSE)
  prov <- prov[prov$phase==phase & prov$selected,,drop=FALSE]
  if(!setequal(prov$key,keys) || anyDuplicated(prov$key)) stop('Control provenance does not cover phase ',phase)
  control_schedule <- stats::setNames(prov$schedule,prov$key)
  sel_dir <- file.path(study,'selection',phase)
  rd <- function(name,cc=c(key='character',fit='character',fit_md5='character')) {f <- file.path(sel_dir,name)
    if(file.exists(f)) utils::read.csv(f,stringsAsFactors=FALSE,colClasses=cc) else NULL}
  sf <- rd('selected-fits.csv');ls <- rd('long-selection.csv')
  mm <- rd('manual-missing.csv',c(key='character'));ms <- rd('manual-selected.csv')
  selected <- if(1 %in% arms) data.frame(sd=1,key=prov$key,fit=prov$control_fit,fit_md5=prov$file_md5,stringsAsFactors=FALSE) else NULL
  first <- NULL;missing <- NULL
  for(sd in setdiff(arms,1)) {
    a <- if(is.null(sf)) sf else sf[sf$role=='new' & sf$sd==sd,,drop=FALSE]
    b <- if(is.null(ms)) ms else ms[ms$sd==sd,,drop=FALSE]
    gone <- if(is.null(mm)) character() else mm$key[mm$sd==sd]
    if(NROW(a) && NROW(b)) stop('SD ',sd,' is in both selected-fits.csv and manual-selected.csv')
    if(NROW(a)) {
      selected <- rbind(selected,data.frame(sd=sd,key=a$key,fit=a$fit,fit_md5=a$fit_md5,stringsAsFactors=FALSE))
      r <- a[a$long_repeat %in% TRUE,,drop=FALSE]
      for(i in seq_len(nrow(r))) {
        l <- ls[ls$role=='new' & ls$sd==sd & ls$key==r$key[i] & ls$schedule==r$first_schedule[i],,drop=FALSE]
        if(nrow(l)!=1L) stop('No recorded first fit of the longer repeat ',r$fit[i])
        first <- rbind(first,data.frame(sd=sd,key=l$key,fit=l$fit,stringsAsFactors=FALSE))
      }
    } else if(NROW(b)) {
      selected <- rbind(selected,data.frame(sd=sd,key=b$key,fit=b$fit,fit_md5=b$fit_md5,stringsAsFactors=FALSE))
      r <- b[b$schedule=='long' & control_schedule[b$key]=='initial',,drop=FALSE]
      if(nrow(r)) first <- rbind(first,data.frame(sd=sd,key=r$key,fit=sprintf('fits/%s/sd%s/initial/%s-fit.rds',phase,format(sd),r$key),
        stringsAsFactors=FALSE))
    }
    both <- intersect(gone,selected$key[selected$sd==sd])
    if(length(both)) stop('Community both selected and recorded as missing for SD ',sd,': ',both[1])
    missing <- rbind(missing,data.frame(sd=rep(sd,length(setdiff(keys,selected$key[selected$sd==sd]))),
      key=setdiff(keys,selected$key[selected$sd==sd]),stringsAsFactors=FALSE))
  }
  flag_rows <- if(dir.exists(sel_dir)) selection_flag_rows(sel_dir) else NULL
  # The flags supplement.R recorded, by flags.R's frozen rule, for an SD select.R could not record.
  mf <- file.path(summary,phase,'manual-flags.csv')
  if(file.exists(mf)) {
    m <- utils::read.csv(mf,stringsAsFactors=FALSE,colClasses=c(role='character',sd='numeric',key='character',fit='character',
      fit_md5='character',reasons='character'))
    flag_rows <- if(is.null(flag_rows)) m else rbind(flag_rows,m[names(flag_rows)])
  }
  list(selected=selected,first=first,missing=missing,control_schedule=control_schedule,flag_rows=flag_rows,
    sources=c(file.path(sel_dir,c('selected-fits.csv','long-selection.csv','repeat-flags.csv','long-fit-flags.csv','control-flags.csv',
      'convergence.csv','manual-missing.csv','manual-selected.csv')),mf))
}

# ---------------------------------------------------------------------------
# Across-community tables and the gate, recomputed from this audit's own values
# ---------------------------------------------------------------------------

VERIFY_GATE_TOL <- 1e-9   # the accepted "at most" allowance (R30), restated

# Criteria 1 to 4 of one new arm (README "Gate", AMENDMENT-1 R21, AMENDMENT-2
# R27): vals holds this audit's primary per-fit groups of both arms (selected
# fits of valid communities only), flags its recomputed flags per selected fit.
audit_gate <- function(phase,vals,flags,sd,with_c4=TRUE) {
  strata <- if(phase=='A1') c('n100','n300') else 'all'
  gated <- if(phase=='A1') 'low' else c('low','rare_1_5pct')
  keys <- phase_keys(phase);ctl <- vals[vals$sd==1,,drop=FALSE];new <- vals[vals$sd==sd,,drop=FALSE]
  out <- list()
  add <- function(criterion,component,stratum,communities,needed,improved,value_new,value_control,statistic,missing_new,pass,decisive)
    out[[length(out)+1L]] <<- data.frame(criterion=criterion,component=component,stratum=stratum,communities=as.integer(communities),
      needed=as.integer(needed),improved=as.integer(improved),value_new=as.numeric(value_new),value_control=as.numeric(value_control),
      statistic=as.numeric(statistic),missing_new=as.integer(missing_new),pass=as.logical(pass),tolerance_decisive=as.logical(decisive),
      stringsAsFactors=FALSE)
  for(g in gated) {
    value <- function(d) {
      d <- d[d$group==g & d$cells>0,,drop=FALSE];u <- unique(d$community)
      v <- vapply(u,function(c) {z <- d[d$community==c,,drop=FALSE]
        if(nrow(z)==length(strata) && setequal(z$stratum,strata)) mean(z$abs_signed_error) else NA_real_},numeric(1))
      stats::setNames(v,u)[!is.na(v)]
    }
    a <- value(new);b <- value(ctl);ids <- intersect(names(a),names(b));n <- length(ids)
    need <- ceiling(2*n/3);imp <- sum(a[ids]<b[ids]);mn <- if(n) mean(a[ids]) else NA;mc <- if(n) mean(b[ids]) else NA
    add('c1',g,if(phase=='A1') 'pooled' else 'all',n,need,imp,mn,mc,mn-mc,NA,if(n==0) NA else imp>=need && mn<mc,NA)
  }
  at_most <- function(criterion,component,h,d,n,a,b,limit) add(criterion,component,h,n,NA,NA,a,b,d,NA,
    if(n==0) NA else d<=limit+VERIFY_GATE_TOL,if(n==0) NA else d>limit && d<=limit+VERIFY_GATE_TOL)
  both <- function(group,h,column) {
    x <- new[new$group==group & new$stratum==h & new$cells>0,,drop=FALSE];y <- ctl[ctl$group==group & ctl$stratum==h & ctl$cells>0,,drop=FALSE]
    k <- intersect(x$key,y$key);list(n=length(k),a=if(length(k)) mean(x[[column]][match(k,x$key)]) else NA,
      b=if(length(k)) mean(y[[column]][match(k,y$key)]) else NA)
  }
  for(h in strata) {
    for(band in c('low','middle','high')) {p <- both(band,h,'abs_signed_error');at_most('c2',band,h,p$a-p$b,p$n,p$a,p$b,1)}
    p <- both('all',h,'mean_abs_cell_error');at_most('c2','mae',h,p$a-p$b,p$n,p$a,p$b,.5)
  }
  for(h in strata) for(band in c('low','middle','high')) {p <- both(band,h,'coverage');at_most('c3',band,h,p$b-p$a,p$n,p$a,p$b,.03)}
  if(with_c4) for(h in strata) {
    expected <- keys[audit_stratum(phase,keys)==h];valid <- unique(new$key[new$stratum==h])
    fn <- if(length(valid)) sum(flags$flagged[flags$sd==sd & flags$key %in% valid]) else NA
    fc <- sum(flags$flagged[flags$sd==1 & flags$stratum==h]);miss <- sum(!expected %in% valid)
    add('c4','flags',h,NA,NA,NA,fn,fc,fn-fc,miss,miss==0 && isTRUE(fn<=fc),NA)
  }
  x <- do.call(rbind,out);rownames(x) <- NULL
  x$result <- if(!with_c4) 'DESCRIPTIVE' else if(any(!x$pass,na.rm=TRUE)) 'FAIL' else if(anyNA(x$pass)) 'UNDETERMINED' else 'PASS'
  x
}

# The documented gate.csv columns of one SD, from a criterion table.
audit_gate_wide <- function(d) {
  v <- list()
  for(i in seq_len(nrow(d))) {
    r <- d[i,]
    p <- switch(r$criterion,c1=paste0('c1_',r$component,'_'),c2=,c3=paste0(r$criterion,'_',r$component,'_',r$stratum,'_'),c4=paste0('c4_',r$stratum,'_'))
    vals <- switch(r$criterion,c1=list(communities=r$communities,improved=r$improved,needed=r$needed,mean_new=r$value_new,mean_control=r$value_control,pass=r$pass),
      c2=list(new=r$value_new,control=r$value_control,diff=r$statistic,pass=r$pass),c3=list(new=r$value_new,control=r$value_control,drop=r$statistic,pass=r$pass),
      c4=list(flagged_new=r$value_new,flagged_control=r$value_control,missing_new=r$missing_new,pass=r$pass))
    for(k in names(vals)) v[[paste0(p,k)]] <- vals[[k]]
  }
  for(k in c('c1','c2','c3','c4')) v[[paste0(k,'_pass')]] <- if(any(d$criterion==k)) all(d$pass[d$criterion==k]) else NA
  v$result <- d$result[1]
  v
}

same_value <- function(a,b,tol=VERIFY_TOL_POINTS) {
  if(length(a)!=1L || length(b)!=1L) return(FALSE)
  if(is.na(a) || is.na(b)) return(is.na(a) && is.na(b))
  if(is.character(a) || is.character(b) || is.logical(a) || is.logical(b)) return(identical(as.character(a),as.character(b)) ||
    (is.logical(a) || is.logical(b)) && identical(as.logical(a),as.logical(b)))
  abs(as.numeric(a)-as.numeric(b))<=tol
}

# One check row per compared value; a missing or extra row is one failing check.
compare_rows <- function(table,mine,theirs,keys,cols) {
  id <- function(x) if(nrow(x)) do.call(paste,c(unname(as.list(x[keys])),sep=' ')) else character()
  m <- id(mine);t <- id(theirs)
  rows <- list(data.frame(table=table,item='row set',own=paste(sort(m),collapse='; '),stored=paste(sort(t),collapse='; '),
    pass=setequal(m,t) && !anyDuplicated(m) && !anyDuplicated(t),stringsAsFactors=FALSE))
  for(i in seq_along(m)) {j <- match(m[i],t);if(is.na(j)) next
    for(k in cols) rows[[length(rows)+1L]] <- data.frame(table=table,item=paste(m[i],k),own=format(mine[[k]][i],digits=17),
      stored=format(theirs[[k]][j],digits=17),pass=same_value(mine[[k]][i],theirs[[k]][j]),stringsAsFactors=FALSE)}
  do.call(rbind,rows)
}

audit_means <- function(g,by,cols,fns) {
  if(!nrow(g)) return(g[0,c(by,'communities',names(cols))])
  s <- split(g,do.call(interaction,c(unname(as.list(g[by])),drop=TRUE)))
  x <- do.call(rbind,lapply(s,function(z) {r <- z[1,by,drop=FALSE];r$communities <- nrow(z)
    for(k in names(cols)) r[[k]] <- mean(z[[cols[[k]]]]);r}))
  rownames(x) <- NULL;x
}

# summary-means.csv, b0-means.csv, convergence.csv, gate-detail.csv, gate.csv
# and schedule-matched.csv against this audit's own per-community values.
audit_gate_tables <- function(phase,arms,own,own_b0,flags,sel,tabs) {
  checks <- list();gates <- list()
  sel_g <- own[own$kind=='selected',,drop=FALSE]
  sm <- audit_means(sel_g,c('sd','stratum','scope','group'),c(mean_signed_error='signed_error',mean_abs_signed_error='abs_signed_error',
    mean_abs_cell_error='mean_abs_cell_error',mean_coverage='coverage'))
  checks[[1]] <- compare_rows('summary-means.csv',sm,tabs$means,c('sd','stratum','scope','group'),
    c('communities','mean_signed_error','mean_abs_signed_error','mean_abs_cell_error','mean_coverage'))
  bm <- audit_means(own_b0[own_b0$kind=='selected',,drop=FALSE],c('sd','stratum'),c(mean_b0_bias='b0_bias',mean_b0_abs_bias='b0_abs_bias',
    mean_b0_coverage='b0_coverage'))
  checks[[2]] <- compare_rows('b0-means.csv',bm,tabs$b0means,c('sd','stratum'),c('communities','mean_b0_bias','mean_b0_abs_bias','mean_b0_coverage'))
  keys <- phase_keys(phase);strata <- unique(audit_stratum(phase,keys))
  cv <- do.call(rbind,lapply(arms,function(sd) do.call(rbind,lapply(strata,function(h) {
    f <- flags[flags$sd==sd & flags$stratum==h,,drop=FALSE];n <- nrow(f)
    firsts <- own[own$kind=='first' & own$sd==sd & own$stratum==h,c('key','flagged'),drop=FALSE];firsts <- unique(firsts)
    first_flag <- ifelse(f$key %in% firsts$key,firsts$flagged[match(f$key,firsts$key)],f$flagged)
    data.frame(sd=sd,stratum=h,selected_fits=n,selected_flagged=if(n) sum(f$flagged) else NA,
      first_flagged=if(sd==1 || !n) NA else sum(first_flag),long_repeats=if(sd==1 || !n) NA else nrow(firsts),
      missing_fits=sum(audit_stratum(phase,keys)==h)-n,stringsAsFactors=FALSE)}))))
  checks[[3]] <- compare_rows('convergence.csv',cv,tabs$convergence,c('sd','stratum'),
    c('selected_fits','selected_flagged','first_flagged','long_repeats','missing_fits'))
  new_sds <- setdiff(arms,1)
  if(1 %in% arms && length(new_sds)) {
    primary <- sel_g[sel_g$scope=='primary',,drop=FALSE]
    for(sd in new_sds) {
      d <- audit_gate(phase,primary[primary$sd %in% c(1,sd),,drop=FALSE],flags,sd);d$sd <- sd;gates[[length(gates)+1L]] <- d
      theirs <- tabs$detail[tabs$detail$sd==sd,,drop=FALSE]
      checks[[length(checks)+1L]] <- compare_rows('gate-detail.csv',d,theirs,c('criterion','component','stratum'),
        c('communities','needed','improved','value_new','value_control','statistic','missing_new','pass','tolerance_decisive'))
      w <- audit_gate_wide(d);row <- tabs$gate[tabs$gate$sd==sd,,drop=FALSE]
      stored_cols <- setdiff(names(tabs$gate),c('sd','subphase'))
      if(nrow(row)!=1L) checks[[length(checks)+1L]] <- data.frame(table='gate.csv',item=paste('row of SD',sd),own='1 row',
        stored=paste(nrow(row),'rows'),pass=FALSE,stringsAsFactors=FALSE) else {
        cs <- data.frame(table='gate.csv',item=paste('SD',sd,'columns'),own=paste(sort(names(w)),collapse=' '),
          stored=paste(sort(stored_cols[!vapply(stored_cols,function(k) is.na(row[[k]]) && !k %in% names(w),logical(1))]),collapse=' '),
          pass=setequal(names(w),stored_cols[!vapply(stored_cols,function(k) is.na(row[[k]]) && !k %in% names(w),logical(1))]),stringsAsFactors=FALSE)
        checks[[length(checks)+1L]] <- rbind(cs,do.call(rbind,lapply(names(w),function(k) data.frame(table='gate.csv',item=paste('SD',sd,k),
          own=format(w[[k]],digits=17),stored=if(k %in% names(row)) format(row[[k]],digits=17) else 'absent',
          pass=k %in% names(row) && same_value(w[[k]],row[[k]]),stringsAsFactors=FALSE))))
      }
      # Schedule-matched sensitivity: first fits in place of their longer repeats.
      fr <- own[own$kind=='first' & own$sd==sd & own$scope=='primary',,drop=FALSE]
      sm_vals <- rbind(primary[primary$sd==1 | (primary$sd==sd & !primary$key %in% fr$key),,drop=FALSE],fr)
      s <- audit_gate(phase,sm_vals[sm_vals$sd %in% c(1,sd),,drop=FALSE],flags,sd,with_c4=FALSE)
      checks[[length(checks)+1L]] <- compare_rows('schedule-matched.csv',s,tabs$matched[tabs$matched$sd==sd,,drop=FALSE],
        c('criterion','component','stratum'),c('communities','needed','improved','value_new','value_control','statistic','pass','tolerance_decisive'))
    }
  }
  list(checks=do.call(rbind,checks),gate=if(length(gates)) do.call(rbind,gates) else NULL)
}

# The package's native kernel basis (the check score.R made against its own
# basis), from the study's fingerprinted library.
native_basis_check <- function(fit,bases) {
  native <- suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid,fit$infos$list_Xs))
  max(vapply(seq_along(bases),function(g) {
    ref <- occJSDM:::KsBproduct(native$Ks_all[,,g],diag(fit$infos$ps),fit$infos$list_Xs$Xs_centers)
    if(!identical(dim(ref),dim(bases[[g]]))) return(Inf)
    max(abs(ref-bases[[g]]))
  },numeric(1)))
}

audit_source <- function(out,phase,files,repo,study) {
  files <- unique(files[file.exists(files)])
  rel <- function(f) {f <- normalizePath(f);for(root in c(repo,study)) {p <- paste0(normalizePath(root),'/')
    if(startsWith(f,p)) return(substring(f,nchar(p)+1L))};f}
  x <- data.frame(file=vapply(files,rel,''),md5=unname(tools::md5sum(files)),stringsAsFactors=FALSE)
  utils::write.csv(x,file.path(out,phase,'audit-source.csv'),row.names=FALSE)
}

# Phase A: every SD passes phase A only if this audit's own A1 and A2 results
# are PASS (R9); the stored summary/A/gate.csv must say the same.
verify_phase_a <- function(repo,study,summary,out) {
  own <- lapply(c('A1','A2'),function(p) {
    f <- file.path(out,p,'gate-audit-arms-1-2-3-5.csv');v <- file.path(out,p,'verify-arms-1-2-3-5.csv')
    if(!file.exists(f) || !file.exists(v)) stop('Run verify.R for phase ',p,' with --arms=1,2,3,5 first')
    if(!all(utils::read.csv(v)$pass) || !all(utils::read.csv(file.path(out,p,'tables-arms-1-2-3-5.csv'))$pass)) stop('The phase ',p,' audit did not pass')
    g <- utils::read.csv(f,stringsAsFactors=FALSE);unique(g[c('sd','result')])})
  sds <- sort(own[[1]]$sd);if(!setequal(sds,own[[2]]$sd)) stop('The A1 and A2 audits cover different SDs')
  r1 <- own[[1]]$result[match(sds,own[[1]]$sd)];r2 <- own[[2]]$result[match(sds,own[[2]]$sd)];pass <- r1=='PASS' & r2=='PASS'
  mine <- data.frame(sd=sds,a1_result=r1,a2_result=r2,result=ifelse(pass,'PASS','FAIL'),
    smallest_passing=pass & sds==(if(any(pass)) min(sds[pass]) else -Inf),stringsAsFactors=FALSE)
  stored <- utils::read.csv(file.path(summary,'A','gate.csv'),stringsAsFactors=FALSE)
  st <- stored[stored$subphase=='A',,drop=FALSE]
  checks <- rbind(compare_rows('A/gate.csv',mine,st,'sd',c('a1_result','a2_result','result','smallest_passing')),
    compare_rows('A/gate.csv sub-phase rows',rbind(data.frame(sd=sds,subphase='A1',result=r1),data.frame(sd=sds,subphase='A2',result=r2)),
      stored[stored$subphase!='A',,drop=FALSE],c('sd','subphase'),'result'))
  dir.create(file.path(out,'A'),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(checks,file.path(out,'A','verify-phase-a.csv'),row.names=FALSE)
  audit_source(out,'A',c(file.path(repo,'dev/simstudy/occupancy-intercept-prior',c('verify.R','jobs.R')),file.path(summary,'A','gate.csv'),
    file.path(out,rep(c('A1','A2'),each=3),c('gate-audit-arms-1-2-3-5.csv','verify-arms-1-2-3-5.csv','tables-arms-1-2-3-5.csv'))),repo,study)
  cat('Phase A audit:',sum(checks$pass),'of',nrow(checks),'checks agree; wrote',file.path(out,'A','verify-phase-a.csv'),'\n')
  if(all(checks$pass)) 0L else 1L
}

verify_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','phase','arms','workers','summary','out'),
    required=c('repo','study','phase'))
  phase <- o$phase;if(!phase %in% c('A1','A2','A')) stop('--phase must be A1, A2 or A')
  repo <- normalizePath(o$repo,mustWork=TRUE);study <- normalizePath(o$study,mustWork=TRUE)
  summary <- o$summary %||% file.path(study,'summary');out <- o$out %||% file.path(study,'verify')
  if(phase=='A') return(verify_phase_a(repo,study,summary,out))
  if(is.null(o$arms) || is.null(o$`inputs-root`)) stop('Missing --arms or --inputs-root')
  arms <- strsplit(o$arms,',',fixed=TRUE)[[1]]
  if(!length(arms) || !all(arms %in% c('1','2','3','5')) || anyDuplicated(arms)) stop('--arms must list arms among 1, 2, 3, 5')
  arms <- sort(as.numeric(arms));workers <- parse_workers(o$workers %||% '1')
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE);inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  cc <- c(key='character',fit='character',fit_md5='character',community='character')
  tf <- function(name) file.path(summary,phase,name)
  read <- function(name,optional=FALSE) {if(optional && !file.exists(tf(name))) return(NULL)
    x <- utils::read.csv(tf(name),stringsAsFactors=FALSE,colClasses=cc[intersect(names(cc),names(utils::read.csv(tf(name),nrows=1)))])
    x[x$sd %in% arms,,drop=FALSE]}
  tabs <- list(band=read('band-error.csv'),coverage=read('coverage.csv'),mae=read('mae.csv'),b0=read('b0.csv'),
    means=read('summary-means.csv'),b0means=read('b0-means.csv'),convergence=read('convergence.csv'))
  new_sds <- setdiff(arms,1)
  if(1 %in% arms && length(new_sds)) {tabs$detail <- read('gate-detail.csv');tabs$gate <- read('gate.csv');tabs$matched <- read('schedule-matched.csv')}
  fits <- unique(tabs$band[c('sd','arm','kind','key','schedule','fit','fit_md5')]);rownames(fits) <- NULL
  sel <- audit_selection(phase,arms,repo,study,summary)
  got <- fits[fits$kind=='selected',];id <- function(x) paste(x$sd,x$key,x$fit,x$fit_md5)
  if(!setequal(id(got),id(sel$selected)) || anyDuplicated(id(got))) stop('The summary does not hold exactly the selected fits of the requested arms')
  gf <- fits[fits$kind=='first',];fid <- function(x) if(is.null(x) || !nrow(x)) character() else paste(x$sd,x$key,x$fit)
  if(!setequal(fid(gf),fid(sel$first))) stop('The summary does not hold exactly the first fits of the longer repeats')
  if(nrow(sel$missing %||% data.frame())) cat('R27: communities without a valid selected fit:',paste(paste0('sd',sel$missing$sd),sel$missing$key,collapse=', '),'\n')
  native <- NULL
  if(phase=='A2' && nrow(fits)) {
    check_library_fingerprint(study,file.path(repo,FINGERPRINT_FILE))
    .libPaths(c(file.path(study,'library'),.libPaths()));suppressPackageStartupMessages(library(occJSDM))
    native <- native_basis_check
  }
  cat('Auditing',nrow(fits),'fits of phase',phase,'\n');flush.console()
  work <- function(i) tryCatch({r <- audit_fit(fits[i,,drop=FALSE],phase,study,archives,inputs_root,tabs,sel$flag_rows,native)
      cat(format(Sys.time()),fits$key[i],paste0('sd',fits$sd[i]),fits$kind[i],if(r$row$pass) 'verified' else 'DISAGREES','\n');flush.console();r},
    error=function(e) structure(conditionMessage(e),class='audit-error'))
  res <- if(workers==1L) lapply(seq_len(nrow(fits)),work) else
    parallel::mclapply(seq_len(nrow(fits)),work,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  bad <- !vapply(res,is.list,logical(1)) | vapply(res,inherits,logical(1),what='audit-error')
  if(any(bad)) {for(i in which(bad)) cat('ERROR',fits$key[i],paste0('sd',fits$sd[i]),':',as.character(res[[i]]),'\n');return(1L)}
  x <- do.call(rbind,lapply(res,`[[`,'row'))
  own <- do.call(rbind,lapply(res,function(r) cbind(r$groups,flagged=r$flagged)))
  own_b0 <- do.call(rbind,lapply(res,`[[`,'b0'))
  flags <- do.call(rbind,lapply(res[fits$kind=='selected'],function(r) data.frame(sd=r$row$sd,key=r$row$key,
    stratum=audit_stratum(phase,r$row$key),flagged=r$flagged,stringsAsFactors=FALSE)))
  ag <- audit_gate_tables(phase,arms,own,own_b0,flags,sel,tabs)
  dir.create(file.path(out,phase),recursive=TRUE,showWarnings=FALSE);tag <- paste(arms,collapse='-')
  dest <- file.path(out,phase,paste0('verify-arms-',tag,'.csv'));utils::write.csv(x,dest,row.names=FALSE)
  tdest <- file.path(out,phase,paste0('tables-arms-',tag,'.csv'));utils::write.csv(ag$checks,tdest,row.names=FALSE)
  if(!is.null(ag$gate)) utils::write.csv(ag$gate,file.path(out,phase,paste0('gate-audit-arms-',tag,'.csv')),row.names=FALSE)
  audit_source(out,phase,c(file.path(repo,'dev/simstudy/occupancy-intercept-prior',c('verify.R','jobs.R','results/control-provenance.csv')),
    file.path(summary,phase,c('band-error.csv','coverage.csv','mae.csv','b0.csv','summary-means.csv','b0-means.csv','convergence.csv',
      'gate-detail.csv','gate.csv','schedule-matched.csv')),sel$sources),repo,study)
  m <- function(k) if(all(is.na(x[[k]]))) 'not applicable' else format(max(x[[k]],na.rm=TRUE),digits=3)
  cat('Verified',sum(x$pass),'of',nrow(x),'fits. Largest differences: signed error',m('max_diff_signed_error'),
    '; mean absolute error',m('max_diff_mean_abs_cell_error'),'; MAE',m('max_diff_mae'),'; coverage',m('max_diff_coverage'),
    '; B0',m('max_diff_b0'),'; selection Rhat',m('sel_group_rhat_diff'),m('sel_element_rhat_diff'),'; archive Rhat',m('archive_rhat_vector_diff'),
    '; archive ESS (relative)',m('archive_ess_relative_diff'),'\n')
  cat('Tables:',sum(ag$checks$pass),'of',nrow(ag$checks),'checks agree across',paste(unique(ag$checks$table),collapse=', '),'\n')
  cat('Wrote',dest,tdest,file.path(out,phase,'audit-source.csv'),sep='\n  ');cat('\n')
  if(!all(x$pass) || !all(ag$checks$pass)) {
    if(!all(x$pass)) print(x[!x$pass,c('sd','key','kind')],row.names=FALSE)
    if(!all(ag$checks$pass)) print(utils::head(ag$checks[!ag$checks$pass,],20),row.names=FALSE)
    return(1L)
  }
  0L
}

if(sys.nframe()==0L) {
  here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
  source(file.path(here,'jobs.R'))
  status <- tryCatch(verify_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
