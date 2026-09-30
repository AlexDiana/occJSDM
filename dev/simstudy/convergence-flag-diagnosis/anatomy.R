# Chain anatomy from saved posterior draws: Task 1 of PLAN.md in this
# directory. It reads saved runOccJSDM fits and their simulation inputs and
# never fits a model; archives are read, never written.
#
# Interfaces
#
#   chain_anatomy(fit_paths, input_path, keep_draws = FALSE)
#     One row per species, chain and quantity (ANATOMY_QUANTITIES), with
#     columns key, species (index 1..S), chain, quantity, chain_mean, chain_sd,
#     split_rhat_within_chain (posterior::rhat of that chain alone, so the
#     rank-normalised split-Rhat of its two halves) and truth. Two further
#     columns, rhat_all_chains and ess_bulk_all_chains (posterior::rhat and
#     posterior::ess_bulk over all chains, constant within a key, species and
#     quantity), let chain_separation() work from the table alone; fit_file and
#     file_chain record which saved file and chain each row came from. Several
#     fit files of one key are read one at a time and their chains concatenated
#     in the given order, renumbered 1..K (ruling R4). Attributes: psi_mean,
#     the reconstructed posterior mean occupancy probability (sites x species,
#     over all files); psi_check, per file the largest absolute difference
#     between the reconstruction and that file's stored psi_output (the call
#     stops if any reaches PSI_TOLERANCE); original_sites; and draws, only when
#     keep_draws is TRUE, one species x iteration x chain array per quantity.
#
#   chain_separation(anatomy)
#     A list of three data frames (ruling R2 asks that the form be recorded):
#     $separation, per key, species and quantity: gap (range of the chain
#       means), pooled_within_sd (root mean square of the chain SDs),
#       separation (gap / pooled_within_sd), max_split_rhat_within_chain, rhat
#       and ess_bulk over all chains, fixed (all draws constant, as for the
#       loadings the package fixes at 0 or 1), and label: 'drifting' if some
#       chain's own split-Rhat is above SPLIT_RHAT_LIMIT, otherwise
#       'separated' if separation is above SEPARATION_LIMIT, otherwise
#       'agrees'.
#     $species, per key and species: label from the non-loading quantities
#       only (ruling R6: 'separated' if any is separated, else 'drifting' if
#       any is drifting, else 'agrees'), the separated and drifting
#       non-loading quantities, and the loadings' own labels, reported but not
#       deciding.
#     $chain_groups, for each key and species with any separated quantity
#       (loadings included, as R2 reads; filter on species_label to keep only
#       species separated on identifiable quantities): the chains sorted by
#       their mean_psi_original_sites chain means and split at the largest gap
#       into chain_group 1 (lower occupancy) and 2 (higher).
#
# Quantities and truths. Truths follow the pr11 scorer, score_fit() in
# nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R (lines 60-89),
# including its checks of the design matrices and the true linear predictor.
#   B0, B_slope1, B_slope2: occupancy intercept and slopes on the standardised
#     covariates; truth jsdmParams_true$B0 and rows of jsdmParams_true$B.
#   theta0: field false-positive probability; truth input$truth$params$theta0.
#   beta_theta_intercept, beta_theta_slope: collection coefficients on the
#     standardised collection covariate; truth the generating intercept plus
#     the mean raw covariate times the generating slope, and the SD of the raw
#     covariate times the generating slope (score_fit lines 78-85).
#   p_primer1/2, q_primer1/2: PCR detection and false-positive probabilities;
#     truth p_true and q_true times the retention probability (a positive PCR
#     keeps at least one read), as score_fit line 86-90 scores them.
#   mean_psi_original_sites: per draw, the mean over the original sites of the
#     occupancy probability; truth the same mean of plogis(true eta). The
#     original sites are 1..input$design$original_sites (100), the pr11
#     scorer's original-site set (make_design_scorer in
#     nonspatial-design-recheck/design_helpers.R; verify.R lines 39-45).
#   L1, L2: latent-factor loadings; truth jsdmParams_true$L when the generating
#     loadings satisfy the package's constraint (L[1,1] = L[2,2] = 1,
#     L[2,1] = 0), otherwise NA.
#
# Occupancy reconstruction (ruling R1): reimplemented from
# current-main-recheck/verify.R lines 29-38, which is a command-line script
# with no function to source. Test (d) in test-anatomy.R is the proof: the
# per-draw probabilities average to the stored psi_output within 1e-10.
#
# Reused by parsing, not copied: extract_functions() from
# current-main-recheck/helpers.R, and ORIGINAL_INPUT_ROOT, remap_input() and
# checked_input() from occupancy-intercept-prior/jobs.R (the input path remap
# and md5 refusal). Their md5s are written to results/anatomy/source-hashes.csv.
#
# Run `Rscript anatomy.R` to write results/anatomy/ for the five flagged and
# three reference fits; tests: `Rscript test-anatomy.R`.

`%||%` <- function(x,y) if(is.null(x)) y else x

ANATOMY_DIR <- local({
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
REPO <- normalizePath(file.path(ANATOMY_DIR,'../../..'))
ARCHIVES <- Sys.getenv('OCCJSDM_ARCHIVES','/Users/douglasyu/src/occJSDM/dev/simstudy/results')
PR11_ARCHIVE <- file.path(ARCHIVES,'pr11-current-20260927')
INPUTS_ROOT <- file.path(ARCHIVES,'intercept-prior-inputs')
RECHECK_RESULTS <- file.path(REPO,'dev/simstudy/current-main-recheck/results')
OUTPUT_DIR <- file.path(ANATOMY_DIR,'results/anatomy')

REUSED_SOURCES <- c(helpers=file.path(REPO,'dev/simstudy/current-main-recheck/helpers.R'),
  jobs=file.path(REPO,'dev/simstudy/occupancy-intercept-prior/jobs.R'))
.reused_helpers <- new.env(parent=globalenv())
sys.source(REUSED_SOURCES[['helpers']],envir=.reused_helpers)
.reused_helpers$extract_functions(REUSED_SOURCES[['jobs']],
  c('ORIGINAL_INPUT_ROOT','remap_input','checked_input'),environment())
rm(.reused_helpers)

ANATOMY_QUANTITIES <- c('B0','theta0','beta_theta_intercept','beta_theta_slope',
  'p_primer1','p_primer2','q_primer1','q_primer2','B_slope1','B_slope2',
  'mean_psi_original_sites','L1','L2')
LOADING_QUANTITIES <- c('L1','L2')
SPLIT_RHAT_LIMIT <- 1.05
SEPARATION_LIMIT <- 3
PSI_TOLERANCE <- 1e-10
FLAGGED_KEYS <- c('design-qfar_K6-sites300-05','design-qnear_K6-sites300-02',
  'design-qfar_K6-sites300-07','design-qfar_K6-sites300-09','design-qfar_K6-field4-09')
REFERENCE_KEYS <- c('design-qfar_K6-sites300-01','design-qfar_K6-sites300-02',
  'design-qfar_K6-sites300-03')

# ---- Selected fits and inputs --------------------------------------------

# The pr11 selection: schedule 'long' when long-selection.csv says run_long,
# else 'initial' (current-main-recheck/select.R line 35), cross-checked against
# selected-manifest.csv.
selection_row <- function(key) {
  plan <- read.csv(file.path(RECHECK_RESULTS,'long-selection.csv'),stringsAsFactors=FALSE)
  manifest <- read.csv(file.path(RECHECK_RESULTS,'selected-manifest.csv'),stringsAsFactors=FALSE)
  p <- plan[plan$key==key,,drop=FALSE];m <- manifest[manifest$key==key,,drop=FALSE]
  if(nrow(p)!=1L || nrow(m)!=1L) stop('Key not found exactly once in the pr11 selection: ',key)
  schedule <- if(isTRUE(p$run_long)) 'long' else 'initial'
  if(!identical(m$schedule,schedule))
    stop('Selected schedule disagrees between long-selection.csv and selected-manifest.csv for ',key)
  m
}

# select.R line 52: a selected fit is flagged if it has a package warning, an
# unresolved Rhat or a scored Rhat above 1.05.
selection_flagged <- function(key) {
  m <- selection_row(key)
  m$warnings>0 || m$unresolved_rhat>0 || m$max_group_rhat>1.05 || m$max_element_rhat>1.05
}

selected_fit <- function(key) {
  m <- selection_row(key)
  f <- file.path(PR11_ARCHIVE,m$schedule,paste0(key,'-fit.rds'))
  if(!identical(basename(sub('-result\\.rds$','-fit.rds',m$result_file)),basename(f)))
    stop('Selected result file does not name the expected fit for ',key)
  if(!file.exists(f)) stop('Missing fit: ',f)
  f
}

selected_input <- function(key) {
  m <- selection_row(key)
  checked_input(remap_input(m$input_file,INPUTS_ROOT),m$input_md5)
}

# ---- Truths ----------------------------------------------------------------

anatomy_truth <- function(fit,input) {
  sim <- input$sim;tp <- sim$true_params;jp <- tp$jsdmParams_true
  S <- length(jp$B0)
  stopifnot(fit$infos$ps==0,identical(fit$infos$model,'two_stage'),
    identical(fit$infos$speciesNames,colnames(sim$data_list$OTU)),S==length(fit$infos$speciesNames))
  info <- sim$data_list$info
  site_info <- info[!duplicated(info$Site),,drop=FALSE]
  site_info <- site_info[order(site_info$Site),,drop=FALSE]
  rawx <- as.matrix(site_info[paste0('X_psi.EnvCov.',seq_len(ncol(fit$X_psi)))])
  stopifnot(max(abs(unname(scale(rawx))-unname(fit$X_psi)))<1e-10)
  true_eta <- sweep(fit$X_psi%*%jp$B+jp$U%*%jp$L,2,jp$B0,'+')
  stopifnot(max(abs(true_eta-jp$eta))<1e-10)
  sample_info <- info[!duplicated(info[c('Site','Sample')]),,drop=FALSE]
  sample_info <- sample_info[order(sample_info$Site,sample_info$Sample),,drop=FALSE]
  rawt <- sample_info$X_theta
  stopifnot(max(abs(unname(cbind(1,scale(rawt)))-unname(fit$X_theta)))<1e-10)
  bt <- tp$beta_theta_true
  bt_fit <- rbind(bt[1,]+mean(rawt)*bt[2,],sd(rawt)*bt[2,])
  stopifnot(max(abs(plogis(fit$X_theta%*%bt_fit)-plogis(cbind(1,rawt)%*%bt)))<1e-10)
  retention <- c(p=pnorm(log(1.5),5,1,lower.tail=FALSE),q=pnorm(log(1.5),1.5,1,lower.tail=FALSE))
  n0 <- input$design$original_sites
  stopifnot(identical(as.integer(n0),100L),nrow(fit$X_psi)>=n0)
  L <- jp$L
  constrained <- nrow(L)==2L && ncol(L)>=2L && L[1,1]==1 && L[2,1]==0 && L[2,2]==1
  loading <- function(k) if(constrained) L[k,] else rep(NA_real_,S)
  truth <- cbind(B0=jp$B0,theta0=input$truth$params$theta0,
    beta_theta_intercept=bt_fit[1,],beta_theta_slope=bt_fit[2,],
    p_primer1=tp$p_true[1,]*retention[['p']],p_primer2=tp$p_true[2,]*retention[['p']],
    q_primer1=tp$q_true[1,]*retention[['q']],q_primer2=tp$q_true[2,]*retention[['q']],
    B_slope1=jp$B[1,],B_slope2=jp$B[2,],
    mean_psi_original_sites=colMeans(plogis(true_eta[seq_len(n0),,drop=FALSE])),
    L1=loading(1L),L2=loading(2L))
  stopifnot(identical(colnames(truth),ANATOMY_QUANTITIES),nrow(truth)==S)
  list(truth=truth,original_sites=seq_len(n0),loadings_constrained=constrained)
}

# ---- Draws of one saved fit ------------------------------------------------

# Species x iteration x chain arrays of every quantity, and the per-site
# occupancy sums used to check the reconstruction against psi_output.
fit_draws <- function(fit,original_sites) {
  ro <- fit$results_output;jo <- ro$jsdm_output
  dd <- dim(jo$B0_output);S <- dd[1];ni <- dd[2];nc <- dd[3]
  n <- nrow(fit$X_psi);ncov <- ncol(fit$X_psi);d <- dim(jo$U_output)[2]
  stopifnot(length(dd)==3L,ncov==2L,d==2L,identical(dim(jo$B_output),c(ncov,S,ni,nc)),
    identical(dim(jo$U_output),c(n,d,ni,nc)),identical(dim(jo$L_output),c(d,S,ni,nc)),
    identical(dim(ro$theta0_output),dd),identical(dim(ro$beta_theta_output),c(2L,S,ni,nc)),
    identical(dim(ro$p_output),c(2L,S,ni,nc)),identical(dim(ro$q_output),c(2L,S,ni,nc)),
    identical(dim(ro$psi_output),c(n,S)))
  row <- function(x,r) array(x[r,,,,drop=FALSE],dd)
  draws <- list(B0=jo$B0_output,theta0=ro$theta0_output,
    beta_theta_intercept=row(ro$beta_theta_output,1L),beta_theta_slope=row(ro$beta_theta_output,2L),
    p_primer1=row(ro$p_output,1L),p_primer2=row(ro$p_output,2L),
    q_primer1=row(ro$q_output,1L),q_primer2=row(ro$q_output,2L),
    B_slope1=row(jo$B_output,1L),B_slope2=row(jo$B_output,2L),
    L1=row(jo$L_output,1L),L2=row(jo$L_output,2L))
  psi_sum <- matrix(0,n,S);mean_psi <- array(NA_real_,dd)
  for(s in seq_len(S)) for(ch in seq_len(nc)) {
    beta <- matrix(jo$B_output[,s,,ch],nrow=ncov,ncol=ni)
    env <- sweep(fit$X_psi%*%beta,2,jo$B0_output[s,,ch],'+')
    hidden <- matrix(0,n,ni)
    for(k in seq_len(d))
      hidden <- hidden+sweep(matrix(jo$U_output[,k,,ch],n,ni),2,jo$L_output[k,s,,ch],'*')
    psi <- plogis(env+hidden)
    psi_sum[,s] <- psi_sum[,s]+rowSums(psi)
    mean_psi[s,,ch] <- colMeans(psi[original_sites,,drop=FALSE])
  }
  draws$mean_psi_original_sites <- mean_psi
  list(draws=draws[ANATOMY_QUANTITIES],psi_sum=psi_sum,ni=ni,nc=nc,
    psi_difference=max(abs(psi_sum/(ni*nc)-ro$psi_output)))
}

# ---- Per-chain rows ---------------------------------------------------------

# Rows for one species and quantity from an iterations x chains matrix.
anatomy_rows <- function(draws,key,species,quantity,truth=NA_real_,
  fit_file=NA_character_,file_chain=NULL) {
  draws <- as.matrix(draws);nc <- ncol(draws)
  within <- vapply(seq_len(nc),function(ch) posterior::rhat(draws[,ch,drop=FALSE]),numeric(1))
  data.frame(key=key,species=as.integer(species),chain=seq_len(nc),quantity=quantity,
    chain_mean=unname(colMeans(draws)),chain_sd=unname(apply(draws,2,sd)),
    split_rhat_within_chain=within,truth=unname(truth),
    rhat_all_chains=posterior::rhat(draws),ess_bulk_all_chains=posterior::ess_bulk(draws),
    fit_file=fit_file,file_chain=file_chain %||% seq_len(nc),
    stringsAsFactors=FALSE,row.names=NULL)
}

chain_anatomy <- function(fit_paths,input_path,keep_draws=FALSE) {
  fit_paths <- as.character(fit_paths)
  if(!length(fit_paths)) stop('No fit files given')
  absent <- fit_paths[!file.exists(fit_paths)]
  if(length(absent)) stop('Missing fit file: ',absent[1])
  if(!file.exists(input_path)) stop('Missing input: ',input_path)
  input_md5 <- unname(tools::md5sum(input_path))
  key <- NULL;input <- NULL;truth <- NULL;reference <- NULL
  parts <- list();checks <- list();psi_sum <- NULL
  for(i in seq_along(fit_paths)) {
    saved <- readRDS(fit_paths[i]);job <- saved$job;fit <- saved$fit;rm(saved)
    if(is.null(job$key)) stop('Fit file has no job key: ',fit_paths[i])
    if(is.null(key)) key <- job$key else if(!identical(job$key,key))
      stop('Fit files belong to different keys: ',key,' and ',job$key)
    if(!identical(job$input_md5,input_md5))
      stop('Input md5 mismatch for ',input_path,': the job of ',basename(fit_paths[i]),
        ' records ',job$input_md5,', found ',input_md5)
    if(!identical(job$family,'design')) stop('Only design-family (two-stage) fits are supported: ',key)
    if(is.null(input)) {
      input <- readRDS(input_path);truth <- anatomy_truth(fit,input)
      reference <- list(X_psi=fit$X_psi,X_theta=fit$X_theta,species=fit$infos$speciesNames)
    } else if(!identical(fit$X_psi,reference$X_psi) || !identical(fit$X_theta,reference$X_theta) ||
      !identical(fit$infos$speciesNames,reference$species))
      stop('Fit files of ',key,' differ in their design matrices or species')
    part <- fit_draws(fit,truth$original_sites);rm(fit);invisible(gc(FALSE))
    if(!is.finite(part$psi_difference) || part$psi_difference>=PSI_TOLERANCE)
      stop('Reconstructed occupancy disagrees with the stored psi_output of ',
        basename(fit_paths[i]),': ',format(part$psi_difference))
    psi_sum <- if(is.null(psi_sum)) part$psi_sum else psi_sum+part$psi_sum
    checks[[i]] <- data.frame(fit_file=basename(fit_paths[i]),n_chains=part$nc,
      n_iterations=part$ni,psi_max_abs_diff=part$psi_difference,stringsAsFactors=FALSE)
    parts[[i]] <- part$draws;rm(part)
  }
  checks <- do.call(rbind,checks)
  ni <- unique(checks$n_iterations)
  if(length(ni)!=1L) stop('Fit files of ',key,' hold different numbers of retained iterations')
  K <- sum(checks$n_chains);S <- nrow(truth$truth)
  # Chains concatenated in the given file order (R4): the chain dimension is
  # the slowest, so joining the arrays' values appends chains.
  draws <- lapply(stats::setNames(ANATOMY_QUANTITIES,ANATOMY_QUANTITIES),function(q)
    array(unlist(lapply(parts,`[[`,q),use.names=FALSE),c(S,ni,K)))
  rm(parts)
  fit_file <- rep(checks$fit_file,checks$n_chains)
  file_chain <- unlist(lapply(checks$n_chains,seq_len))
  rows <- list()
  for(s in seq_len(S)) for(q in ANATOMY_QUANTITIES)
    rows[[length(rows)+1L]] <- anatomy_rows(matrix(draws[[q]][s,,],ni,K),key,s,q,
      truth$truth[s,q],fit_file,file_chain)
  out <- do.call(rbind,rows);rownames(out) <- NULL
  attr(out,'psi_mean') <- psi_sum/(ni*K)
  attr(out,'psi_check') <- checks
  attr(out,'original_sites') <- truth$original_sites
  attr(out,'loadings_constrained') <- truth$loadings_constrained
  if(keep_draws) attr(out,'draws') <- draws
  out
}

# ---- Separation labels ------------------------------------------------------

quantity_order <- function(q) match(q,c(ANATOMY_QUANTITIES,setdiff(unique(q),ANATOMY_QUANTITIES)))
join_labels <- function(q) paste(q[order(quantity_order(q))],collapse=';')
safe_max <- function(x) if(!length(x) || all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)

chain_separation <- function(anatomy) {
  need <- c('key','species','chain','quantity','chain_mean','chain_sd',
    'split_rhat_within_chain','truth','rhat_all_chains','ess_bulk_all_chains')
  absent <- setdiff(need,names(anatomy))
  if(length(absent)) stop('Anatomy lacks columns: ',paste(absent,collapse=', '))
  cells <- split(anatomy,list(anatomy$key,anatomy$species,anatomy$quantity),drop=TRUE,sep='\r')
  separation <- do.call(rbind,lapply(cells,function(g) {
    g <- g[order(g$chain),,drop=FALSE]
    if(anyDuplicated(g$chain) || length(unique(g$rhat_all_chains))!=1L)
      stop('Anatomy rows are not one per chain for ',g$key[1],' species ',g$species[1],' ',g$quantity[1])
    gap <- diff(range(g$chain_mean));pooled <- sqrt(mean(g$chain_sd^2))
    separation <- if(gap==0) 0 else if(pooled==0) Inf else gap/pooled
    within <- g$split_rhat_within_chain
    label <- if(any(within>SPLIT_RHAT_LIMIT,na.rm=TRUE)) 'drifting' else
      if(separation>SEPARATION_LIMIT) 'separated' else 'agrees'
    data.frame(key=g$key[1],species=g$species[1],quantity=g$quantity[1],n_chains=nrow(g),
      pooled_mean=mean(g$chain_mean),truth=g$truth[1],gap=gap,pooled_within_sd=pooled,
      separation=separation,max_split_rhat_within_chain=safe_max(within),
      rhat=g$rhat_all_chains[1],ess_bulk=g$ess_bulk_all_chains[1],
      fixed=all(g$chain_sd==0) && gap==0,label=label,
      lowest_chain=g$chain[which.min(g$chain_mean)],highest_chain=g$chain[which.max(g$chain_mean)],
      stringsAsFactors=FALSE)
  }))
  separation <- separation[order(separation$key,separation$species,quantity_order(separation$quantity)),,drop=FALSE]
  rownames(separation) <- NULL
  species <- do.call(rbind,lapply(split(separation,list(separation$key,separation$species),drop=TRUE,sep='\r'),function(g) {
    core <- g[!g$quantity %in% LOADING_QUANTITIES,,drop=FALSE]
    load <- g[g$quantity %in% LOADING_QUANTITIES,,drop=FALSE]
    label <- if(any(core$label=='separated')) 'separated' else
      if(any(core$label=='drifting')) 'drifting' else 'agrees'
    data.frame(key=g$key[1],species=g$species[1],label=label,
      separated_quantities=join_labels(core$quantity[core$label=='separated']),
      drifting_quantities=join_labels(core$quantity[core$label=='drifting']),
      separated_loadings=join_labels(load$quantity[load$label=='separated']),
      drifting_loadings=join_labels(load$quantity[load$label=='drifting']),
      max_separation=safe_max(core$separation[is.finite(core$separation)]),
      max_rhat=safe_max(core$rhat),stringsAsFactors=FALSE)
  }))
  species <- species[order(species$key,species$species),,drop=FALSE];rownames(species) <- NULL
  chain_groups <- lapply(split(separation,list(separation$key,separation$species),drop=TRUE,sep='\r'),function(g) {
    if(!any(g$label=='separated')) return(NULL)
    k <- g$key[1];s <- g$species[1]
    psi <- anatomy[anatomy$key==k & anatomy$species==s & anatomy$quantity=='mean_psi_original_sites',,drop=FALSE]
    if(nrow(psi)<2L) stop('Chain grouping needs mean_psi_original_sites chain means for ',k,' species ',s)
    psi <- psi[order(psi$chain),,drop=FALSE]
    o <- order(psi$chain_mean);cut <- which.max(diff(psi$chain_mean[o]))
    group <- integer(nrow(psi));group[o[seq_len(cut)]] <- 1L;group[o[-seq_len(cut)]] <- 2L
    data.frame(key=k,species=s,chain=psi$chain,chain_group=group,chain_mean_psi=psi$chain_mean,
      species_label=species$label[species$key==k & species$species==s],
      separated_any=join_labels(g$quantity[g$label=='separated']),stringsAsFactors=FALSE)
  })
  chain_groups <- do.call(rbind,chain_groups)
  if(is.null(chain_groups)) chain_groups <- data.frame(key=character(),species=integer(),
    chain=integer(),chain_group=integer(),chain_mean_psi=numeric(),species_label=character(),
    separated_any=character(),stringsAsFactors=FALSE)
  rownames(chain_groups) <- NULL
  list(separation=separation,species=species,chain_groups=chain_groups)
}

# ---- Trace figures ------------------------------------------------------------

# Categorical slots 1-8 of the dataviz reference palette, in its fixed order
# (validated: slots 1-4 pass the adjacent CVD and normal-vision checks on the
# light surface; the contrast warning is relieved by the legend and the CSVs).
CHAIN_COLOURS <- c('#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#4a3aa7','#e34948')
INK <- c(primary='#0b0b0b',secondary='#52514e',muted='#8a8984',grid='#e4e3df',surface='#fcfcfb')

# One figure for one key: a trace panel per species and quantity labelled
# separated or drifting (all chains overlaid, every kth retained draw), with
# the generating value dashed. With no labelled quantity, the four highest
# all-chain Rhats are shown instead and the title says so.
plot_traces <- function(anatomy,separation,file,max_panels=24L) {
  draws <- attr(anatomy,'draws')
  if(is.null(draws)) stop('plot_traces needs chain_anatomy(..., keep_draws = TRUE)')
  key <- unique(anatomy$key);stopifnot(length(key)==1L)
  sep <- separation$separation[separation$separation$key==key,,drop=FALSE]
  chosen <- sep[sep$label %in% c('separated','drifting'),,drop=FALSE]
  note <- NULL
  if(!nrow(chosen)) {
    chosen <- head(sep[order(-sep$rhat),,drop=FALSE],4L)
    note <- 'No quantity is labelled separated or drifting; the four highest all-chain Rhats are shown.'
  }
  chosen <- chosen[order(chosen$label!='separated',chosen$species,quantity_order(chosen$quantity)),,drop=FALSE]
  if(nrow(chosen)>max_panels) {
    note <- sprintf('%d labelled quantities; the %d with the largest separation are shown.',nrow(chosen),max_panels)
    chosen <- head(chosen[order(-chosen$separation),,drop=FALSE],max_panels)
  }
  K <- dim(draws[[1]])[3];ni <- dim(draws[[1]])[2]
  if(K>length(CHAIN_COLOURS)) stop('The trace figure distinguishes at most ',length(CHAIN_COLOURS),' chains')
  colours <- grDevices::adjustcolor(CHAIN_COLOURS[seq_len(K)],alpha.f=.85)
  keep <- seq(1L,ni,by=max(1L,ni%/%1500L))
  n <- nrow(chosen);columns <- min(3L,n);rows <- ceiling(n/columns)
  grDevices::png(file,width=560*columns,height=300*rows+170,res=110,bg=INK[['surface']])
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::par(mfrow=c(rows,columns),mar=c(3.2,4.2,3.2,1),oma=c(3.2,0,3.6,0),
    col.axis=INK[['secondary']],col.lab=INK[['secondary']],fg=INK[['muted']],las=1,cex.axis=.85)
  for(i in seq_len(n)) {
    r <- chosen[i,];x <- matrix(draws[[r$quantity]][r$species,,],ni,K)[keep,,drop=FALSE]
    ylim <- range(c(x,r$truth),na.rm=TRUE)
    graphics::plot(NA,xlim=range(keep),ylim=ylim,xlab='',ylab='',bty='n')
    graphics::abline(h=pretty(ylim),col=INK[['grid']],lwd=.8)
    graphics::matlines(keep,x,lty=1,lwd=1,col=colours)
    if(is.finite(r$truth)) graphics::abline(h=r$truth,lty=2,lwd=1.6,col=INK[['primary']])
    graphics::mtext(sprintf('Species %d, %s',r$species,r$quantity),side=3,line=1.6,adj=0,cex=.9,
      font=2,col=INK[['primary']])
    graphics::mtext(sprintf('%s; gap %.1f SD; Rhat %.2f; max chain split-Rhat %.3f',r$label,
      r$separation,r$rhat,r$max_split_rhat_within_chain),side=3,line=.5,adj=0,cex=.72,col=INK[['secondary']])
    graphics::mtext('Retained iteration',side=1,line=2.1,cex=.7,col=INK[['secondary']])
  }
  graphics::mtext(sprintf('%s: per-chain traces of retained draws (every %dth draw)',key,max(1L,ni%/%1500L)),
    outer=TRUE,side=3,line=2,adj=0,cex=1,font=2,col=INK[['primary']])
  if(!is.null(note)) graphics::mtext(note,outer=TRUE,side=3,line=.8,adj=0,cex=.8,col=INK[['secondary']])
  graphics::par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  graphics::plot(0,0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
  graphics::legend('bottom',legend=c(paste('Chain',seq_len(K)),'Generating value'),
    col=c(CHAIN_COLOURS[seq_len(K)],INK[['primary']]),lty=c(rep(1,K),2),lwd=c(rep(2,K),1.6),
    horiz=TRUE,bty='n',cex=.85,text.col=INK[['secondary']],inset=c(0,.005),xpd=NA)
  invisible(file)
}

# ---- Runner ---------------------------------------------------------------------

relative_to <- function(path,root) {
  prefix <- paste0(normalizePath(root),'/');path <- normalizePath(path)
  if(startsWith(path,prefix)) substring(path,nchar(prefix)+1L) else path
}

# Numeric columns to seven significant digits, so the CSVs stay compact.
write_compact <- function(x,file) {
  for(nm in names(x)) if(is.double(x[[nm]])) x[[nm]] <- signif(x[[nm]],7L)
  utils::write.csv(x,file,row.names=FALSE)
}

run_anatomy <- function(keys=c(FLAGGED_KEYS,REFERENCE_KEYS),out=OUTPUT_DIR) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(k in intersect(keys,FLAGGED_KEYS)) if(!selection_flagged(k)) stop(k,' is not flagged in the pr11 selection')
  for(k in intersect(keys,REFERENCE_KEYS)) if(selection_flagged(k)) stop(k,' is flagged in the pr11 selection')
  anatomy <- list();separation <- list();species <- list();groups <- list();provenance <- list()
  for(key in keys) {
    started <- Sys.time();flagged <- key %in% FLAGGED_KEYS
    fit <- selected_fit(key);input <- selected_input(key)
    a <- chain_anatomy(fit,input,keep_draws=flagged)
    s <- chain_separation(a)
    figure <- ''
    if(flagged) {figure <- paste0('trace-',key,'.png');plot_traces(a,s,file.path(out,figure))}
    check <- attr(a,'psi_check')
    provenance[[key]] <- data.frame(key=key,role=if(flagged) 'flagged' else 'reference',
      schedule=selection_row(key)$schedule,fit_file=relative_to(fit,ARCHIVES),
      fit_md5=unname(tools::md5sum(fit)),input_file=relative_to(input,ARCHIVES),
      input_md5=unname(tools::md5sum(input)),n_chains=sum(check$n_chains),
      n_iterations=check$n_iterations[1],psi_max_abs_diff=max(check$psi_max_abs_diff),
      loadings_truth_constrained=attr(a,'loadings_constrained'),trace_figure=figure,stringsAsFactors=FALSE)
    attr(a,'draws') <- NULL
    anatomy[[key]] <- a;separation[[key]] <- s$separation;species[[key]] <- s$species
    groups[[key]] <- s$chain_groups
    labelled <- s$species[s$species$label!='agrees',,drop=FALSE]
    cat(format(Sys.time()),key,sprintf('done in %.0f s;',as.numeric(difftime(Sys.time(),started,units='secs'))),
      sprintf('psi check %.2e;',max(check$psi_max_abs_diff)),
      if(nrow(labelled)) paste0('species ',labelled$species,' ',labelled$label,collapse='; ') else 'all species agree',
      '\n')
    flush.console();rm(a,s);invisible(gc(FALSE))
  }
  bind <- function(x) {y <- do.call(rbind,unname(x));rownames(y) <- NULL;y}
  write_compact(bind(anatomy),file.path(out,'chain-summary.csv'))
  write_compact(bind(separation),file.path(out,'separation.csv'))
  write_compact(bind(species),file.path(out,'species-labels.csv'))
  write_compact(bind(groups),file.path(out,'chain-groups.csv'))
  write_compact(bind(provenance),file.path(out,'provenance.csv'))
  sources <- c(file.path(ANATOMY_DIR,c('anatomy.R','test-anatomy.R')),unname(REUSED_SOURCES))
  utils::write.csv(data.frame(file=vapply(sources,relative_to,'',root=REPO),
    md5=unname(tools::md5sum(sources))),file.path(out,'source-hashes.csv'),row.names=FALSE)
  cat('Wrote',length(keys),'fits to',out,'\n')
  invisible(TRUE)
}

if(sys.nframe()==0L) run_anatomy()
