# Community 5 data features and truth: Task 2 of PLAN.md in this directory.
# It reads the saved simulation input and the saved pr11 fit of
# design-qfar_K6-sites300-05 and never fits a model; archives are read, never
# written.
#
# Question (Q3 of PLAN.md): what in the data allows both explanations of
# species 6 that the chains found (Task 1: chains 1 and 3 near the generating
# values, chains 2 and 4 with a high field false-positive probability theta0,
# a low collection intercept and flipped environmental slopes)?
#
# Interfaces
#
#   site_features(input, species, threshold = READ_THRESHOLD)
#     One row per site for one species (index into the OTU columns): the true
#     occupancy state (true_occupied, from z_true), the true occupancy
#     probability (true_psi, plogis of the generating linear predictor), the
#     true collection state of each field sample (collected_sample<m>, from
#     w_true: whether the sample holds the species' DNA, by true occupancy or by
#     field contamination), the number of PCRs that are positive in total, per
#     primer and per sample and primer, and positive_samples, the number of
#     field samples with any positive PCR. A PCR is positive when it has at
#     least `threshold` reads, which is how runOccJSDM turns reads into
#     detections. original_site marks the first input$design$original_sites
#     sites (the 100 that the pr11 scorer evaluates).
#
#   positive_routes(features, scope = c('all','original'))
#     One row summarising one species' site features: positive PCRs at truly
#     unoccupied sites (explainable only by field or laboratory false
#     positives) and at truly occupied sites, each split by whether the field
#     sample truly holds DNA (collected) or not (uncollected); occupied sites
#     with no positive PCR at all; and the naive rates that ignore all false
#     positives. Also the collected samples by the true state of their site,
#     the positive PCR rate in collected and uncollected samples, and
#     sample_positive_count_auc: how well the count of positive PCRs in a
#     sample separates collected from uncollected samples (0.5 no better than
#     chance, 1 perfectly).
#
#   choose_comparison_species(prevalence, species_labels, target, n)
#     The n species labelled 'agrees' by Task 1 (species_labels as in
#     results/anatomy/species-labels.csv) whose prevalence is closest to the
#     target's, closest first, ties to the lower index. Prevalence is the mean
#     true occupancy probability over the fitted sites, as the pr11 scorer
#     defines it (colMeans(true_psi)).
#
#   subset_fit_chains(fit, chains), chain_reconstruction(fit, original_sites)
#     A fit reduced to some of its chains, and, per chain, anatomy.R's
#     fit_draws() applied to that one chain: the draws of every anatomy
#     quantity and the per-site sum of the reconstructed occupancy
#     probabilities over the chain's draws. fit_draws() is reused, not copied,
#     so the occupancy reconstruction is exactly Task 1's; the test sums the
#     chains and checks them against the fit's stored psi_output.
#
#   group_parameter_rows(recon, chains, species, truth, key, chain_group)
#     For one chain group: the group's pooled posterior mean of each quantity
#     (draws of the chains in the group pooled), its 95 percent central
#     interval, the generating value, the signed and absolute error, and
#     whether the interval contains the generating value.
#
#   site_occupancy_errors(psi, truth, sites)
#     MAE, signed error (estimate minus truth), RMSE and correlation of a
#     posterior-mean occupancy probability vector against the true one, over
#     the given sites.
#
#   psi_by_true_state(sites)
#     The mean posterior occupancy probability of each chain group among the
#     truly occupied and the truly unoccupied sites, from community5_modes()$sites.
#
#   community5_modes(fit, input, groups, key, species)
#     modes-vs-truth.csv rows and per-site occupancy rows for each Task 1
#     chain group of a species, plus the pooled chains for reference.
#
# Truths and the quantity list are anatomy.R's (anatomy_truth(),
# ANATOMY_QUANTITIES): they follow the pr11 scorer.
#
# Read threshold. The fit was made by current-main-recheck/run.R, which calls
# runOccJSDM(..., threshold = 1) for design-family jobs, so a PCR is positive
# with at least one read. The test checks that the input's reads are exactly
# the reads stored in the fit.
#
# Run `Rscript data-features.R` to write results/community5/; tests:
# `Rscript test-data-features.R`.

FEATURES_DIR <- local({
  here <- NULL
  for(i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if(!is.null(f)) {here <- dirname(normalizePath(f));break}
  }
  if(is.null(here)) {
    a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
    if(length(a)) here <- dirname(normalizePath(sub('^--file=','',a[1])))
  }
  if(is.null(here)) getwd() else here
})
source(file.path(FEATURES_DIR,'anatomy.R'))

COMMUNITY5_KEY <- 'design-qfar_K6-sites300-05'
TARGET_SPECIES <- 6L
N_COMPARISON <- 2L
READ_THRESHOLD <- 1
ANATOMY_RESULTS <- file.path(FEATURES_DIR,'results/anatomy')
COMMUNITY5_DIR <- file.path(FEATURES_DIR,'results/community5')

# ---- Site features ------------------------------------------------------------

site_features <- function(input,species,threshold=READ_THRESHOLD) {
  dl <- input$sim$data_list;tp <- input$sim$true_params
  info <- dl$info;Y <- dl$OTU;species <- as.integer(species)
  stopifnot(length(species)==1L,!is.na(species),species>=1L,species<=ncol(Y),nrow(info)==nrow(Y),
    is.numeric(threshold),length(threshold)==1L,threshold>0)
  if(anyNA(Y[,species])) stop('Reads of species ',species,' contain NA; the counts assume none')
  z <- tp$z_true;w <- tp$w_true;eta <- tp$jsdmParams_true$eta
  sites <- sort(unique(info$Site));n <- length(sites)
  stopifnot(identical(as.integer(sites),seq_len(n)),nrow(z)==n,nrow(eta)==n,ncol(z)==ncol(Y),ncol(eta)==ncol(Y))
  # Field samples: one Sample id belongs to one site, ids run 1..n*M in site
  # order, and w_true has one row per sample id (as the pr11 scorer indexes it).
  smp <- unique(info[c('Site','Sample')]);smp <- smp[order(smp$Site,smp$Sample),,drop=FALSE]
  stopifnot(!anyDuplicated(smp$Sample),identical(as.integer(smp$Sample),seq_len(nrow(smp))),nrow(w)==nrow(smp),
    ncol(w)==ncol(Y))
  smp$slot <- as.integer(ave(smp$Sample,smp$Site,FUN=seq_along))
  M <- max(smp$slot);stopifnot(all(table(smp$Site)==M))
  primers <- sort(unique(info$Primer));P <- length(primers)
  positive <- Y[,species]>=threshold
  cells <- list(factor(info$Sample,smp$Sample),factor(info$Primer,primers))
  counts <- tapply(positive,cells,sum);pcrs <- tapply(positive,cells,length)
  stopifnot(!anyNA(counts),length(unique(as.vector(pcrs)))==1L)
  K <- as.integer(pcrs[1])
  sample_of <- function(m) smp$Sample[smp$slot==m]
  out <- data.frame(site=as.integer(sites),
    original_site=sites<=(input$design$original_sites %||% n),species=species,
    species_name=(colnames(Y)[species] %||% as.character(species)),
    true_occupied=as.integer(z[,species]),true_psi=unname(plogis(eta[,species])),stringsAsFactors=FALSE)
  for(m in seq_len(M)) out[[paste0('collected_sample',m)]] <- as.integer(w[sample_of(m),species])
  out$n_collected_samples_true <- as.integer(rowSums(out[paste0('collected_sample',seq_len(M))]))
  out$field_samples <- M;out$pcrs_per_site <- as.integer(M*P*K)
  for(p in seq_len(P)) out[[paste0('positive_pcr_primer',p)]] <-
    as.integer(Reduce(`+`,lapply(seq_len(M),function(m) counts[sample_of(m),p])))
  out$positive_pcr_total <- as.integer(rowSums(out[paste0('positive_pcr_primer',seq_len(P))]))
  out$positive_samples <- as.integer(Reduce(`+`,lapply(seq_len(M),function(m)
    rowSums(counts[sample_of(m),,drop=FALSE])>0)))
  out$naive_detected <- out$positive_pcr_total>0L
  for(m in seq_len(M)) for(p in seq_len(P))
    out[[paste0('positives_sample',m,'_primer',p)]] <- as.integer(counts[sample_of(m),p])
  rownames(out) <- NULL
  out
}

# ---- Positive routes ------------------------------------------------------------

rate <- function(x,n) if(n>0) x/n else NaN

# Probability that a random collected sample has more positive PCRs than a
# random uncollected one, ties counted half (the Wilcoxon statistic scaled).
count_auc <- function(collected,uncollected) {
  if(!length(collected) || !length(uncollected)) return(NA_real_)
  r <- rank(c(collected,uncollected))
  (sum(r[seq_along(collected)])-length(collected)*(length(collected)+1)/2)/
    (length(collected)*length(uncollected))
}

positive_routes <- function(features,scope=c('all','original')) {
  scope <- match.arg(scope)
  f <- if(scope=='original') features[features$original_site,,drop=FALSE] else features
  stopifnot(nrow(f)>0L,length(unique(f$species))==1L,length(unique(f$field_samples))==1L)
  M <- f$field_samples[1];per_sample <- f$pcrs_per_site[1]/M
  sample_pos <- sapply(seq_len(M),function(m)
    rowSums(f[grep(paste0('^positives_sample',m,'_primer'),names(f))]))
  sample_col <- sapply(seq_len(M),function(m) f[[paste0('collected_sample',m)]])
  sample_pos <- matrix(sample_pos,nrow(f),M);sample_col <- matrix(sample_col,nrow(f),M)
  occ <- f$true_occupied==1L;site_positive <- f$positive_pcr_total>0L
  route <- function(site_mask,collected) as.integer(sum(sample_pos[site_mask,,drop=FALSE][
    sample_col[site_mask,,drop=FALSE]==collected]))
  n_pcr <- as.integer(sum(f$pcrs_per_site));n_pos <- as.integer(sum(f$positive_pcr_total))
  n_collected <- sum(sample_col==1L);n_uncollected <- sum(sample_col==0L)
  occ_matrix <- matrix(occ,nrow(f),M)
  pos_collected <- sum(sample_pos[sample_col==1L]);pos_uncollected <- sum(sample_pos[sample_col==0L])
  data.frame(species=f$species[1],species_name=f$species_name[1],
    scope=if(scope=='original') 'original_sites' else 'all_sites',
    n_sites=nrow(f),n_occupied_sites=sum(occ),n_unoccupied_sites=sum(!occ),
    true_prevalence=mean(occ),n_pcr=n_pcr,n_positive_pcr=n_pos,
    positive_pcr_at_occupied_sites=as.integer(sum(f$positive_pcr_total[occ])),
    positive_pcr_at_unoccupied_sites=as.integer(sum(f$positive_pcr_total[!occ])),
    share_positive_pcr_at_unoccupied_sites=rate(sum(f$positive_pcr_total[!occ]),n_pos),
    positive_pcr_occupied_collected=route(occ,1L),positive_pcr_occupied_uncollected=route(occ,0L),
    positive_pcr_unoccupied_collected=route(!occ,1L),positive_pcr_unoccupied_uncollected=route(!occ,0L),
    occupied_sites_without_positive=sum(occ & !site_positive),
    unoccupied_sites_with_positive=sum(!occ & site_positive),
    unoccupied_sites_with_collected_sample=sum(!occ & rowSums(sample_col)>0L),
    occupied_sites_without_collected_sample=sum(occ & rowSums(sample_col)==0L),
    naive_detection_rate_sites=rate(sum(occ & site_positive),sum(occ)),
    naive_false_positive_rate_sites=rate(sum(!occ & site_positive),sum(!occ)),
    naive_positive_rate_sites=mean(site_positive),
    positive_pcr_rate=rate(n_pos,n_pcr),
    n_collected_samples=as.integer(n_collected),n_uncollected_samples=as.integer(n_uncollected),
    collected_samples_at_occupied_sites=as.integer(sum(sample_col==1L & occ_matrix)),
    collected_samples_at_unoccupied_sites=as.integer(sum(sample_col==1L & !occ_matrix)),
    sample_positive_count_auc=count_auc(sample_pos[sample_col==1L],sample_pos[sample_col==0L]),
    positive_pcr_rate_collected=rate(pos_collected,n_collected*per_sample),
    positive_pcr_rate_uncollected=rate(pos_uncollected,n_uncollected*per_sample),
    stringsAsFactors=FALSE,row.names=NULL)
}

# ---- Comparison species -----------------------------------------------------------

choose_comparison_species <- function(prevalence,species_labels,target,n=N_COMPARISON) {
  if(!target %in% species_labels$species) stop('The target species is not in the species labels')
  stopifnot(length(prevalence)>=max(species_labels$species))
  candidates <- species_labels$species[species_labels$label=='agrees' & species_labels$species!=target]
  if(length(candidates)<n) stop('Only ',length(candidates),' candidate species are labelled agrees; ',n,' wanted')
  distance <- abs(prevalence[candidates]-prevalence[target])
  candidates[order(distance,candidates)][seq_len(n)]
}

# ---- Occupancy errors -------------------------------------------------------------

site_occupancy_errors <- function(psi,truth,sites=seq_along(psi)) {
  stopifnot(length(psi)==length(truth))
  p <- psi[sites];t <- truth[sites]
  data.frame(n_sites=length(sites),mean_estimate=mean(p),mean_truth=mean(t),
    mae=mean(abs(p-t)),signed_error=mean(p-t),rmse=sqrt(mean((p-t)^2)),
    correlation=if(length(sites)>2L && sd(p)>0 && sd(t)>0) cor(p,t) else NA_real_)
}

# ---- Chain groups against the generating values -------------------------------------

# A fit reduced to some of its chains. Only the arrays fit_draws() reads are
# reduced; the stored posterior means (psi_output and the others) still cover
# all chains and must not be used on a subset.
subset_fit_chains <- function(fit,chains) {
  ro <- fit$results_output;nc <- dim(ro$jsdm_output$B0_output)[3]
  if(!length(chains) || anyNA(chains) || any(chains<1L | chains>nc) || anyDuplicated(chains))
    stop('Invalid chain selection ',paste(chains,collapse=','),': the fit has ',nc,' chains')
  keep <- function(x) {
    d <- dim(x);if(d[length(d)]!=nc) stop('Array does not end in the chain dimension')
    idx <- rep(list(TRUE),length(d));idx[[length(d)]] <- as.integer(chains)
    do.call(`[`,c(list(x),idx,list(drop=FALSE)))
  }
  for(nm in c('B0_output','B_output','U_output','L_output')) ro$jsdm_output[[nm]] <- keep(ro$jsdm_output[[nm]])
  for(nm in c('theta0_output','beta_theta_output','p_output','q_output')) ro[[nm]] <- keep(ro[[nm]])
  fit$results_output <- ro
  fit
}

chain_reconstruction <- function(fit,original_sites) {
  nc <- dim(fit$results_output$jsdm_output$B0_output)[3]
  lapply(seq_len(nc),function(ch) {
    part <- fit_draws(subset_fit_chains(fit,ch),original_sites)
    list(chain=ch,draws=part$draws,psi_sum=part$psi_sum,ni=part$ni)
  })
}

group_psi <- function(recon,chains) {
  ni <- unique(vapply(recon[chains],`[[`,integer(1),'ni'));stopifnot(length(ni)==1L)
  Reduce(`+`,lapply(recon[chains],`[[`,'psi_sum'))/(ni*length(chains))
}

group_parameter_rows <- function(recon,chains,species,truth,key,chain_group) {
  rows <- lapply(names(truth),function(q) {
    x <- unlist(lapply(chains,function(ch) recon[[ch]]$draws[[q]][species,,]),use.names=FALSE)
    m <- mean(x);ci <- unname(quantile(x,c(.025,.975)));t <- unname(truth[[q]])
    data.frame(key=key,species=as.integer(species),chain_group=as.character(chain_group),
      chains=paste(chains,collapse=';'),n_chains=length(chains),kind='parameter',scope=NA_character_,
      quantity=q,n_draws=length(x),n_sites=NA_integer_,group_mean=m,group_sd=sd(x),q025=ci[1],q975=ci[2],
      truth=t,error=m-t,abs_error=abs(m-t),truth_in_interval=t>=ci[1] & t<=ci[2],stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

# Occupancy summary rows in the same layout: the value is the group's summary
# and the truth is what perfect recovery would give (0 for the errors, 1 for
# the correlation), so error = value - truth.
group_occupancy_rows <- function(psi,true_psi,original_sites,key,species,chains,chain_group) {
  scopes <- list(original_sites=original_sites,all_sites=seq_along(psi))
  rows <- list()
  for(sc in names(scopes)) {
    e <- site_occupancy_errors(psi,true_psi,scopes[[sc]])
    value <- c(occupancy_mae=e$mae,occupancy_signed_error=e$signed_error,occupancy_rmse=e$rmse,
      occupancy_correlation=e$correlation)
    target <- c(occupancy_mae=0,occupancy_signed_error=0,occupancy_rmse=0,occupancy_correlation=1)
    rows[[sc]] <- data.frame(key=key,species=as.integer(species),chain_group=as.character(chain_group),
      chains=paste(chains,collapse=';'),n_chains=length(chains),kind='site_occupancy',scope=sc,
      quantity=names(value),n_draws=NA_integer_,n_sites=e$n_sites,group_mean=unname(value),group_sd=NA_real_,
      q025=NA_real_,q975=NA_real_,truth=unname(target),error=unname(value-target),
      abs_error=unname(abs(value-target)),truth_in_interval=NA,stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}

# Mean posterior occupancy probability of each chain group among the truly
# occupied and the truly unoccupied sites (of the original sites and of all).
psi_by_true_state <- function(sites) {
  scopes <- list(original_sites=sites[sites$original_site,,drop=FALSE],all_sites=sites)
  rows <- list()
  for(sc in names(scopes)) {
    d <- scopes[[sc]]
    for(g in unique(d$chain_group)) for(z in c(1L,0L)) {
      x <- d[d$chain_group==g & d$true_occupied==z,,drop=FALSE]
      if(!nrow(x)) next
      rows[[length(rows)+1L]] <- data.frame(key=x$key[1],species=x$species[1],chain_group=g,
        chains=x$chains[1],scope=sc,true_occupied=z,n_sites=nrow(x),
        mean_posterior_psi=mean(x$posterior_mean_psi),mean_true_psi=mean(x$true_psi),stringsAsFactors=FALSE)
    }
  }
  out <- do.call(rbind,rows);rownames(out) <- NULL;out
}

# Task 1's chain groups for one species of a key, from its chain_separation()
# output (the rows of results/anatomy/chain-groups.csv), keeping only a species
# labelled separated on quantities that are identified (loadings do not count).
separated_groups <- function(groups,key,species) {
  g <- groups[groups$key==key & groups$species==species,,drop=FALSE]
  if(!nrow(g)) stop('No chain groups for ',key,' species ',species)
  if(!all(g$species_label=='separated')) stop('Species ',species,' of ',key,' is not labelled separated')
  g[order(g$chain),,drop=FALSE]
}

community5_modes <- function(fit,input,groups,key=COMMUNITY5_KEY,species=TARGET_SPECIES) {
  truth <- anatomy_truth(fit,input)
  nc <- dim(fit$results_output$jsdm_output$B0_output)[3]
  g <- separated_groups(groups,key,species)
  if(!identical(g$chain,seq_len(nc))) stop('Chain groups do not cover the fit\'s ',nc,' chains')
  recon <- chain_reconstruction(fit,truth$original_sites)
  psi_all <- group_psi(recon,seq_len(nc))
  psi_difference <- max(abs(psi_all-fit$results_output$psi_output))
  if(!is.finite(psi_difference) || psi_difference>=PSI_TOLERANCE)
    stop('Chain reconstructions disagree with the stored psi_output: ',format(psi_difference))
  true_psi <- plogis(input$sim$true_params$jsdmParams_true$eta[,species])
  z <- input$sim$true_params$z_true[,species]
  truth_row <- truth$truth[species,]
  modes <- list();sites <- list()
  add <- function(chains,label) {
    psi <- group_psi(recon,chains)[,species]
    modes[[length(modes)+1L]] <<- group_parameter_rows(recon,chains,species,truth_row,key,label)
    modes[[length(modes)+1L]] <<- group_occupancy_rows(psi,true_psi,truth$original_sites,key,species,chains,label)
    sites[[length(sites)+1L]] <<- data.frame(key=key,species=as.integer(species),chain_group=as.character(label),
      chains=paste(chains,collapse=';'),site=seq_along(psi),original_site=seq_along(psi) %in% truth$original_sites,
      true_occupied=as.integer(z),true_psi=true_psi,posterior_mean_psi=psi,error=psi-true_psi,
      stringsAsFactors=FALSE)
  }
  for(grp in sort(unique(g$chain_group))) add(g$chain[g$chain_group==grp],grp)
  add(seq_len(nc),'all')
  bind <- function(x) {y <- do.call(rbind,x);rownames(y) <- NULL;y}
  list(modes=bind(modes),sites=bind(sites),psi_max_abs_diff=psi_difference,truth=truth)
}

# ---- Runner --------------------------------------------------------------------------

run_data_features <- function(key=COMMUNITY5_KEY,species=TARGET_SPECIES,out=COMMUNITY5_DIR) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  input_path <- selected_input(key);fit_path <- selected_fit(key)
  input <- readRDS(input_path)
  saved <- readRDS(fit_path);fit <- saved$fit;job <- saved$job;rm(saved)
  stopifnot(identical(job$key,key),identical(job$input_md5,unname(tools::md5sum(input_path))),
    identical(unname(fit$infos$OTU),unname(input$sim$data_list$OTU)))
  tp <- input$sim$true_params
  labels <- read.csv(file.path(ANATOMY_RESULTS,'species-labels.csv'),stringsAsFactors=FALSE)
  labels <- labels[labels$key==key,,drop=FALSE]
  stopifnot(nrow(labels)==ncol(input$sim$data_list$OTU),labels$label[labels$species==species]=='separated')
  groups <- read.csv(file.path(ANATOMY_RESULTS,'chain-groups.csv'),stringsAsFactors=FALSE)
  prevalence <- colMeans(plogis(tp$jsdmParams_true$eta))
  comparison <- choose_comparison_species(prevalence,labels,species,N_COMPARISON)
  S <- ncol(input$sim$data_list$OTU)
  role <- rep('other',S);role[comparison] <- 'comparison';role[species] <- 'target'
  # Per-site data of the target and the comparison species.
  shown <- c(species,comparison)
  site_rows <- do.call(rbind,lapply(shown,function(s) {
    f <- site_features(input,s);cbind(key=key,role=role[s],f,stringsAsFactors=FALSE)
  }))
  write_compact(site_rows,file.path(out,'species6-data.csv'))
  # Positive routes of every species, so the target can be read against all.
  truth <- anatomy_truth(fit,input)$truth
  routes <- do.call(rbind,lapply(seq_len(S),function(s) {
    f <- site_features(input,s)
    do.call(rbind,lapply(c('all','original'),function(sc) {
      r <- positive_routes(f,sc)
      cbind(data.frame(key=key,role=role[s],species=s,species_name=f$species_name[1],
        chain_label=labels$label[labels$species==s],prevalence_psi=unname(prevalence[s]),
        theta0_true=truth[s,'theta0'],
        p_effective_primer1=truth[s,'p_primer1'],p_effective_primer2=truth[s,'p_primer2'],
        q_effective_primer1=truth[s,'q_primer1'],q_effective_primer2=truth[s,'q_primer2'],
        stringsAsFactors=FALSE),r[-(1:2)])
    }))
  }))
  write_compact(routes,file.path(out,'positive-routes.csv'))
  # Chain groups against the generating values.
  m <- community5_modes(fit,input,groups,key,species)
  write_compact(m$modes,file.path(out,'modes-vs-truth.csv'))
  write_compact(m$sites,file.path(out,'modes-site-occupancy.csv'))
  write_compact(psi_by_true_state(m$sites),file.path(out,'modes-psi-by-true-state.csv'))
  provenance <- data.frame(key=key,input_file=relative_to(input_path,ARCHIVES),
    input_md5=unname(tools::md5sum(input_path)),fit_file=relative_to(fit_path,ARCHIVES),
    fit_md5=unname(tools::md5sum(fit_path)),read_threshold=READ_THRESHOLD,target_species=species,
    comparison_species=paste(comparison,collapse=';'),
    comparison_rule=paste0('species labelled agrees in results/anatomy/species-labels.csv with the ',
      'closest mean true occupancy probability over the 300 fitted sites to the target, closest first'),
    n_chains=dim(fit$results_output$jsdm_output$B0_output)[3],
    n_iterations=dim(fit$results_output$jsdm_output$B0_output)[2],
    psi_max_abs_diff=m$psi_max_abs_diff,stringsAsFactors=FALSE)
  write_compact(provenance,file.path(out,'provenance.csv'))
  sources <- c(file.path(FEATURES_DIR,c('data-features.R','test-data-features.R','anatomy.R')),unname(REUSED_SOURCES),
    file.path(ANATOMY_RESULTS,c('species-labels.csv','chain-groups.csv','chain-summary.csv')))
  utils::write.csv(data.frame(file=vapply(sources,relative_to,'',root=REPO),
    md5=unname(tools::md5sum(sources))),file.path(out,'source-hashes.csv'),row.names=FALSE)
  cat('Wrote',out,'; comparison species',paste(comparison,collapse=', '),
    sprintf('; psi check %.2e\n',m$psi_max_abs_diff))
  invisible(list(site_rows=site_rows,routes=routes,modes=m$modes,sites=m$sites))
}

if(sys.nframe()==0L) run_data_features()
