#!/usr/bin/env Rscript
# Analysis of one stage of the mirror-frequency check (README.md in this
# directory): which anchored region each chain of each repeat settled in.
#
#   Rscript analyse.R --repo=REPO --study=ARCHIVE --inputs-root=DIR
#     --stage=a|b|pilot [--archives=DIR] [--repeats=K1-K2] [--out=DIR]
#
# It reads the stage's fit files in repeat order through the convergence-flag
# diagnosis's chain_anatomy() (chains concatenated, four per file), takes
# species 6's draws of the anchor quantities, labels every draw with the
# diagnosis's frozen anchored classifier (AMENDMENT-1, rulings R17 and R18:
# near-truth or mirror by Gaussian density, unknown region when more than 5%
# of a chain's draws are far from both anchors) and gives each chain its
# region. modes.R, anatomy.R and the anchor file are refused if their md5s
# differ from the ones the amendment froze. Writes, under
# results/stage-<stage>/ (or --out): chain-regions.csv (one row per chain),
# repeat-summary.csv (one row per repeat, with the rank-normalised Rhat of
# species 6's quantities over that repeat's four chains), summary.csv (one
# row: counts, the mirror share with its Clopper-Pearson 95% interval, mirror
# counts by chain position) and decision.txt ("run_stage_b: yes" when any
# chain is outside the near-truth region, else "no").

`%||%` <- function(x,y) if(is.null(x)) y else x
here <- local({
  a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  f <- NULL
  for(i in rev(seq_len(sys.nframe()))) {x <- sys.frame(i)$ofile;if(!is.null(x)) {f <- x;break}}
  dirname(normalizePath(f %||% sub('^--file=','',a[1])))
})
source(file.path(here,'run.R'))

DIAG_REL <- 'dev/simstudy/convergence-flag-diagnosis'
# Frozen by AMENDMENT-1 of the diagnosis ("Files Task 3c must check").
DIAG_FROZEN_MD5 <- c('modes.R'='bab7c1c3934ba2e0c83b8f26232a0aaf','anatomy.R'='bd98ff3f40cad972f96bdeca5a07e206')
RHAT_QUANTITIES <- c('theta0','B0','B_slope1','B_slope2','beta_theta_intercept','mean_psi_original_sites')
REGIONS <- c('near-truth','mirror','both','unknown')

load_diagnosis_code <- function(repo) {
  dir <- file.path(repo,DIAG_REL)
  for(f in names(DIAG_FROZEN_MD5)) {
    found <- unname(tools::md5sum(file.path(dir,f)))
    if(!identical(found,DIAG_FROZEN_MD5[[f]])) stop('Diagnosis file changed since its amendment froze it: ',f,' (md5 ',found,')')
  }
  modes <- new.env(parent=globalenv());sys.source(file.path(dir,'modes.R'),envir=modes)
  anatomy <- new.env(parent=globalenv());sys.source(file.path(dir,'anatomy.R'),envir=anatomy)
  list(modes=modes,anatomy=anatomy,dir=dir)
}

stage_fit_paths <- function(study,stage,repeats) {
  paths <- vapply(repeats,function(k) fit_file(study,list(stage=stage,repeat_index=k)),character(1))
  absent <- paths[!file.exists(paths)]
  if(length(absent)) stop('Missing fit file(s): ',paste(basename(absent),collapse=', '))
  paths
}

# The small part of each fit file the tables need; the fit itself is dropped.
fit_records <- function(paths) {
  rows <- lapply(paths,function(p) {
    s <- readRDS(p)
    m <- data.frame(repeat_index=as.integer(s$repeat_index),seed=as.integer(s$seed),source_revision=s$source_revision,
      nchain=as.integer(s$mcmc$nchain),nburn=as.integer(s$mcmc$nburn),niter=as.integer(s$mcmc$niter),nthin=as.integer(s$mcmc$nthin),
      elapsed_seconds=as.numeric(s$elapsed_seconds),n_warnings=length(s$warnings),stringsAsFactors=FALSE)
    rm(s);invisible(gc(FALSE));m
  })
  do.call(rbind,rows)
}

clopper_pearson <- function(x,n) if(n>0) stats::binom.test(x,n)$conf.int[1:2] else c(NA_real_,NA_real_)

analyse_stage <- function(repo,study,inputs_root,archives,stage,repeats=NULL,out=NULL) {
  if(!stage %in% STAGES) stop('stage must be one of ',paste(STAGES,collapse=', '))
  code <- load_diagnosis_code(repo);modes <- code$modes
  repeats <- as.integer(repeats %||% REPEATS[[stage]])
  paths <- stage_fit_paths(study,stage,repeats)
  source <- pr11_job(archives,inputs_root)
  input_path <- checked_input(source$input_file,source$input_md5)
  meta <- fit_records(paths)
  if(!identical(meta$repeat_index,repeats)) stop('Fit files record repeat indices other than the requested ones')
  revision <- unique(meta$source_revision)
  if(length(revision)!=1L) stop('Fits come from different library revisions: ',paste(revision,collapse=', '))
  nchain <- unique(meta$nchain)
  if(length(nchain)!=1L) stop('Fits hold different numbers of chains')
  a <- code$anatomy$chain_anatomy(paths,input_path,keep_draws=TRUE)
  checks <- attr(a,'psi_check')
  if(!identical(checks$fit_file,basename(paths)) || !all(checks$n_chains==nchain)) stop('Fit files were read out of order or with unexpected chain counts')
  quantities <- unique(c(modes$ANCHOR_QUANTITIES,RHAT_QUANTITIES))
  d <- modes$species_mode_draws(a,TARGET_SPECIES,quantities)
  rm(a);invisible(gc(FALSE))
  nc <- ncol(d$theta0)
  if(nc!=length(repeats)*nchain) stop('Chain count ',nc,' is not ',length(repeats),' repeats of ',nchain,' chains')
  anchor <- modes$read_anchor(file.path(code$dir,modes$ANCHOR_FILE))
  r <- modes$assign_anchored(d,anchor)
  reg <- modes$anchored_regions(r)
  if(!identical(reg$chain,seq_len(nc)) || !all(reg$region %in% REGIONS)) stop('Unexpected region table')
  reg$repeat_index <- repeats[ceiling(reg$chain/nchain)]
  reg$position <- (reg$chain-1L)%%nchain+1L
  reg$seed <- meta$seed[match(reg$repeat_index,meta$repeat_index)]
  reg$share_theta0_above_cut <- colMeans(d$theta0>modes$THETA0_CUT)
  for(q in RHAT_QUANTITIES) reg[[paste0('mean_',q)]] <- colMeans(d[[q]])
  reg <- reg[c('repeat_index','position','seed','chain','region','share_near_truth','share_mirror','share_far','visits_both',
    'share_theta0_above_cut',paste0('mean_',RHAT_QUANTITIES))]
  per_repeat <- do.call(rbind,lapply(seq_along(repeats),function(i) {
    ch <- which(ceiling(seq_len(nc)/nchain)==i)
    rh <- vapply(RHAT_QUANTITIES,function(q) posterior::rhat(d[[q]][,ch,drop=FALSE]),numeric(1))
    regs <- reg$region[reg$chain %in% ch]
    data.frame(repeat_index=repeats[i],seed=meta$seed[i],n_near_truth=sum(regs=='near-truth'),n_mirror=sum(regs=='mirror'),
      n_both=sum(regs=='both'),n_unknown=sum(regs=='unknown'),regions=paste(regs,collapse=';'),
      as.list(stats::setNames(rh,paste0('rhat_',RHAT_QUANTITIES))),elapsed_seconds=meta$elapsed_seconds[i],
      n_warnings=meta$n_warnings[i],stringsAsFactors=FALSE)
  }))
  n <- nrow(reg);counts <- vapply(REGIONS,function(x) sum(reg$region==x),integer(1))
  ci <- clopper_pearson(counts[['mirror']],n)
  by_position <- vapply(seq_len(nchain),function(p) sum(reg$region=='mirror' & reg$position==p),integer(1))
  outside <- reg$region!='near-truth'
  summary <- data.frame(stage=stage,library_revision=revision,n_repeats=length(repeats),chains_per_repeat=nchain,n_chains=n,
    n_near_truth=counts[['near-truth']],n_mirror=counts[['mirror']],n_both=counts[['both']],n_unknown=counts[['unknown']],
    n_outside_near_truth=sum(outside),mirror_share=counts[['mirror']]/n,mirror_ci_low=ci[1],mirror_ci_high=ci[2],
    repeats_with_mirror=length(unique(reg$repeat_index[reg$region=='mirror'])),
    repeats_with_any_outside_near_truth=length(unique(reg$repeat_index[outside])),
    mirror_by_position=paste(by_position,collapse=';'),nburn=paste(unique(meta$nburn),collapse=';'),
    niter=paste(unique(meta$niter),collapse=';'),nthin=paste(unique(meta$nthin),collapse=';'),
    mean_elapsed_seconds=mean(meta$elapsed_seconds),anchor_md5=unname(tools::md5sum(file.path(code$dir,modes$ANCHOR_FILE))),
    analysed=format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),stringsAsFactors=FALSE)
  out <- out %||% file.path(repo,STUDY_REL,'results',paste0('stage-',stage))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(reg,file.path(out,'chain-regions.csv'),row.names=FALSE)
  utils::write.csv(per_repeat,file.path(out,'repeat-summary.csv'),row.names=FALSE)
  utils::write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
  writeLines(paste0('run_stage_b: ',if(sum(outside)>0) 'yes' else 'no'),file.path(out,'decision.txt'))
  cat(sprintf('Stage %s: %d repeats, %d chains; near-truth %d, mirror %d, both %d, unknown %d; mirror share %.3f (95%% CI %.3f to %.3f); mirror by position %s; results in %s\n',
    stage,length(repeats),n,counts[['near-truth']],counts[['mirror']],counts[['both']],counts[['unknown']],
    counts[['mirror']]/n,ci[1],ci[2],summary$mirror_by_position,out))
  invisible(list(regions=reg,per_repeat=per_repeat,summary=summary))
}

analyse_main <- function(args) {
  o <- parse_options(args,known=c('repo','study','inputs-root','archives','stage','repeats','out'),
    required=c('repo','study','inputs-root','stage'))
  repo <- normalizePath(o$repo,mustWork=TRUE)
  if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this analyse.R')
  study <- normalizePath(o$study,mustWork=TRUE);inputs_root <- normalizePath(o$`inputs-root`,mustWork=TRUE)
  archives <- normalizePath(o$archives %||% dirname(study),mustWork=TRUE)
  repeats <- NULL
  if(!is.null(o$repeats)) {
    if(!grepl('^[0-9]+(-[0-9]+)?$',o$repeats)) stop('--repeats must be K or K1-K2')
    b <- as.integer(strsplit(o$repeats,'-',fixed=TRUE)[[1]]);repeats <- if(length(b)==2L) seq(b[1],b[2]) else b
  }
  analyse_stage(repo,study,inputs_root,archives,o$stage,repeats,o$out)
  0L
}

if(sys.nframe()==0L) {
  status <- tryCatch(analyse_main(commandArgs(trailingOnly=TRUE)),error=function(e) {cat('ERROR:',conditionMessage(e),'\n');1L})
  quit(save='no',status=status)
}
