#!/usr/bin/env Rscript
# Post-fit sensitivity only. No additional fitting and no excluded communities.
args <- commandArgs(trailingOnly=TRUE)
option <- function(name,default=NULL) {
  hit <- grep(paste0('^--',name,'='),args,value=TRUE)
  if(length(hit)>1L)stop('Repeated option')
  if(!length(hit)){if(is.null(default))stop('Missing --',name);return(default)}
  substring(hit,nchar(name)+4L)
}
repo <- normalizePath(option('repo','.'))
summary_dir <- normalizePath(option('summary'))
out <- option('out',file.path(summary_dir,'chain-review'))
if(dir.exists(out) && length(list.files(out)))stop('Refusing to replace chain review: ',out)
source(file.path(repo,'dev/simstudy/spatial-targeted-recheck/chain-sensitivity.R'))
selection <- read.csv(file.path(summary_dir,'selected-fits.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(selection)==81L,!anyDuplicated(selection$key),
  all(selection$phase[selection$needs_long]=='long'),
  identical(unname(tools::md5sum(selection$selected_file)),selection$selected_md5))
chains <- pooled <- ranges <- list()
for(i in seq_len(nrow(selection))) {
  a <- readRDS(selection$selected_file[i])
  stopifnot(identical(a$job$key,selection$key[i]),ncol(a$chain_probability)==a$mcmc$nchain)
  prefix <- selection[i,c('key','community','grid_index','replicate','arm','knots')]
  add_prefix <- function(d)cbind(prefix[rep(1,nrow(d)),,drop=FALSE],d)
  chains[[i]] <- add_prefix(score_observed_probability_chains(a))
  pooled[[i]] <- add_prefix(a$groups[a$groups$metric=='occupancy',c('group','bias','mae','rmse')])
  ranges[[i]] <- add_prefix(diagnose_range_chains(a$range))
}
chains <- do.call(rbind,chains);pooled <- do.call(rbind,pooled)
aggregate <- summarise_observed_chains(chains,pooled)
paired <- paired_observed_chains(chains,pooled)
check_reference <- function(observed,reference,keys) {
  joined <- merge(observed,reference,by=keys,all=TRUE)
  stopifnot(nrow(joined)==nrow(observed),nrow(joined)==nrow(reference),
    all(is.finite(joined$pooled)),all(is.finite(joined$mean)),
    max(abs(joined$pooled-joined$mean))<1e-12)
}
reference <- read.csv(file.path(summary_dir,'aggregate.csv'),stringsAsFactors=FALSE)
reference <- reference[reference$metric=='occupancy' & reference$quantity %in% c('bias','mae','rmse'),]
check_reference(aggregate,reference,c('arm','knots','group','quantity'))
reference <- read.csv(file.path(summary_dir,'paired.csv'),stringsAsFactors=FALSE)
reference <- reference[reference$quantity %in% c('bias','mae','rmse'),]
check_reference(paired,reference,c('comparison','group','quantity'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(chains,file.path(out,'observed-chain-scores.csv'),row.names=FALSE)
write.csv(aggregate,file.path(out,'observed-chain-aggregate.csv'),row.names=FALSE)
write.csv(paired,file.path(out,'observed-chain-paired.csv'),row.names=FALSE)
write.csv(do.call(rbind,ranges),file.path(out,'range-chain-diagnostics.csv'),row.names=FALSE)
writeLines(c('All nine communities retained. These are observed-chain sensitivity summaries, not confidence intervals or bounds on MCMC error.',
  'Each per-chain estimate averages posterior probabilities within that chain before calculating bias, MAE and RMSE.',
  'Chain-min/max averages choose each community minimum/maximum, then average equally across the nine communities.',
  'Paired minima use minimum A minus maximum B for the same community; maxima use maximum A minus minimum B.',
  'The pooled MAE can be smaller than every chain-only MAE because opposite chain errors can cancel.',
  'Agreement between observed chains cannot exclude an unvisited mode. Constant range chains require separate examination.'),file.path(out,'README.txt'))
cat('Inspected all 81 selected fits; pooled estimates reproduce final aggregate and paired results.\n')
