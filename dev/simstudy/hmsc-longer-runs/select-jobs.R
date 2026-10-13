#!/usr/bin/env Rscript
# Write results/jobs.csv: the Hmsc fits of the Lesson 6 extension that stayed
# flagged after their longer attempt with at least one failing coefficient
# quantity (README.md, "Which fits"). Refuses to overwrite an existing file.
#
#   Rscript select-jobs.R --repo=REPO --root=EXTENSION_ARCHIVE
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(here,'common.R'))
o <- parse_options(commandArgs(trailingOnly=TRUE),known=c('repo','root'),required=c('repo','root'))
repo <- normalizePath(o$repo,mustWork=TRUE);root <- normalizePath(o$root,mustWork=TRUE)
if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this script')
out <- file.path(repo,JOBS_FILE)
if(file.exists(out)) stop('Job selection already written: ',out)
status <- utils::read.csv(file.path(root,'status.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(status)==640L,all(status$completed))
flagged <- status[status$package=='Hmsc' & status$fit_ok & !status$diagnostic_pass,,drop=FALSE]
rows <- lapply(seq_len(nrow(flagged)),function(i) {
  j <- flagged[i,];d <- file.path(root,'jobs',j$job)
  if(!file.exists(file.path(d,'attempt-2-diagnostics.csv'))) stop('No longer attempt for ',j$job)
  t <- utils::read.csv(file.path(d,'attempt-2-diagnostics.csv'),stringsAsFactors=FALSE)
  s <- jsonlite::read_json(file.path(d,'status.json'))
  fails <- function(block) sum(!t$pass & t$block==block)
  data.frame(j[c('job','job_number','dataset','input','replicate','scenario','n_sites','n_species','response','use_traits','input_md5')],
    n_quantities=nrow(t),n_fail=sum(!t$pass),coef_fail=fails('coefficient'),cov_fail=fails('covariance'),
    prob_fail=fails('new_site_probability'),trait_fail=fails('trait'),
    max_coef_rhat=max(t$rhat[t$block=='coefficient']),min_coef_ess=min(pmin(t$bulk_ess,t$tail_ess)[t$block=='coefficient']),
    badly_unmixed=sum(!t$pass)>20L,saved_elapsed_min=as.numeric(s$elapsed_seconds)/60,stringsAsFactors=FALSE)
})
tally <- do.call(rbind,rows)
selected <- tally[tally$coef_fail>0,,drop=FALSE]
selected <- selected[order(selected$job_number),];rownames(selected) <- NULL
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
utils::write.csv(selected,out,row.names=FALSE)
cat(sprintf('Flagged Hmsc fits %d; with coefficient failures %d (badly unmixed %d; max coefficient Rhat %.3f); written to %s\n',
  nrow(tally),nrow(selected),sum(selected$badly_unmixed),max(selected$max_coef_rhat),out))
