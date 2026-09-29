#!/usr/bin/env Rscript
# Provenance of every reused control fit (initial and longer, where present):
# which installed library produced it, read from the fit's own saved hashes.
#
#   Rscript controls.R --repo=REPO --study=STUDY [--archives=DIR] [--inputs-root=DIR] [--workers=4]
#
# Reads archived fits only. Writes STUDY/controls/control-provenance.csv and a
# compact copy to dev/simstudy/occupancy-intercept-prior/results/.
args <- commandArgs(trailingOnly=TRUE)
repo_arg <- sub('^--repo=','',grep('^--repo=',args,value=TRUE))
if(length(repo_arg)!=1L) stop('Missing or repeated --repo')
scripts <- file.path(normalizePath(repo_arg),'dev/simstudy/occupancy-intercept-prior')
source(file.path(scripts,'jobs.R'));source(file.path(scripts,'verify-helpers.R'))
o <- parse_options(args,known=c('repo','study','archives','inputs-root','workers'),required=c('repo','study'))
repo <- normalizePath(o$repo);study <- normalizePath(o$study)
archives <- normalizePath(o$archives %||% dirname(study))
inputs_root <- normalizePath(o$`inputs-root` %||% file.path(archives,'intercept-prior-inputs'))
workers <- parse_workers(o$workers %||% '4')
out <- file.path(study,'controls');dir.create(out,showWarnings=FALSE)

libraries <- do.call(rbind,lapply(names(ARCHIVES),function(name) {
  a <- archive_dir(archives,name)
  so <- file.path(a,'library/occJSDM/libs/occJSDM.so')
  data.frame(library_archive=name,library=file.path(a,'library'),so_file=so,
    so_md5=unname(tools::md5sum(so)),revision=readLines(file.path(a,'source-revision.txt')),stringsAsFactors=FALSE)
}))

candidates <- do.call(rbind,lapply(PHASES,function(phase) {
  selected <- control_schedules(phase,archives)
  keys <- phase_keys(phase)
  files <- if(phase=='A2') {
    rbind(data.frame(key=keys,schedule='initial',control_fit=file.path(archive_dir(archives,'targeted'),'initial',paste0(keys,'-fit.rds'))),
      data.frame(key=keys,schedule='long',control_fit=ifelse(
        file.exists(file.path(archive_dir(archives,'amplitude'),'inverse_gamma/long',paste0(keys,'-fit.rds'))),
        file.path(archive_dir(archives,'amplitude'),'inverse_gamma/long',paste0(keys,'-fit.rds')),
        file.path(archive_dir(archives,'targeted'),'long',paste0(keys,'-fit.rds')))))
  } else {
    rbind(data.frame(key=keys,schedule='initial',control_fit=file.path(archive_dir(archives,'pr11'),'initial',paste0(keys,'-fit.rds'))),
      data.frame(key=keys,schedule='long',control_fit=file.path(archive_dir(archives,'pr11'),'long',paste0(keys,'-fit.rds'))))
  }
  files <- files[file.exists(files$control_fit),]
  files$phase <- phase
  files$selected <- files$control_fit %in% selected$control_fit
  stopifnot(all(selected$control_fit %in% files$control_fit),sum(files$selected)==length(keys))
  files$recorded_md5 <- selected$control_fit_md5[match(files$control_fit,selected$control_fit)]
  files
}))
key_order <- unlist(lapply(PHASES,phase_keys))
candidates <- candidates[order(match(candidates$phase,PHASES),match(candidates$key,key_order),candidates$schedule),]
rownames(candidates) <- NULL

jobs <- unlist(lapply(PHASES,function(p) phase_jobs(p,archives,inputs_root)),recursive=FALSE)

inspect <- function(i) {
  row <- candidates[i,];saved <- readRDS(row$control_fit)
  hashes <- saved$fit_hashes %||% saved$source_hashes
  so <- hashes[endsWith(names(hashes),'/libs/occJSDM.so')];stopifnot(length(so)==1L)
  lib <- libraries[libraries$so_file==names(so),]
  # Production entries (exported source and installed package) must still match.
  production <- hashes[grepl('/(source|source-main|library)/',names(hashes))]
  current <- tools::md5sum(names(production))
  spec <- jobs[[row$key]]
  job_same <- if(row$phase=='A2') {
    setequal(names(saved$job),names(spec$job)) && identical(saved$job[names(spec$job)],spec$job)
  } else identical(saved$job,spec$job)
  b0 <- saved$fit$results_output$jsdm_output$B0_output
  m <- saved$mcmc
  data.frame(row,file_md5=unname(tools::md5sum(row$control_fit)),size_mb=round(file.size(row$control_fit)/2^20,1),
    job_identical=job_same,
    mcmc=paste(unlist(m),collapse='/'),mcmc_matches_schedule=isTRUE(all(unlist(m)==PROTOCOL_MCMC[[row$schedule]])),
    recorded_so=names(so),recorded_so_md5=unname(so),
    library_archive=if(nrow(lib)==1L) lib$library_archive else NA_character_,
    library_revision=if(nrow(lib)==1L) lib$revision else NA_character_,
    so_unchanged=nrow(lib)==1L && identical(lib$so_md5,unname(so)),
    production_entries=length(production),production_unchanged=identical(unname(current),unname(production)),
    warnings=length(saved$warnings),
    elapsed_seconds=round(as.numeric(difftime(saved$finished,saved$started,units='secs')),1),
    started=format(saved$started),finished=format(saved$finished),
    b0_dim=paste(dim(b0),collapse='x'),
    b0_matches_mcmc=identical(unname(dim(b0))[2:3],as.integer(c(m$niter/m$nthin,m$nchain))),
    posterior_finite=all(is.finite(b0)),
    intercept_prior_recorded=!is.null(saved$fit$infos$intercept_prior),
    stringsAsFactors=FALSE)
}
rows <- parallel::mclapply(seq_len(nrow(candidates)),function(i)
  tryCatch(inspect(i),error=function(e) stop(candidates$control_fit[i],': ',conditionMessage(e))),
  mc.cores=workers,mc.preschedule=FALSE)
failed <- !vapply(rows,is.data.frame,logical(1))
if(any(failed)) stop(paste(vapply(rows[failed],as.character,character(1)),collapse='\n'))
provenance <- do.call(rbind,rows)
provenance$recorded_md5_matches <- is.na(provenance$recorded_md5) | provenance$recorded_md5==provenance$file_md5
stopifnot(all(provenance$job_identical),all(provenance$mcmc_matches_schedule),all(!is.na(provenance$library_archive)),
  all(provenance$so_unchanged),all(provenance$production_unchanged),all(provenance$b0_matches_mcmc),
  all(provenance$posterior_finite),all(provenance$recorded_md5_matches),!any(provenance$intercept_prior_recorded))
write.csv(provenance,file.path(out,'control-provenance.csv'),row.names=FALSE)
write.csv(libraries,file.path(out,'control-libraries.csv'),row.names=FALSE)
compact <- provenance[,c('phase','key','schedule','selected','library_archive','library_revision','recorded_so_md5',
  'so_unchanged','production_unchanged','job_identical','mcmc','warnings','elapsed_seconds','b0_dim','file_md5')]
compact$control_fit <- sub(paste0('^',archives,'/'),'',provenance$control_fit)
results <- file.path(scripts,'results');dir.create(results,showWarnings=FALSE)
write.csv(compact,file.path(results,'control-provenance.csv'),row.names=FALSE)
cat('Control fits inspected:',nrow(provenance),'\n')
print(table(provenance$phase,provenance$library_archive,provenance$selected,dnn=c('phase','library','selected')))
cat('All control fits trace to unchanged archive libraries.\n')
