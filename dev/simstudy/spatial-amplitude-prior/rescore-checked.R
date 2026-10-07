#!/usr/bin/env Rscript
# Bind every robust cache to the unchanged selection/scoring dependencies.
# This complements rescore.R's existing hashes without modifying its frozen code.
args<-commandArgs(TRUE);stopifnot(length(args) %in% c(2L,3L))
dependency_repo<-normalizePath(args[1]);dependency_study<-normalizePath(args[2])
dependency_scripts<-file.path(dependency_repo,'dev/simstudy/spatial-amplitude-prior')
dependency_manifest<-file.path(dependency_scripts,'scoring-dependencies.csv')
dependency_rows<-read.csv(dependency_manifest,stringsAsFactors=FALSE)
dependency_paths<-file.path(dependency_repo,dependency_rows$file)
dependency_hashes<-tools::md5sum(dependency_paths)
stopifnot(identical(unname(dependency_hashes),dependency_rows$md5))
dependency_self<-normalizePath(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))]))
dependency_control_hashes<-tools::md5sum(c(dependency_manifest,dependency_self))
sys.source(file.path(dependency_scripts,'rescore.R'),envir=new.env(parent=globalenv()))
stopifnot(identical(unname(tools::md5sum(dependency_paths)),dependency_rows$md5),
  identical(tools::md5sum(names(dependency_control_hashes)),dependency_control_hashes))
dependency_root<-file.path(dependency_study,'robust-v1')
dependency_results<-unlist(lapply(c('inverse_gamma','half_cauchy'),function(prior)
  list.files(file.path(dependency_root,prior),pattern='-result[.]rds$',recursive=TRUE,full.names=TRUE)))
stopifnot(length(dependency_results)>0L)
dependency_record<-list(checked_at=Sys.time(),dependency_hashes=dependency_hashes,
  control_hashes=dependency_control_hashes,result_hashes=tools::md5sum(dependency_results))
dependency_dest<-file.path(dependency_root,'dependency-provenance.rds')
saveRDS(dependency_record,paste0(dependency_dest,'.tmp'))
stopifnot(file.rename(paste0(dependency_dest,'.tmp'),dependency_dest))
cat('Bound',length(dependency_results),'robust results to verified scoring/selection dependencies.\n')
