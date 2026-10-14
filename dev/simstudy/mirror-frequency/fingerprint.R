#!/usr/bin/env Rscript
# Write results/library-fingerprint.csv for the build under ARCHIVE (README.md,
# "Library"): the revision in ARCHIVE/source-revision.txt and the md5 of every
# exported source file and installed library file. Refuses to overwrite an
# existing fingerprint. Called by build.sh.
#
#   Rscript fingerprint.R --repo=REPO --study=ARCHIVE
`%||%` <- function(x,y) if(is.null(x)) y else x
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(here,'run.R'))
PR11$extract_functions(file.path(SCRIPT_REPO,'dev/simstudy/occupancy-intercept-prior/verify-helpers.R'),
  'write_library_fingerprint',environment())
o <- parse_options(commandArgs(trailingOnly=TRUE),known=c('repo','study'),required=c('repo','study'))
repo <- normalizePath(o$repo,mustWork=TRUE)
if(!identical(repo,SCRIPT_REPO)) stop('--repo is not the checkout holding this fingerprint.R')
study <- normalizePath(o$study,mustWork=TRUE)
csv <- file.path(repo,FINGERPRINT_FILE)
if(file.exists(csv)) stop('Fingerprint already written: ',csv)
dir.create(dirname(csv),recursive=TRUE,showWarnings=FALSE)
write_library_fingerprint(study,csv)
revision <- check_library_fingerprint(study,csv)
cat('Fingerprint written for revision',revision,'to',csv,'\n')
