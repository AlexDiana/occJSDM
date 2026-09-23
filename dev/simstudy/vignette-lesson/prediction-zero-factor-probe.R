# Shape-level wrapper reproducer only. This does not fit a zero-factor model.
# Rscript prediction-zero-factor-probe.R ORIGINAL_ARCHIVE
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
archive <- normalizePath(args[1L],mustWork=TRUE)
.libPaths(c(file.path(archive,"library"),.libPaths()))
suppressPackageStartupMessages(library(occJSDM))
f <- readRDS(file.path(archive,"default-fit.rds"))$fit
f$infos$n_factors <- 0L
f$results_output$jsdm_output$L_output <-
  f$results_output$jsdm_output$L_output[FALSE,,,,drop=FALSE]
raw <- f$infos$data_info[!duplicated(f$infos$data_info$Site),
                       f$infos$list_X_psi_mat$names_df,drop=FALSE]
message <- tryCatch({
  predictNewSites(f,X_psi=as.data.frame(raw[1:3,,drop=FALSE]),useSpatial=FALSE,
                  useBiotic=FALSE,verbose=FALSE)
  NA_character_
},error=function(e)conditionMessage(e))
stopifnot(!is.na(message),grepl("perm.*wrong length",message))
cat("Verified empty-loading-array wrapper failure (not a full zero-factor fit):\n")
cat(message,"\n")
