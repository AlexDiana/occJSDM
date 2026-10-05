# Coordinator runs after each stdio acceptance with a cutoff recorded before the run.
# Arguments: project root, prior directory-name JSON array or Unix cutoff, R project, optional report.
a <- commandArgs(trailingOnly=TRUE);root <- normalizePath(a[1]);cutoff <- suppressWarnings(as.numeric(a[2]))
out <- file.path(root,'reports/verification/quickstart')
Sys.setenv(P2A_R_PROJECT=normalizePath(a[3]));source(file.path(a[3],'activate.R'))
suppressPackageStartupMessages(library(occJSDM));suppressPackageStartupMessages(library(jsonlite))
stopifnot(startsWith(normalizePath(find.package('occJSDM')),paste0(normalizePath(a[3]),'/')))
sha <- unique(as.character(stats::na.omit(read.dcf(file.path(find.package('occJSDM'),'DESCRIPTION'),fields='RemoteSha')[,'RemoteSha'])))
stopifnot(identical(sha,'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'))
prior <- if(file.exists(a[2]))basename(unlist(read_json(a[2]),use.names=FALSE)) else character()
if(!file.exists(a[2]) && !is.finite(cutoff))stop('Second argument must be a prior call-directory JSON array or finite Unix cutoff')
specs <- read_json(file.path(out,'independent_inputs.json'),simplifyVector=FALSE)
exercise <- read_json(file.path(out,'mcp_exercised_results.json'),simplifyVector=FALSE)
lookup <- list()
for(nm in names(specs)) {
 if(!is.null(specs[[nm]]$error))next
 arts <- exercise[[nm]]$fit_model$artifacts
 path <- arts[[which(vapply(arts,function(x)x$kind=='fit',logical(1)))]]$path
 lookup[[path]] <- nm
}
check <- function(x,y,label) {
 z <- all.equal(x,y,tolerance=0,check.attributes=TRUE)
 if(!isTRUE(z))stop(paste(label,paste(z,collapse='; ')))
}
summary_native <- function(f) {
 x <- returnOccupancyRates(f)
 data.frame(species=colnames(x),mean=apply(x,2,mean),median=apply(x,2,quantile,.5),
 q2.5=apply(x,2,quantile,.025),q97.5=apply(x,2,quantile,.975),row.names=NULL)
}
compare_csv <- function(table,path,label) {
 tmp <- tempfile(tmpdir=out);on.exit(unlink(tmp));write.csv(table,tmp,row.names=FALSE,na='NA')
 stopifnot(identical(readBin(tmp,'raw',file.info(tmp)$size),readBin(path,'raw',file.info(path)$size)))
}
calls <- list.dirs(file.path(root,'artifacts/calls'),recursive=FALSE)
fit_checks <- list();table_checks <- list()
for(dir in calls) {
 request <- file.path(dir,'request.json')
 if(!file.exists(request))next
 if(length(prior)) {if(basename(dir) %in% prior)next} else if(as.numeric(file.info(request)$mtime)<cutoff)next
 req <- read_json(request,simplifyVector=FALSE)
 if(req$operation=='fit_model' && file.exists(file.path(dir,'manifest.json'))) {
   m <- read_json(file.path(dir,'manifest.json'),simplifyVector=FALSE)
   found <- names(specs)[vapply(specs,function(s)identical(m$input_path,s$data_path)&&identical(m$seed,s$seed),logical(1))]
   stopifnot(length(found)==1L);nm <- found[[1]]
   check(m$requested_settings,specs[[nm]]$options,paste(nm,'requested settings'))
   actual <- readRDS(m$fit_path);oracle <- readRDS(file.path(out,paste0(nm,'_oracle_fit.rds')))
   for(k in c('results_output','X_psi','X_theta','Xs','Tr','infos'))check(actual[[k]],oracle[[k]],paste(nm,k))
   fit_checks[[m$fit_id]] <- list(case=nm,fit_path=m$fit_path,full_native_match=TRUE)
 } else if(req$operation %in% c('diagnostics','summarise_fit') && file.exists(file.path(dir,'diagnostics.csv'))) {
   nm <- lookup[[req$fit_path]]
   if(is.null(nm)) {
     # The explicit constant-chain registry fixture is independently source-backed.
     f <- readRDS(req$fit_path)
     ref <- readRDS(file.path(root,'notebooks/quickstart/data/diagnostic_constant_species_fit.rds'))
     check(f,ref,'constant acceptance fit');oracle <- ref;nm <- 'constant_species'
   } else oracle <- readRDS(file.path(out,paste0(nm,'_oracle_fit.rds')))
   compare_csv(returnConvergenceDiagnostics(oracle),file.path(dir,'diagnostics.csv'),paste(nm,'diagnostics'))
   if(req$operation=='summarise_fit')compare_csv(summary_native(oracle),file.path(dir,'summary.csv'),paste(nm,'summary'))
   table_checks[[basename(dir)]] <- list(case=nm,operation=req$operation,full_native_match=TRUE)
 }
}
stopifnot(length(fit_checks)>=4L,length(table_checks)>=5L)
report <- list(success=TRUE,tolerance=0,cutoff=cutoff,R_project=Sys.getenv('P2A_R_PROJECT'),fits=fit_checks,tables=table_checks)
write_json(report,if(length(a)>=4L)a[4] else stdout(),auto_unbox=TRUE,pretty=TRUE)
