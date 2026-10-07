# One isolated activation batches all reference and fresh native comparisons.
a <- commandArgs(trailingOnly=TRUE); root <- normalizePath(a[1]);out <- file.path(root,'reports/verification/quickstart')
Sys.setenv(P2A_R_PROJECT=file.path(root,'r-runtime'));source(file.path(root,'r-runtime','activate.R'))
suppressPackageStartupMessages(library(occJSDM));suppressPackageStartupMessages(library(jsonlite))
check <- function(x,y,label) {
 c <- all.equal(x,y,tolerance=0,check.attributes=TRUE)
 if(!isTRUE(c)) stop(paste(label,paste(c,collapse='; ')))
}
equal_fit <- function(x,y,label) {
 for(k in c('results_output','X_psi','X_theta','Xs','Tr','infos'))check(x[[k]],y[[k]],paste(label,k))
}
summary_native <- function(f) {
 r <- returnOccupancyRates(f)
 stopifnot(is.matrix(r),ncol(r)==length(f$infos$speciesNames),identical(colnames(r),f$infos$speciesNames))
 data.frame(species=colnames(r),mean=apply(r,2,mean),median=apply(r,2,quantile,.5),
            q2.5=apply(r,2,quantile,.025),q97.5=apply(r,2,quantile,.975),row.names=NULL)
}
csv_compare <- function(table,reference,label) {
 file <- file.path(out,paste0(label,'_recomputed.csv'));write.csv(table,file,row.names=FALSE,na='NA')
 stopifnot(identical(readBin(file,'raw',file.info(file)$size),readBin(reference,'raw',file.info(reference)$size)))
}
# Test installed function bodies/formals against the exact pinned scientific source.
env <- new.env(parent=asNamespace('occJSDM'))
for(file in c('runOccJSDM.R','jsdmfun.R','diagnostics.R','output.R'))sys.source(file.path(root,'repo/occJSDM/R',file),env)
for(nm in c('inferDataModel','process_covariates','create_covariates_matrix','runOccJSDM','returnConvergenceDiagnostics','returnOccupancyRates')) {
 native <- get(nm,asNamespace('occJSDM'));pinned <- get(nm,env)
 check(body(native),body(pinned),paste(nm,'body'));check(formals(native),formals(pinned),paste(nm,'formals'))
}
saved <- read_json(file.path(out,'saved_manifest_selection.json'),simplifyVector=FALSE)
passed <- list()
for(nm in names(saved)) {
 manifest <- saved[[nm]];actual <- readRDS(manifest$fit_path)
 ref <- readRDS(file.path(root,'notebooks/quickstart/reference',paste0(nm,'_fit.rds')))
 call <- readRDS(file.path(root,'notebooks/quickstart/data',paste0(nm,'_call.rds')))
 set.seed(call$seed);direct <- do.call(runOccJSDM,call$args)
 equal_fit(actual,ref,paste(nm,'wrapper/reference'))
 equal_fit(direct,ref,paste(nm,'fresh direct/reference'))
 dg <- returnConvergenceDiagnostics(direct);rates <- returnOccupancyRates(direct)
 check(dg,readRDS(file.path(root,'notebooks/quickstart/reference',paste0(nm,'_diagnostics.rds'))),paste(nm,'all diagnostic values'))
 check(rates,readRDS(file.path(root,'notebooks/quickstart/reference',paste0(nm,'_occupancy_draws.rds'))),paste(nm,'all occupancy draws'))
 csv_compare(summary_native(direct),file.path(root,'notebooks/quickstart/reference',paste0(nm,'_occupancy_summary.csv')),paste0(nm,'_saved_summary'))
 passed[[nm]] <- list(full_arrays=TRUE,dimensions_identifiers_classes=TRUE,diagnostics=TRUE,baseline=TRUE)
}
changed <- read_json(file.path(out,'mcp_exercised_results.json'),simplifyVector=FALSE)
constant <- readRDS(file.path(root,'notebooks/quickstart/data/diagnostic_constant_species_fit.rds'))
check(returnConvergenceDiagnostics(constant),readRDS(file.path(root,'notebooks/quickstart/reference/constant_species_diagnostics.rds')),'constant native full diagnostic table')
bundled <- readRDS(file.path(root,'notebooks/quickstart/reference/bundled_sampleresults_fit.rds'))
check(returnConvergenceDiagnostics(bundled),readRDS(file.path(root,'notebooks/quickstart/reference/bundled_sampleresults_diagnostics.rds')),'bundled native diagnostics')
csv_compare(summary_native(bundled),file.path(root,'notebooks/quickstart/reference/bundled_sampleresults_occupancy_summary.csv'),'bundled_native_summary')
fresh <- list()
for(nm in names(changed)) {
 rec <- changed[[nm]];if(!is.null(rec$native_error))next
 artifacts <- rec$fit_model$artifacts
 fp <- artifacts[[which(vapply(artifacts,function(x)x$kind=='fit',logical(1)))]]$path
 actual <- readRDS(fp);oracle <- readRDS(file.path(out,paste0(nm,'_oracle_fit.rds')))
 equal_fit(actual,oracle,paste(nm,'changed input direct oracle'))
 fresh[[nm]] <- list(full_arrays=TRUE,dimensions_identifiers_classes=TRUE,
   B0_dimensions=as.list(dim(actual$results_output$jsdm_output$B0_output)),
   baseline_dimensions=as.list(dim(returnOccupancyRates(actual))))
}
write_json(list(all_passed=TRUE,installed_source_bodies_identical=TRUE,saved_reference_cases=passed,changed_cases=fresh,
 tolerance=0,nonfinite='Native R all.equal preserves NA/NaN locations and attributes; diagnostic CSV NA remains null in JSON'),
 file.path(out,'numerical_audit.json'),auto_unbox=TRUE,pretty=TRUE)
cat('All saved references and fresh changed native artifacts agree exactly\n')
