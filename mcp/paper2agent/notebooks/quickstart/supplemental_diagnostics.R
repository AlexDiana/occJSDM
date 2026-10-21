args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1],mustWork=TRUE)
Sys.setenv(P2A_R_PROJECT=file.path(root,"r-runtime"))
source(file.path(root,"r-runtime","activate.R"))
library(occJSDM)
library(jsonlite)
setTimeLimit(elapsed=600,transient=FALSE)
out <- file.path(root,"notebooks","quickstart")
capture <- function(expr) tryCatch(list(status="ok",value=expr),error=function(e)list(status="error",error=conditionMessage(e)))
fit <- readRDS(file.path(out,"reference","binary_fit.rds"))
# Diagnostic boundary inputs are derived from retained native arrays, analogous
# to source test-diagnostics.R's synthetic-array contract tests. No sampler repair.
one <- fit
one$results_output$jsdm_output$B0_output <- fit$results_output$jsdm_output$B0_output[,1,,drop=FALSE]
one$results_output$jsdm_output$B_output <- fit$results_output$jsdm_output$B_output[,,1,,drop=FALSE]
saveRDS(one,file.path(out,"data","diagnostic_one_draw_fit.rds"))
one_result <- capture(returnConvergenceDiagnostics(one))
constant <- fit
constant$results_output$jsdm_output$B0_output[1,,] <- 3.7
saveRDS(constant,file.path(out,"data","diagnostic_constant_species_fit.rds"))
constant_result <- capture(returnConvergenceDiagnostics(constant))
if(constant_result$status=="ok") {
 write.csv(constant_result$value,file.path(out,"reference","constant_species_diagnostics.csv"),row.names=FALSE,na="NA")
 saveRDS(constant_result$value,file.path(out,"reference","constant_species_diagnostics.rds"))
 constant_result$rhat_first_species <- constant_result$value$rhat[constant_result$value$param=="beta0_psi" & constant_result$value$idx1==1]
 constant_result$ess_first_species <- constant_result$value$ess[constant_result$value$param=="beta0_psi" & constant_result$value$idx1==1]
 constant_result$pass_rhat_unavailable <- all(is.na(constant_result$rhat_first_species))
 constant_result$value <- NULL
}
saveRDS(list(one_draw=one_result,constant_species=constant_result),file.path(out,"reference","supplemental_diagnostics_results.rds"))
write_json(list(one_draw=one_result,constant_species=constant_result),file.path(out,"reference","supplemental_diagnostics_results.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null",digits=17)
