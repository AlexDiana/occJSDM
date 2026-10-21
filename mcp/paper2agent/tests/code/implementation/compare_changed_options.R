a <- commandArgs(trailingOnly=TRUE)
Sys.setenv(P2A_R_PROJECT=file.path(a[1],'r-runtime'))
source(file.path(a[1],'r-runtime','activate.R'))
library(occJSDM)
request <- jsonlite::read_json(a[3],simplifyVector=FALSE)
saved <- readRDS(file.path(a[1],'notebooks/quickstart/data/two_stage_call.rds'))
saved$args$listParams <- request$options$listParams
saved$args$threshold <- request$options$threshold
saved$args$MCMCparams <- request$options$MCMCparams
set.seed(request$seed)
expected <- do.call(runOccJSDM,saved$args)
actual <- readRDS(a[2])
for(k in c('results_output','X_psi','X_theta','Tr','infos')) {
 cmp <- all.equal(actual[[k]],expected[[k]],tolerance=0,check.attributes=TRUE)
 if(!isTRUE(cmp)) stop(paste(k,paste(cmp,collapse='; ')))
}
cat('Changed options match direct native R exactly\n')
