args <- commandArgs(trailingOnly=TRUE)
Sys.setenv(P2A_R_PROJECT=file.path(args[1],'r-runtime'))
source(file.path(args[1],'r-runtime','activate.R'))
a <- readRDS(args[2]); b <- readRDS(args[3])
# infos$execution records wall time, worker pids and BLAS, which differ between runs of the same seed (see ?runOccJSDM); everything else must match exactly.
a$infos$execution <- NULL; b$infos$execution <- NULL
for(k in c('results_output','X_psi','X_theta','Tr','infos')) {
  cmp <- all.equal(a[[k]],b[[k]],tolerance=0,check.attributes=TRUE)
  if(!isTRUE(cmp)) stop(paste(k,paste(cmp,collapse='; ')))
}
cat('Exact native equivalence verified\n')
