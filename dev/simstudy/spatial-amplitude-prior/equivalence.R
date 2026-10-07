# Run once with each isolated library in separate R sessions.
a <- commandArgs(TRUE);stopifnot(length(a)==3L)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
.libPaths(c(normalizePath(a[1]),.libPaths()))
library(occJSDM);RcppParallel::setThreadOptions(numThreads=1)
input <- readRDS(a[2]);assign('.Random.seed',input$fit_rng,.GlobalEnv)
warnings <- character();started <- Sys.time()
fit <- withCallingHandlers(suppressMessages(runOccJSDM(input$data$binary,
  listParams=list(n_factors=0L,n_lattrait=0L,n_supportpoints=100L),
  threshold=1,occCovariates='environment',spatCovariates=c('longitude','latitude'),
  MCMCparams=list(nchain=2L,nburn=40L,niter=60L,nthin=1L))),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
saveRDS(list(results=fit$results_output,warnings=warnings,rng=.Random.seed,
  elapsed=as.numeric(difftime(Sys.time(),started,units='secs')),library=find.package('occJSDM'),
  input_md5=tools::md5sum(a[2])),a[3])
