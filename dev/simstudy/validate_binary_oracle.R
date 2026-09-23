#!/usr/bin/env Rscript
# Independent known-parameter reference and historical draw summarization.
# No occJSDM package fit is called. The sibling C++ kernel is unchanged.
args <- commandArgs(TRUE)
option <- function(name,default=NULL) {
  hit <- args[startsWith(args,paste0("--",name,"="))]
  if(length(hit)>1L) stop("Duplicate option: ",name)
  if(length(hit)) substring(hit,nchar(name)+4L) else if(!is.null(default)) default else stop("Missing --",name)
}
if("--help" %in% args) {
  cat("New reference: --input=data-truth.rds --out=NEW_DIRECTORY --burn=N --iter=N --thin=N [--chains=4 --seed-base=202609130 --initial-scales=0,1,2,3]\nSummarize existing draws without sampling: --input=data-truth.rds --out=NEW_DIRECTORY --draws=EXISTING_DRAW_DIRECTORY\nOptional: --validation=PACKAGE_VALIDATION_RDS\n")
  quit(status=0)
}
integer_option <- function(name,default=NULL,minimum=1L) {
  x <- suppressWarnings(as.numeric(option(name,default)))
  if(length(x)!=1L||!is.finite(x)||x!=floor(x)||x<minimum) stop("Invalid integer: ",name)
  as.integer(x)
}
input <- normalizePath(option("input"),mustWork=TRUE)
if(dir.exists(input)) input <- normalizePath(file.path(input,"data-truth.rds"),mustWork=TRUE)
out <- option("out")
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Output directory must be empty: ",out)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")])
stopifnot(length(script)==1L)
script <- normalizePath(script,mustWork=TRUE)
kernel <- file.path(dirname(script),"binary_oracle.cpp")
stopifnot(file.exists(kernel))
existing <- option("draws","")
validation_path <- option("validation","")
d <- readRDS(input);truth <- d$truth
n <- nrow(d$data$OTU);S <- ncol(d$data$OTU)
loc <- truth$location
first <- which(!duplicated(loc))
x <- truth$standardized_coordinates[first,,drop=FALSE]
offset <- matrix(truth$B0,n,S,byrow=TRUE)+truth$standardized_environment %o% truth$B
stopifnot(max(abs(offset+truth$field-truth$eta))<1e-12,
          identical(loc[first],seq_along(first)))
provenance <- list(input=input,input_hash=tools::md5sum(input),
  script_hash=tools::md5sum(script),kernel_hash=tools::md5sum(kernel),
  command=commandArgs(),started=Sys.time(),session=sessionInfo(),
  diagnostic_method="Classical univariate coda Gelman-Rubin Rhat with autoburnin=FALSE; not rank-normalized. Stratum MCSE combines per-chain variance/ESS.")
if(nzchar(existing)) {
  existing <- normalizePath(existing,mustWork=TRUE)
  files <- list.files(existing,pattern="^true-covariance-chain[0-9]+[.]rds$",full.names=TRUE)
  if(!length(files)) stop("No existing oracle chain files in ",existing)
  id <- as.integer(sub("^true-covariance-chain([0-9]+)[.]rds$","\\1",basename(files)))
  files <- files[order(id)]
  chains <- lapply(files,readRDS)
  origin <- file.path(existing,"oracle-provenance.rds")
  origin <- if(file.exists(origin))readRDS(origin) else NULL
  if(!is.null(origin$input_hash)) stopifnot(identical(unname(origin$input_hash),unname(tools::md5sum(input))))
  provenance$mode <- "Read-only summary of existing saved draws; no sampling"
  provenance$draw_hashes <- tools::md5sum(files)
  provenance$original_provenance <- origin
  provenance$settings <- origin$settings
} else {
  nchains <- integer_option("chains","4",minimum=2L)
  settings <- list(burn=integer_option("burn",minimum=0L),kept=integer_option("iter"),
    thin=integer_option("thin"),chains=nchains,seed_base=integer_option("seed-base","202609130",0L))
  settings$initial_scale <- as.numeric(strsplit(option("initial-scales",paste(seq_len(nchains)-1L,collapse=",")),",",fixed=TRUE)[[1]])
  stopifnot(length(settings$initial_scale)==nchains,all(is.finite(settings$initial_scale)))
  settings$seeds <- settings$seed_base+seq_len(nchains)
  provenance$settings <- settings
  provenance$mode <- "Known true mean, range, spatial SD and full generating covariance"
  saveRDS(provenance,file.path(out,"oracle-provenance.rds"))
  Rcpp::sourceCpp(kernel)
  K <- exp(-as.matrix(dist(x))^2/(2*truth$range^2))+diag(1e-10,length(first))
  chains <- vector("list",nchains)
  for(ch in seq_len(nchains)) {
    set.seed(settings$seeds[ch])
    elapsed <- system.time(chains[[ch]] <- binary_oracle(t(chol(K)),d$data$OTU,offset,loc-1L,
      settings$burn,settings$kept,settings$thin,settings$initial_scale[ch]))
    saveRDS(chains[[ch]],file.path(out,paste0("true-covariance-chain",ch,".rds")))
    cat("Chain",ch,"elapsed seconds",elapsed[[3]],"\n")
  }
  files <- file.path(out,paste0("true-covariance-chain",seq_len(nchains),".rds"))
  provenance$draw_hashes <- tools::md5sum(files)
}
stopifnot(length(chains)>=2L,all(vapply(chains,function(z)
  identical(dim(z$field_draws)[1:2],c(length(first),S))&&identical(dim(z$probability_mean),c(n,S)),logical(1))))
nchains <- length(chains)
q <- Reduce("+",lapply(chains,`[[`,"probability_mean"))/nchains
f <- Reduce("+",lapply(chains,function(z)apply(z$field_draws,c(1,2),mean)))/nchains
score <- function(q) c(bias=mean(q-truth$psi),rmse=sqrt(mean((q-truth$psi)^2)),
  lower_bias=mean((q-truth$psi)[truth$psi<.2]),upper_bias=mean((q-truth$psi)[truth$psi>.8]))
validation <- if(nzchar(validation_path)) readRDS(normalizePath(validation_path,mustWork=TRUE)) else NULL
metrics <- if(is.null(validation))rbind(oracle=score(q)) else rbind(package=score(validation$psi_mean),oracle=score(q))
field_metrics <- c(field_rmse=sqrt(mean((f-truth$field_locations)^2)),
  field_slope=unname(coef(lm(as.vector(f)~as.vector(truth$field_locations)))[2]),
  package_oracle_rmse=if(is.null(validation))NA_real_ else sqrt(mean((q-validation$psi_mean)^2)))
ml <- do.call(coda::mcmc.list,lapply(chains,function(z)coda::mcmc(t(matrix(z$field_draws,nrow=length(f))))))
diagnostics <- c(max_rhat=max(coda::gelman.diag(ml,multivariate=FALSE,autoburnin=FALSE)$psrf[,1]),
  min_ess=min(coda::effectiveSize(ml)))
qtrue <- truth$psi
groups <- list(low=qtrue<.2,medium=qtrue>=.2&qtrue<=.8,high=qtrue>.8,
  bin_0_0.2=qtrue<=.2,bin_0.2_0.4=qtrue>.2&qtrue<=.4,
  bin_0.4_0.6=qtrue>.4&qtrue<=.6,bin_0.6_0.8=qtrue>.6&qtrue<=.8,bin_0.8_1=qtrue>.8,
  rare_below_0.05=qtrue<.05,rare_below_0.1=qtrue<.1,rare_above_0.9=qtrue>.9,rare_above_0.95=qtrue>.95)
for(s in seq_len(S)) {
  mask <- matrix(FALSE,n,S);mask[,s] <- TRUE
  groups[[paste0("species",s)]] <- mask
}
nn <- vapply(groups,sum,numeric(1));groups <- groups[nn>0];nn <- nn[nn>0]
traces <- lapply(chains,function(z) {
  ndraw <- dim(z$field_draws)[3]
  values <- matrix(0,ndraw,length(groups),dimnames=list(NULL,names(groups)))
  for(s in seq_len(S)) {
    qs <- plogis(sweep(z$field_draws[loc,s,],1,truth$eta[,s]-truth$field[,s],"+"))
    for(g in seq_along(groups)) {
      mask <- groups[[g]][,s]
      if(any(mask)) values[,g] <- values[,g]+colSums(qs[mask,,drop=FALSE])/nn[g]
    }
  }
  values
})
meantrue <- vapply(groups,function(i)mean(qtrue[i]),numeric(1))
meanestimate <- vapply(groups,function(i)mean(q[i]),numeric(1))
chainvars <- do.call(rbind,lapply(traces,function(z)apply(z,2,var)))
ess <- do.call(rbind,lapply(traces,function(z)coda::effectiveSize(coda::mcmc(z))))
mcse <- sqrt(colSums(chainvars/ess))/nchains
ml <- do.call(coda::mcmc.list,lapply(traces,coda::mcmc))
rhat <- coda::gelman.diag(ml,multivariate=FALSE,autoburnin=FALSE)$psrf[,1]
strata <- data.frame(stratum=names(groups),n=nn,true_probability=meantrue,
  estimated_probability=meanestimate,error=meanestimate-meantrue,mcse=mcse,ess=colSums(ess),rhat=rhat)
stopifnot(max(abs(colMeans(do.call(rbind,traces))-meanestimate))<1e-12)
saveRDS(list(q=q,f=f,metrics=metrics,diagnostics=diagnostics,field_metrics=field_metrics),file.path(out,"oracle-summary.rds"))
write.csv(metrics,file.path(out,"oracle-metrics.csv"))
write.csv(data.frame(metric=names(field_metrics),value=field_metrics),file.path(out,"field-metrics.csv"),row.names=FALSE)
write.csv(strata,file.path(out,"oracle-strata.csv"),row.names=FALSE)
saveRDS(traces,file.path(out,"probability-strata-traces.rds"))
write.csv(data.frame(species=seq_len(S),true_prevalence=colMeans(qtrue),observed_prevalence=colMeans(d$data$OTU)),file.path(out,"species-prevalence.csv"),row.names=FALSE)
provenance$finished <- Sys.time()
saveRDS(provenance,file.path(out,"oracle-provenance.rds"))
writeLines(c("# Independent binary oracle", "",provenance$mode,"",
  "This estimates the posterior under the known generating parameters. It is not a complete package fit or independent repeated-dataset bias validation.","",
  provenance$diagnostic_method,"",
  "The sampler uses elliptical slice sampling: https://proceedings.mlr.press/v9/murray10a.html. All paths, hashes, settings and available originating provenance are retained in oracle-provenance.rds.","",
  "See oracle-metrics.csv and oracle-strata.csv. Monte Carlo uncertainty is separate from uncertainty across independent simulated datasets."),file.path(out,"README.md"))
print(metrics);print(diagnostics);print(strata)
