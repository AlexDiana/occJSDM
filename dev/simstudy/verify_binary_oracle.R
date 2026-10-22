#!/usr/bin/env Rscript
# Original independent one-dimensional quadrature verification, made portable.
args <- commandArgs(TRUE)
arg <- args[startsWith(args,"--out=")]
if(length(arg)!=1L)stop("Supply --out=NEW_DIRECTORY")
out <- sub("^--out=","",arg)
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Output directory must be empty: ",out)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")])
stopifnot(length(script)==1L)
script <- normalizePath(script,mustWork=TRUE)
kernel <- file.path(dirname(script),"binary_oracle.cpp")
Rcpp::sourceCpp(kernel)
y <- c(1,0,0,1,0,0)
offset <- c(-.9,-.3,0,.1,.3,.5)
ll <- function(f)vapply(f,function(z)sum(dbinom(y,1,plogis(offset+z),log=TRUE)),numeric(1))
unnorm <- function(f)dnorm(f)*exp(ll(f))
normalizer <- integrate(unnorm,-Inf,Inf,rel.tol=1e-12)$value
expected_f <- integrate(function(f)f*unnorm(f),-Inf,Inf,rel.tol=1e-12)$value/normalizer
expected_q <- vapply(offset,function(off)integrate(function(f)plogis(off+f)*unnorm(f),-Inf,Inf,rel.tol=1e-12)$value/normalizer,numeric(1))
settings <- list(seed=2026091340L,burn=1000L,kept=40000L,thin=1L,initial_scale=3)
set.seed(settings$seed)
draw <- binary_oracle(matrix(1),matrix(y),matrix(offset),rep(0L,6),settings$burn,settings$kept,settings$thin,settings$initial_scale)
fs <- as.vector(draw$field_draws);q1 <- plogis(offset[1]+fs)
mcse <- function(z)sd(z)/sqrt(coda::effectiveSize(z))
result <- data.frame(target=c("field_mean","probability_first_row"),
  quadrature=c(expected_f,expected_q[1]),mcmc=c(mean(fs),mean(q1)),mcse=c(mcse(fs),mcse(q1)))
result$standardized_discrepancy <- with(result,(mcmc-quadrature)/mcse)
stopifnot(all(abs(result$standardized_discrepancy)<4))
write.csv(result,file.path(out,"oracle-quadrature-verification.csv"),row.names=FALSE)
saveRDS(list(settings=settings,y=y,offset=offset,results=result,draw=draw,
  hashes=tools::md5sum(c(script,kernel)),session=sessionInfo(),
  label="Original one-dimensional quadrature check; not a new spatial oracle or package fit"),file.path(out,"verification-provenance.rds"))
print(result)
