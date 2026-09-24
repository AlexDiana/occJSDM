# Independent observation-space integration for the repeated continuous design.
# This reads saved fits; it does not call occJSDM or run MCMC.
args <- commandArgs(TRUE)
option <- function(name) {
  hit <- args[startsWith(args,paste0("--",name,"="))]
  if (length(hit)!=1L) stop("Supply --",name,"=...")
  substring(hit,nchar(name)+4L)
}
input <- normalizePath(option("input"),mustWork=TRUE)
out <- option("out")
if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Output directory must be empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files <- file.path(input,c("fit.rds","data-truth.rds","provenance.rds"))
fit <- readRDS(files[1]); data <- readRDS(files[2]); js <- fit$results_output$jsdm_output
stopifnot(fit$infos$model=="continuous",fit$infos$n_factors==0L,
          ncol(fit$Tr)==0L,dim(js$A_output)[2]==0L)
xs <- fit$infos$list_Xs
index <- xs$Xs_index
first <- match(seq_len(nrow(xs$X_s)),index)
stopifnot(!anyNA(first),max(abs(xs$X_tilde-xs$X_s))<1e-12)
X <- fit$X_psi[first,,drop=FALSE]
stopifnot(max(abs(fit$X_psi-X[index,,drop=FALSE]))<1e-12)
counts <- tabulate(index)
y <- rowsum(data$data$OTU,index,reorder=TRUE)/counts
n <- nrow(y); S <- ncol(y)
grid <- fit$infos$l_s_grid
distance2 <- as.matrix(dist(xs$X_s))^2
spatial_covariance <- lapply(grid,function(l) {
  K <- exp(-distance2/(2*l^2))
  K %*% solve(K+diag(1e-5,n),K)
})
# Three dispersed variance states per chain. The mean parameters are
# integrated out, rather than held at the sampled or generating values.
iterations <- unique(as.integer(round(seq(1,nrow(js$idx_ls_output),length.out=3))))
results <- list()
for (ch in seq_len(ncol(js$idx_ls_output))) for (it in iterations) {
  scores <- vapply(seq_along(grid),function(j) {
    base <- matrix(1,n,n)+js$sigmab_output[it,ch]^2*tcrossprod(X)+
      js$sigmabs_output[it,ch]^2*spatial_covariance[[j]]
    sum(vapply(seq_len(S),function(s) {
      V <- base+diag(js$tau_output[s,it,ch]^2/counts,n)
      C <- chol(V)
      residual <- forwardsolve(t(C),y[,s])
      -.5*sum(residual^2)-sum(log(diag(C)))
    },numeric(1)))+dgamma(grid[j],1,1,log=TRUE)
  },numeric(1))
  weights <- exp(scores-max(scores)); weights <- weights/sum(weights)
  true_index <- which.min(abs(grid-data$truth$range))
  results[[length(results)+1L]] <- data.frame(chain=ch,iteration=it,
    true_range=data$truth$range,conditional_true_probability=weights[true_index],
    best_range=grid[which.max(scores)],
    true_vs_best_alternative_logodds=scores[true_index]-max(scores[-true_index]))
}
results <- do.call(rbind,results)
write.csv(results,file.path(out,"range-reference.csv"),row.names=FALSE)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")])
saveRDS(list(results=results,input=input,hashes=tools::md5sum(c(files,script)),
  session=sessionInfo(),method="Independent Gaussian marginal covariance of repeated-location means; current default coefficient and range priors"),
  file.path(out,"reference-provenance.rds"))
print(results)
