S <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
suppressWarnings(library(posterior))
for (job in c("t1-r0.08-n100-rep1-uniform","t1-r0.08-n300-rep1-uniform","t1-r0.08-n900-rep1-uniform")) {
  r <- readRDS(file.path(S,"jobs",paste0(job,".rds"))); cat("\n==", job, " phases:", r$phases, " proposals:", round(r$mean_proposals,2), "\n")
  tr <- r$traces; ni <- dim(tr)[1]
  # drift: Rhat of first vs second half means of each chain (chain-level mean)
  half <- function(j) { m <- tr[,,j]; c(first=mean(m[1:(ni/2),]), second=mean(m[(ni/2+1):ni,])) }
  hs <- t(sapply(1:36, half)); cat("trace mean drift first->second half: max |diff| =", round(max(abs(hs[,1]-hs[,2])),4), " (trace SD median:", round(median(apply(tr,3,sd)),4), ")\n")
  # between-chain spread of chain means relative to within SD for each trace
  cm <- sapply(1:36, function(j) { m <- tr[,,j]; sd(colMeans(m))/sd(as.vector(m)) }); cat("SD of chain means / pooled SD: median", round(median(cm),3), " max", round(max(cm),3), "\n")
  # IAT per trace via ESS on pooled draws
  ess <- sapply(1:36, function(j) ess_bulk(tr[,,j])); cat("IAT (iterations) = total draws/ESS: median", round(median(ni*4/ess)), " range", round(range(ni*4/ess)), "\n")
  # per-block: thinned draws; ESS and Rhat of block mean and block SD, 400 per chain x 4
  nk <- ncol(r$site_draws[[1]])/4
  bl <- t(sapply(1:14, function(b) { d <- r$site_draws[[b]]; mm <- array(colMeans(d), c(nk,4)); ss <- array(apply(d,2,sd), c(nk,4))
     p <- array(colSums(d * rowMeans(d)) / sum(rowMeans(d)^2), c(nk,4))   # projection on block posterior mean
     c(ess_mean=ess_bulk(mm), ess_sd=ess_bulk(ss), ess_proj=ess_bulk(p), rhat_proj=rhat(p), rhat_sd=rhat(ss)) }))
  rownames(bl) <- c("shared1","shared2",paste0("indep",1:12)); print(round(bl,1))
  cat("thinned draws per chain:", nk, "(thin of stored =", ni/nk, "iterations)\n")
}
