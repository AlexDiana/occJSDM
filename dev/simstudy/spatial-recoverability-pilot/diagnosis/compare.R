args <- commandArgs(TRUE); job <- args[1]; nburn <- as.integer(args[2]); niter <- as.integer(args[3])
for (f in c("kernel.R","seeds.R","generator.R","oracle.R")) source(file.path("dev/simstudy/spatial-recoverability-pilot", f))
suppressWarnings(library(posterior))
compile_pilot_samplers(".")
Rcpp::sourceCpp("/private/tmp/claude-501/-Users-douglasyu-src-occJSDM/09b04f53-d80f-4fb3-b6aa-96ede0ede2a5/scratchpad/diag/block-ellipse.cpp")
S0 <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
cm <- readRDS(file.path(S0,"inputs",paste0(sub("-(uniform|clustered)$","",job),".rds"))); arr <- sub(".*-(uniform|clustered)$","\\1",job)
idx <- cm$index[[arr]]; xy <- cm$points[idx,,drop=FALSE]; n <- length(idx); Sp <- nrow(cm$species)
f <- stable_factor(sq_exp(xy,xy,cm$range)); offset <- training_offset(cm, idx); w <- joint_weights(cm); ci <- sqrt(1-cm$shared_fraction)
truth <- cm$spatial[idx,,drop=FALSE]; tc <- sweep(truth,2,colMeans(truth),"-"); blocks <- 2L+Sp; nchain <- 4L
set.seed(1); initial <- cbind(rep(0,blocks*n), do.call(rbind, replicate(blocks, f$lower %*% matrix(rnorm(n*3),n,3), simplify=FALSE)))
Y <- cm$z[idx,,drop=FALSE]
summ <- function(tr) { e <- sapply(1:36, function(j) ess_bulk(tr[,,j])); r <- sapply(1:36, function(j) rhat(tr[,,j])); c(ess_med=median(e), ess_min=min(e), rhat_med=median(r), rhat_max=max(r), frac_pass=mean(e>=400 & r<=1.01)) }
set.seed(2); t0 <- Sys.time(); a <- joint_ellipse_draws_cpp(f$lower, Y, offset, w, 2L, ci, tc, initial, nburn, niter, 10L); ta <- as.numeric(Sys.time()-t0, units="secs")
set.seed(2); t0 <- Sys.time(); b <- block_ellipse_draws_cpp(f$lower, Y, offset, w, 2L, ci, tc, initial, nburn, niter, 10L); tb <- as.numeric(Sys.time()-t0, units="secs")
cat(job, "nburn", nburn, "niter", niter, "\n")
cat("joint : secs", round(ta), " proposals/step", round(mean(a$mean_proposals),2), "\n"); print(round(summ(a$traces),3))
cat("block : secs", round(tb), " proposals/block-step", round(mean(b$mean_proposals),2), "\n"); print(round(summ(b$traces),3))
cat("ESS per second  joint:", round(median(sapply(1:36,function(j) ess_bulk(a$traces[,,j])))/ta,2), " block:", round(median(sapply(1:36,function(j) ess_bulk(b$traces[,,j])))/tb,2), "\n")
