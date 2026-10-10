for (f in c("kernel.R","seeds.R","generator.R","oracle.R")) source(file.path("dev/simstudy/spatial-recoverability-pilot", f))
compile_pilot_samplers("."); Rcpp::sourceCpp("/private/tmp/claude-501/-Users-douglasyu-src-occJSDM/09b04f53-d80f-4fb3-b6aa-96ede0ede2a5/scratchpad/diag/block-ellipse.cpp")
S0 <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
cm <- readRDS(file.path(S0,"inputs","t1-r0.08-n100-rep1.rds")); idx <- cm$index$uniform; xy <- cm$points[idx,,drop=FALSE]; n <- length(idx)
f <- stable_factor(sq_exp(xy,xy,cm$range)); offset <- training_offset(cm, idx); w <- joint_weights(cm); ci <- sqrt(1-cm$shared_fraction)
tc <- sweep(cm$spatial[idx,,drop=FALSE],2,colMeans(cm$spatial[idx,,drop=FALSE]),"-"); blocks <- 14L
set.seed(1); init <- cbind(rep(0,blocks*n), do.call(rbind, replicate(blocks, f$lower %*% matrix(rnorm(n*3),n,3), simplify=FALSE)))
Y <- cm$z[idx,,drop=FALSE]
set.seed(3); a <- joint_ellipse_draws_cpp(f$lower,Y,offset,w,2L,ci,tc,init,3000L,40000L,40L)
set.seed(4); b <- block_ellipse_draws_cpp(f$lower,Y,offset,w,2L,ci,tc,init,3000L,40000L,40L)
ma <- apply(a$draws,1,mean); mb <- apply(b$draws,1,mean); sa <- apply(a$draws,1,sd); sb <- apply(b$draws,1,sd)
cat("n=100, 40000 iters x4 chains each. site-field posterior mean: cor", round(cor(ma,mb),4), " max|diff|", round(max(abs(ma-mb)),3), " (median posterior SD", round(median(sa),3), ")\n")
cat("posterior SD ratio block/joint: median", round(median(sb/sa),3), "\n")
cat("trace means (36) corr:", round(cor(apply(a$traces,3,mean), apply(b$traces,3,mean)),5), " max|diff|", round(max(abs(apply(a$traces,3,mean)-apply(b$traces,3,mean))),4), "\n")
