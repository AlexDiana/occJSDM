S <- "/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010"
sc <- read.csv(file.path(S,"scores/screen.csv")); j <- read.csv(file.path(S,"scores/jobs.csv"))
sc <- merge(sc, j[c("job","cell","arrangement","mean_proposals")], by="job")
sc$truth <- substr(sc$cell,1,2)
sc$rhat_ok <- sc$rhat <= 1.01; sc$bulk_ok <- sc$ess_bulk >= 400; sc$tail_ok <- sc$ess_tail >= 400
agg <- function(d) c(n=nrow(d), pass=mean(d$pass), rhat_ok=mean(d$rhat_ok), bulk_ok=mean(d$bulk_ok), tail_ok=mean(d$tail_ok),
   rhat_med=median(d$rhat), bulk_med=median(d$ess_bulk), tail_med=median(d$ess_tail))
t1 <- sc[sc$truth=="t1",]
print(round(do.call(rbind, lapply(split(t1, t1$cell), agg)),3))
print(round(do.call(rbind, lapply(split(t1, t1$quantity), agg)),3))
print(round(do.call(rbind, lapply(split(t1, t1$arrangement), agg)),3))
cat("proposals by cell\n"); print(round(tapply(j$mean_proposals, j$cell, mean),2))
cat("t2 by cell/class not available; t2 by cell\n")
t2 <- sc[sc$truth=="t2",]; print(round(do.call(rbind, lapply(split(t2, t2$cell), agg)),3))
