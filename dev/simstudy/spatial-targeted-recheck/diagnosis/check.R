i <- readRDS(file.path(study,'inputs/range6-rep01.rds'))
# Independent enumeration of z,w1,w2 for each of first three sites.
for(s in c(1,3,8)) for(arm in c('low','high'))for(site in 1:3) {
  L <- site_likelihood(i,s,arm);d <- i$data[[arm]];sm <- ((site-1)*2+1):(site*2)
  for(z in 0:1) {
    direct <- 0
    for(w1 in 0:1)for(w2 in 0:1) {
      w <- c(w1,w2);term <- 1
      for(m in 1:2) {
        ids <- which(d$info$Sample==sm[m]);theta <- if(z) i$truth$theta[sm[m],s] else i$truth$theta0[s]
        rate <- if(w[m]) i$truth$p[d$info$Primer[ids],s] else i$truth$q[[arm]][d$info$Primer[ids],s]
        term <- term*dbinom(w[m],1,theta)*prod(dbinom(d$OTU[ids,s],1,rate))
      }
      direct <- direct+term
    }
    stopifnot(abs(log(direct)-L[site,z+1])<1e-10)
  }
}
L <- site_likelihood(i,1,'low');b <- c(-10,-5,0,2);off <- i$truth$X*.4+i$truth$field[,1]
cpp <- collapsed_grid(b,off,matrix(L[,1]),matrix(L[,2]))[,1]
ref <- vapply(b,function(a)sum(logadd(plogis(a+off,lower.tail=FALSE,log.p=TRUE)+L[,1],plogis(a+off,log.p=TRUE)+L[,2])),numeric(1))
stopifnot(max(abs(cpp-ref))<1e-10)
# Exact binomial integral gives the same binary calculation up to a likelihood constant.
for(k in c(0,1,5))for(sd in c(1,2.5)) {
  L <- cbind(ifelse(seq_len(100)<=k,-Inf,0),ifelse(seq_len(100)<=k,0,-Inf))
  q <- intercept_posterior(L,rep(0,100),sd)
  f <- function(b) exp(dbinom(k,100,plogis(b),log=TRUE)+dnorm(b,0,sd,log=TRUE))
  exact <- integrate(function(b)plogis(b)*f(b),-20,12,rel.tol=1e-10)$value/integrate(f,-20,12,rel.tol=1e-10)$value
  stopifnot(abs(q$mean-exact)<1e-8,q$boundary<1e-8)
}
cat('Independent enumeration, native grid/R agreement, and binomial integrals passed.\n')
