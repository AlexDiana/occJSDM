# Independent, deterministic likelihood calculations for the saved study.
# Oracle controls fix nuisance parameters at simulation truth; no model fits changed.
logadd <- function(a,b) {
  m <- pmax(a,b); ans <- m+log(exp(a-m)+exp(b-m));ans[is.infinite(m)&m<0] <- -Inf;ans
}
# Return log Pr(data_i | z_i=0/1), with field samples conditionally independent.
site_likelihood <- function(input,s,arm,collection_intercept=NULL) {
  t <- input$truth;n <- input$settings$n; M <- input$settings$M
  theta <- if(is.null(collection_intercept)) t$theta[,s] else
    plogis(collection_intercept+t$Xt[,2]*t$beta_theta[2,s])
  if(arm=='binary') return(cbind(ifelse(t$z[,s]==0,0,-Inf),ifelse(t$z[,s]==1,0,-Inf)))
  if(arm=='field') {
    w <- t$w[,s];l0 <- dbinom(w,1,t$theta0[s],log=TRUE);l1 <- dbinom(w,1,theta,log=TRUE)
  } else {
    d <- input$data[[arm]]; y <- d$OTU[,s];info <- d$info
    lp <- rowsum(dbinom(y,1,t$p[info$Primer,s],log=TRUE),info$Sample,reorder=FALSE)[,1]
    lq <- rowsum(dbinom(y,1,t$q[[arm]][info$Primer,s],log=TRUE),info$Sample,reorder=FALSE)[,1]
    l0 <- logadd(log(t$theta0[s])+lp,log1p(-t$theta0[s])+lq)
    l1 <- logadd(log(theta)+lp,log1p(-theta)+lq)
  }
  cbind(colSums(matrix(l0,M,n)),colSums(matrix(l1,M,n)))
}
intercept_posterior <- function(site_loglik,offset,prior_sd=1,step=.025) {
  b <- seq(-24,16,by=step)
  eta <- outer(offset,b,'+'); pr <- plogis(eta)
  logp <- plogis(eta,log.p=TRUE);log1p <- plogis(eta,lower.tail=FALSE,log.p=TRUE)
  ll <- colSums(logadd(log1p+site_loglik[,1],logp+site_loglik[,2]))
  weights <- exp(ll+dnorm(b,0,prior_sd,log=TRUE)-max(ll+dnorm(b,0,prior_sd,log=TRUE)))
  weights[c(1,length(b))] <- weights[c(1,length(b))]/2
  weights <- weights/sum(weights)
  list(mean=sum(colMeans(pr)*weights),mean_b=sum(b*weights),
    lower=approx(cumsum(weights),colMeans(pr),xout=.025,ties='ordered')$y,
    upper=approx(cumsum(weights),colMeans(pr),xout=.975,ties='ordered')$y,
    boundary=sum(weights[b<=-23|b>=15]),ll=ll,b=b,weights=weights)
}
# Gauss-Legendre integration is unnecessary here: a regular two-intercept grid
# allows a direct resolution check, including the full tails of both priors.
collection_marginal <- function(input,s,arm,step=.1) {
  offset <- input$truth$X*input$truth$B[s]+input$truth$field[,s]
  b <- seq(-24,16,by=step);ct <- seq(-12,10,by=step)
  sites <- lapply(ct,function(c)site_likelihood(input,s,arm,c))
  l0 <- vapply(sites,function(x)x[,1],numeric(input$settings$n))
  l1 <- vapply(sites,function(x)x[,2],numeric(input$settings$n))
  ll <- collapsed_grid(b,offset,l0,l1)
  logw <- sweep(ll,2,dnorm(ct,0,1,log=TRUE),'+')
  measures <- lapply(c(1,2.5),function(sd) {
    v <- sweep(logw,1,dnorm(b,0,sd,log=TRUE),'+');w <- exp(v-max(v))
    w[c(1,nrow(w)),] <- w[c(1,nrow(w)),]/2;w[,c(1,ncol(w))] <- w[,c(1,ncol(w))]/2;w <- w/sum(w)
    pr <- rowMeans(plogis(outer(b,offset,'+')))
    data.frame(prior_sd=sd,estimate=sum(rowSums(w)*pr),
      collection_mean=sum(colSums(w)*vapply(ct,function(c)mean(plogis(c+input$truth$Xt[,2]*input$truth$beta_theta[2,s])),numeric(1))),
      boundary=sum(w[b<=-23|b>=15,])+sum(w[,ct<=-11|ct>=9]))
  })
  do.call(rbind,measures)
}
