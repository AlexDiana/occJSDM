# Observation-space integration of pseudo-data kappa/Omega under
# Normal(Z*m, diag(1/Omega) + Z*D^-1*Z'). It is independent of the
# coefficient-space block crossproduct formulas used by both helpers.
grouped_range_reference <- function(a) {
  n <- nrow(a$X);p <- ncol(a$X);d <- ncol(a$U);ps <- nrow(a$M_Bs);S <- ncol(a$Omega)
  vv <- c(1,rep(a$sigma_b^2,p),rep(1,d),rep(a$sigma_bs^2,ps))
  mm <- rbind(rep(0,S),a$M_B,matrix(0,d,S),a$M_Bs)
  vapply(seq_along(a$list_SoRSummaries$l_s_grid),function(j) {
    H <- matrix(0,n,ps)
    Ks <- matrix(a$list_SoRSummaries$Ks_all[,,j],nrow=n)
    for(i in seq_len(n)) if(ps) H[i,a$Xs_centers[i,]] <- Ks[i,]
    Z <- cbind(1,a$X,a$U,H)
    ans <- 0
    for(s in seq_len(S)) {
      V <- diag(1/a$Omega[,s],n)+tcrossprod(sweep(Z,2,sqrt(vv),"*"))
      residual <- a$kappa[,s]/a$Omega[,s]-drop(Z%*%mm[,s])
      C <- chol(V)
      z <- forwardsolve(t(C),residual)
      ans <- ans-sum(log(diag(C)))-.5*sum(z^2)
    }
    ans+dgamma(a$list_SoRSummaries$l_s_grid[j],a$a_l_s,a$b_l_s,log=TRUE)
  },numeric(1))
}

grouped_range_fixture <- function(p,d,ps=3L,S=3L,permuted=FALSE,constant=FALSE) {
  n <- 17L; m <- 5L; grid <- c(.03,.09,.16,.27)
  group <- sample(rep(seq_len(m),length.out=n))
  ids <- c("location_Q","location_7","location_alpha","location_99","location_2")[group]
  centers <- matrix(rep(seq_len(ps),each=n),n,ps)
  if(permuted&&ps) for(i in seq_len(n)) centers[i,] <- sample(seq_len(ps))
  bases <- array(0,c(n,ps,length(grid)))
  for(j in seq_along(grid)) {
    hh <- matrix(rnorm(m*ps,sd=j/2),m,ps)
    if(ps) for(i in seq_len(n)) bases[i,,j] <- hh[group[i],centers[i,]]
  }
  Omega <- matrix(exp(rnorm(n*S,sd=.8)),n,S)
  if(constant) Omega[,] <- rep(seq(.2,1,length.out=S),each=n)
  list(X=matrix(rnorm(n*p),n,p),U=matrix(rnorm(n*d),n,d),
       M_B=matrix(rnorm(p*S,mean=.4),p,S),M_Bs=matrix(rnorm(ps*S,mean=-.3),ps,S),
       sigma_b=.6,sigma_bs=1.2,kappa=matrix(rnorm(n*S),n,S),Omega=Omega,
       Xs_centers=centers,list_SoRSummaries=list(Ks_all=bases,l_s_grid=grid),
       a_l_s=2.3,b_l_s=3.4,location=ids)
}

test_that("repeated locations preserve the independently integrated range probabilities", {
  set.seed(2026091361)
  cases <- list(grouped_range_fixture(2,3),grouped_range_fixture(1,0),
    grouped_range_fixture(0,1),grouped_range_fixture(0,0),grouped_range_fixture(1,1),
    grouped_range_fixture(1,1,1,1),grouped_range_fixture(1,1,0),
    grouped_range_fixture(2,1,3,3,TRUE),grouped_range_fixture(1,1,constant=TRUE))
  for (a in cases) {
    expected <- grouped_range_reference(a)
    actual <- do.call(spatial_range_logweights,a)
    expect_equal(actual-actual[1],expected-expected[1],tolerance=1e-10)
    # Compare the repeated-location path with the ungrouped arithmetic too.
    dense <- do.call(spatial_range_logweights,a[names(a)!="location"])
    expect_equal(actual,dense,tolerance=1e-10)
    perm <- sample(seq_len(nrow(a$X)))
    for (term in c("X","U","kappa","Omega","Xs_centers"))
      a[[term]] <- a[[term]][perm,,drop=FALSE]
    a$location <- a$location[perm]
    a$list_SoRSummaries$Ks_all <- a$list_SoRSummaries$Ks_all[perm,,,drop=FALSE]
    expect_equal(do.call(spatial_range_logweights,a),actual,tolerance=1e-10)
  }
})

test_that("grouping cannot silently combine different spatial designs", {
  set.seed(792)
  a <- grouped_range_fixture(1,0)
  different <- which(a$location!=a$location[1])[1]
  a$location[1] <- a$location[different]
  expect_error(do.call(spatial_range_logweights,a),"spatial design rows differ")
  a$location[1] <- NA_character_
  expect_error(do.call(spatial_range_logweights,a),"location")
  a$location <- "one label"
  expect_error(do.call(spatial_range_logweights,a),"location")
})
