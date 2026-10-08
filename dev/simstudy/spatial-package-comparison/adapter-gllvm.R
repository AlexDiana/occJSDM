conditional_spatial_moments <- function(train_mean,train_cov,Ktt,Knt,Knn) {
  A <- t(solve(Ktt,t(Knt)))
  covariance <- Knn-A %*% t(Knt)+A %*% train_cov %*% t(A)
  covariance <- (covariance+t(covariance))/2
  if(min(eigen(covariance,symmetric=TRUE,only.values=TRUE)$values)< -1e-8) stop("Conditional covariance is not positive semidefinite")
  list(mean=as.vector(A %*% train_mean),covariance=covariance,A=A,
       covariance_with_training=A %*% train_cov)
}
fit_gllvm <- function(training,basis,settings) {
  validate_training(training); started <- Sys.time(); warnings <- character()
  set.seed(settings$seed); rng_start <- .Random.seed
  X <- occ_design(training,basis)
  args <- list(y=training$y,X=X,family=binomial("logit"),link="logit",
    num.lv=if(basis$arm=="spatial2") 2L else 0L,method="VA",row.eff=FALSE,
    scale.X=FALSE,sd.errors=FALSE,seed=settings$seed,
    control=list(maxit=settings$gllvm$maxit,max.iter=settings$gllvm$maxit),
    control.start=list(starting.val="res",n.init=1L,jitter.var=.1))
  if(basis$arm=="spatial2") {
    args$lvCor <- ~corExp(1|site)
    args$studyDesign <- data.frame(site=factor(training$site_id,levels=training$site_id))
    args$distLV <- training$xy
  }
  fit <- withCallingHandlers(do.call(gllvm::gllvm,args),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  })
  stopifnot(identical(fit$y,training$y),max(abs(fit$X-as.matrix(X)))<1e-12)
  param <- fit$TMBfn$env$last.par.best
  gradient <- fit$TMBfn$gr(param); objective <- fit$TMBfn$fn(param)
  if(basis$arm=="spatial2") {
    stopifnot(fit$num.lv==2L,fit$corP$cstruclv=="corExp",fit$Lambda.struc=="unstructured",
      identical(dim(fit$A),c(nrow(training$y),nrow(training$y),2L)),
      max(abs(as.matrix(fit$TMBfn$env$data$dLV)-diag(nrow(training$y))))<1e-12,
      max(abs(fit$TMBfn$env$data$dc_lv-xy_distances(training$xy,training$xy)))<1e-12)
  }
  # Retain native fit plus numerical moments; no predictor depends on a live TMB pointer.
  native <- list(fit=fit,params=fit$params,lvs=fit$lvs,A=fit$A,
    objective=objective,gradient=gradient,convergence=fit$convergence,
    covariance_status=if(basis$arm=="spatial2") "variational covariance available; fixed parameter uncertainty omitted" else "fixed parameter uncertainty not estimated",
    kernel="exp(-distance/rho.lv); raw unit-square Euclidean distance; unit marginal latent variance before sigma.lv")
  record_fit("gllvm",training,basis,settings,native,started,warnings,
    list(rng_start=rng_start,training_environment=training$environment,training_xy=training$xy))
}
gllvm_field_moments <- function(fit_record,prediction) {
  native <- fit_record$native; xy <- fit_record$training_xy; n <- nrow(xy)
  L <- sweep(native$params$theta,2,native$params$sigma.lv,"*")
  N <- nrow(prediction$xy); S <- nrow(L)
  mean <- variance <- centered_variance <- matrix(0,N,S)
  center <- center_var <- numeric(S)
  for(q in 1:2) {
    range <- native$params$rho.lv[q]
    Ktt <- exp(-xy_distances(xy,xy)/range)
    Knt <- exp(-xy_distances(prediction$xy,xy)/range)
    Knn <- exp(-xy_distances(prediction$xy,prediction$xy)/range)
    C <- native$A[,,q]; mu <- native$lvs[,q]
    conditional <- conditional_spatial_moments(mu,C,Ktt,Knt,Knn)
    vc <- mean(C); cross <- rowMeans(conditional$covariance_with_training)
    cv <- pmax(diag(conditional$covariance)+vc-2*cross,0)
    mean <- mean+conditional$mean %o% L[,q]
    variance <- variance+diag(conditional$covariance) %o% (L[,q]^2)
    centered_variance <- centered_variance+cv %o% (L[,q]^2)
    center <- center+mean(mu)*L[,q]; center_var <- center_var+vc*L[,q]^2
  }
  list(mean=mean,variance=pmax(variance,0),center=center,center_variance=pmax(center_var,0),centered_variance=centered_variance)
}
predict_gllvm <- function(fit_record,prediction,basis) {
  validate_prediction_data(prediction,fit_record)
  X <- as.matrix(occ_design(prediction,basis)); B <- t(fit_record$native$params$Xcoef)
  intercept <- fit_record$native$params$beta0; S <- length(intercept); n <- nrow(X)
  E <- X[,1,drop=FALSE] %*% B[1,,drop=FALSE]
  sp <- matrix(0,n,S); Sc <- rep(0,S); ec <- rep(0,S)
  if(ncol(X)>1L) sp <- X[,-1,drop=FALSE] %*% B[-1,,drop=FALSE]
  native_moments <- NULL
  if(basis$arm=="spatial2") {
    # No lattice-sized covariance allocation: use <=25-site blocks.
    blocks <- split(seq_len(n),ceiling(seq_len(n)/25L))
    native_moments <- lapply(blocks,function(idx) gllvm_field_moments(fit_record,subset_prediction(prediction,idx)))
    sp <- do.call(rbind,lapply(native_moments,`[[`,"mean")); Sc <- native_moments[[1]]$center
  }
  terms <- center_terms(intercept,E,sp,ec,Sc)
  out <- point_prediction(prediction,terms,plogis,"logit","Fitted coefficients; variational and conditional spatial uncertainty integrated when spatial factors are present")
  if(!is.null(native_moments)) {
    variance <- do.call(rbind,lapply(native_moments,`[[`,"variance"))
    sd <- sqrt(variance); center_sd <- matrix(rep(sqrt(native_moments[[1]]$center_variance),each=n),n,S)
    out$p_full <- normal_integral(terms$intercept+terms$environment+terms$spatial,sd)
    out$p_env_zero <- normal_integral(terms$intercept+terms$spatial,sd)
    out$p_space_zero <- normal_integral(terms$intercept+terms$environment,center_sd)
    out$p_both_zero <- normal_integral(terms$intercept,center_sd)
    # Compare a higher-order quadrature before reporting supported probabilities.
    checks <- c(full=max(abs(out$p_full-normal_integral(terms$intercept+terms$environment+terms$spatial,sd,n=100L))),
      zero=max(abs(out$p_space_zero-normal_integral(terms$intercept+terms$environment,center_sd,n=100L))))
    if(max(checks)>1e-6) stop("gllvm Gaussian quadrature refinement failed")
    out$integration_checks <- checks
    field_sd <- sqrt(do.call(rbind,lapply(native_moments,`[[`,"centered_variance")))
    mu <- terms$intercept+terms$environment+terms$spatial
    out$intervals <- list(probability_lower=plogis(mu-qnorm(.975)*sd),probability_upper=plogis(mu+qnorm(.975)*sd),
      field_lower=terms$spatial-qnorm(.975)*field_sd,field_upper=terms$spatial+qnorm(.975)*field_sd)
    out$uncertainty <- list(status="available",reason="Variational training latent covariance and residual conditional new-site covariance; fixed coefficient, loading and range uncertainty omitted")
    out$delta_p <- out$p_full-out$p_space_zero
  }
  validate_prediction(out,prediction); out
}
