fit_hmsc <- function(training,basis,settings) {
  validate_training(training); started <- Sys.time(); warnings <- character()
  set.seed(settings$seed); rng_start <- .Random.seed
  X <- occ_design(training,basis)
  args <- list(Y=training$y,XData=X,XFormula=~.,XScale=FALSE,distr="probit")
  if(basis$arm=="spatial2") {
    xy <- training$xy; rownames(xy) <- training$site_id
    level <- Hmsc::setPriors(Hmsc::HmscRandomLevel(sData=xy,sMethod="Full"),nfMin=2L,nfMax=2L)
    args$studyDesign <- data.frame(site=factor(training$site_id,levels=training$site_id))
    args$ranLevels <- list(site=level)
  }
  native <- withCallingHandlers({
    model <- do.call(Hmsc::Hmsc,args)
    Hmsc::sampleMcmc(model,samples=settings$mcmc[["niter"]],transient=settings$mcmc[["nburn"]],
      thin=1L,nChains=settings$mcmc[["nchain"]],nParallel=1L,verbose=0)
  },warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
  stopifnot(identical(native$Y,training$y),max(abs(native$X-cbind(1,as.matrix(X))))<1e-12,
    max(abs(native$X-native$XScaled))<1e-12)
  if(basis$arm=="spatial2") stopifnot(all(vapply(native$postList,function(ch)
    all(vapply(ch,function(d) nrow(d$Lambda[[1]])==2L,logical(1))),logical(1))))
  record_fit("Hmsc",training,basis,settings,native,started,warnings,
    list(rng_start=rng_start,training_environment=training$environment,training_xy=training$xy))
}
hmsc_prediction_level <- function(fit_record,prediction) {
  old <- fit_record$site_id; xy <- fit_record$training_xy
  ids <- prediction$site_id; newxy <- prediction$xy
  # Coordinate matches refer to existing latent variables even if IDs differ.
  d <- xy_distances(newxy,xy)
  for(i in seq_along(ids)) if(any(d[i,]==0)) ids[i] <- old[which(d[i,]==0)[1]]
  if(anyDuplicated(ids)) stop("Repeated prediction coordinates require separate rows after latent prediction")
  add <- !ids %in% old
  unionxy <- rbind(xy,newxy[add,,drop=FALSE]); rownames(unionxy) <- c(old,ids[add])
  level <- fit_record$native$ranLevels$site
  level$s <- unionxy; level$N <- nrow(unionxy); level$units <- rownames(unionxy)
  list(level=level,ids=ids,union_xy=unionxy)
}
# Same Full conditional distribution as Hmsc::predictLatentFactor, with cached
# conditioning matrices for the discrete range grid and explicit Cholesky failure.
hmsc_conditional_draws <- function(unitsPred,units,postEta,postAlpha,level) {
  old <- unitsPred %in% units; fresh <- which(!old); n <- length(unitsPred)
  cache <- new.env(parent=emptyenv())
  if(length(fresh)) {
    dt <- xy_distances(level$s[units,,drop=FALSE],level$s[units,,drop=FALSE])
    dn <- xy_distances(level$s[unitsPred[fresh],,drop=FALSE],level$s[units,,drop=FALSE])
    dnn <- xy_distances(level$s[unitsPred[fresh],,drop=FALSE],level$s[unitsPred[fresh],,drop=FALSE])
  }
  lapply(seq_along(postEta),function(it) {
    eta <- postEta[[it]]; ans <- matrix(NA_real_,n,ncol(eta),dimnames=list(unitsPred,NULL))
    ans[old,] <- eta[match(unitsPred[old],units),,drop=FALSE]
    if(length(fresh)) for(q in seq_len(ncol(eta))) {
      a <- level$alphapw[postAlpha[[it]][q],1]; key <- as.character(postAlpha[[it]][q])
      if(a==0) ans[fresh,q] <- rnorm(length(fresh)) else {
        if(!exists(key,cache,inherits=FALSE)) {
          A <- t(solve(exp(-dt/a),t(exp(-dn/a))))
          W <- exp(-dnn/a)-A %*% t(exp(-dn/a)); W <- (W+t(W))/2
          cache[[key]] <- list(A=A,L=t(chol(W)))
        }
        k <- cache[[key]]
        ans[fresh,q] <- k$A %*% eta[,q]+k$L %*% rnorm(length(fresh))
      }
    }
    ans
  })
}
hmsc_draw_extractor <- function(fit_record,prediction,basis) {
  validate_prediction_data(prediction,fit_record)
  fit <- fit_record$native; X <- as.matrix(occ_design(prediction,basis)); S <- length(fit_record$species_id)
  ni <- length(fit$postList[[1]]); nc <- length(fit$postList)
  trainX <- fit$X[,-1,drop=FALSE]; envcol <- match("environment",colnames(X))
  spcol <- setdiff(seq_len(ncol(X)),envcol)
  if(basis$arm=="spatial2") level <- hmsc_prediction_level(fit_record,prediction)
  function(idx) {
    stopifnot(length(idx)<=25L)
    dims <- c(length(idx),S,ni,nc)
    out <- setNames(lapply(1:5,function(x) array(NA_real_,dims)),c("p_full","p_space_zero","p_env_zero","p_both_zero","field"))
    for(ch in seq_len(nc)) {
      post <- fit$postList[[ch]]
      if(basis$arm=="spatial2") {
        etas <- hmsc_conditional_draws(unitsPred=level$ids[idx],units=fit_record$site_id,
          postEta=lapply(post,function(d) d$Eta[[1]]),postAlpha=lapply(post,function(d) d$Alpha[[1]]),
          level=level$level)
      }
      for(it in seq_len(ni)) {
        draw <- post[[it]]; B <- draw$Beta[-1,,drop=FALSE]
        env <- X[idx,envcol,drop=FALSE] %*% B[envcol,,drop=FALSE]
        sp <- matrix(0,length(idx),S); center <- rep(0,S)
        if(length(spcol)) {
          sp <- X[idx,spcol,drop=FALSE] %*% B[spcol,,drop=FALSE]
          center <- as.vector(colMeans(trainX)[spcol] %*% B[spcol,,drop=FALSE])
        }
        if(basis$arm=="spatial2") {
          sp <- etas[[it]] %*% draw$Lambda[[1]]
          center <- colMeans(draw$Eta[[1]] %*% draw$Lambda[[1]])
        }
        terms <- center_terms(draw$Beta[1,],env,sp,as.vector(mean(trainX[,envcol])*B[envcol,]),center)
        prob <- term_probabilities(terms,pnorm)
        for(nm in names(prob)) out[[nm]][,,it,ch] <- prob[[nm]]
        out$field[,,it,ch] <- terms$spatial
      }
    }
    out
  }
}
predict_hmsc <- function(fit_record,prediction,basis) {
  # Canonical site order fixes block composition and RNG consumption under row permutations.
  ordering <- order(prediction$site_id,method="radix"); canonical <- subset_prediction(prediction,ordering)
  seed <- fit_record$settings$prediction_seed
  if(is.null(seed)) seed <- fit_record$settings$seed
  set.seed(seed)
  out <- summarize_draw_blocks(canonical,hmsc_draw_extractor(fit_record,canonical,basis),"probit",
    "Hmsc posterior and certified Full conditional latent fields; blocks are joint internally but independent conditional simulations across blocks")
  undo <- order(ordering)
  out$site_id <- prediction$site_id
  for(nm in c("p_full","p_space_zero","p_env_zero","p_both_zero","delta_p","field_mean","field_median")) out[[nm]] <- out[[nm]][undo,,drop=FALSE]
  out$intervals <- lapply(out$intervals,function(x) x[undo,,drop=FALSE])
  out$mcse <- lapply(out$mcse,function(x) x[undo,,drop=FALSE])
  out$prediction_seed <- seed
  validate_prediction(out,prediction); out
}
