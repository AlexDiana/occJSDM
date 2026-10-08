fit_oracle <- function(training,input,arrangement,species,settings,source_root,cache) {
  validate_training(training); started <- Sys.time(); set.seed(settings$seed); rng_start <- .Random.seed
  # Original sampler and dependencies are frozen byte-for-byte, never rewritten here.
  env <- new.env(parent=globalenv())
  sys.source(file.path(source_root,"dev/simstudy/spatial-design-sweep/generator.R"),env)
  sys.source(file.path(source_root,"dev/simstudy/spatial-design-sweep/oracle.R"),env)
  env$compile_oracle(source_root,cache)
  native <- env$run_oracle(input,arrangement,species,settings$seed,
    nburn=settings$mcmc[["nburn"]],niter=settings$mcmc[["niter"]],nchain=settings$mcmc[["nchain"]])
  native$B0 <- input$landscape$B0[species]; native$B <- input$landscape$B[species]; native$range <- input$landscape$range
  # The oracle knows nuisance parameters but receives no realised field or held-out outcomes.
  list(schema_version=1L,package="oracle",arm="oracle",dataset_id=training$dataset_id,
    site_id=training$site_id,species_id=training$species_id,input_hashes=training$input_hashes,
    settings=settings,source_hashes=settings$source_hashes,native=native,warnings=character(),
    started=started,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),native_status="returned",
    runtime=list(R=R.version.string),rng_start=rng_start,rng_state=.Random.seed,
    training_xy=training$xy,training_environment=training$environment)
}
oracle_conditional <- function(fit_record,prediction) {
  f <- fit_record; kernel <- function(a,b) exp(-xy_distances(a,b)^2/(2*f$native$range^2))
  Ktt <- kernel(f$training_xy,f$training_xy)+diag(1e-8,nrow(f$training_xy))
  Knt <- kernel(prediction$xy,f$training_xy)
  A <- t(solve(Ktt,t(Knt)))
  V <- kernel(prediction$xy,prediction$xy)-A %*% t(Knt)+diag(1e-8,nrow(prediction$xy))
  # Exact training coordinates denote the already sampled latent variables.
  matches <- match(apply(prediction$xy,1,paste,collapse="|"),apply(f$training_xy,1,paste,collapse="|"))
  same <- which(!is.na(matches))
  if(length(same)) { A[same,] <- diag(nrow(Ktt))[matches[same],,drop=FALSE]; V[same,] <- 0; V[,same] <- 0 }
  list(A=A,covariance=(V+t(V))/2,same=same)
}
oracle_draw_extractor <- function(fit_record,prediction,basis=NULL) {
  f <- fit_record; validate_prediction_data(prediction,f); ni <- dim(f$native$draws)[2]; nc <- dim(f$native$draws)[3]
  train_center <- apply(f$native$draws,c(2,3),mean)
  ec <- mean(f$training_environment)*f$native$B
  function(idx) {
    stopifnot(length(idx)<=25L); p <- subset_prediction(prediction,idx); k <- oracle_conditional(f,p)
    fresh <- setdiff(seq_along(idx),k$same)
    if(length(fresh)) L <- t(chol(k$covariance[fresh,fresh,drop=FALSE]))
    dims <- c(length(idx),1L,ni,nc)
    out <- setNames(lapply(1:5,function(x) array(NA_real_,dims)),c("p_full","p_space_zero","p_env_zero","p_both_zero","field"))
    if(isTRUE(getOption("spatial.retain_raw_field"))) out$field_raw <- array(NA_real_,dims)
    for(ch in seq_len(nc)) for(it in seq_len(ni)) {
      field <- as.vector(k$A %*% f$native$draws[,it,ch])
      if(length(fresh)) field[fresh] <- field[fresh]+as.vector(L %*% rnorm(length(fresh)))
      terms <- center_terms(f$native$B0,p$environment*f$native$B,matrix(field,ncol=1),ec,train_center[it,ch])
      prob <- term_probabilities(terms,plogis)
      for(nm in names(prob)) out[[nm]][,,it,ch] <- prob[[nm]]
      out$field[,,it,ch] <- terms$spatial
      if(!is.null(out$field_raw)) out$field_raw[,,it,ch] <- field
    }
    out
  }
}
predict_oracle <- function(fit_record,prediction,basis=NULL) {
  ordering <- order(prediction$site_id,method="radix"); canonical <- subset_prediction(prediction,ordering)
  set.seed(fit_record$settings$prediction_seed)
  out <- summarize_draw_blocks(canonical,oracle_draw_extractor(fit_record,canonical),"logit",
    "Known nuisance parameters; posterior training field and true-GP conditional residual variance; independent conditional simulations across blocks")
  undo <- order(ordering); out$site_id <- prediction$site_id
  for(nm in c("p_full","p_space_zero","p_env_zero","p_both_zero","delta_p","field_mean","field_median")) out[[nm]] <- out[[nm]][undo,,drop=FALSE]
  out$intervals <- lapply(out$intervals,function(x) x[undo,,drop=FALSE])
  out$mcse <- lapply(out$mcse,function(x) x[undo,,drop=FALSE])
  if(!is.null(out$historical_field_median_centered)) out$historical_field_median_centered <- out$historical_field_median_centered[undo,,drop=FALSE]
  out$prediction_seed <- fit_record$settings$prediction_seed
  validate_prediction(out,prediction); out
}
