occ_design <- function(data,basis) {
  spatial <- predict_basis(basis,data)
  out <- data.frame(environment=as.numeric(apply_transform(data$environment,basis$environment_transform)))
  if (basis$arm %in% c("trend","mem")) {
    colnames(spatial) <- paste0("sp_",colnames(spatial)); out <- cbind(out,spatial)
  }
  out
}
fit_occjsdm <- function(training,basis,settings) {
  validate_training(training); stopifnot(identical(training$site_id,basis$training_ids))
  started <- Sys.time(); warnings <- character(); set.seed(settings$seed); rng_start <- .Random.seed
  X <- occ_design(training,basis)
  input <- list(info=cbind(data.frame(Site=seq_len(nrow(X))),X),OTU=training$y,traits=NULL)
  if (basis$arm=="gp") input$info <- cbind(input$info,longitude=training$xy[,1],latitude=training$xy[,2])
  args <- list(data=input,listParams=list(n_factors=0L,n_lattrait=0L,n_supportpoints=as.integer(settings$support_points)),
    occCovariates=colnames(X),MCMCparams=as.list(settings$mcmc))
  if (basis$arm=="gp") args$spatCovariates <- c("longitude","latitude")
  native <- withCallingHandlers(do.call(occJSDM::runOccJSDM,args),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  })
  stopifnot(identical(native$infos$speciesNames,training$species_id),native$infos$n_factors==0L,
    nrow(native$X_psi)==nrow(training$y),native$infos$model=="binary",
    identical(native$infos$OTU,training$y))
  if(basis$arm=="gp") stopifnot(nrow(native$infos$list_Xs$X_tilde)==as.integer(settings$support_points))
  record_fit("occJSDM",training,basis,settings,native,started,warnings,
    list(rng_start=rng_start,coefficient_columns=colnames(native$X_psi),
         training_environment=training$environment,training_xy=training$xy,support_points=settings$support_points))
}
occ_gp_bases <- function(fit,xy=NULL) {
  xs <- fit$infos$list_Xs
  kernel <- function(a,b,l) exp(-xy_distances(a,b)^2/(2*l^2))
  if(is.null(xy)) locations <- xs$X_s[xs$Xs_index,,drop=FALSE] else {
    locations <- as.matrix(occJSDM:::transform_new_covariates(
      data.frame(longitude=xy[,1],latitude=xy[,2]),fit$infos$list_Xs_mat,remove_intercept=TRUE))
  }
  lapply(fit$infos$l_s_grid,function(l) {
    L <- t(chol(kernel(xs$X_tilde,xs$X_tilde,l)+diag(1e-5,nrow(xs$X_tilde))))
    t(solve(L,t(kernel(locations,xs$X_tilde,l))))
  })
}
occ_draw_extractor <- function(fit_record,prediction,basis) {
  validate_prediction_data(prediction,fit_record)
  fit <- fit_record$native; js <- fit$results_output$jsdm_output
  X <- as.matrix(occJSDM:::transform_new_covariates(occ_design(prediction,basis),fit$infos$list_X_psi_mat,remove_intercept=TRUE))
  envcol <- match("environment",colnames(X)); spcol <- setdiff(seq_len(ncol(X)),envcol)
  stopifnot(!is.na(envcol),identical(colnames(X),colnames(fit$X_psi)))
  S <- length(fit_record$species_id); ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  env_center <- mean(fit$X_psi[,envcol]); fixed_center <- colMeans(fit$X_psi)[spcol]
  if (basis$arm=="gp") {
    H <- occ_gp_bases(fit,prediction$xy)
    Hcenter <- lapply(occ_gp_bases(fit),colMeans)
  }
  function(idx) {
    if(length(idx)>25L) stop("Prediction draw block exceeds 25 sites")
    dims <- c(length(idx),S,ni,nc)
    result <- setNames(lapply(1:5,function(x) array(NA_real_,dims)),c("p_full","p_space_zero","p_env_zero","p_both_zero","field"))
    if(isTRUE(getOption("spatial.retain_raw_field"))) result$field_raw <- array(NA_real_,dims)
    for (ch in seq_len(nc)) for (it in seq_len(ni)) {
      B <- matrix(js$B_output[,,it,ch],ncol(X),S)
      env <- X[idx,envcol,drop=FALSE] %*% B[envcol,,drop=FALSE]
      sp <- matrix(0,length(idx),S); center <- rep(0,S)
      if(length(spcol)) {
        sp <- X[idx,spcol,drop=FALSE] %*% B[spcol,,drop=FALSE]
        center <- as.vector(fixed_center %*% B[spcol,,drop=FALSE])
      }
      if(basis$arm=="gp") {
        g <- js$idx_ls_output[it,ch]; bs <- matrix(js$Bs_output[,,it,ch],fit$infos$ps,S)
        sp <- H[[g]][idx,,drop=FALSE] %*% bs; center <- as.vector(Hcenter[[g]] %*% bs)
      }
      terms <- center_terms(js$B0_output[,it,ch],env,sp,as.vector(env_center*B[envcol,]),center)
      probs <- term_probabilities(terms,plogis)
      for(nm in names(probs)) result[[nm]][,,it,ch] <- probs[[nm]]
      result$field[,,it,ch] <- terms$spatial
      if(!is.null(result$field_raw)) result$field_raw[,,it,ch] <- sp
    }
    result
  }
}
predict_occjsdm <- function(fit_record,prediction,basis) {
  out <- summarize_draw_blocks(prediction,occ_draw_extractor(fit_record,prediction,basis),"logit",
    "Posterior coefficients and fitted SoR field; no added full-GP residual variance")
  validate_prediction(out,prediction); out
}
