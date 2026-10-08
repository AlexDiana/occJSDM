.sjsdm_code_dir <- local({
  paths <- Filter(function(x) is.character(x) && length(x)==1L, lapply(sys.frames(),function(e) e$ofile))
  if (!length(paths)) stop("Source adapter-sjsdm.R with source() so its support module can be located")
  dirname(normalizePath(tail(paths,1L)[[1L]]))
})
sjsdm_marginal <- function(mu,n=60L) normal_integral(mu,1,
  function(x) .999999*plogis(1.70169*x)+.0000005,n)
network_forward <- function(x,layers) {
  x <- as.matrix(x)
  for (layer in layers) {
    if (layer$kind=="linear") {
      x <- x %*% t(layer$weight)
      if(!is.null(layer$bias)) x <- sweep(x,2,as.vector(layer$bias),"+")
    } else if(layer$kind=="selu") {
      x <- 1.0507009873554805*ifelse(x>0,x,1.6732632423543772*expm1(pmin(x,0)))
    } else stop("Uncertified saved network layer")
  }
  unname(x)
}
fit_sjsdm <- function(training,basis,settings) {
  validate_training(training); stopifnot(identical(training$site_id,basis$training_ids))
  started <- Sys.time(); warnings <- character(); set.seed(settings$seed); rng_start <- .Random.seed
  library(sjSDM)
  support <- reticulate::import_from_path("sjsdm_support",path=.sjsdm_code_dir,convert=TRUE)
  control <- settings$sjsdm
  env <- as.data.frame(apply_transform(training$environment,basis$environment_transform)); names(env) <- "environment"
  sp <- predict_basis(basis,training); spatial <- NULL
  if (basis$arm %in% c("trend","mem")) spatial <- sjSDM::linear(as.data.frame(sp),~0+.,lambda=control$lambda,alpha=control$alpha)
  if (basis$arm=="dnn") spatial <- sjSDM::DNN(as.data.frame(sp),~0+.,hidden=control$hidden,
    activation=control$activation,bias=control$bias,dropout=control$dropout,lambda=control$lambda,alpha=control$alpha)
  fit <- withCallingHandlers(sjSDM::sjSDM(Y=training$y,env=sjSDM::linear(env,lambda=0),
    spatial=spatial,biotic=sjSDM::bioticStruct(diag=TRUE,lambda=0),family=binomial("probit"),
    device="cpu",dtype="float64",iter=control$epochs[1],sampling=control$sampling,
    step_size=nrow(training$y),learning_rate=control$lr[1],parallel=0L,
    control=sjSDM::sjSDMControl(optimizer=sjSDM::RMSprop(weight_decay=control$weight_decay),
      scheduler=0,early_stopping_training=0),seed=settings$seed,se=FALSE,verbose=FALSE),
    warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
  stopifnot(identical(fit$data$Y,training$y),length(fit$history)==control$epochs[1])
  model <- fit$model; support$check_diagonal(model)
  initial <- support$snapshot(model)
  l1 <- if(is.null(spatial)) 0 else spatial$l1_coef
  l2 <- if(is.null(spatial)) 0 else spatial$l2_coef
  SP <- if(is.null(spatial)) NULL else spatial$X
  before <- support$exact_objective(model,fit$data$X,training$y,SP,l1,l2,control$weight_decay)
  history <- fit$history
  support$restart_optimizer(model,control$lr[2],control$weight_decay)
  model$fit(fit$data$X,training$y,SP=SP,batch_size=as.integer(nrow(training$y)),
    epochs=as.integer(control$epochs[2]),sampling=as.integer(control$sampling),parallel=0L,
    early_stopping_training=-1L,verbose=FALSE)
  history <- c(history,as.vector(sjSDM:::force_r(model$history)))
  saved <- support$snapshot(model)
  check60 <- support$exact_objective(model,fit$data$X,training$y,SP,l1,l2,control$weight_decay,60L)
  check100 <- support$exact_objective(model,fit$data$X,training$y,SP,l1,l2,control$weight_decay,100L)
  integration <- max(abs(check60$probabilities-check100$probabilities))
  if(integration>1e-6 || abs(check60$objective-check100$objective)>1e-6) stop("sjSDM quadrature refinement failed")
  native_raw <- sjSDM:::force_r(model$predict(fit$data$X,SP=SP,link=FALSE,batch_size=as.integer(nrow(training$y))))
  raw <- network_forward(fit$data$X,saved$environment)
  if(length(saved$spatial)) raw <- raw+network_forward(SP,saved$spatial)
  if(max(abs(raw*.999999+.0000005-native_raw))>1e-6) stop("Saved network does not reproduce native affine raw predictor")
  native_probability <- sjSDM:::force_r(model$predict(fit$data$X,SP=SP,link=TRUE,batch_size=as.integer(nrow(training$y))))
  if(max(abs(native_probability-(pnorm(raw)*.999999+.0000005)))>1e-6) stop("Ordinary native probit prediction mismatch")
  backend_rng <- support$rng_state()
  calibration <- support$certify_training_marginal(model,as.integer(settings$seed))
  calibration$reference <- sjsdm_marginal(calibration$mu)
  calibration$passed <- abs(calibration$mean-calibration$reference)<=5*calibration$se+1e-6
  if(!all(calibration$passed)) stop("Native Monte Carlo marginal check failed")
  native <- list(network=saved,first_phase=initial,history=history,l1=l1,l2=l2,backend_rng=backend_rng,
    objective=check100,before=before,integration_error=integration,calibration=calibration,
    ordinary_native_probability=native_probability,
    raw_prediction_convention="Native raw prediction also applies 0.999999*mu+0.0000005",
    native_raw=native_raw,diagonal_verified=TRUE)
  record_fit("sjSDM",training,basis,settings,native,started,warnings,
    list(rng_start=rng_start,training_environment=training$environment,training_xy=training$xy,
         python=reticulate::py_config()$python))
}
predict_sjsdm <- function(fit_record,prediction,basis) {
  validate_prediction_data(prediction,fit_record)
  net <- fit_record$native$network
  env <- apply_transform(prediction$environment,basis$environment_transform)
  train_env <- apply_transform(fit_record$training_environment,basis$environment_transform)
  intercept <- as.vector(network_forward(matrix(c(1,0),1),net$environment))
  fullenv <- network_forward(cbind(1,env),net$environment)
  E <- sweep(fullenv,2,intercept,"-")
  Ec <- colMeans(network_forward(cbind(1,train_env),net$environment))-intercept
  S <- matrix(0,nrow(env),length(fit_record$species_id)); Sc <- rep(0,ncol(S))
  if(length(net$spatial)) {
    S <- network_forward(predict_basis(basis,prediction),net$spatial)
    Sc <- colMeans(network_forward(basis$training_matrix,net$spatial))
  }
  terms <- center_terms(intercept,E,S,Ec,Sc)
  ans <- point_prediction(prediction,terms,sjsdm_marginal,"probit-labelled bounded logistic-normal",
    "Independent standard-normal noise integrated under the pinned training likelihood; fitted weights fixed")
  ans$ordinary_native_probability <- pnorm(fullenv+S)*.999999+.0000005
  validate_prediction(ans,prediction); ans
}
