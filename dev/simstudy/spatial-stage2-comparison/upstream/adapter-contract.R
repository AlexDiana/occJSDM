validate_training <- function(x) {
  allowed <- c("schema_version","dataset_id","site_id","species_id","y","environment","xy","input_hashes")
  if (!setequal(names(x),allowed)) stop("Unexpected or missing TrainingData fields; truth is forbidden")
  if (!identical(x$schema_version,1L) || anyDuplicated(x$site_id) || anyDuplicated(x$species_id) || anyNA(c(x$site_id,x$species_id))) stop("Invalid schema or IDs")
  stopifnot(is.matrix(x$y),identical(dim(x$y),c(length(x$site_id),length(x$species_id))),
    identical(colnames(x$y),x$species_id),all(x$y %in% 0:1),
    identical(dim(x$environment),c(nrow(x$y),1L)),identical(dim(x$xy),c(nrow(x$y),2L)),
    all(is.finite(x$environment)),all(is.finite(x$xy)),length(x$input_hashes)>0L)
  invisible(TRUE)
}
validate_prediction_data <- function(x,record) {
  if (any(c("y","z","psi","truth","field","landscape","B","B0") %in% names(x))) stop("PredictionData must not contain responses or truth")
  if (!identical(x$species_id,record$species_id)) stop("Prediction species order mismatch")
  if (length(x$site_id)!=nrow(x$xy) || anyDuplicated(x$site_id)) stop("Prediction site IDs mismatch")
  invisible(TRUE)
}
validate_prediction <- function(record,expected_ids) {
  if (!identical(record$site_id,expected_ids$site_id) || !identical(record$species_id,expected_ids$species_id)) stop("Prediction ID order mismatch")
  for (nm in c("p_full","p_space_zero","p_env_zero","p_both_zero")) {
    x <- record[[nm]]
    if (!identical(dim(x),c(length(record$site_id),length(record$species_id))) || any(!is.finite(x)) || any(x<0|x>1)) stop("Invalid probability matrix: ",nm)
  }
  if (max(abs(record$delta_p-(record$p_full-record$p_space_zero)))>1e-12) stop("Probability contrast mismatch")
  if (is.null(record$uncertainty$status)) stop("Uncertainty capability is required")
  if (record$uncertainty$status=="unavailable" && !is.null(record$intervals)) stop("Unavailable intervals cannot be fabricated")
  invisible(TRUE)
}
center_terms <- function(intercept,environment,spatial,environment_center,spatial_center) {
  n <- nrow(environment); S <- ncol(environment)
  stopifnot(identical(dim(environment),dim(spatial)),length(intercept)==S,
            length(environment_center)==S,length(spatial_center)==S)
  list(intercept=matrix(rep(intercept+environment_center+spatial_center,each=n),n,S),
       environment=sweep(environment,2,environment_center,"-"),
       spatial=sweep(spatial,2,spatial_center,"-"))
}
term_probabilities <- function(terms,inverse_link) {
  i <- terms$intercept; e <- terms$environment; s <- terms$spatial
  list(p_full=inverse_link(i+e+s),p_space_zero=inverse_link(i+e),
       p_env_zero=inverse_link(i+s),p_both_zero=inverse_link(i))
}
point_prediction <- function(ids,terms,inverse_link,link,conditioning) {
  p <- term_probabilities(terms,inverse_link)
  c(list(schema_version=1L,site_id=ids$site_id,species_id=ids$species_id,
    target="realised-landscape occurrence probability",link=link,conditioning=conditioning),p,
    list(delta_p=p$p_full-p$p_space_zero,field_mean=terms$spatial,
      field_median=terms$spatial,intervals=NULL,
      uncertainty=list(status="unavailable",reason="Point fit; no certified probability interval")))
}
normal_quadrature <- local({
  cache <- new.env()
  function(n=60L) {
    key <- as.character(n); if (exists(key,cache,inherits=FALSE)) return(cache[[key]])
    J <- matrix(0,n,n); J[cbind(1:(n-1),2:n)] <- sqrt(1:(n-1)); J <- J+t(J)
    e <- eigen(J,symmetric=TRUE); o <- order(e$values)
    cache[[key]] <- list(nodes=e$values[o],weights=e$vectors[1,o]^2); cache[[key]]
  }
})
normal_integral <- function(mu,sd=1,inverse_link=plogis,n=60L) {
  q <- normal_quadrature(n); ans <- mu*0
  for (k in seq_along(q$nodes)) ans <- ans+q$weights[k]*inverse_link(mu+sd*q$nodes[k])
  ans
}
record_fit <- function(package,training,basis,settings,native,started,warnings=character(),extra=list()) {
  for (nm in c("seed","job_id","attempt_id","source_hashes")) if (is.null(settings[[nm]])) stop("Missing fit identity: ",nm)
  c(list(schema_version=1L,package=package,arm=basis$arm,dataset_id=training$dataset_id,
    site_id=training$site_id,species_id=training$species_id,input_hashes=training$input_hashes,
    settings=settings,source_hashes=settings$source_hashes,native=native,warnings=unique(warnings),
    started=started,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
    native_status="returned",runtime=list(R=R.version.string,package_version=as.character(utils::packageVersion(package)),
    package_path=find.package(package)),rng_state=.Random.seed),extra)
}
# Each raw draw block is <=25 sites. Native posterior parameters are retained losslessly.
summarize_draw_blocks <- function(prediction,extractor,link,conditioning) {
  n <- length(prediction$site_id); S <- length(prediction$species_id)
  fields <- c("p_full","p_space_zero","p_env_zero","p_both_zero","field_mean","field_median")
  out <- setNames(lapply(fields,function(x) matrix(NA_real_,n,S)),fields)
  intervals <- setNames(lapply(c("probability_lower","probability_upper","field_lower","field_upper"),function(x) matrix(NA_real_,n,S)),c("probability_lower","probability_upper","field_lower","field_upper"))
  mcse <- setNames(lapply(c(fields[1:4],"delta_p","field"),function(x) matrix(NA_real_,n,S)),c(fields[1:4],"delta_p","field"))
  for (idx in split(seq_len(n),ceiling(seq_len(n)/25L))) {
    rng_before <- if(exists(".Random.seed",globalenv(),inherits=FALSE)) get(".Random.seed",globalenv()) else NULL
    draws <- extractor(idx)
    attr(draws,"rng_before") <- rng_before
    sink <- getOption("spatial.draw_sink")
    if(is.function(sink)) sink(prediction$site_id[idx],draws)
    ni <- dim(draws$p_full)[3]; nc <- dim(draws$p_full)[4]
    for(nm in names(mcse)) {
      values <- if(nm=="delta_p") draws$p_full-draws$p_space_zero else draws[[nm]]
      mcse[[nm]][idx,] <- apply(values,c(1,2),function(x) {
        if(all(x==x[1])) 0 else suppressWarnings(posterior::mcse_mean(matrix(x,ni,nc)))
      })
    }
    for (nm in fields[1:4]) out[[nm]][idx,] <- apply(draws[[nm]],c(1,2),mean)
    out$field_mean[idx,] <- apply(draws$field,c(1,2),mean)
    out$field_median[idx,] <- apply(draws$field,c(1,2),median)
    if(!is.null(draws$field_raw)) {
      if(is.null(out$historical_field_median_centered)) out$historical_field_median_centered <- matrix(NA_real_,n,S)
      out$historical_field_median_centered[idx,] <- apply(draws$field_raw,c(1,2),median)
    }
    for (pair in list(c("p_full","probability"),c("field","field"))) {
      intervals[[paste0(pair[2],"_lower")]][idx,] <- apply(draws[[pair[1]]],c(1,2),quantile,probs=.025)
      intervals[[paste0(pair[2],"_upper")]][idx,] <- apply(draws[[pair[1]]],c(1,2),quantile,probs=.975)
    }
  }
  if(!is.null(out$historical_field_median_centered)) out$historical_field_median_centered <- sweep(out$historical_field_median_centered,2,colMeans(out$historical_field_median_centered),"-")
  c(list(schema_version=1L,site_id=prediction$site_id,species_id=prediction$species_id,
    target="realised-landscape occurrence probability",link=link,conditioning=conditioning),out,
    list(delta_p=out$p_full-out$p_space_zero,intervals=intervals,mcse=mcse,
    uncertainty=list(status="available",reason="Unthinned posterior draws; Bernoulli state uncertainty excluded")))
}
subset_prediction <- function(prediction,idx) {
  out <- prediction; out$site_id <- prediction$site_id[idx]
  out$xy <- prediction$xy[idx,,drop=FALSE]; out$environment <- prediction$environment[idx,,drop=FALSE]; out
}
certified_normal_integral <- function(mu,sd,inverse_link=plogis) {
  low <- normal_integral(mu,sd,inverse_link,n=60L); high <- normal_integral(mu,sd,inverse_link,n=100L)
  error <- abs(low-high); bad <- which(!is.finite(error)|error>1e-6)
  adaptive_error <- numeric(length(bad))
  for(k in seq_along(bad)) {
    i <- bad[k]; m <- mu[i]; s <- sd[i]
    integrand <- function(z) inverse_link(m+s*z)*dnorm(z)
    a <- integrate(integrand,-Inf,Inf,rel.tol=1e-9,abs.tol=1e-10,subdivisions=200L)
    b <- integrate(integrand,-Inf,Inf,rel.tol=1e-11,abs.tol=1e-12,subdivisions=400L)
    adaptive_error[k] <- max(abs(a$value-b$value),b$abs.error)
    if(!is.finite(adaptive_error[k])||adaptive_error[k]>1e-6) stop("Adaptive normal quadrature failed certification")
    high[i] <- b$value
  }
  list(value=high,adaptive_cells=length(bad),max_refinement_error=max(c(if(length(bad)) error[-bad] else error,adaptive_error,0)),
    method="60/100-node Gaussian quadrature with independently refined adaptive integration for discrepant cells")
}
