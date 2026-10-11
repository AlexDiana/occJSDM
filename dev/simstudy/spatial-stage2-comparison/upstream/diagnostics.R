trace_diagnostic <- function(values,label,structural=FALSE,config=study_config()) {
  x <- as.matrix(values)
  constant <- all(is.finite(x)) && length(unique(as.vector(x)))==1L
  if(structural && !constant) stop("Declared structural constant varies: ",label)
  stats <- if(structural) c(rhat=NA_real_,ess_bulk=NA_real_,ess_tail=NA_real_) else
    suppressWarnings(c(rhat=posterior::rhat(x),ess_bulk=posterior::ess_bulk(x),ess_tail=posterior::ess_tail(x)))
  failed <- !structural && (constant || any(!is.finite(stats)) || stats[1]>config$diagnostics$rhat || any(stats[2:3]<config$diagnostics$ess))
  data.frame(quantity=label,rhat=stats[1],ess_bulk=stats[2],ess_tail=stats[3],structural=structural,
    failed=failed,reason=if(structural) "structural constant" else if(constant) "degenerate sampled quantity" else if(failed) "convergence threshold" else "passed",row.names=NULL)
}
diagnose_fit <- function(fit_record,training,basis,config=study_config()) {
  f <- fit_record; rows <- list()
  add <- function(x,label,structural=FALSE) rows[[length(rows)+1L]] <<- trace_diagnostic(x,label,structural,config)
  if(f$package %in% c("sjSDM","gllvm")) {
    if(f$package=="sjSDM") {
      objective <- f$native$objective$objective; gradient <- f$native$objective$max_gradient
      native_ok <- isTRUE(f$native$diagonal_verified); integration <- f$native$integration_error
    } else {
      objective <- f$native$objective; gradient <- max(abs(f$native$gradient)); native_ok <- isTRUE(f$native$convergence)
      integration <- 0
    }
    ok <- native_ok && is.finite(objective) && is.finite(gradient) && gradient<config$diagnostics$gradient && integration<=1e-6
    return(data.frame(quantity="objective/stationarity",objective=objective,stationarity=gradient,native_ok=native_ok,
      integration_error=integration,failed=!ok,reason=if(ok) "passed" else "native convergence, finite objective or stationarity failed"))
  }
  pred <- training; pred$y <- NULL
  extractor <- switch(f$package,occJSDM=occ_draw_extractor(f,pred,basis),Hmsc=hmsc_draw_extractor(f,pred,basis),oracle=oracle_draw_extractor(f,pred,basis))
  sums <- squares <- probability <- NULL
  for(idx in split(seq_along(training$site_id),ceiling(seq_along(training$site_id)/25L))) {
    d <- extractor(idx); dims <- dim(d$field); ni <- dims[3]; nc <- dims[4]
    for(j in seq_along(f$species_id)) for(k in seq_along(idx))
      add(matrix(d$field[k,j,,],ni,nc),paste0("field/",f$species_id[j],"/",training$site_id[idx[k]]),f$arm=="env")
    fs <- apply(d$field,c(2,3,4),sum); fq <- apply(d$field^2,c(2,3,4),sum); ps <- apply(d$p_full,c(2,3,4),sum)
    if(is.null(sums)) { sums <- fs; squares <- fq; probability <- ps } else { sums <- sums+fs; squares <- squares+fq; probability <- probability+ps }
  }
  n <- length(training$site_id); S <- length(f$species_id)
  for(j in seq_len(S)) {
    add(matrix(probability[j,,]/n,ni,nc),paste0("probability_mean/",f$species_id[j]))
    sd <- sqrt(pmax((squares[j,,]-sums[j,,]^2/n)/(n-1),0))
    add(matrix(sd,ni,nc),paste0("realised_field_sd/",f$species_id[j]),f$arm=="env")
  }
  add_array <- function(x,label) {
    dims <- dim(x); if(is.null(x)||length(x)==0) return(invisible(NULL))
    k <- length(x)/(ni*nc); x <- array(x,c(k,ni,nc))
    for(i in seq_len(k)) add(matrix(x[i,,],ni,nc),paste0(label,"/",i))
  }
  if(f$package=="occJSDM") {
    js <- f$native$results_output$jsdm_output
    add_array(js$B0_output,"intercept"); add_array(js$B_output,"coefficient")
    add_array(js$sigmab_output,"coefficient_prior_sd")
    if(f$arm=="gp") {
      add_array(js$sigmabs_output,"spatial_prior_sd")
      range <- array(f$native$infos$l_s_grid[js$idx_ls_output],dim(js$idx_ls_output))
      add_array(range,"range")
    }
  } else if(f$package=="Hmsc") {
    posts <- f$native$postList; beta <- unlist(lapply(posts,function(ch) lapply(ch,`[[`,"Beta")),use.names=FALSE)
    add_array(array(beta,c(length(posts[[1]][[1]]$Beta),ni,nc)),"coefficient")
    if(f$arm=="spatial2") {
      # Species covariance at a fixed raw-coordinate distance is invariant to factor signs/permutations.
      corhalf <- array(NA_real_,c(S,ni,nc))
      for(ch in seq_len(nc)) for(it in seq_len(ni)) {
        p <- posts[[ch]][[it]]; a <- f$native$ranLevels$site$alphapw[p$Alpha[[1]],1]
        rho <- ifelse(a==0,0,exp(-.5/a)); L2 <- p$Lambda[[1]]^2
        corhalf[,it,ch] <- colSums(L2*rho)/colSums(L2)
      }
      add_array(corhalf,"spatial_correlation_at_distance_0.5")
    }
  }
  out <- do.call(rbind,rows)
  convergence_warnings <- f$warnings[grepl("converg|mixing|divergen|positive.definite|singular|NaN|non.finite",f$warnings,ignore.case=TRUE)]
  if(length(convergence_warnings)) out <- rbind(out,data.frame(quantity="native_convergence_warnings",rhat=NA_real_,ess_bulk=NA_real_,ess_tail=NA_real_,structural=FALSE,failed=TRUE,reason=paste(convergence_warnings,collapse="; ")))
  out
}
