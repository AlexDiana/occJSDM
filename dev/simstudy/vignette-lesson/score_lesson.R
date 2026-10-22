# Summarize fits without exposing simulation truth to the fitter.
validate_lesson_fit_identity <- function(fit,bundle) {
  i <- fit$infos
  if(i$model=="binary") {
    z <- bundle$sim$true_params$z_true
    sites <- match(as.character(i$siteNames),rownames(z))
    species <- match(i$speciesNames,colnames(z))
    if(anyNA(c(sites,species))) stop("Unmatched fitted identity")
    expected <- z[sites,species,drop=FALSE]
    observed <- i$OTU
  } else {
    dat <- bundle$sim$data_list
    rows <- order(dat$info$Site,dat$info$Sample,dat$info$Primer)
    species <- match(i$speciesNames,colnames(dat$OTU))
    if(anyNA(species)) stop("Unmatched fitted species identity")
    # runOccJSDM adds SiteSample internally; compare all original columns.
    if(!all(names(dat$info) %in% names(i$data_info)) ||
       !isTRUE(all.equal(i$data_info[,names(dat$info),drop=FALSE],
                         dat$info[rows,,drop=FALSE],check.attributes=FALSE)))
      stop("Unmatched fitted observation identities")
    expected <- dat$OTU[rows,species,drop=FALSE]
    observed <- i$OTU
  }
  if(!isTRUE(all.equal(unname(observed),unname(expected),check.attributes=FALSE)))
    stop("Stored fitted observations do not match the declared experiment")
  invisible(TRUE)
}

lesson_probability_draws <- function(fit) {
  j <- fit$results_output$jsdm_output
  n <- nrow(fit$X_psi); S <- dim(j$B0_output)[1]
  ni <- dim(j$B0_output)[2]; nc <- dim(j$B0_output)[3]
  out <- array(NA_real_,c(n*S,ni,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    B <- matrix(j$B_output[,,it,ch],nrow=ncol(fit$X_psi),ncol=S)
    U <- matrix(j$U_output[,,it,ch],nrow=n)
    L <- matrix(j$L_output[,,it,ch],nrow=ncol(U),ncol=S)
    eta <- sweep(fit$X_psi %*% B + U %*% L,2,j$B0_output[,it,ch],"+")
    out[,it,ch] <- plogis(eta)
  }
  out
}

lesson_draw_summary <- function(a,diagnostics=FALSE) {
  dd <- dim(a); ni <- dd[length(dd)-1L]; nc <- tail(dd,1)
  a <- array(a,c(prod(head(dd,-2)),ni,nc))
  rows <- lapply(seq_len(dim(a)[1]),function(i) {
    x <- matrix(a[i,,],ni,nc)
    q <- quantile(x,c(.025,.975),names=FALSE)
    result <- data.frame(estimate=mean(x),lower=q[1],upper=q[2])
    if(diagnostics) {
      result$rhat <- posterior::rhat(x)
      result$ess <- posterior::ess_mean(x)
      result$mcse <- posterior::mcse_mean(x)
    }
    result
  })
  do.call(rbind,rows)
}

score_lesson_fit <- function(result,bundle) {
  f <- result$fit; arm <- result$arm
  validate_lesson_fit_identity(f,bundle)
  stopifnot(f$infos$ps==0)
  p <- lesson_probability_draws(f)
  cells <- expand.grid(Site=as.character(f$infos$siteNames),
                       species=f$infos$speciesNames,stringsAsFactors=FALSE)
  truth <- bundle$sim$true_params
  ix <- cbind(match(cells$Site,rownames(truth$z_true)),
              match(cells$species,colnames(truth$z_true)))
  stopifnot(!anyNA(ix))
  cells$truth <- plogis(truth$jsdmParams_true$eta[ix])
  cells$z <- truth$z_true[ix]
  cells <- cbind(arm=arm,cells,lesson_draw_summary(p,TRUE))
  saved <- f$results_output$psi_output
  reconstruction_difference <- if(is.null(saved)) NA_real_ else max(abs(cells$estimate-c(saved)))
  if(!is.null(saved)) stopifnot(reconstruction_difference<1e-10)
  cells$conditional <- if(arm=="perfect") NA_real_ else c(f$results_output$z_output)
  bands <- list(All=rep(TRUE,nrow(cells)),Low=cells$truth<.2,
                Middle=cells$truth>=.2 & cells$truth<=.8,High=cells$truth>.8)
  groups <- do.call(rbind,lapply(names(bands),function(band) {
    keep <- bands[[band]]; err <- cells$estimate[keep]-cells$truth[keep]
    tr <- apply(p[keep,,,drop=FALSE],c(2,3),mean)
    data.frame(arm=arm,band=band,cells=sum(keep),truth=mean(cells$truth[keep]),
               estimate=mean(cells$estimate[keep]),signed_error=mean(err),mae=mean(abs(err)),
               rhat=posterior::rhat(tr),ess=posterior::ess_mean(tr),mcse=posterior::mcse_mean(tr))
  }))
  rm(p); invisible(gc())
  rates <- samples <- NULL
  if(arm!="perfect") {
    sp <- f$infos$speciesNames
    pr <- as.character(f$infos$primerNames)
    ixsp <- match(sp,colnames(bundle$sim$data_list$OTU))
    ixpr <- match(pr,as.character(seq_len(bundle$settings$P)))
    stopifnot(!anyNA(c(ixsp,ixpr)))
    pt <- effective_detection_rate(bundle$params$p,bundle$params$mu1,bundle$params$sigma1)
    qt <- effective_detection_rate(bundle$params$q,bundle$params$mu0,bundle$params$sigma0)
    rates <- do.call(rbind,lapply(c("p","q","theta0"),function(param) {
      a <- f$results_output[[paste0(param,"_output")]]
      labels <- if(param=="theta0") data.frame(Primer="Field",species=sp) else
                  expand.grid(Primer=pr,species=sp,stringsAsFactors=FALSE)
      tv <- switch(param,p=c(pt[ixpr,ixsp,drop=FALSE]),q=c(qt[ixpr,ixsp,drop=FALSE]),
                    theta0=bundle$params$theta0[ixsp])
      cbind(arm=arm,param=param,labels,truth=tv,lesson_draw_summary(a,TRUE))
    }))
    si <- f$infos$data_info[!duplicated(f$infos$data_info[,c("Site","Sample")]),c("Site","Sample")]
    samples <- do.call(rbind,lapply(seq_along(sp),function(s) {
      wi <- match(as.character(si$Sample),rownames(truth$w_true))
      zs <- match(as.character(si$Site),rownames(truth$z_true))
      stopifnot(!anyNA(c(wi,zs)))
      data.frame(arm=arm,species=sp[s],si,z=truth$z_true[zs,ixsp[s]],
                 w=truth$w_true[wi,ixsp[s]],
                 sample_probability=f$results_output$w_output[,s],
                 collection_probability=f$results_output$theta_output[,s])
    }))
    key <- paste(cells$species,cells$Site)
    ix <- match(paste(samples$species,samples$Site),key)
    stopifnot(!anyNA(ix))
    samples$site_probability <- cells$conditional[ix]
    samples$occupancy_estimate <- cells$estimate[ix]
    samples$occupancy_truth <- cells$truth[ix]
  }
  diagnostics <- as.data.frame(occJSDM::returnConvergenceDiagnostics(f))
  list(cells=cells,groups=groups,rates=rates,samples=samples,
       diagnostics=diagnostics,warnings=result$warnings,
       reconstruction_difference=reconstruction_difference)
}
