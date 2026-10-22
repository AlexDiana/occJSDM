verify_saved_probability <- function(estimate,stored) {
  if(is.null(stored)) return(NA_real_)
  if(length(dim(stored))==4L) stored <- apply(stored,c(1,2),mean)
  stopifnot(identical(dim(estimate),dim(stored)),length(stored)>0L,
    all(is.finite(estimate)),all(is.finite(stored)))
  difference <- max(abs(estimate-stored))
  stopifnot(is.finite(difference),difference<1e-10)
  difference
}

trace_diagnostics <- function(x) {
  x <- as.matrix(x)
  constant <- any(apply(x,2,sd)==0)
  diagnostic <- function(f) if(constant || nrow(x)<4L) NA_real_ else
    suppressWarnings(as.numeric(f(x)))
  c(rhat=diagnostic(posterior::rhat),ess_bulk=diagnostic(posterior::ess_bulk),
    ess_mean=diagnostic(posterior::ess_mean),mcse=diagnostic(posterior::mcse_mean),
    chain_gap=diff(range(colMeans(x))),constant_chain=as.numeric(constant))
}

independent_bases <- function(fit) {
  xs <- fit$infos$list_Xs
  kernel <- function(a,b,l) exp(-(outer(a[,1],b[,1],'-')^2+
    outer(a[,2],b[,2],'-')^2)/(2*l^2))
  lapply(fit$infos$l_s_grid,function(l) {
    lower <- t(chol(kernel(xs$X_tilde,xs$X_tilde,l)+diag(1e-5,fit$infos$ps)))
    t(solve(lower,t(kernel(xs$X_s,xs$X_tilde,l))))[xs$Xs_index,,drop=FALSE]
  })
}

reconstruct_spatial_draws <- function(fit) {
  js <- fit$results_output$jsdm_output
  n <- nrow(fit$Xs); S <- dim(js$B0_output)[1]
  ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  H <- independent_bases(fit)
  probability <- array(NA_real_,c(n*S,ni,nc))
  field_chains <- array(0,c(n,S,nc))
  for(ch in seq_len(nc)) for(it in seq_len(ni)) {
    field <- H[[js$idx_ls_output[it,ch]]] %*% matrix(js$Bs_output[,,it,ch],fit$infos$ps,S)
    eta <- sweep(field+fit$X_psi %*% matrix(js$B_output[,,it,ch],ncol(fit$X_psi),S),
      2,js$B0_output[,it,ch],'+')
    stopifnot(fit$infos$n_factors==0L)
    probability[,it,ch] <- as.vector(plogis(eta))
    field_chains[,,ch] <- field_chains[,,ch]+field/ni
  }
  list(probability=probability,field_mean=apply(field_chains,c(1,2),mean),
    field_chains=field_chains)
}

spatial_groups <- function(psi,target_prevalence) {
  g <- list(all=matrix(TRUE,nrow(psi),ncol(psi)),low=psi<.2,
    medium=psi>=.2 & psi<=.8,high=psi>.8)
  for(p in c(.01,.05,.25,.75)) g[[paste0('prevalence_',p*100,'pct')]] <-
    matrix(rep(abs(target_prevalence-p)<1e-10,each=nrow(psi)),nrow(psi),ncol(psi))
  g
}

score_draw_block <- function(draws,truth,metric,group='all') {
  ni <- dim(draws)[length(dim(draws))-1L]; nc <- tail(dim(draws),1)
  truth <- as.vector(truth); x <- array(draws,c(length(truth),ni,nc))
  means <- rowMeans(matrix(x,nrow=length(truth)))
  quantiles <- t(vapply(seq_along(truth),function(i)
    unname(quantile(x[i,,],c(.025,.975))),numeric(2)))
  diagnostics <- t(vapply(seq_along(truth),function(i)
    trace_diagnostics(matrix(x[i,,],ni,nc)),numeric(6)))
  covered <- quantiles[,1]<=truth & truth<=quantiles[,2]
  elements <- data.frame(metric=metric,element=seq_along(truth),truth=truth,
    estimate=means,bias=means-truth,lower=quantiles[,1],upper=quantiles[,2],
    covered=covered,diagnostics)
  trace <- apply(x,c(2,3),mean)
  summary <- data.frame(metric=metric,group=group,n=length(truth),truth=mean(truth),
    estimate=mean(means),bias=mean(means-truth),mae=mean(abs(means-truth)),
    rmse=sqrt(mean((means-truth)^2)),coverage=mean(covered),
    interval_width=mean(quantiles[,2]-quantiles[,1]),
    as.list(trace_diagnostics(trace)))
  list(elements=elements,summary=summary,trace=trace)
}

score_spatial_fit <- function(fit,input,arm,knots) {
  t <- input$truth; ro <- fit$results_output; js <- ro$jsdm_output
  n <- input$settings$n; S <- input$settings$S
  ni <- dim(js$B0_output)[2]; nc <- dim(js$B0_output)[3]
  stopifnot(fit$infos$ps==knots,fit$infos$n_factors==0L,
    identical(fit$infos$speciesNames,colnames(input$data[[arm]]$OTU)),
    fit$infos$model==if(arm=='binary')'binary' else 'two_stage',
    max(abs(fit$Xs-t$Xs))<1e-12,max(abs(fit$X_psi-t$X))<1e-12)
  if(arm!='binary') stopifnot(max(abs(fit$X_theta-t$Xt))<1e-12)
  reconstructed <- reconstruct_spatial_draws(fit)
  x <- reconstructed$probability
  estimate <- matrix(rowMeans(matrix(x,n*S)),n,S)
  saved_probability_difference <- verify_saved_probability(estimate,ro$psi_output)
  # Binary fits do not store psi_output. Verify every independently built
  # spatial basis against the package's own kernel/whitening/C++ projection.
  H <- independent_bases(fit)
  native <- suppressMessages(occJSDM:::precomputeSORmatrices(fit$infos$l_s_grid,fit$infos$list_Xs))
  basis_difference <- max(vapply(seq_along(H),function(g) {
    reference <- occJSDM:::KsBproduct(native$Ks_all[,,g],diag(knots),fit$infos$list_Xs$Xs_centers)
    stopifnot(identical(dim(reference),dim(H[[g]])),all(is.finite(reference)))
    max(abs(H[[g]]-reference))
  },numeric(1)))
  stopifnot(is.finite(basis_difference),basis_difference<1e-10)
  draw_difference <- 0
  for(ch in seq_len(nc)) for(it in unique(c(1L,ceiling(ni/2),ni))) {
    g <- js$idx_ls_output[it,ch]
    reference <- occJSDM:::computePsiCoef(fit$X_psi,native$Ks_all[,,g],
      fit$infos$list_Xs$Xs_centers,Tr=matrix(0,S,0),B0=js$B0_output[,it,ch],
      G=matrix(0,0,ncol(fit$X_psi)),A=matrix(0,S,0),C=matrix(0,0,ncol(fit$X_psi)),
      Bt=t(matrix(js$B_output[,,it,ch],ncol(fit$X_psi),S)),
      Gs=matrix(0,0,knots),As=matrix(0,S,0),Cs=matrix(0,0,knots),
      Bst=t(matrix(js$Bs_output[,,it,ch],knots,S)),U=matrix(0,n,0),L=matrix(0,0,S))$eta
    stopifnot(identical(dim(reference),c(n,S)),all(is.finite(reference)))
    draw_difference <- max(draw_difference,max(abs(plogis(reference)-matrix(x[,it,ch],n,S))))
  }
  stopifnot(is.finite(draw_difference),draw_difference<1e-10)
  difference <- max(c(basis_difference,draw_difference,saved_probability_difference),na.rm=TRUE)
  main <- score_draw_block(x,t$psi,'occupancy')
  elements <- main$elements; groups <- list(main$summary)
  traces <- list(occupancy_all=main$trace)
  masks <- spatial_groups(t$psi,t$target_prevalence)
  for(g in setdiff(names(masks),'all')) {
    idx <- which(masks[[g]])
    if(!length(idx)) next
    # Reuse per-element summaries, but compute group-trace diagnostics directly.
    el <- elements[idx,]; trace <- apply(x[idx,,,drop=FALSE],c(2,3),mean)
    groups[[length(groups)+1L]] <- data.frame(metric='occupancy',group=g,n=length(idx),
      truth=mean(el$truth),estimate=mean(el$estimate),bias=mean(el$bias),
      mae=mean(abs(el$bias)),rmse=sqrt(mean(el$bias^2)),coverage=mean(el$covered),
      interval_width=mean(el$upper-el$lower),as.list(trace_diagnostics(trace)))
    traces[[paste0('occupancy_',g)]] <- trace
  }
  blocks <- list(intercept=list(js$B0_output,t$B0),environment_slope=list(js$B_output,t$B),
    range=list(matrix(fit$infos$l_s_grid[js$idx_ls_output],ni,nc),t$range),
    spatial_sd=list(js$sigmabs_output,t$field_sd))
  if(arm!='binary') blocks <- c(blocks,list(
    collection_coefficient=list(ro$beta_theta_output,t$beta_theta),
    theta0=list(ro$theta0_output,t$theta0),p=list(ro$p_output,t$p_positive),
    q=list(ro$q_output,t$q_positive[[arm]])))
  for(nm in names(blocks)) {
    b <- score_draw_block(blocks[[nm]][[1]],blocks[[nm]][[2]],nm)
    elements <- rbind(elements,b$elements);groups[[length(groups)+1L]] <- b$summary
    traces[[nm]] <- b$trace
  }
  species <- do.call(rbind,lapply(seq_len(S),function(s) {
    idx <- ((s-1L)*n+1L):(s*n); el <- elements[elements$metric=='occupancy',][idx,]
    data.frame(species=colnames(input$data[[arm]]$OTU)[s],target=t$target_prevalence[s],
      occupied=sum(t$z[,s]),detections=sum(input$data[[arm]]$OTU[,s]),
      truth=mean(el$truth),estimate=mean(el$estimate),bias=mean(el$bias),
      relative_bias=mean(el$bias)/mean(el$truth),mae=mean(abs(el$bias)),
      rmse=sqrt(mean(el$bias^2)),coverage=mean(el$covered),
      as.list(trace_diagnostics(apply(x[idx,,,drop=FALSE],c(2,3),mean))))
  }))
  field <- reconstructed$field_mean
  Htrue <- independent_bases(fit)[[t$grid_index]]
  projection <- qr.fitted(qr(Htrue),t$field)
  field_score <- data.frame(bias=mean(field-t$field),mae=mean(abs(field-t$field)),
    rmse=sqrt(mean((field-t$field)^2)),
    slope=sum((field-mean(field))*(t$field-mean(t$field)))/sum((t$field-mean(t$field))^2),
    projection_rmse=sqrt(mean((projection-t$field)^2)))
  range_rows <- do.call(rbind,lapply(seq_len(nc),function(ch)
    data.frame(chain=ch,range=fit$infos$l_s_grid,
      frequency=tabulate(js$idx_ls_output[,ch],nbins=length(fit$infos$l_s_grid))/ni)))
  list(groups=do.call(rbind,groups),elements=elements,species=species,
    traces=traces,probability=estimate,field=field_score,range=range_rows,
    chain_probability=apply(x,c(1,3),mean),reconstruction=difference,
    verification=list(basis=basis_difference,selected_draws=draw_difference,
      saved_probability=saved_probability_difference,
      method=if(is.na(saved_probability_difference)) 'Independent bases and selected probabilities versus native calculations' else
        'Independent bases and selected probabilities versus native calculations; saved probability means'))
}
