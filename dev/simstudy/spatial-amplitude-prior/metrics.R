# Reconstruct in blocks to avoid retaining a second full field-draw array.
summarise_field_draws <- function(coefficients,index,bases,truth,amplitude) {
  d <- dim(coefficients); stopifnot(length(d)==4L,length(dim(truth))==2L)
  ps <- d[1]; S <- d[2]; ni <- d[3]; nc <- d[4]; n <- nrow(truth)
  stopifnot(ncol(truth)==S,identical(dim(index),c(ni,nc)),
    identical(dim(amplitude),c(ni,nc)),all(is.finite(coefficients)),
    all(index %in% seq_along(bases)),all(vapply(bases,function(b)
      identical(dim(b),c(n,ps)),logical(1))))
  tc <- sweep(truth,2,colMeans(truth),'-'); denom <- colSums(tc^2)
  stopifnot(all(denom>0))
  labels <- c('amplitude',paste0('field_mean_',seq_len(S)),
    paste0('field_rms_',seq_len(S)),paste0('field_projection_',seq_len(S)))
  traces <- array(NA_real_,c(length(labels),ni,nc),dimnames=list(labels,NULL,NULL))
  traces[1,,] <- amplitude
  total <- matrix(0,n,S)
  for(ch in seq_len(nc)) for(g in sort(unique(index[,ch]))) {
    ii <- which(index[,ch]==g)
    for(first in seq.int(1L,length(ii),by=100L)) {
      its <- ii[first:min(length(ii),first+99L)]; k <- length(its)
      field <- bases[[g]] %*% matrix(coefficients[,,its,ch,drop=FALSE],ps,S*k)
      stopifnot(all(is.finite(field)))
      means <- matrix(colMeans(field),S,k)
      rms <- matrix(sqrt(colMeans(field^2)),S,k)
      projection <- matrix(colSums(field*tc[,rep(seq_len(S),k),drop=FALSE]),S,k)/denom
      traces[1L+seq_len(S),its,ch] <- means
      traces[1L+S+seq_len(S),its,ch] <- rms
      traces[1L+2L*S+seq_len(S),its,ch] <- projection
      for(s in seq_len(S)) total[,s] <- total[,s]+rowSums(field[,seq.int(s,S*k,by=S),drop=FALSE])
    }
  }
  field <- total/(ni*nc); fc <- sweep(field,2,colMeans(field),'-')
  score <- data.frame(raw_rmse=sqrt(mean((field-truth)^2)),raw_mae=mean(abs(field-truth)),
    centred_rmse=sqrt(mean((fc-tc)^2)),centred_mae=mean(abs(fc-tc)),
    centred_correlation=cor(as.vector(fc),as.vector(tc)),
    centred_slope=sum(fc*tc)/sum(tc^2),posterior_mean_rms=sqrt(mean(field^2)),
    true_rms=sqrt(mean(truth^2)))
  diagnostics <- do.call(rbind,lapply(seq_along(labels),function(i)
    data.frame(quantity=labels[i],as.list(trace_diagnostics(matrix(traces[i,,],ni,nc))))))
  list(field_mean=field,score=score,traces=traces,diagnostics=diagnostics)
}

spatial_flags <- function(d) {
  reasons <- character()
  for(i in seq_len(nrow(d))) {
    if(!is.finite(d$rhat[i]) || !is.finite(d$ess_mean[i]))
      reasons <- c(reasons,paste(d$quantity[i],'diagnostics unavailable'))
    if(is.finite(d$rhat[i]) && d$rhat[i]>1.05)
      reasons <- c(reasons,paste(d$quantity[i],'Rhat > 1.05'))
    if(is.finite(d$ess_mean[i]) && d$ess_mean[i]<100)
      reasons <- c(reasons,paste(d$quantity[i],'ESS < 100'))
  }
  reasons
}

with_spatial_starts <- function(fun,starts) {
  stopifnot(is.function(fun),is.numeric(starts),length(starts)>0L,
    all(is.finite(starts)),all(starts>0))
  count <- 0L
  walk <- function(x) {
    if(!is.call(x)) return(x)
    if(identical(x[[1]],as.name('<-')) && identical(x[[2]],as.name('sigma_bs'))) {
      count <<- count+1L
      x[[3]] <- substitute(STARTS[chain],list(STARTS=starts));return(x)
    }
    for(i in seq_along(x)[-1L]) if(!identical(x[[i]],quote(expr=))) x[i] <- list(walk(x[[i]]))
    x
  }
  b <- walk(body(fun))
  if(count!=1L) stop('Expected exactly one sigma_bs initializer; found ',count)
  body(fun) <- b
  fun
}

score_amplitude_fit <- function(fit,input) {
  js <- fit$results_output$jsdm_output
  summarise_field_draws(js$Bs_output,js$idx_ls_output,independent_bases(fit),
    input$truth$field,js$sigmabs_output)
}
