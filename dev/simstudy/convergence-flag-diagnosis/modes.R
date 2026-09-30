# Mode assignment for the Task 3 diagnostic fits (PLAN.md in this directory).
# It reads draws already in memory and never fits a model.
#
# Interfaces
#
#   assign_modes(draws, quantities = MODE_QUANTITIES,
#                transforms = MODE_TRANSFORMS[quantities])
#     `draws` is a named list of iterations x chains matrices, one per
#     quantity, for one species (species_mode_draws() builds it from
#     chain_anatomy(..., keep_draws = TRUE)). Returns a list:
#       labels: iterations x chains integer matrix, the mode of each draw;
#       prob_mode2: the same shape, the fitted probability of mode 2 (NA when
#         one mode);
#       n_modes: 1L or 2L;
#       share: per chain and mode, draws and share (chain, mode, draws, share);
#       overall: per mode, draws and share over all chains;
#       components: the two fitted components (weight, and mean and SD of each
#         quantity on its own scale), whether or not they form two modes;
#       quantities_used, dropped_constant, transforms, ridgeline_maxima,
#       loglik, em_iterations, em_converged, and em_starts (one row per EM
#       start: its log-likelihood, iterations and convergence).
#
#   chain_regions(result, min_share = VISIT_SHARE)
#     Per chain: the share of each mode, visits_both (each mode holds at least
#     min_share of the chain's draws) and region ('1', '2' or 'both').
#
#   mode_mass_table(result, run)
#     The share table joined with chain_regions(), for results/.../mode-mass.csv.
#
#   species_mode_draws(anatomy, species, quantities = MODE_QUANTITIES)
#     The draws list for one species from chain_anatomy(keep_draws = TRUE).
#
#   plot_chain_strips(draws, file, quantities = names(draws), truth = NULL,
#                     modes = NULL, title = NULL)
#     A figure for any number of chains (the Task 1 trace figure takes at most
#     eight): one panel per quantity, one row per chain, with the chain's
#     central 50% and 90% intervals and its mean; chains coloured by their
#     majority mode when `modes` is given, with each chain's share of mode 2
#     at the right; the generating value dashed.
#
# Method (frozen in README.md before any diagnostic fit).
#   1. Each quantity is used on its own scale (MODE_TRANSFORMS: identity for
#      all three; a logit option exists but is not used, see below), centred
#      and scaled by its pooled mean and SD over all draws. Probabilities
#      (theta0, mean_psi_original_sites) outside [0, 1] are refused. A quantity
#      with no variation (variant (b)'s fixed theta0) is dropped and reported.
#   2. A two-component Gaussian mixture with full covariances is fitted to the
#      pooled draws by EM from several deterministic starts (no random numbers
#      are used): k-means with two centres initialised at the means of the
#      draws below and above the median of the first principal component, and
#      for each quantity a split at its median. The fit with the highest
#      log-likelihood is kept.
#   3. The fitted mixture has two modes only if its density has two local
#      maxima along the ridgeline of the two components (Ray and Lindsay
#      2005, Annals of Statistics 33:2042-2065: every mode of a two-component
#      Gaussian mixture lies on that curve) and the smaller component's weight
#      is at least MIN_MODE_WEIGHT. Otherwise every draw is mode 1.
#   4. With two modes, each draw takes the component of higher fitted
#      probability; mode 1 is the component with the lower mean of the first
#      quantity used (theta0 unless it was dropped).
# Calibration on the saved pr11 fit of community 5 (Task 1 data, 4 chains of
# 12,000 draws, done before this method was frozen): on the logit scale with
# the principal-component start only, the low-theta0 mode's long left tail
# (theta0 down to 6e-6) widened its component and 170 of chain 2's draws
# (1.4%), with theta0 0.12 to 0.39 and so in the high region, were labelled
# mode 1; on the probability scale the principal-component start converged to
# one overlapping pair (one mode), while the per-quantity starts gave the two
# modes with at most 26 draws per chain (0.22%) labelled against the chain's
# region. test-modes.R holds both cases as synthetic regression tests.

`%||%` <- function(x,y) if(is.null(x)) y else x

MODE_QUANTITIES <- c('theta0','B0','mean_psi_original_sites')
MODE_TRANSFORMS <- c(theta0='identity',B0='identity',mean_psi_original_sites='identity')
PROBABILITY_QUANTITIES <- c('theta0','mean_psi_original_sites')
MIN_MODE_WEIGHT <- .01
VISIT_SHARE <- .01
EM_MAX_ITERATIONS <- 1000L
EM_TOLERANCE <- 1e-10
COVARIANCE_RIDGE <- 1e-8
RIDGELINE_POINTS <- 2001L

# ---- Mixture fitting ----------------------------------------------------------

log_normal_density <- function(X,mean,cov) {
  R <- chol(cov)
  Z <- backsolve(R,t(X)-mean,transpose=TRUE)
  -.5*colSums(Z^2)-sum(log(diag(R)))-.5*ncol(X)*log(2*pi)
}
log_sum_exp2 <- function(a,b) {m <- pmax(a,b);m+log(exp(a-m)+exp(b-m))}

# EM for two full-covariance Gaussian components from initial hard labels.
fit_two_gaussians <- function(X,start) {
  n <- nrow(X);d <- ncol(X)
  resp <- cbind(start==1L,start==2L)*1
  loglik <- -Inf;converged <- FALSE
  for(it in seq_len(EM_MAX_ITERATIONS)) {
    nk <- colSums(resp)
    if(any(nk<d+1)) break
    weight <- nk/n
    mean <- lapply(1:2,function(k) colSums(X*resp[,k])/nk[k])
    cov <- lapply(1:2,function(k) {
      Xc <- sweep(X,2,mean[[k]])
      crossprod(Xc*sqrt(resp[,k]))/nk[k]+diag(COVARIANCE_RIDGE,d)
    })
    lp1 <- log(weight[1])+log_normal_density(X,mean[[1]],cov[[1]])
    lp2 <- log(weight[2])+log_normal_density(X,mean[[2]],cov[[2]])
    total <- log_sum_exp2(lp1,lp2)
    new <- sum(total)
    resp <- cbind(exp(lp1-total),exp(lp2-total))
    if(is.finite(loglik) && abs(new-loglik)<=EM_TOLERANCE*abs(new)) {loglik <- new;converged <- TRUE;break}
    loglik <- new
  }
  nk <- colSums(resp)
  if(any(nk<d+1)) return(NULL)
  list(weight=weight,mean=mean,cov=cov,resp=resp,loglik=loglik,iterations=it,converged=converged)
}

# Number of local maxima of the mixture density along the ridgeline
# x(a) = [(1-a) S1^-1 + a S2^-1]^-1 [(1-a) S1^-1 m1 + a S2^-1 m2], a in [0, 1].
ridgeline_maxima <- function(fit,points=RIDGELINE_POINTS) {
  P1 <- solve(fit$cov[[1]]);P2 <- solve(fit$cov[[2]])
  m1 <- fit$mean[[1]];m2 <- fit$mean[[2]]
  a <- seq(0,1,length.out=points)
  x <- t(vapply(a,function(t) solve((1-t)*P1+t*P2,(1-t)*P1%*%m1+t*P2%*%m2),numeric(length(m1))))
  if(length(m1)==1L) x <- matrix(x,ncol=1L)
  h <- log_sum_exp2(log(fit$weight[1])+log_normal_density(x,m1,fit$cov[[1]]),
    log(fit$weight[2])+log_normal_density(x,m2,fit$cov[[2]]))
  up <- diff(h)
  interior <- sum(up[-length(up)]>0 & up[-1]<0)
  interior+(up[1]<0)+(up[length(up)]>0)
}

# Deterministic k-means start: centres at the means of the draws below and
# above the median of the first principal component.
kmeans_start <- function(X) {
  pc <- X%*%prcomp(X,center=FALSE,scale.=FALSE)$rotation[,1]
  low <- pc<=stats::median(pc)
  centres <- rbind(colMeans(X[low,,drop=FALSE]),colMeans(X[!low,,drop=FALSE]))
  stats::kmeans(X,centers=centres,iter.max=100L)$cluster
}

# All EM starts: the k-means start and a median split of each quantity.
em_starts <- function(X,names) {
  s <- list(principal_component=tryCatch(kmeans_start(X),error=function(e) NULL))
  for(j in seq_len(ncol(X))) s[[paste0('median_',names[j])]] <- ifelse(X[,j]>stats::median(X[,j]),2L,1L)
  s
}

# ---- Mode assignment ----------------------------------------------------------

transform_quantity <- function(x,how,name) {
  if(name %in% PROBABILITY_QUANTITIES && any(x<0 | x>1)) stop('Quantity ',name,' has draws outside [0, 1]')
  switch(how,identity=x,
    logit={
      if(any(x<=0 | x>=1)) stop('Quantity ',name,' has draws outside (0, 1); cannot take its logit')
      stats::qlogis(x)
    },
    stop('Unknown transform ',how,' for ',name))
}

assign_modes <- function(draws,quantities=MODE_QUANTITIES,transforms=MODE_TRANSFORMS[quantities]) {
  quantities <- as.character(quantities)
  absent <- setdiff(quantities,names(draws))
  if(length(absent)) stop('Quantities not in draws: ',paste(absent,collapse=', '))
  if(!length(quantities)) stop('No quantities given')
  transforms <- stats::setNames(as.character(transforms),quantities)
  transforms[is.na(transforms)] <- 'identity'
  mats <- lapply(draws[quantities],as.matrix)
  shape <- dim(mats[[1]])
  if(!all(vapply(mats,function(m) identical(dim(m),shape),logical(1))))
    stop('Every quantity must have the same iterations x chains shape')
  if(any(vapply(mats,function(m) any(!is.finite(m)),logical(1)))) stop('Draws must be finite')
  ni <- shape[1];nc <- shape[2]
  z <- lapply(quantities,function(q) transform_quantity(as.vector(mats[[q]]),transforms[[q]],q))
  names(z) <- quantities
  sds <- vapply(z,stats::sd,numeric(1))
  constant <- quantities[!(sds>0)]
  used <- setdiff(quantities,constant)
  one <- function(extra=list()) {
    labels <- matrix(1L,ni,nc)
    c(list(labels=labels,prob_mode2=matrix(NA_real_,ni,nc),n_modes=1L,
      share=chain_share(labels,1L),overall=overall_share(labels,1L),
      quantities_used=used,dropped_constant=constant,transforms=transforms[used]),extra)
  }
  if(!length(used)) return(one(list(components=NULL,ridgeline_maxima=NA_integer_,loglik=NA_real_,
    em_iterations=NA_integer_,em_converged=NA)))
  centre <- vapply(z[used],mean,numeric(1));scale <- sds[used]
  X <- vapply(used,function(q) (z[[q]]-centre[[q]])/scale[[q]],numeric(ni*nc))
  X <- matrix(X,ncol=length(used))
  starts <- em_starts(X,used)
  fits <- lapply(starts,function(st) if(is.null(st)) NULL else fit_two_gaussians(X,st))
  field <- function(f,k,na) if(is.null(f)) na else f[[k]]
  start_table <- data.frame(start=names(starts),loglik=vapply(fits,field,numeric(1),k='loglik',na=NA_real_),
    iterations=vapply(fits,function(f) as.integer(field(f,'iterations',NA_integer_)),integer(1)),
    converged=vapply(fits,field,logical(1),k='converged',na=NA),row.names=NULL,stringsAsFactors=FALSE)
  if(all(is.na(start_table$loglik))) return(one(list(components=NULL,ridgeline_maxima=1L,loglik=NA_real_,
    em_iterations=NA_integer_,em_converged=FALSE,em_starts=start_table)))
  fit <- fits[[which.max(start_table$loglik)]]
  # Order components so that mode 1 has the lower mean of the first quantity.
  o <- order(vapply(fit$mean,`[`,numeric(1),1L))
  fit$weight <- fit$weight[o];fit$mean <- fit$mean[o];fit$cov <- fit$cov[o];fit$resp <- fit$resp[,o]
  components <- data.frame(mode=1:2,weight=fit$weight)
  back <- function(v,q) if(transforms[[q]]=='logit') stats::plogis(v) else v
  for(j in seq_along(used)) {
    q <- used[j]
    # Mean and SD on the quantity's own scale, from the draws each component claims.
    w <- fit$resp
    components[[paste0('mean_',q)]] <- vapply(1:2,function(k) sum(w[,k]*mats[[q]])/sum(w[,k]),numeric(1))
    components[[paste0('sd_',q)]] <- vapply(1:2,function(k) {
      m <- sum(w[,k]*mats[[q]])/sum(w[,k]);sqrt(sum(w[,k]*(mats[[q]]-m)^2)/sum(w[,k]))
    },numeric(1))
    components[[paste0('centre_',q)]] <- vapply(1:2,function(k) back(fit$mean[[k]][j]*scale[[q]]+centre[[q]],q),numeric(1))
  }
  maxima <- ridgeline_maxima(fit)
  extra <- list(components=components,ridgeline_maxima=as.integer(maxima),loglik=fit$loglik,
    em_iterations=as.integer(fit$iterations),em_converged=fit$converged,em_starts=start_table)
  if(maxima<2L || min(fit$weight)<MIN_MODE_WEIGHT) return(one(extra))
  labels <- matrix(ifelse(fit$resp[,2]>fit$resp[,1],2L,1L),ni,nc)
  c(list(labels=labels,prob_mode2=matrix(fit$resp[,2],ni,nc),n_modes=2L,
    share=chain_share(labels,2L),overall=overall_share(labels,2L),
    quantities_used=used,dropped_constant=constant,transforms=transforms[used]),extra)
}

chain_share <- function(labels,n_modes) {
  nc <- ncol(labels);modes <- seq_len(n_modes)
  out <- data.frame(chain=rep(seq_len(nc),each=n_modes),mode=rep(modes,nc))
  out$draws <- as.integer(mapply(function(ch,m) sum(labels[,ch]==m),out$chain,out$mode))
  out$share <- out$draws/nrow(labels)
  out
}
overall_share <- function(labels,n_modes) {
  modes <- seq_len(n_modes)
  data.frame(mode=modes,draws=vapply(modes,function(m) sum(labels==m),integer(1)),
    share=vapply(modes,function(m) mean(labels==m),numeric(1)))
}

chain_regions <- function(result,min_share=VISIT_SHARE) {
  s <- result$share;chains <- sort(unique(s$chain))
  share_of <- function(ch,m) {x <- s$share[s$chain==ch & s$mode==m];if(length(x)) x else 0}
  s1 <- vapply(chains,share_of,numeric(1),m=1L);s2 <- vapply(chains,share_of,numeric(1),m=2L)
  both <- s1>=min_share & s2>=min_share
  data.frame(chain=chains,share_mode1=s1,share_mode2=s2,visits_both=both,
    region=ifelse(both,'both',ifelse(s1>=s2,'1','2')),stringsAsFactors=FALSE)
}

mode_mass_table <- function(result,run) {
  r <- chain_regions(result)
  m <- merge(result$share,r[c('chain','visits_both','region')],by='chain',sort=FALSE)
  m <- m[order(m$chain,m$mode),,drop=FALSE];rownames(m) <- NULL
  data.frame(run=run,m[c('chain','mode','draws','share')],n_modes=result$n_modes,
    m[c('visits_both','region')],stringsAsFactors=FALSE)
}

species_mode_draws <- function(anatomy,species,quantities=MODE_QUANTITIES) {
  draws <- attr(anatomy,'draws')
  if(is.null(draws)) stop('species_mode_draws needs chain_anatomy(..., keep_draws = TRUE)')
  absent <- setdiff(quantities,names(draws))
  if(length(absent)) stop('Anatomy draws lack ',paste(absent,collapse=', '))
  lapply(stats::setNames(quantities,quantities),function(q) {
    x <- draws[[q]];matrix(x[species,,],dim(x)[2],dim(x)[3])
  })
}

# ---- Figure for many chains -----------------------------------------------------

# Categorical slots 1 and 2 of the dataviz reference palette (the pair passes
# the lightness, chroma, CVD, normal-vision and contrast checks on this
# surface), and its ink tokens.
MODE_COLOURS <- c('#2a78d6','#eb6834')
STRIP_INK <- c(primary='#0b0b0b',secondary='#52514e',muted='#8a8984',grid='#e4e3df',surface='#fcfcfb')

plot_chain_strips <- function(draws,file,quantities=names(draws),truth=NULL,modes=NULL,title=NULL) {
  absent <- setdiff(quantities,names(draws))
  if(length(absent)) stop('Quantities not in draws: ',paste(absent,collapse=', '))
  mats <- lapply(draws[quantities],as.matrix);K <- ncol(mats[[1]])
  majority <- rep(1L,K);share2 <- rep(NA_real_,K)
  if(!is.null(modes)) {
    r <- chain_regions(modes)
    if(nrow(r)!=K) stop('Mode result has ',nrow(r),' chains; draws have ',K)
    majority <- ifelse(r$share_mode2>r$share_mode1,2L,1L);share2 <- r$share_mode2
  }
  n <- length(quantities);columns <- min(3L,n);rows <- ceiling(n/columns)
  grDevices::png(file,width=520*columns+60,height=(90+22*K)*rows+150,res=110,bg=STRIP_INK[['surface']],
    type=if(isTRUE(capabilities('cairo'))) 'cairo' else getOption('bitmapType'))
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::par(mfrow=c(rows,columns),mar=c(3.2,4.2,2.6,if(is.null(modes)) 1 else 4.2),oma=c(2.8,.5,3.2,0),
    col.axis=STRIP_INK[['secondary']],col.lab=STRIP_INK[['secondary']],fg=STRIP_INK[['muted']],las=1,cex.axis=.8)
  for(q in quantities) {
    x <- mats[[q]]
    qs <- apply(x,2,stats::quantile,probs=c(.05,.25,.75,.95),names=FALSE)
    mu <- colMeans(x);tv <- if(!is.null(truth) && q %in% names(truth)) truth[[q]] else NA_real_
    xlim <- range(c(qs,tv),na.rm=TRUE);if(diff(xlim)==0) xlim <- xlim+c(-.5,.5)*max(abs(xlim[1]),1)*.1
    y <- rev(seq_len(K))
    graphics::plot(NA,xlim=xlim,ylim=c(.5,K+.5),xlab='',ylab='',yaxt='n',bty='n')
    graphics::abline(v=pretty(xlim),col=STRIP_INK[['grid']],lwd=.8)
    graphics::axis(2,at=y,labels=paste('Chain',seq_len(K)),tick=FALSE,cex.axis=.7)
    col <- MODE_COLOURS[majority]
    graphics::segments(qs[1,],y,qs[4,],y,col=col,lwd=1.2,lend=1)
    graphics::segments(qs[2,],y,qs[3,],y,col=col,lwd=4,lend=1)
    graphics::points(mu,y,pch=21,bg=col,col=STRIP_INK[['surface']],cex=1.1,lwd=1.2)
    if(is.finite(tv)) graphics::abline(v=tv,lty=2,lwd=1.4,col=STRIP_INK[['primary']])
    if(!is.null(modes)) graphics::mtext(sprintf('%.0f%%',100*share2),side=4,at=y,line=.3,cex=.6,
      col=STRIP_INK[['secondary']],las=1)
    graphics::mtext(q,side=3,line=.9,adj=0,cex=.85,font=2,col=STRIP_INK[['primary']])
  }
  main <- title %||% sprintf('Per-chain distributions, %d chains (thin line 90%%, thick 50%%, dot mean)',K)
  graphics::mtext(main,outer=TRUE,side=3,line=1.8,adj=0,at=.01,cex=.95,font=2,col=STRIP_INK[['primary']])
  if(!is.null(modes)) graphics::mtext(sprintf('%d mode(s) found; right margin: share of the chain\'s draws in mode 2',
    modes$n_modes),outer=TRUE,side=3,line=.6,adj=0,at=.01,cex=.75,col=STRIP_INK[['secondary']])
  graphics::par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  graphics::plot(0,0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
  keys <- if(is.null(modes)) character() else c('Mostly mode 1 (lower theta0)','Mostly mode 2')
  if(length(keys) || !is.null(truth)) graphics::legend('bottom',legend=c(keys,if(!is.null(truth)) 'Generating value'),
    col=c(MODE_COLOURS[seq_along(keys)],if(!is.null(truth)) STRIP_INK[['primary']]),
    lty=c(rep(1,length(keys)),if(!is.null(truth)) 2),lwd=c(rep(4,length(keys)),if(!is.null(truth)) 1.4),
    horiz=TRUE,bty='n',cex=.8,text.col=STRIP_INK[['secondary']],inset=c(0,.005),xpd=NA)
  invisible(file)
}
