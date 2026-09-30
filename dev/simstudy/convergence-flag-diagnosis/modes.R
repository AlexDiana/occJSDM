# Mode assignment for the Task 3 diagnostic fits (PLAN.md in this directory).
# It reads draws already in memory and never fits a model.
#
# AMENDMENT-1.md (ruling R17) makes the anchored classifier below the primary
# assignment and the refitted mixture (assign_modes) and the theta0 cut
# secondary checks reported beside it. README.md's method text describes the
# refitted mixture as first frozen; AMENDMENT-1 supersedes it.
#
# Interfaces: primary (anchored) assignment
#
#   calibrate_anchor(draws, near_chains, mirror_chains,
#                    quantities = ANCHOR_QUANTITIES)
#     A two-component classifier: component 1 ('near-truth') the mean and
#     covariance of the given quantities over near_chains, component 2
#     ('mirror') the same over mirror_chains, equal weights. The frozen one
#     (results/modes/anchor-classifier.csv, md5 ANCHOR_MD5) is calibrated on
#     the saved pr11 fit of design-qfar_K6-sites300-05 species 6, chains 1 and
#     3 near-truth and 2 and 4 mirror (verify.R --mode=anchor-calibration).
#   write_anchor(anchor, file), read_anchor(file, md5 = ANCHOR_MD5)
#     The classifier as a CSV of weights, means and covariances, each number
#     both to 17 significant digits (`value`, for reading) and as a
#     hexadecimal float (`hex`, which R reads back exactly; its decimal
#     conversion is not exact for every 17-digit value); read_anchor refuses
#     a file whose md5 differs (md5 = NULL skips the check).
#   assign_anchored(draws, anchor)
#     Each draw takes the component of higher density (equal weights). Returns
#     method 'anchored', labels (1 near-truth, 2 mirror), log_ratio (log
#     density of mirror minus near-truth), mode_names, share, overall and
#     atypical: per chain, the share of draws whose squared Mahalanobis
#     distance to both components exceeds the ATYPICAL_LEVEL quantile of the
#     chi-squared distribution (reported, never deciding). theta0 and B0 are
#     not used, so the assignment is the same whether they vary or are fixed.
#   anchored_regions(result, min_share = VISIT_SHARE)
#     Per chain: share_near_truth, share_mirror, visits_both (each at least
#     min_share) and region ('near-truth', 'mirror' or 'both').
#   region_pattern(regions)
#     'every chain visits both', 'each chain in one region', 'one mode only:
#     near-truth', 'one mode only: mirror' or 'other' (AMENDMENT-1, R15).
#
# Interfaces: secondary checks
#
#   theta0_cut_assignment(draws, cut = THETA0_CUT)
#     Mirror when theta0 exceeds 0.135, the valley of the pooled theta0
#     density of the pr11 fit; only where theta0 is free (the extended run and
#     variant a). Same result structure as assign_anchored().
#   compare_assignments(primary, cut = NULL, refit = NULL)
#     Per chain: the primary region and share of mirror draws, the same for
#     the theta0 cut and the refitted mixture (its modes named by majority
#     overlap with the primary labels), and whether each agrees with the
#     primary; attribute draw_agreement gives the share of draws labelled the
#     same.
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
#     The refitted mixture, a secondary check under AMENDMENT-1.
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
#      for each quantity a split at its median; and (added by AMENDMENT-1)
#      chain partitions: each chain against the rest, and for each quantity
#      the chains split at the largest gap in their chain means. The fit with
#      the highest log-likelihood is kept.
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
# Review of that method (AMENDMENT-1): on pseudo-chains cut from the same fit
# it failed when one mode held 1 or 2 of 16 chains or 1 of 8. The
# chain-partition starts fix the cases where EM stuck (14 and 2, 2 and 14, 7
# and 1, 1 and 7), but with 15 and 1 or 1 and 15 the mixture's likelihood
# itself prefers splitting the majority mode's non-Gaussian shape, which no
# start fixes. Hence the anchored classifier as the primary assignment.

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
# Anchored classifier (AMENDMENT-1, ruling R17): quantities that stay free in
# every variant (theta0 is fixed in (b) and squeezed by the prior in (a)).
ANCHOR_QUANTITIES <- c('B_slope1','B_slope2','beta_theta_intercept','mean_psi_original_sites')
ANCHOR_MODES <- c('near-truth','mirror')
ANCHOR_FILE <- 'results/modes/anchor-classifier.csv'
ANCHOR_MD5 <- 'e235c2fa641eb05c36232bec0d513bc8'
ATYPICAL_LEVEL <- .999
THETA0_CUT <- .135

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

# All EM starts: the k-means start and a median split of each quantity, and
# with ni draws in each of nc chains (nc at least 2) the chain partitions:
# each chain against the rest, and for each quantity the chains split at the
# largest gap in their chain means.
em_starts <- function(X,names,ni=NULL,nc=NULL) {
  s <- list(principal_component=tryCatch(kmeans_start(X),error=function(e) NULL))
  for(j in seq_len(ncol(X))) s[[paste0('median_',names[j])]] <- ifelse(X[,j]>stats::median(X[,j]),2L,1L)
  if(!is.null(nc) && nc>=2L) {
    for(k in seq_len(nc)) s[[paste0('chain_',k)]] <- rep(ifelse(seq_len(nc)==k,2L,1L),each=ni)
    for(j in seq_len(ncol(X))) {
      mu <- colMeans(matrix(X[,j],ni,nc));o <- order(mu);cut <- which.max(diff(mu[o]))
      g <- rep(1L,nc);g[o[-seq_len(cut)]] <- 2L
      s[[paste0('chain_gap_',names[j])]] <- rep(g,each=ni)
    }
  }
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
    c(list(method='refit',mode_names='single mode',labels=labels,prob_mode2=matrix(NA_real_,ni,nc),n_modes=1L,
      share=chain_share(labels,1L),overall=overall_share(labels,1L),
      quantities_used=used,dropped_constant=constant,transforms=transforms[used]),extra)
  }
  if(!length(used)) return(one(list(components=NULL,ridgeline_maxima=NA_integer_,loglik=NA_real_,
    em_iterations=NA_integer_,em_converged=NA)))
  centre <- vapply(z[used],mean,numeric(1));scale <- sds[used]
  X <- vapply(used,function(q) (z[[q]]-centre[[q]])/scale[[q]],numeric(ni*nc))
  X <- matrix(X,ncol=length(used))
  starts <- em_starts(X,used,ni,nc)
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
  c(list(method='refit',mode_names=paste0(c('lower ','higher '),used[1]),labels=labels,
    prob_mode2=matrix(fit$resp[,2],ni,nc),n_modes=2L,
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

# One row per chain and mode. For the anchored and theta0-cut assignments the
# modes are named (near-truth, mirror) and n_modes is NA; for the refitted
# mixture n_modes is the number of modes it found.
mode_mass_table <- function(result,run) {
  method <- result$method %||% 'refit'
  labels <- result$mode_names %||% as.character(seq_len(result$n_modes))
  r <- if(method=='refit') chain_regions(result) else {
    a <- anchored_regions(result);data.frame(chain=a$chain,visits_both=a$visits_both,region=a$region)
  }
  m <- merge(result$share,r[c('chain','visits_both','region')],by='chain',sort=FALSE)
  m <- m[order(m$chain,m$mode),,drop=FALSE];rownames(m) <- NULL
  data.frame(run=run,method=method,m[c('chain','mode')],mode_name=labels[m$mode],m[c('draws','share')],
    n_modes=if(method=='refit') result$n_modes else NA_integer_,m[c('visits_both','region')],stringsAsFactors=FALSE)
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

# ---- Anchored assignment (AMENDMENT-1, ruling R17) -----------------------------

draw_matrix <- function(draws,quantities) {
  absent <- setdiff(quantities,names(draws))
  if(length(absent)) stop('Quantities not in draws: ',paste(absent,collapse=', '))
  mats <- lapply(draws[quantities],as.matrix);shape <- dim(mats[[1]])
  if(!all(vapply(mats,function(m) identical(dim(m),shape),logical(1))))
    stop('Every quantity must have the same iterations x chains shape')
  if(any(vapply(mats,function(m) any(!is.finite(m)),logical(1)))) stop('Draws must be finite')
  Z <- matrix(vapply(mats,as.vector,numeric(prod(shape))),ncol=length(quantities),dimnames=list(NULL,quantities))
  list(Z=Z,ni=shape[1],nc=shape[2])
}

calibrate_anchor <- function(draws,near_chains,mirror_chains,quantities=ANCHOR_QUANTITIES) {
  if(length(intersect(near_chains,mirror_chains))) stop('A chain cannot calibrate both components')
  pick <- function(ch) draw_matrix(lapply(draws[quantities],function(m) as.matrix(m)[,ch,drop=FALSE]),quantities)$Z
  near <- pick(near_chains);mirror <- pick(mirror_chains)
  list(quantities=quantities,modes=ANCHOR_MODES,weight=c(.5,.5),
    mean=stats::setNames(list(colMeans(near),colMeans(mirror)),ANCHOR_MODES),
    cov=stats::setNames(list(stats::cov(near),stats::cov(mirror)),ANCHOR_MODES),
    n_draws=c(nrow(near),nrow(mirror)),near_chains=as.integer(near_chains),mirror_chains=as.integer(mirror_chains))
}

number_text <- function(x) sprintf('%.17g',x)
hex_text <- function(x) sprintf('%a',x)

write_anchor <- function(anchor,file,source=character()) {
  q <- anchor$quantities;rows <- list()
  add <- function(component,kind,row,col,value,hex='') rows[[length(rows)+1L]] <<-
    data.frame(component=component,kind=kind,row=row,col=col,value=value,hex=hex,stringsAsFactors=FALSE)
  num <- function(component,kind,row,col,x) add(component,kind,row,col,number_text(x),hex_text(x))
  for(k in seq_along(anchor$modes)) {
    m <- anchor$modes[k]
    num(m,'weight','','',anchor$weight[k])
    add(m,'n_draws','','',as.character(anchor$n_draws[k]))
    add(m,'chains','','',paste(if(k==1L) anchor$near_chains else anchor$mirror_chains,collapse=';'))
    num(m,'mean',q,'',unname(anchor$mean[[m]][q]))
    S <- anchor$cov[[m]]
    for(i in q) num(m,'cov',i,q,unname(S[i,q]))
  }
  for(n in names(source)) add('','source',n,'',source[[n]])
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(do.call(rbind,rows),file,row.names=FALSE)
  invisible(file)
}

read_anchor <- function(file,md5=ANCHOR_MD5) {
  if(!file.exists(file)) stop('Missing anchor classifier: ',file)
  if(!is.null(md5)) {
    found <- unname(tools::md5sum(file))
    if(!identical(found,md5)) stop('Anchor classifier md5 mismatch for ',file,': expected ',md5,', found ',found)
  }
  x <- utils::read.csv(file,colClasses='character')
  modes <- unique(x$component[x$kind=='mean']);if(!identical(modes,ANCHOR_MODES)) stop('Malformed anchor classifier: ',file)
  q <- x$row[x$kind=='mean' & x$component==modes[1]]
  get <- function(m,kind) x[x$component==m & x$kind==kind,,drop=FALSE]
  mean <- lapply(modes,function(m) {g <- get(m,'mean');stats::setNames(as.numeric(g$hex[match(q,g$row)]),q)})
  cov <- lapply(modes,function(m) {
    g <- get(m,'cov');S <- matrix(NA_real_,length(q),length(q),dimnames=list(q,q))
    S[cbind(match(g$row,q),match(g$col,q))] <- as.numeric(g$hex)
    if(anyNA(S)) stop('Incomplete covariance in ',file);S
  })
  chains <- function(m) as.integer(strsplit(get(m,'chains')$value,';',fixed=TRUE)[[1]])
  list(quantities=q,modes=modes,weight=vapply(modes,function(m) as.numeric(get(m,'weight')$hex),numeric(1),USE.NAMES=FALSE),
    mean=stats::setNames(mean,modes),cov=stats::setNames(cov,modes),
    n_draws=vapply(modes,function(m) as.integer(get(m,'n_draws')$value),integer(1),USE.NAMES=FALSE),
    near_chains=chains(modes[1]),mirror_chains=chains(modes[2]),
    source=stats::setNames(x$value[x$kind=='source'],x$row[x$kind=='source']))
}

named_result <- function(method,labels,mode_names=ANCHOR_MODES,extra=list())
  c(list(method=method,labels=labels,n_modes=2L,mode_names=mode_names,
    share=chain_share(labels,2L),overall=overall_share(labels,2L)),extra)

assign_anchored <- function(draws,anchor) {
  dm <- draw_matrix(draws,anchor$quantities);Z <- dm$Z;p <- ncol(Z)
  l <- lapply(1:2,function(k) log(anchor$weight[k])+log_normal_density(Z,anchor$mean[[k]],anchor$cov[[k]]))
  labels <- matrix(ifelse(l[[2]]>l[[1]],2L,1L),dm$ni,dm$nc)
  d2 <- vapply(1:2,function(k) stats::mahalanobis(Z,anchor$mean[[k]],anchor$cov[[k]]),numeric(nrow(Z)))
  far <- matrix(pmin(d2[,1],d2[,2])>stats::qchisq(ATYPICAL_LEVEL,p),dm$ni,dm$nc)
  named_result('anchored',labels,anchor$modes,list(log_ratio=matrix(l[[2]]-l[[1]],dm$ni,dm$nc),
    atypical=data.frame(chain=seq_len(dm$nc),share=colMeans(far)),quantities_used=anchor$quantities))
}

theta0_cut_assignment <- function(draws,cut=THETA0_CUT) {
  th <- draw_matrix(draws,'theta0')
  if(stats::sd(th$Z[,1])==0) stop('theta0 is constant; the theta0 cut applies only where theta0 is free')
  named_result('theta0 cut',matrix(ifelse(th$Z[,1]>cut,2L,1L),th$ni,th$nc),extra=list(cut=cut))
}

anchored_regions <- function(result,min_share=VISIT_SHARE) {
  names <- result$mode_names %||% ANCHOR_MODES
  r <- chain_regions(result,min_share)
  data.frame(chain=r$chain,share_near_truth=r$share_mode1,share_mirror=r$share_mode2,visits_both=r$visits_both,
    region=ifelse(r$region=='both','both',names[match(r$region,c('1','2'))]),stringsAsFactors=FALSE)
}

region_pattern <- function(regions) {
  r <- regions$region
  if(all(r=='both')) return('every chain visits both')
  if(all(r %in% ANCHOR_MODES)) return(if(length(unique(r))==1L) paste('one mode only:',r[1]) else 'each chain in one region')
  'other'
}

# The refitted mixture's labels renamed by the primary mode most of each
# refitted mode's draws carry, coded 1 (near-truth) and 2 (mirror).
refit_in_anchor_names <- function(refit,primary) {
  map <- vapply(sort(unique(as.vector(refit$labels))),function(m) {
    t <- tabulate(primary$labels[refit$labels==m],2L);which.max(t)
  },integer(1))
  names(map) <- sort(unique(as.vector(refit$labels)))
  matrix(unname(map[as.character(refit$labels)]),nrow(refit$labels),ncol(refit$labels))
}

compare_assignments <- function(primary,cut=NULL,refit=NULL) {
  pr <- anchored_regions(primary)
  out <- data.frame(chain=pr$chain,primary_region=pr$region,primary_share_mirror=pr$share_mirror,stringsAsFactors=FALSE)
  agreement <- c()
  if(!is.null(cut)) {
    cr <- anchored_regions(cut)
    out$cut_region <- cr$region;out$cut_share_mirror <- cr$share_mirror;out$agree_cut <- cr$region==pr$region
    agreement['cut'] <- mean(cut$labels==primary$labels)
  }
  if(!is.null(refit)) {
    coded <- refit_in_anchor_names(refit,primary)
    rr <- anchored_regions(named_result('refit',coded))
    out$refit_n_modes <- refit$n_modes;out$refit_region <- rr$region;out$refit_share_mirror <- rr$share_mirror
    out$agree_refit <- rr$region==pr$region
    agreement['refit'] <- mean(coded==primary$labels)
  }
  attr(out,'draw_agreement') <- agreement
  out
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
  names2 <- if(is.null(modes)) NULL else if(identical(modes$method %||% 'refit','refit'))
    c('mode 1 (lower theta0)','mode 2') else modes$mode_names
  if(!is.null(modes)) graphics::mtext(if(identical(modes$method %||% 'refit','refit'))
    sprintf('Refitted mixture: %d mode(s) found; right margin: share of the chain\'s draws in mode 2',modes$n_modes) else
    sprintf('%s assignment; right margin: share of the chain\'s draws in the %s mode',
      modes$method,names2[2]),outer=TRUE,side=3,line=.6,adj=0,at=.01,cex=.75,col=STRIP_INK[['secondary']])
  graphics::par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  graphics::plot(0,0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
  keys <- if(is.null(modes)) character() else paste('Mostly',names2)
  if(length(keys) || !is.null(truth)) graphics::legend('bottom',legend=c(keys,if(!is.null(truth)) 'Generating value'),
    col=c(MODE_COLOURS[seq_along(keys)],if(!is.null(truth)) STRIP_INK[['primary']]),
    lty=c(rep(1,length(keys)),if(!is.null(truth)) 2),lwd=c(rep(4,length(keys)),if(!is.null(truth)) 1.4),
    horiz=TRUE,bty='n',cex=.8,text.col=STRIP_INK[['secondary']],inset=c(0,.005),xpd=NA)
  invisible(file)
}
