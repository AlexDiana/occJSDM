# These summaries describe differences among observed chains. They are not
# confidence intervals, bounds on MCMC error, or evidence about unvisited modes.
score_observed_probability_chains <- function(result) {
  x <- result$chain_probability
  cells <- result$elements[result$elements$metric=='occupancy',]
  n <- nrow(result$probability);s <- ncol(result$probability)
  stopifnot(is.matrix(x),nrow(x)==n*s,ncol(x)>=2L,
    nrow(cells)==n*s,all(cells$element==seq_len(n*s)),all(is.finite(x)),
    max(abs(rowMeans(x)-as.vector(result$probability)))<1e-12,
    nrow(result$species)==s)
  truth <- cells$truth
  masks <- list(all=rep(TRUE,n*s),low=truth<.2,
    medium=truth>=.2 & truth<=.8,high=truth>.8)
  for(p in c(.01,.05,.25,.75)) masks[[paste0('prevalence_',100*p,'pct')]] <-
    rep(abs(result$species$target-p)<1e-10,each=n)
  rows <- list()
  for(g in names(masks)) {
    idx <- which(masks[[g]])
    stopifnot(length(idx)>0L)
    for(ch in seq_len(ncol(x))) {
      error <- x[idx,ch]-truth[idx]
      rows[[length(rows)+1L]] <- data.frame(group=g,chain=ch,
        bias=mean(error),mae=mean(abs(error)),rmse=sqrt(mean(error^2)))
    }
  }
  do.call(rbind,rows)
}

collapse_observed_chains <- function(chains,pooled) {
  keys <- c('community','grid_index','replicate','arm','knots','group')
  values <- c('bias','mae','rmse')
  stopifnot(!anyDuplicated(pooled[keys]),!anyDuplicated(chains[c(keys,'chain')]))
  cells <- split(chains,interaction(chains$community,chains$arm,chains$knots,chains$group,drop=TRUE))
  rows <- lapply(cells,function(d) {
    stopifnot(nrow(d)>=2L,nrow(unique(d[keys]))==1L)
    do.call(rbind,lapply(values,function(q) {
      stopifnot(all(is.finite(d[[q]])))
      data.frame(d[1,keys],quantity=q,chain_min=min(d[[q]]),chain_max=max(d[[q]]))
    }))
  })
  ranges <- do.call(rbind,rows)
  point <- do.call(rbind,lapply(values,function(q)data.frame(pooled[keys],quantity=q,pooled=pooled[[q]])))
  joined <- merge(ranges,point,by=c(keys,'quantity'),all=TRUE)
  stopifnot(nrow(joined)==nrow(ranges),nrow(joined)==nrow(point),
    all(is.finite(as.matrix(joined[c('chain_min','chain_max','pooled')]))))
  joined
}

observed_chain_average <- function(d) {
  stopifnot(nrow(d)==9L,!anyDuplicated(d$community),
    setequal(unique(d$grid_index),c(4,6,8)),all(table(d$grid_index)==3L))
  c(pooled=mean(d$pooled),chain_min_mean=mean(d$chain_min),
    chain_max_mean=mean(d$chain_max),n=nrow(d))
}

summarise_observed_chains <- function(chains,pooled) {
  d <- collapse_observed_chains(chains,pooled)
  cells <- split(d,interaction(d$arm,d$knots,d$group,d$quantity,drop=TRUE))
  do.call(rbind,lapply(cells,function(x)data.frame(
    x[1,c('arm','knots','group','quantity')],as.list(observed_chain_average(x)))))
}

paired_observed_chains <- function(chains,pooled) {
  d <- collapse_observed_chains(chains,pooled)
  comparisons <- list()
  for(arm in c('binary','low','high')) for(k in c(50L,100L))
    comparisons[[paste0(arm,'_k',k,'_minus_k20')]] <- list(
      d[d$arm==arm & d$knots==k,],d[d$arm==arm & d$knots==20L,])
  for(k in c(20L,50L,100L)) for(arm in c('low','high'))
    comparisons[[paste0(arm,'_minus_binary_k',k)]] <- list(
      d[d$arm==arm & d$knots==k,],d[d$arm=='binary' & d$knots==k,])
  rows <- list()
  for(name in names(comparisons)) {
    pair <- comparisons[[name]]
    x <- merge(pair[[1]],pair[[2]],by=c('community','grid_index','replicate','group','quantity'),
      all=TRUE,suffixes=c('_a','_b'))
    stopifnot(nrow(x)>0L,nrow(x)==nrow(pair[[1]]),nrow(x)==nrow(pair[[2]]))
    x$pooled <- x$pooled_a-x$pooled_b
    x$chain_min <- x$chain_min_a-x$chain_max_b
    x$chain_max <- x$chain_max_a-x$chain_min_b
    stopifnot(all(is.finite(as.matrix(x[c('pooled','chain_min','chain_max')]))))
    cells <- split(x,interaction(x$group,x$quantity,drop=TRUE))
    rows[[name]] <- do.call(rbind,lapply(cells,function(a)data.frame(
      comparison=name,a[1,c('group','quantity')],as.list(observed_chain_average(a)))))
  }
  do.call(rbind,rows)
}

diagnose_range_chains <- function(d) {
  stopifnot(all(is.finite(as.matrix(d))),!anyDuplicated(d[c('chain','range')]),
    all(d$frequency>=0),all(d$frequency<=1))
  cells <- split(d,d$chain)
  grid <- sort(unique(d$range))
  frequency <- vapply(cells,function(x) {
    stopifnot(setequal(x$range,grid),abs(sum(x$frequency)-1)<1e-12)
    x$frequency[match(grid,x$range)]
  },numeric(length(grid)))
  stopifnot(ncol(frequency)>=2L)
  visited <- colSums(frequency>0)
  constant <- visited==1L
  state <- if(all(constant)) {
    if(length(unique(max.col(t(frequency))))==1L)'same_point_mass' else 'distinct_point_masses'
  } else if(any(constant)) 'some_constant_chains' else 'all_variable'
  pairs <- combn(ncol(frequency),2)
  tv <- apply(pairs,2,function(pair)sum(abs(frequency[,pair[1]]-frequency[,pair[2]]))/2)
  data.frame(state=state,chains=ncol(frequency),constant_chains=sum(constant),
    grid_points_visited=sum(rowSums(frequency)>0),
    chain_mean_gap=diff(range(colSums(frequency*grid))),max_total_variation=max(tv))
}
