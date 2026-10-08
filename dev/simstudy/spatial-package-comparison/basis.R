# A BasisSpec never retains responses. All transforms are trained exactly once.
fit_transform <- function(x) {
  x <- as.matrix(x); center <- colMeans(x); scale <- apply(x,2,sd)
  if (any(!is.finite(scale) | scale<=0)) stop("Cannot standardise a constant predictor")
  list(center=center,scale=scale)
}
apply_transform <- function(x,transform) {
  sweep(sweep(as.matrix(x),2,transform$center,"-"),2,transform$scale,"/")
}
trend_columns <- function(xy) {
  cbind(x=unname(xy[,1]),y=unname(xy[,2]),x2=unname(xy[,1]^2),xy=unname(xy[,1]*xy[,2]),y2=unname(xy[,2]^2))
}
xy_distances <- function(a,b) {
  sqrt(outer(a[,1],b[,1],"-")^2+outer(a[,2],b[,2],"-")^2)
}
graph_components <- function(W) {
  unseen <- seq_len(nrow(W)); count <- 0L
  while(length(unseen)) {
    count <- count+1L; frontier <- unseen[1]; visited <- integer()
    while(length(frontier)) {
      visited <- union(visited,frontier)
      frontier <- setdiff(which(colSums(W[frontier,,drop=FALSE])>0),visited)
    }
    unseen <- setdiff(unseen,visited)
  }
  count
}
fit_neighbor_graph <- function(xy,site_id,k=8L) {
  n <- nrow(xy); k <- min(k,n-1L)
  if (k<1L) stop("At least two training sites required")
  d <- xy_distances(xy,xy); diag(d) <- Inf
  W <- matrix(0,n,n); radius <- numeric(n); boundary_id <- character(n)
  for (i in seq_len(n)) {
    chosen <- order(d[i,],site_id,method="radix")[seq_len(k)]
    W[i,chosen] <- 1
    radius[i] <- d[i,chosen[k]]; boundary_id[i] <- site_id[chosen[k]]
  }
  W <- 1*((W+t(W))>0)
  list(W=W,k=k,radius=radius,boundary_id=boundary_id,components=graph_components(W))
}
fit_basis <- function(training,arm) {
  if (!arm %in% c("env","trend","mem","dnn","gp","spatial2")) stop("Unknown spatial arm")
  xy <- as.matrix(training$xy); ids <- as.character(training$site_id); n <- nrow(xy)
  if (ncol(xy)!=2L || any(!is.finite(xy)) || length(ids)!=n || anyNA(ids) || anyDuplicated(ids) || anyDuplicated(xy)) stop("Invalid training coordinates or site IDs")
  b <- list(schema_version=1L,arm=arm,status="available",reason=NA_character_,
    dataset_id=training$dataset_id,training_ids=ids,training_xy=xy,
    environment_transform=fit_transform(training$environment))
  if (arm=="env") b$training_matrix <- matrix(numeric(),n,0L)
  if (arm %in% c("gp","spatial2")) b$training_matrix <- xy
  if (arm %in% c("trend","dnn")) {
    b$coordinate_transform <- fit_transform(xy)
    standardized <- apply_transform(xy,b$coordinate_transform)
    if (arm=="dnn") b$training_matrix <- standardized else {
      raw <- trend_columns(standardized); b$column_transform <- fit_transform(raw)
      b$training_matrix <- apply_transform(raw,b$column_transform)
    }
  }
  if (arm=="mem") {
    b$graph <- fit_neighbor_graph(xy,ids)
    W <- b$graph$W; centered <- sweep(sweep(W,1,rowMeans(W),"-"),2,colMeans(W),"-")+mean(W)
    eig <- eigen(centered,symmetric=TRUE)
    b$positive_tolerance <- 1e-10*max(1,max(abs(eig$values)))
    keep <- head(which(eig$values>b$positive_tolerance),10L)
    b$eigenvalues <- eig$values[keep]; b$eigenvectors <- eig$vectors[,keep,drop=FALSE]
    for (j in seq_along(keep)) {
      anchor <- which.max(abs(b$eigenvectors[,j]))
      b$eigenvectors[,j] <- b$eigenvectors[,j]*sign(b$eigenvectors[anchor,j])
    }
    colnames(b$eigenvectors) <- sprintf("mem%02d",seq_along(keep))
    if (!length(keep)) {
      b$status <- "unavailable"; b$reason <- "no positive Moran eigenvalues"
      b$training_matrix <- matrix(numeric(),n,0L)
    } else {
      b$column_transform <- fit_transform(b$eigenvectors)
      b$training_matrix <- apply_transform(b$eigenvectors,b$column_transform)
    }
  }
  b$column_names <- colnames(b$training_matrix); b
}
mem_extend_adjacency <- function(basis,w) {
  W <- basis$graph$W
  centered <- sweep(sweep(as.matrix(w),1,rowMeans(w),"-"),2,colMeans(W),"-")+mean(W)
  raw <- sweep(centered %*% basis$eigenvectors,2,basis$eigenvalues,"/")
  apply_transform(raw,basis$column_transform)
}
predict_basis <- function(basis,prediction) {
  if (basis$status!="available") stop(basis$reason)
  xy <- as.matrix(prediction$xy); n <- nrow(xy)
  if (ncol(xy)!=2L || any(!is.finite(xy))) stop("Invalid prediction coordinates")
  if (basis$arm=="env") return(matrix(numeric(),n,0L))
  if (basis$arm %in% c("gp","spatial2")) return(xy)
  if (basis$arm %in% c("trend","dnn")) {
    raw <- apply_transform(xy,basis$coordinate_transform)
    return(if (basis$arm=="trend") apply_transform(trend_columns(raw),basis$column_transform) else raw)
  }
  ids <- prediction$site_id
  if (length(ids)!=n || anyNA(ids)) stop("Prediction site IDs are required for stable neighbour ties")
  d <- xy_distances(xy,basis$training_xy); w <- matrix(0,n,nrow(basis$training_xy))
  for (i in seq_len(n)) {
    nearest <- order(d[i,],basis$training_ids,method="radix")[seq_len(basis$graph$k)]
    reverse <- d[i,]<basis$graph$radius | (d[i,]==basis$graph$radius & ids[i]<=basis$graph$boundary_id)
    w[i,union(nearest,which(reverse))] <- 1
  }
  result <- mem_extend_adjacency(basis,w)
  for (i in seq_len(n)) {
    matched <- which(d[i,]==0)
    if (length(matched)) result[i,] <- basis$training_matrix[matched[1],]
  }
  result
}
