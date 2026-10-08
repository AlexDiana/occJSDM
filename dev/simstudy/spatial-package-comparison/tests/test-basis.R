library(testthat)
if (file.exists("../basis.R")) source("../basis.R")
fixture <- function(n=36L) {
  xy <- as.matrix(expand.grid(x=seq(0,1,length.out=6),y=seq(0,1,length.out=6)))[seq_len(n),,drop=FALSE]
  list(schema_version=1L,dataset_id="fixture",site_id=sprintf("s%03d",seq_len(n)),
       xy=xy,environment=matrix(seq_len(n),ncol=1),y=matrix(0,n,2))
}
test_that("trend transformations are trained once and survive serialisation", {
  tr <- fixture(); b <- fit_basis(tr,"trend")
  expect_identical(colnames(b$training_matrix),c("x","y","x2","xy","y2"))
  expect_equal(unname(colMeans(b$training_matrix)),rep(0,5),tolerance=1e-12)
  expect_equal(unname(apply(b$training_matrix,2,sd)),rep(1,5),tolerance=1e-12)
  p <- tempfile(); saveRDS(b,p); b2 <- readRDS(p)
  expect_identical(predict_basis(b2,tr),b$training_matrix)
  pred <- tr; pred$xy <- pred$xy+4
  expect_identical(predict_basis(b,tr),b$training_matrix)
  expect_false(isTRUE(all.equal(predict_basis(b,pred),b$training_matrix)))
  expect_equal(predict_basis(b,pred)[1,,drop=FALSE],predict_basis(b,lapply(pred,function(x) if(is.matrix(x)) x[1,,drop=FALSE] else if(length(x)==36) x[1] else x)))
})
test_that("MEM uses no outcomes and reproduces training rows without its cache", {
  tr <- fixture(); b <- fit_basis(tr,"mem"); tr2 <- tr; tr2$y[,] <- 1
  expect_identical(fit_basis(tr2,"mem"),b)
  expect_equal(predict_basis(b,tr),b$training_matrix,tolerance=1e-10)
  expect_equal(mem_extend_adjacency(b,b$graph$W),b$training_matrix,tolerance=1e-10)
  expect_equal(b$graph$W,t(b$graph$W))
  expect_equal(diag(b$graph$W),rep(0,36))
  expect_equal(colMeans(b$training_matrix),rep(0,ncol(b$training_matrix)),ignore_attr=TRUE,tolerance=1e-12)
  expect_lte(ncol(b$training_matrix),10L)
  H <- diag(36)-matrix(1/36,36,36)
  ev <- eigen(H %*% b$graph$W %*% H,symmetric=TRUE)
  # Eigenspace projection comparison is invariant to sign and repeated-root rotations.
  V <- b$eigenvectors; ref <- ev$vectors[,seq_len(ncol(V)),drop=FALSE]
  expect_equal(V %*% t(V),ref %*% t(ref),tolerance=1e-10)
  expect_true(all(b$eigenvalues > b$positive_tolerance))
})
test_that("ties use stable IDs and prediction permutations only reorder rows", {
  tr <- fixture(); b <- fit_basis(tr,"mem")
  ids <- rev(seq_len(nrow(tr$xy))); tr2 <- tr
  tr2$xy <- tr$xy[ids,,drop=FALSE]; tr2$site_id <- tr$site_id[ids]
  b2 <- fit_basis(tr2,"mem")
  expect_equal(b2$graph$W[order(ids),order(ids)],b$graph$W)
  pred <- list(site_id=c("newA","newB","newC"),xy=rbind(c(.3,.3),c(100,100),c(0,0)))
  p <- predict_basis(b,pred)
  expect_true(all(is.finite(p)))
  pred$xy <- pred$xy[3:1,,drop=FALSE]; pred$site_id <- pred$site_id[3:1]
  expect_equal(predict_basis(b,pred),p[3:1,,drop=FALSE])
})
test_that("Nyström cross-row centring matches a direct calculation", {
  b <- fit_basis(fixture(),"mem"); w <- matrix(c(rep(1,8),rep(0,28)),1)
  c <- w-mean(w)-matrix(colMeans(b$graph$W),1)+mean(b$graph$W)
  raw <- sweep(c %*% b$eigenvectors,2,b$eigenvalues,"/")
  expected <- sweep(sweep(raw,2,b$column_transform$center,"-"),2,b$column_transform$scale,"/")
  expect_equal(unname(mem_extend_adjacency(b,w)),unname(expected),tolerance=1e-12)
})
test_that("disconnected graphs and scarce positive eigenvalues remain explicit", {
  tr <- fixture(20); tr$xy[11:20,] <- tr$xy[11:20,]+100
  b <- fit_basis(tr,"mem")
  expect_equal(b$graph$components,2L)
  expect_lt(ncol(b$training_matrix),10L)
  small <- fixture(8); b0 <- fit_basis(small,"mem")
  expect_identical(b0$status,"unavailable")
  expect_identical(ncol(b0$training_matrix),0L)
  expect_error(predict_basis(b0,small),"no positive")
})
test_that("environment and DNN use training-only transforms and native arms raw coordinates", {
  tr <- fixture(); b <- fit_basis(tr,"dnn")
  expect_equal(unname(colMeans(b$training_matrix)),c(0,0),tolerance=1e-12)
  expect_equal(unname(apply(b$training_matrix,2,sd)),c(1,1))
  expect_equal(dim(fit_basis(tr,"env")$training_matrix),c(36L,0L))
  expect_identical(fit_basis(tr,"gp")$training_matrix,tr$xy)
  expect_identical(fit_basis(tr,"spatial2")$training_matrix,tr$xy)
})
