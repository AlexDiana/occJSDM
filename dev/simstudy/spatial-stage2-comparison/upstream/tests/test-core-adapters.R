library(testthat)
source("../../config.R"); for (f in c("basis.R","adapter-contract.R","adapter-occjsdm.R")) source(file.path("..",f))
small_training <- function(n=20L) {
  set.seed(341)
  list(schema_version=1L,dataset_id="synthetic",site_id=sprintf("s%03d",1:n),species_id=c("absent","present"),
       y=cbind(absent=rep(0,n),present=rbinom(n,1,.5)),environment=matrix(rnorm(n),ncol=1,dimnames=list(NULL,"environment")),
       xy=matrix(runif(n*2),n,2,dimnames=list(NULL,c("x","y"))),input_hashes=c(synthetic="fixture"))
}
test_that("fit inputs cannot carry truth or unordered duplicate IDs", {
  tr <- small_training(); expect_silent(validate_training(tr))
  tr$psi <- tr$y; expect_error(validate_training(tr),"Unexpected")
  tr$psi <- NULL; tr$site_id[2] <- tr$site_id[1]; expect_error(validate_training(tr),"IDs")
})
test_that("gauge centring preserves full probabilities and retains structural zero space", {
  E <- matrix(c(1,2,3,4),2); S <- E/2; B <- c(-1,1)
  terms <- center_terms(B,E,S,colMeans(E),colMeans(S))
  expect_equal(terms$intercept+terms$environment+terms$spatial,sweep(E+S,2,B,"+"))
  rec <- point_prediction(list(site_id=c("a","b"),species_id=c("u","v")),terms,plogis,"logit", "fixture")
  expect_silent(validate_prediction(rec,list(site_id=c("a","b"),species_id=c("u","v"))))
  expect_equal(rec$delta_p,rec$p_full-rec$p_space_zero)
  expect_identical(rec$uncertainty$status,"unavailable")
  expect_null(rec$intervals)
  zero <- center_terms(B,E,S*0,colMeans(E),c(0,0))
  no <- point_prediction(list(site_id=c("a","b"),species_id=c("u","v")),zero,plogis,"logit","fixture")
  expect_equal(no$p_full,no$p_space_zero)
  bad <- rec; bad$species_id <- rev(bad$species_id)
  expect_error(validate_prediction(bad,list(site_id=c("a","b"),species_id=c("u","v"))),"order")
})
test_that("native adapters preserve zeros, predictor order, and fresh-process prediction", {
  launcher <- normalizePath("../../run-r"); code <- normalizePath("..")
  root <- Sys.getenv("SPATIAL_SMOKE_ROOT",unset=tempfile()); dir.create(root,recursive=TRUE,showWarnings=FALSE); saveRDS(small_training(),file.path(root,"training.rds"))
  fixture <- file.path(code,"tests/native-core-fixture.R")
  expect_true(file.exists(fixture))
  if (!file.exists(fixture)) return(invisible(NULL))
  for (package in c("occJSDM")) {
    out <- suppressWarnings(system2(launcher,c("-f",shQuote(fixture),"--args",shQuote(code),shQuote(root),package,"fit"),stdout=TRUE,stderr=TRUE))
    expect_null(attr(out,"status"),info=paste(tail(out,20),collapse="\n"))
    out <- suppressWarnings(system2(launcher,c("-f",shQuote(fixture),"--args",shQuote(code),shQuote(root),package,"predict"),stdout=TRUE,stderr=TRUE))
    expect_null(attr(out,"status"),info=paste(tail(out,20),collapse="\n"))
  }
})
