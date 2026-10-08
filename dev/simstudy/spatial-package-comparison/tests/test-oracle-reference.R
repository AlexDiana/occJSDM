library(testthat)
for(f in c("config.R","archive.R","basis.R","adapter-contract.R","adapter-gllvm.R","oracle-reference.R")) if(file.exists(file.path("..",f))) source(file.path("..",f))
test_that("oracle prediction uses full conditional variance and the common gauge", {
  tr <- list(schema_version=1L,dataset_id="oracle-fixture",site_id=c("a","b"),species_id="s",xy=rbind(c(0,0),c(.5,.5)),environment=matrix(c(-1,1),2),y=matrix(c(0,1),2,dimnames=list(NULL,"s")),input_hashes=c(fixture="a"))
  pred <- tr; pred$y <- NULL
  draws <- array(c(-1,.7,-.8,.6,-1.2,.8,-1,.5),c(2,2,2))
  f <- list(package="oracle",arm="oracle",site_id=tr$site_id,species_id=tr$species_id,
    training_xy=tr$xy,training_environment=tr$environment,native=list(draws=draws,B0=-1,B=.4,range=.1),settings=list(prediction_seed=8L))
  ex <- oracle_draw_extractor(f,pred,NULL); d <- ex(1:2)
  expect_equal(d$p_full[,1,,],plogis(draws+c(-1.4,-.6)),tolerance=1e-12)
  expect_lt(max(abs(apply(d$field,c(2,3,4),sum))),1e-12)
  expect_gt(max(abs(apply(d$p_full,c(1,2),mean)-plogis(apply(draws,1,mean)+c(-1.4,-.6)))),1e-5)
  far <- subset_prediction(pred,1); far$xy[,] <- 10; far$site_id <- "far"
  moments <- oracle_conditional(f,far)
  expect_equal(moments$A,matrix(0,1,2)); expect_gt(moments$covariance[1,1],.99)
  expect_equal(moments$covariance,matrix(1+1e-8,1,1),tolerance=1e-12)
  p <- predict_oracle(f,pred,NULL); expect_equal(p$p_full,apply(d$p_full,c(1,2),mean))
})
test_that("oracle layer reproduces original lattice probability means for the same draws and RNG", {
  env <- new.env(parent=globalenv()); source("../../spatial-design-sweep/generator.R",local=env); source("../../spatial-design-sweep/oracle.R",local=env)
  set.seed(842); trxy <- matrix(runif(12),6,2); xy <- matrix(runif(8),4,2); e <- matrix(rnorm(4),4,1)
  f <- list(package="oracle",arm="oracle",site_id=paste0("t",1:6),species_id="s",training_xy=trxy,
    training_environment=matrix(rnorm(6),6,1),native=list(draws=array(rnorm(6*100*4),c(6,100,4)),B0=-1,B=.4,range=.1),settings=list(prediction_seed=143L))
  p <- list(site_id=paste0("n",1:4),species_id="s",xy=xy,environment=e)
  old <- env$oracle_lattice(f$native$draws,trxy,xy,.1,as.vector(-1+.4*e),thin=1L,seed=143L)
  current <- predict_oracle(f,p)
  expect_equal(as.vector(current$p_full),old$probability_mean,tolerance=1e-10)
})
