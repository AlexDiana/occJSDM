library(testthat)
source('../spatial-targeted-recheck/score.R')
source('metrics.R')

test_that('field reconstruction respects species, iterations, ranges and centring', {
  truth <- matrix(c(-1,0,2,3,1,-2),3,2)
  index <- matrix(rep(c(1L,2L),8),8,2)
  coefficients <- array(0,c(3,2,8,2))
  for(ch in 1:2)for(it in 1:8)
    coefficients[,,it,ch] <- truth/index[it,ch]
  result <- summarise_field_draws(coefficients,index,list(diag(3),2*diag(3)),truth,
                                  matrix(1,8,2))
  expect_equal(result$field_mean,truth)
  expect_equal(result$score$raw_rmse,0,tolerance=1e-12)
  expect_equal(result$score$centred_rmse,0,tolerance=1e-12)
  expect_equal(result$score$centred_correlation,1,tolerance=1e-12)
  expect_equal(result$score$centred_slope,1,tolerance=1e-12)
  for(ch in 1:2)for(it in 1:8)
    coefficients[,,it,ch] <- sweep(truth,2,c(4,-3),'+')/index[it,ch]
  offset <- summarise_field_draws(coefficients,index,list(diag(3),2*diag(3)),truth,
                                  matrix(1,8,2))
  expect_equal(offset$score$raw_rmse,sqrt(12.5))
  expect_equal(offset$score$centred_rmse,0,tolerance=1e-12)
  expect_equal(offset$score$centred_correlation,1,tolerance=1e-12)
  expect_equal(as.numeric(offset$traces['field_mean_1',,]),rep(mean(truth[,1])+4,16))
  expect_equal(as.numeric(offset$traces['field_projection_2',,]),rep(1,16))
})

test_that('amplitude and field mixing checks cannot pass low ESS or missing diagnostics', {
  d <- data.frame(quantity=c('amplitude','field_mean_1'),rhat=c(1.01,1.02),ess_mean=c(200,300))
  expect_length(spatial_flags(d),0)
  d$ess_mean[1] <- 50
  expect_match(paste(spatial_flags(d),collapse=';'),'amplitude.*ESS')
  d$ess_mean[1] <- 200;d$rhat[2] <- 1.1
  expect_match(paste(spatial_flags(d),collapse=';'),'field_mean_1.*Rhat')
  d$rhat[2] <- NA
  expect_match(paste(spatial_flags(d),collapse=';'),'unavailable')
})

test_that('start sensitivity changes exactly the amplitude initializer', {
  f <- function(nchain=4) {
    x <- 8;out <- numeric(nchain)
    for(chain in 1:nchain) {sigma_bs <- .001;out[chain] <- sigma_bs+x}
    out
  }
  altered <- with_spatial_starts(f,c(.1,.3,1,3))
  expect_equal(f(),rep(8.001,4))
  expect_equal(altered(),8+c(.1,.3,1,3))
  expect_identical(environment(altered),environment(f))
  expect_error(with_spatial_starts(function() 1,1),'exactly one')
})

test_that('the real fitting function can be cloned without changing its other expressions', {
  e <- new.env();sys.source('../../../R/runOccJSDM.R',e)
  original <- e$runOccJSDM
  altered <- with_spatial_starts(original,c(.1,.3,1,3))
  expect_identical(formals(altered),formals(original))
  expect_identical(environment(altered),environment(original))
  # Restoring the sole initializer recovers the complete original body.
  restore <- function(x) {
    if(!is.call(x))return(x)
    if(identical(x[[1]],as.name('<-')) && identical(x[[2]],as.name('sigma_bs'))) {
      x[[3]] <- quote(if (identical(spatial_sd_prior$type,"half_cauchy")) spatial_sd_prior$scale else .001)
      return(x)
    }
    for(i in seq_along(x)[-1L])if(!identical(x[[i]],quote(expr=)))x[i] <- list(restore(x[[i]]))
    x
  }
  expect_identical(restore(body(altered)),body(original))
})

test_that('varying field draws agree with a direct iteration-wise oracle', {
  set.seed(948)
  b <- array(rnorm(2*3*12*2),c(2,3,12,2));index <- matrix(sample(1:2,24,TRUE),12,2)
  bases <- list(matrix(rnorm(10),5,2),matrix(rnorm(10),5,2))
  truth <- matrix(rnorm(15),5,3);amplitude <- matrix(runif(24),12,2)
  fields <- array(0,c(5,3,12,2));projection <- array(0,c(3,12,2))
  tc <- sweep(truth,2,colMeans(truth),'-')
  for(ch in 1:2)for(it in 1:12) {
    fields[,,it,ch] <- bases[[index[it,ch]]]%*%b[,,it,ch]
    projection[,it,ch] <- colSums(fields[,,it,ch]*tc)/colSums(tc^2)
  }
  actual <- summarise_field_draws(b,index,bases,truth,amplitude)
  expect_equal(actual$field_mean,apply(fields,c(1,2),mean),tolerance=1e-12)
  expect_equal(unname(actual$traces[paste0('field_mean_',1:3),,]),apply(fields,c(2,3,4),mean))
  expect_equal(unname(actual$traces[paste0('field_rms_',1:3),,]),sqrt(apply(fields^2,c(2,3,4),mean)))
  expect_equal(unname(actual$traces[paste0('field_projection_',1:3),,]),projection,tolerance=1e-12)
})
