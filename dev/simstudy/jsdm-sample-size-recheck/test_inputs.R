here <- 'work/jsdm-sample-size-20260919'
source(file.path(here, 'helpers.R'))
root <- 'work/nonspatial-bias-recheck-20260914'
design <- 'work/nonspatial-design-recheck-20260919'
for (r in 1:10) {
 a <- readRDS(file.path(root, 'operational-k6', sprintf('qnear_K6-data-%02d.rds', r)))
 b <- readRDS(file.path(design, 'inputs', sprintf('qnear_K6-sites300-%02d.rds', r)))
 far <- readRDS(file.path(root, 'operational-k6', sprintf('qfar_K6-data-%02d.rds', r)))
 stopifnot(identical(a$sim$true_params$z_true, far$sim$true_params$z_true),
           identical(a$sim$true_params$jsdmParams_true, far$sim$true_params$jsdmParams_true))
 inputs <- lapply(c(100L,300L,1000L), function(n) make_input(a,b,n))
 for (j in 1:3) {
  x <- inputs[[j]]; n <- c(100L,300L,1000L)[j]
  stopifnot(nrow(x$data$OTU)==n, ncol(x$data$OTU)==10, !anyDuplicated(x$data$info$Site),
            !any(c('Sample','Primer','X_theta') %in% names(x$data$info)),
            all(x$data$OTU %in% c(0,1)), identical(x$data$traits,a$sim$data_list$traits),
            identical(unname(x$data$OTU[1:100,]), unname(a$sim$true_params$z_true)),
            max(abs(x$truth$eta[1:100,]-a$sim$true_params$jsdmParams_true$eta)) < 1e-12,
            identical(x, make_input(a,b,n)))
  xx <- as.matrix(x$data$info[,x$covariates])
  rawB <- sweep(x$truth$B,1,apply(xx,2,sd),'/')
  rawB0 <- x$truth$B0-as.vector(colMeans(xx)%*%rawB)
  if(j==1) { B <- rawB; B0 <- rawB0 } else {
   stopifnot(max(abs(rawB-B))<1e-12,max(abs(rawB0-B0))<1e-12)
  }
 }
 stopifnot(identical(inputs[[2]]$data$OTU,inputs[[3]]$data$OTU[1:300,]),
           identical(inputs[[2]]$data$info,inputs[[3]]$data$info[1:300,]),
           max(abs(inputs[[2]]$truth$eta-inputs[[3]]$truth$eta[1:300,]))<1e-12)
}
cat('PASS: ten nested communities; exact z, raw relationships, traits, site identity and probabilities preserved.\n')
