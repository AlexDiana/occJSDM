args <- commandArgs(trailingOnly = TRUE)
code <- if (length(args)) args[1] else "dev/simstudy/jsdm-package-comparison/extension"
suppressPackageStartupMessages(library(sjSDM))
sjSDM:::check_module()
reticulate::source_python(file.path(code, "native-sjsdm.py"))
model <- sjSDM:::pkg.env$fa$Model_sjSDM(device = "cpu", dtype = "float64", seed = 18L)
model$add_env(2L, 2L, intercept = TRUE, l1 = 0, l2 = 0)
model$build(df = 2L, link = "logit", alpha = 1, scheduler = FALSE)
b <- matrix(c(.3, -.2, -.4, .5), 2, 2)
model$set_env_weights(list(t(b))); model$set_sigma(matrix(0, 2, 2))
x <- cbind(1, seq(-2, 2, length.out = 25))
y <- cbind(rep(c(0, 1), length.out = 25), rep(c(1, 0, 1), length.out = 25))
z <- native_covariance(model, x, y, 100L, 123L)
for (s in 1:2) {
  # Zero loadings reduce to independent logistic regressions. sjSDM's tiny
  # probability shrinkage changes the exact Hessian by less than this tolerance.
  p <- plogis(x %*% b[,s])
  h <- crossprod(x, as.vector(p*(1-p))*x)
  stopifnot(max(abs(z$regularized_hessian[s,,] - (h + .001))) < 1e-4)
}
reticulate::py$torch$manual_seed(123L)
direct <- sjSDM:::force_r(model$se(x, y, batch_size=25L, sampling=100L, parallel=0L, verbose=FALSE))
stopifnot(max(abs(do.call(rbind, direct) - z$se)) < 1e-12)
cat("Captured native sjSDM covariance matches native SEs and independent logistic information.\n")
