args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "pilot-math.R"))
source(file.path(args[1], "sjsdm-gradient.R"))
x <- data.frame(environment_1 = c(-1, .2, 1), environment_2 = c(.4, -.2, .8))
y <- matrix(c(1, 0, 1, 0, 1, 0), 3, 2)
p <- list(beta = matrix(c(.2, -.3, .6, -.4, .8, -.7), 3, 2),
          loading = matrix(c(.5, -.2, .3, .7), 2, 2))
f <- make_joint_gradient(x, y, nodes = 31L)
theta <- c(p$beta, p$loading)
a <- f(theta)
reference <- joint_integration(p, x, y, 61L, clamp = TRUE)
stopifnot(abs(a$value + reference$loglik) < 1e-9)
finite <- vapply(seq_along(theta), function(k) {
  delta <- theta * 0
  delta[k] <- 1e-5
  (f(theta + delta)$value - f(theta - delta)$value) / 2e-5
}, numeric(1))
stopifnot(max(abs(a$gradient - finite)) < 1e-7)
# Check independent species with no latent factor effect against a logistic GLM score.
p$loading[,] <- 0
a <- f(c(p$beta, p$loading))
eta <- cbind(1, as.matrix(x)) %*% p$beta
v <- plogis(eta)
prob <- .999999 * v + .0000005
score <- (y / prob - (1-y) / (1-prob)) * .999999 * v * (1-v)
stopifnot(max(abs(a$gradient[1:6] + as.vector(crossprod(cbind(1, as.matrix(x)), score)))) < 1e-10,
          max(abs(a$gradient[7:10])) < 1e-10)
cat("Checked likelihood, all gradient components, and zero-factor limit.\n")
