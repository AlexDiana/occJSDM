args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "pilot-math.R"))
r <- normal_rule(31L)
stopifnot(abs(sum(r$weights) - 1) < 1e-12,
          abs(sum(r$weights * r$nodes^2) - 1) < 1e-12,
          max(abs(logistic_normal(c(-1, 0, 1), 0) - plogis(c(-1, 0, 1)))) < 1e-12,
          abs(logistic_normal(0, 2) - 0.5) < 1e-12)
for (mu in c(-3, 0.2, 4)) for (sd in c(.1, .8, 2)) {
  reference <- integrate(function(z) plogis(mu + sd*z) * dnorm(z),
                          -Inf, Inf, rel.tol = 1e-11)$value
  stopifnot(abs(logistic_normal_adaptive(mu, sd) - reference) < 1e-8,
            abs(logistic_normal(mu, sd, 121L) - reference) < 1e-8)
}
stopifnot(abs(logistic_normal_adaptive(10000, 100000) - pnorm(.1)) < 1e-7)
x <- data.frame(environment_1 = c(-1, 1), environment_2 = c(0, 0))
b <- matrix(c(-.3, .8, -.2, .7, -.4, .3), nrow = 3)
y <- matrix(c(1, 0, 0, 1), 2)
p <- list(beta = b, loading = matrix(0, 2, 2))
a <- joint_integration(p, x, y)
expected <- plogis(cbind(1, as.matrix(x)) %*% b)
stopifnot(max(abs(a$conditional - expected)) < 1e-12,
          abs(a$loglik - sum(dbinom(y, 1, expected, log = TRUE))) < 1e-12)
p$loading <- matrix(c(.7, -.3, .2, .5), 2)
a <- joint_integration(p, x, y, 31L)
better <- joint_integration(p, x, y, 61L)
rotation <- matrix(c(0, 1, -1, 0), 2)
p$loading <- rotation %*% p$loading
rotated <- joint_integration(p, x, y, 61L)
stopifnot(max(abs(a$conditional - better$conditional)) < 1e-8,
          abs(a$loglik - better$loglik) < 1e-8,
          abs(rotated$loglik - better$loglik) < 1e-8)
cat("Probability identities, independent integration, zero-factor likelihood and rotation checks passed.\n")
