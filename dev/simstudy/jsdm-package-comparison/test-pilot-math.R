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
# At very large SD, the near-step transition can sit between quadrature nodes.
# Its large-SD limit is accurate here to far better than the asserted tolerance.
stopifnot(max(abs(logistic_normal_adaptive(c(-100, -30, 30, 100), 100000) -
                   pnorm(c(-100, -30, 30, 100) / 100000))) < 1e-9,
          max(abs(logistic_normal_adaptive(c(-1e6, 1e6), .1) - c(0, 1))) < 1e-12)
# Splitting the normal integral at a very sharp transition alone can also fail.
# These references were independently checked with the probability identity.
stopifnot(abs(logistic_normal_adaptive(5000, 10000) - 0.691462458378) < 1e-10,
          max(abs(logistic_normal_adaptive(c(-5000, 5000), 1000) -
                    c(2.86663799967e-7, 1 - 2.86663799967e-7))) < 1e-10)
# Job 39: a rare prediction in the curved sjSDM fit caused the probability-
# domain integral to report divergence. Reference checked with 121-node normal
# quadrature; the reflected case checks the complementary probability as well.
rare_mu <- -11.608873808283207
rare_sd <- 0.97504942700673214
rare_expected <- 1.4613702701051533e-5
stopifnot(max(abs(logistic_normal_adaptive(c(rare_mu, -rare_mu), rare_sd) -
                   c(rare_expected, 1 - rare_expected))) < 1e-12)
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
