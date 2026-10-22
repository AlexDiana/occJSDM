args <- commandArgs(trailingOnly = TRUE)
code <- if (length(args)) args[1] else "dev/simstudy/jsdm-package-comparison/extension"
source(file.path(code, "calibration-math.R"))
if (file.exists(file.path(code, "native-interval-math.R"))) source(file.path(code, "native-interval-math.R"))
stopifnot(exists("native_coefficient_intervals", mode = "function"))
# A correlated intercept/slope needs their covariance when undoing centring.
b <- matrix(c(10, 4, 6), 3, 1)
V <- matrix(c(4, 1, 0, 1, 9, 2, 0, 2, 16), 3)
z <- native_coefficient_intervals(b, list(V), c(1, 2), c(2, 3))
a <- c(1, -.5, -2/3)
expected <- c(as.numeric(t(a) %*% V %*% a), 9/4, 16/9)
stopifnot(max(abs(z$estimate - c(4, 2, 2))) < 1e-12,
  max(abs(z$se^2 - expected)) < 1e-12,
  max(abs(z$upper - z$estimate - qnorm(.975)*sqrt(expected))) < 1e-12)
# An indefinite covariance can have positive diagonals; never repair it silently.
bad <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), 3)
z <- native_coefficient_intervals(b, list(bad), c(1, 2), c(2, 3))
stopifnot(all(is.na(z$lower)), all(is.na(z$upper)), all(is.finite(z$estimate)))
cat("Native interval covariance transformation and invalid covariance checks passed.\n")
