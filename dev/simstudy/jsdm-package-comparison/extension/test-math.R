args <- commandArgs(trailingOnly = TRUE)
code <- args[1]
source(file.path(dirname(code), "pilot-math.R"))
if (file.exists(file.path(code, "math.R"))) source(file.path(code, "math.R"))
stopifnot(exists("posterior_marginal", mode = "function"))
# Any hard-coded original dimensions (100 sites, 10 species, 2 slopes) break this.
p <- list(beta = array(0, c(4, 3, 5, 2)), loading = array(0, c(2, 3, 5, 2)), link = "logit")
p$beta[1, , , ] <- qlogis(.2)
p$beta[4, , , ] <- 1
x <- data.frame(a = c(0, 0), b = c(0, 0), square = c(0, 1))
pred <- posterior_marginal(p, x, retain_draws = TRUE)
stopifnot(identical(dim(pred), c(2L, 3L, 5L, 2L)),
          max(abs(pred[1, , , ] - .2)) < 1e-12,
          max(abs(pred[2, , , ] - plogis(qlogis(.2) + 1))) < 1e-12)
p$link <- "probit"; p$beta[] <- 0; p$beta[1, , , ] <- 1
p$loading[1, , , ] <- 2
pred <- posterior_marginal(p, x)
stopifnot(max(abs(pred - pnorm(1 / sqrt(5)))) < 1e-12)
cat("General dimensions and logit/probit marginal probability identities passed.\n")
