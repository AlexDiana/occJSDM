args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
source(file.path(args[2], "pilot-math.R"))
input <- readRDS(file.path(root, "inputs/training.rds"))
records <- list.files(file.path(root, "fits"), "record[.]rds$", full.names = TRUE)
hash <- unname(tools::md5sum(file.path(root, "inputs/training.rds")))
stopifnot(all(vapply(records, function(f) identical(readRDS(f)$input_hash, hash), logical(1))))

suppressPackageStartupMessages(library(gllvm))
path <- file.path(root, "fits/gllvm-VA-attempt-1-start-1.rds")
g <- readRDS(path)
p <- point_parameters(g, "gllvm")
stopifnot(identical(g$y, input$y), all(g$link == "logit"), g$method == "VA",
          g$num.lv == 2, g$num.lv.c == 0, g$num.RR == 0,
          max(abs(tcrossprod(t(p$loading)) - getResidualCov(g, adjust = 0)$cov)) < 1e-9,
          max(abs(plogis(cbind(1, as.matrix(input$x)) %*% p$beta) -
                    predict(g, newX = input$x, level = 0, type = "response"))) < 1e-9)

suppressPackageStartupMessages(library(Hmsc))
f <- readRDS(file.path(root, "fits/Hmsc-attempt-1-start-1.rds"))
p <- bayesian_parameters(f, "Hmsc")
source(file.path(Sys.getenv("LESSON_N_ROOT"), "hmsc-serial.R"))
d <- f$postList[[1]][[1]]
set.seed(26092299)
native <- predict_hmsc_serial(f, post = list(d, d), XData = input$x, expected = TRUE)[[1]]
manual <- pnorm(cbind(1, as.matrix(input$x)) %*% p$beta[, , 1, 1] +
                 p$score[, , 1, 1] %*% p$loading[, , 1, 1])
stopifnot(max(abs(native - manual)) < 1e-10)
# An independent numerical check of the probit-normal identity.
for (mu in c(-2, 0.4, 3)) for (sd in c(.2, .8, 2)) {
  numeric <- integrate(function(z) pnorm(mu + sd*z)*dnorm(z), -Inf, Inf)$value
  stopifnot(abs(numeric - pnorm(mu / sqrt(1+sd^2))) < 1e-8)
}

f <- readRDS(file.path(root, "fits/occJSDM-attempt-1-start-1.rds"))
p <- bayesian_parameters(f, "occJSDM")
j <- f$results_output$jsdm_output
stopifnot(identical(f$infos$OTU, input$y),
          identical(f$infos$speciesNames, colnames(input$y)),
          dim(j$A_output)[2] == 0, f$infos$ps == 0)
manual <- plogis(sweep(f$X_psi %*% j$B_output[, , 1, 1],
                        2, j$B0_output[, 1, 1], "+") +
                    j$U_output[, , 1, 1] %*% j$L_output[, , 1, 1])
reconstructed <- plogis(cbind(1, as.matrix(input$x)) %*% p$beta[, , 1, 1] +
                         p$score[, , 1, 1] %*% p$loading[, , 1, 1])
stopifnot(max(abs(manual - reconstructed)) < 1e-12)

cat("All ordered training hashes match; gllvm covariance and native zero-factor output, Hmsc paired-draw native output, probit integration and occJSDM paired-draw reconstruction passed.\n")
