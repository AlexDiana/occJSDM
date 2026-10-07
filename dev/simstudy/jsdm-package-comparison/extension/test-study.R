# Catch broken pairing, leakage through preprocessing, and mislabelled truth.
code <- commandArgs(trailingOnly = TRUE)[1]
if (file.exists(file.path(code, "study.R"))) source(file.path(code, "study.R"))
stopifnot(exists("simulate_community", mode = "function"))
base <- simulate_community(1, "baseline")
a <- make_input(base, 100, 10, "linear", FALSE)
b <- make_input(base, 300, 10, "linear", FALSE)
stopifnot(identical(a$y, b$y[1:100, , drop = FALSE]),
          identical(rownames(a$test_x), rownames(b$test_x)),
          !any(rownames(a$x) %in% rownames(a$test_x)),
          max(abs(colMeans(a$x))) < 1e-12,
          max(abs(apply(a$x, 2, sd) - 1)) < 1e-12)
# Altering held-out conditions cannot alter fitted inputs or their scaling.
changed <- base
changed$raw_x[301:600, ] <- changed$raw_x[301:600, ] + 100
c <- make_input(changed, 100, 10, "linear", FALSE)
stopifnot(identical(a$x, c$x), identical(a$y, c$y), identical(a$centre, c$centre))
rare <- simulate_community(1, "rare")
stopifnot(identical(base$raw_x, rare$raw_x), identical(base$hidden, rare$hidden),
          identical(base$probability[, 4:10], rare$probability[, 4:10]),
          all(rare$probability[, 1:3] < base$probability[, 1:3]))
curved <- simulate_community(1, "curved")
cx <- make_input(curved, 100, 10, "quadratic", FALSE)
beta <- scaled_truth(curved, cx)
eta <- cbind(1, as.matrix(cx$x)) %*% beta + curved$hidden[1:100, ] %*% curved$loading
stopifnot(max(abs(plogis(eta) - curved$probability[1:100, ])) < 1e-12,
          max(abs(curved$probability - base$probability)) > .1,
          identical(curved$uniforms, base$uniforms))
traits <- simulate_community(1, "traits")
t10 <- make_input(traits, 100, 10, "linear", TRUE)
t30 <- make_input(traits, 100, 30, "linear", TRUE)
stopifnot(identical(t10$y, t30$y[, 1:10]), identical(t10$traits, t30$traits[1:10, ]),
          all(traits$trait_effect[2, ] == 0),
          any(traits$trait_effect[1, ] != 0),
          !identical(simulate_community(1, "baseline")$y, simulate_community(2, "baseline")$y))
# Two independent signed errors can cancel; absolute errors must not.
m <- error_summary(c(.2, .8), c(.3, .7))
stopifnot(abs(m$bias_pp) < 1e-10, abs(m$mae_pp - 10) < 1e-10)
# Independent integration identity: no hidden variation gives inverse-link truth.
tiny <- base
tiny$loading[] <- 0
truth <- marginal_truth(tiny, tiny$raw_x[1:3, ])
expected <- plogis(cbind(1, tiny$raw_x[1:3, ]) %*% tiny$beta)
stopifnot(max(abs(truth - expected)) < 1e-10)
cat("Paired sites/species, independent replicates, no test leakage, curved truth, null trait and error summaries passed.\n")
