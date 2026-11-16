# Reproduction of the site-order defect in one-row-per-site fits (16 November 2026).
# Usage: Rscript --vanilla reproduce.R <library holding the occJSDM version to test>
# 150 sites, 6 species, true slope 1.5 on one covariate, one chain of 300 + 300.
# A version with the defect gives a slope near 0 for every label set that is
# not in sorted order; a fixed version gives the same slope for all of them.
lib <- commandArgs(TRUE)[1]
suppressMessages(library(occJSDM, lib.loc = lib))
cat("occJSDM", as.character(packageVersion("occJSDM", lib.loc = lib)), "dplyr", as.character(packageVersion("dplyr")), "\n")
set.seed(1)
n <- 150; x <- rnorm(n)
Y <- sapply(1:6, function(s) rbinom(n, 1, plogis(-0.5 + 1.5 * x)))
colnames(Y) <- paste0("sp", 1:6)
run <- function(site) {
  set.seed(2)
  fit <- suppressMessages(runOccJSDM(
    data = list(info = data.frame(Site = site, x = x), OTU = Y),
    occCovariates = "x", listParams = list(n_factors = 0, n_lattrait = 0),
    MCMCparams = list(nchain = 1, nburn = 300, niter = 300, nthin = 1, ncores = 1)))
  c(slope = mean(fit$results_output$jsdm_output$B_output),
    cor_X_x = cor(fit$X_psi[, 1], x),
    X_is_scaled_x = isTRUE(all.equal(unname(fit$X_psi[, 1]), as.numeric(scale(x)))))
}
set.seed(3); shuf <- sample(sprintf("P%03d", seq_len(n)))
res <- rbind(
  integer = run(seq_len(n)),
  padded = run(sprintf("P%03d", seq_len(n))),
  FHT = run(paste0("FHT", seq_len(n))),
  shuffled = run(shuf),
  factor_FHT = run(factor(paste0("FHT", seq_len(n)))),
  factor_FHT_levels_in_order = run(factor(paste0("FHT", seq_len(n)), levels = paste0("FHT", seq_len(n)))),
  integer_reversed = run(rev(seq_len(n))))
print(round(res, 3))
