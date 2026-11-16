# Check that fits whose sites were already in sorted order are unchanged by the fix.
# Usage: Rscript --vanilla sorted-fits-identity.R <library> <output.rds>
# Run once against a library holding main and once against one holding the fix,
# then compare: identical(readRDS("main.rds"), readRDS("fix.rds")).
args <- commandArgs(TRUE); lib <- args[1]; out <- args[2]
suppressMessages(library(occJSDM, lib.loc = lib))
mcmc <- list(nchain = 2, nburn = 30, niter = 30, nthin = 1)
sim <- function(model, spatial, seed = 1L, n = 40L, S = 4L) {
  set.seed(seed)
  M <- if (model %in% c("occupancy", "two_stage")) rep(2L, n) else rep(1L, n)
  N <- sum(M); K <- if (model == "two_stage") rep(2L, N) else rep(1L, N)
  simulateOccJSDMData(
    list(n = n, S = S, g = 2L, M = M, P = 1L, K = K, ncov_psi = 2L, ncov_theta = 1L),
    list(p = matrix(runif(S, .3, .6), 1, S), q = matrix(runif(S, .01, .05), 1, S),
         theta0 = runif(S, .02, .1), theta_baseline = runif(S, .2, .5)),
    list(gt = 2L, d = 2L, ds = 0, sigma_b = .5, sigma_bs = .5, sigma_ts = .5, sigma_h = 1,
         sigma_s = .5, l_s = .3, tau = rep(1, S), useSpatField = spatial),
    model = model)$data_list
}
fit <- function(data, occ, coll = NULL, spat = NULL, params = list(n_factors = 2)) {
  set.seed(99)
  f <- suppressMessages(suppressWarnings(runOccJSDM(data, occCovariates = occ, collCovariates = coll,
    spatCovariates = spat, listParams = params, MCMCparams = mcmc)))
  f$infos$execution <- NULL   # process ids, timings, library paths
  f
}
fits <- list()
for (model in c("binary", "continuous", "occupancy", "two_stage")) for (spatial in c(FALSE, TRUE)) {
  d <- sim(model, spatial)
  occ <- grep("^X_psi", names(d$info), value = TRUE)
  coll <- if (model %in% c("occupancy", "two_stage")) grep("^X_theta", names(d$info), value = TRUE) else NULL
  spat <- if (spatial) c("Xs.1", "Xs.2") else NULL
  params <- list(n_factors = 2); if (spatial) params$n_supportpoints <- 8
  nm <- paste(model, if (spatial) "spatial" else "nonspatial")
  fits[[nm]] <- fit(d, occ, coll, spat, params)
  if (model %in% c("binary", "continuous")) {
    # the same data with an explicit Site column: integer, sorted character, sorted factor
    for (kind in c("integer", "character", "factor")) {
      d2 <- d; n <- nrow(d2$info)
      d2$info$Site <- switch(kind, integer = seq_len(n), character = sprintf("P%03d", seq_len(n)),
                             factor = factor(sprintf("P%03d", seq_len(n))))
      fits[[paste(nm, "Site", kind)]] <- fit(d2, occ, coll, spat, params)
    }
  }
}
# The shipped example data (two-stage, integer sites in sorted order), short chains.
e <- new.env(); data("sampledata", package = "occJSDM", lib.loc = lib, envir = e)
data("quickstartdata", package = "occJSDM", lib.loc = lib, envir = e)
sd <- e$sampledata; qd <- e$quickstartdata
fits[["sampledata"]] <- fit(sd, grep("^X_psi", names(sd$info), value = TRUE),
                            grep("^X_theta", names(sd$info), value = TRUE))
qocc <- setdiff(names(qd$info), c("Site", "Sample", "Primer"))
fits[["quickstartdata"]] <- tryCatch(fit(qd, qocc[1:2]), error = function(e) conditionMessage(e))
saveRDS(fits, out)
cat(length(fits), "fits saved\n")
