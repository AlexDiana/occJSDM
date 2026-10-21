# Time the quickstart chain under a development load and under an installed build.
# Usage: Rscript compare-load.R <loadall|installed> <nburn> <niter> [lib]
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; nburn <- as.integer(args[2]); niter <- as.integer(args[3])
cat("JIT level:", compiler::enableJIT(-1), "\n")
if (mode == "loadall") {
  suppressMessages(pkgload::load_all(Sys.getenv("OCCJSDM_PKG", unset = "."), quiet = TRUE))
} else {
  suppressMessages(library(occJSDM, lib.loc = args[4]))
}
f <- occJSDM:::spatial_range_logweights
cat("is bytecode:", grepl("bytecode", paste(capture.output(print(f)), collapse = " ")), "\n")
set.seed(20261014)
t0 <- proc.time()
fit <- suppressMessages(runOccJSDM(data = sampledata, listParams = list(n_factors = 2), threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"), collCovariates = c("X_theta.1", "X_theta.2"),
  spatCovariates = c("Xs.1", "Xs.2"),
  MCMCparams = list(nchain = 1, nburn = nburn, niter = niter, nthin = 1, ncores = 1),
  summarisedLatentPresences = TRUE))
el <- proc.time() - t0
cat("mode:", mode, " elapsed:", el[["elapsed"]], " chain elapsed:", fit$infos$execution$elapsed,
    " ms per iteration:", 1000 * fit$infos$execution$elapsed / (nburn + niter), "\n")
