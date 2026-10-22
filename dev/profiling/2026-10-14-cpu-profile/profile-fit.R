# Profile or time one chain of a runOccJSDM() fit.
# Usage: Rscript profile-fit.R <quickstart|binary900> <nburn> <niter> <outprefix> [data.rds]
# Env: OCCJSDM_PKG is the package checkout to load (default: the working directory);
# NOPROF=1 skips Rprof (timing only); PROTO=1 installs the prototype range step
# from proto-logweights.R, which must sit beside <outprefix>.
args <- commandArgs(trailingOnly = TRUE)
which_fit <- args[1]; nburn <- as.integer(args[2]); niter <- as.integer(args[3]); out <- args[4]
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
suppressMessages(pkgload::load_all(Sys.getenv("OCCJSDM_PKG", unset = "."), quiet = TRUE))
if (Sys.getenv("PROTO") == "1") source(file.path(dirname(out), "proto-logweights.R"))
if (which_fit == "quickstart") {
  data <- sampledata
  call_fit <- function() runOccJSDM(data = data, listParams = list(n_factors = 2), threshold = 1,
    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
    collCovariates = c("X_theta.1", "X_theta.2"),
    spatCovariates = c("Xs.1", "Xs.2"),
    MCMCparams = list(nchain = 1, nburn = nburn, niter = niter, nthin = 1, ncores = 1),
    summarisedLatentPresences = TRUE)
} else {
  data <- readRDS(args[5])
  call_fit <- function() runOccJSDM(data = data, listParams = list(n_factors = 2), threshold = 1,
    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
    spatCovariates = c("Xs.1", "Xs.2"),
    MCMCparams = list(nchain = 1, nburn = nburn, niter = niter, nthin = 1, ncores = 1),
    summarisedLatentPresences = TRUE)
}
prof <- Sys.getenv("NOPROF") != "1"
set.seed(20261014)
if (prof) Rprof(paste0(out, ".Rprof"), interval = 0.005, line.profiling = FALSE)
t0 <- proc.time()
fit <- suppressMessages(call_fit())
el <- proc.time() - t0
if (prof) Rprof(NULL)
cat("fit:", which_fit, " elapsed:", el[["elapsed"]], "user:", el[["user.self"]], "sys:", el[["sys.self"]],
    " chain elapsed:", fit$infos$execution$elapsed, "\n")
cat("iterations:", nburn + niter, " ms per iteration (whole fit):", 1000 * el[["elapsed"]] / (nburn + niter), "\n")
if (!prof) quit(save = "no")
s <- summaryRprof(paste0(out, ".Rprof"))
cat("sampling.time:", s$sampling.time, "\n")
cat("\n== by.self (top 40) ==\n"); print(head(s$by.self, 40))
cat("\n== by.total (top 60) ==\n"); print(head(s$by.total, 60))
saveRDS(list(summary = s, elapsed = el, nburn = nburn, niter = niter, fit_infos = fit$infos$execution), paste0(out, "-summary.rds"))
