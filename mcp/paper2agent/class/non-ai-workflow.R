# Ordinary R fallback for the pilot. Run from the student's project root.
# Edit the input path and covariates to match the reviewed data.
# The 20/20 schedule below tests the workflow only. It cannot support ecological inference.
library(occJSDM)

data_path <- "inputs/my_data.rds"
fit_path <- "results/pilot-fit.rds"
if (file.exists(fit_path)) stop("Choose a new fit_path to preserve the existing result")
dat <- readRDS(data_path)
set.seed(1702)
fit <- runOccJSDM(
  data = dat,
  listParams = list(n_factors = 2),
  threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = "X_theta",
  spatCovariates = character(),
  MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1),
  summarisedLatentPresences = TRUE
)

dir.create(dirname(fit_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, fit_path)
diag <- returnConvergenceDiagnostics(fit)
print(diag, n = Inf)
print(c(
  rows = nrow(diag),
  rhat_above_1_01 = sum(is.finite(diag$rhat) & diag$rhat > 1.01),
  ess_below_400 = sum(is.finite(diag$ess) & diag$ess < 400),
  unavailable_rhat = sum(!is.finite(diag$rhat)),
  unavailable_ess = sum(!is.finite(diag$ess))
))
# These diagnostics cover selected coefficients, not every latent variable or assumption.
# ESS may exceed retained draw count. Favourable diagnostics do not prove convergence.

baseline_draws <- returnOccupancyRates(fit)
baseline_summary <- data.frame(
  species = colnames(baseline_draws),
  mean = colMeans(baseline_draws),
  q2.5 = apply(baseline_draws, 2, quantile, probs = 0.025),
  q97.5 = apply(baseline_draws, 2, quantile, probs = 0.975),
  row.names = NULL
)
print(baseline_summary)
# Intercept-based baseline occupancy is not occupancy at a particular site.
# Investigate flagged/unavailable diagnostics before interpreting a substantive fit.
# Faithful execution does not establish unbiased estimates or calibrated intervals.
