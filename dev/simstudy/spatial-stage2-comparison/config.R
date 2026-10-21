# Frozen protocol for stage 2. Scientific changes need a dated spec amendment first.
stage2_config <- function() {
  list(schema_version = 1L,
       arms = data.frame(package = c("occJSDM", "occJSDM", "Hmsc", "Hmsc", "gllvm"),
                         arm = c("env", "gp", "spatial2", "spatial6", "spatial2"), stringsAsFactors = FALSE),
       cells = c("t1-r0.08-n300", "t1-r0.08-n900", "t2-n300", "t2-n900"), replicates = 1:10, arrangement = "uniform",
       gate_occupancy = c(0.25, 0.50, 0.75),
       mcmc = list(initial = c(nchain = 4L, nburn = 3000L, niter = 4000L, nthin = 1L),
                   extended = c(nchain = 4L, nburn = 8000L, niter = 16000L, nthin = 1L)),
       extension_factor = 4.43,
       diagnostics = list(rhat = 1.01, ess = 400, gradient = 0.01, objective_agreement = 0.1,
                          probability_agreement = 0.01, initial_starts = 3L, max_starts = 6L),
       gllvm = list(maxit = 12000L, method = "VA", family = "binomial", link = "logit"),
       support_points_300 = 300L,
       # Frozen from the timing pilot before dispatch (run.R --mode=freeze); NA until then.
       frozen = list(support_points_900 = NA_integer_, hmsc_method = NA_character_, hmsc_knots = NA_integer_, prediction_thin = NA_integer_),
       ladder = list(support_points_900 = c(900L, 600L, 450L, 300L), hmsc = list(c("Full", NA), c("GPP", 300L))),
       delta = 0.10, min_paired = 6L,
       comparisons = list(c("occJSDM/gp", "Hmsc/spatial2"), c("occJSDM/gp", "Hmsc/spatial6"), c("occJSDM/gp", "gllvm/spatial2"),
                          c("Hmsc/spatial6", "Hmsc/spatial2"), c("Hmsc/spatial2", "gllvm/spatial2")),
       limits = list(attempts = 400L, jobs = 200L), wall_hours = 48, workers = 8L, prediction_block = 25L)
}

# study_config() is the name the upstream diagnostics and selection code call.
study_config <- function() {
  cfg <- stage2_config()
  list(diagnostics = cfg$diagnostics, gllvm = cfg$gllvm, mcmc = cfg$mcmc, prediction_block = cfg$prediction_block)
}

arm_key <- function(package, arm) paste(package, arm, sep = "/")
