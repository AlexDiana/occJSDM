# A rare species that one column separates (evidence step 2 of
# dev/superpowers/specs/2026-11-25-axis-slab-design.md).
#
# The data set, built here and not by simulateOccJSDMData(), with seed 3001:
#   - 150 sites and 10 embedding columns, drawn as independent standard
#     Normals (one 150 by 10 matrix, filled column by column) and then
#     standardised with scale().
#   - 6 species. Species 1 is present at exactly the 6 sites with the largest
#     values of column 1 and absent elsewhere. Species 2 to 6 are present
#     independently at each site with probability 0.3 (one 150 by 5 matrix of
#     rbinom() draws, filled column by column, drawn after the columns) and
#     respond to nothing.
# The fit: model type binary (one row per site), all 10 columns in
# embeddingVars and in occCovariates, no measured covariate, no traits,
# no spatial covariates, n_lattrait = 2, n_factors = 2, 4 chains of 5,000
# burn-in and 5,000 kept sweeps, fit seed 3101 (set.seed and setOccJSDMSeed).
# The data are passed to runOccJSDM() as tests/testthat/helper-grouped-prior.R
# passes a binary fit with embedding columns: data$info with the covariate
# columns and data$OTU with one row per site.
#
# --slab default fits with listPriors = list(), so the cap is the package's
# default (checked to be 2); --slab Inf fits with listPriors =
# list(slab_scale = Inf). Any other value is passed as slab_scale.
#
# The number read is the 99% point, over the kept draws of all chains, of the
# absolute value of species 1's coefficient on column 1; the largest absolute
# draw is reported beside it. Reading rule, fixed in the design before the
# run (limit 14): with the cap, under 14 passes and otherwise the cap fails;
# without the cap, over 14 reads "the tail did harm here" and under 14 reads
# "no harm seen in this data set".
#
# Usage, from anywhere:
#   Rscript --vanilla separation.R --lib PATH --slab default|Inf [--ncores 1]
#     [--commit SHA] [--nburn 5000] [--niter 5000] [--label NAME]
# --nburn, --niter and --label exist for the short check of the script only.
# Writes dev/simstudy/results/axis-slab-20261125/separation-<label>.rds, where
# the label is the --slab value unless --label is given.

Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing --", name)
    return(default)
  }
  args[i + 1L]
}
lib <- normalizePath(arg("lib"), mustWork = TRUE)
slab_arg <- arg("slab")
label <- arg("label", slab_arg)
ncores <- as.integer(arg("ncores", "1"))
nburn <- as.integer(arg("nburn", "5000"))
niter <- as.integer(arg("niter", "5000"))
commit <- arg("commit", NA_character_)
priors <- if (slab_arg == "default") list() else list(slab_scale = as.numeric(slab_arg))
if (length(priors) && (length(priors$slab_scale) != 1L || is.na(priors$slab_scale)))
  stop("--slab must be 'default' or a number")

library(occJSDM, lib.loc = lib)
if (normalizePath(dirname(system.file(package = "occJSDM"))) != lib)
  stop("occJSDM was not loaded from ", lib)
RcppParallel::setThreadOptions(numThreads = 1)

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
study_dir <- dirname(normalizePath(script))
out_dir <- file.path(study_dir, "..", "results", "axis-slab-20261125")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
out_file <- file.path(out_dir, sprintf("separation-%s.rds", label))

power_state <- function() {
  paste(c(system("pmset -g batt", intern = TRUE)[1],
          grep("lowpowermode", system("pmset -g", intern = TRUE), value = TRUE)),
        collapse = "; ")
}

LIMIT <- 14
N_SITES <- 150L
N_COLS <- 10L
N_SPECIES <- 6L

# The data set.
set.seed(3001L)
X <- scale(matrix(rnorm(N_SITES * N_COLS), N_SITES, N_COLS))
emb <- paste0("X_emb.", seq_len(N_COLS))
colnames(X) <- emb
Y <- matrix(0L, N_SITES, N_SPECIES)
Y[order(X[, 1], decreasing = TRUE)[1:6], 1] <- 1L
Y[, 2:N_SPECIES] <- matrix(stats::rbinom(N_SITES * (N_SPECIES - 1L), 1L, 0.3),
                           N_SITES, N_SPECIES - 1L)
colnames(Y) <- paste0("OTU_", seq_len(N_SPECIES))
dat <- list(info = as.data.frame(unclass(X)[, emb, drop = FALSE]), OTU = Y)
stopifnot(sum(Y[, 1]) == 6L, all(Y[Y[, 1] == 1L, 1] == 1L),
          min(X[Y[, 1] == 1L, 1]) > max(X[Y[, 1] == 0L, 1]))

MCMCparams <- list(nchain = 4L, nburn = nburn, niter = niter, nthin = 1L,
                   ncores = ncores)
cat("occJSDM", as.character(utils::packageVersion("occJSDM", lib.loc = lib)), "from", lib, "\n")
cat("arm", label, ": listPriors =", deparse(priors), "\n")
power_start <- power_state()
cat("power at start:", power_start, "\n")
started <- Sys.time()
cat("started", format(started, tz = "GMT", usetz = TRUE), "\n")

set.seed(3101L)
occJSDM:::setOccJSDMSeed(3101L)
fit <- runOccJSDM(dat,
  listParams = list(n_factors = 2L, n_lattrait = 2L, embeddingVars = emb),
  listPriors = priors,
  occCovariates = emb, spatCovariates = NULL, MCMCparams = MCMCparams)
finished <- Sys.time()
minutes <- as.numeric(difftime(finished, started, units = "mins"))
power_end <- power_state()

jo <- fit$results_output$jsdm_output
cp <- fit$infos$coef_prior
if (!identical(cp$type, "grouped") || !all(cp$is_emb)) stop("the fit is not grouped on all 10 columns")
if (slab_arg == "default" && !identical(cp$slab_scale, 2)) stop("the default cap is not 2")
if (!identical(colnames(fit$X_psi)[1], "X_emb.1")) stop("column 1 of X_psi is not X_emb.1")
if (!identical(as.character(fit$infos$speciesNames[1]), "OTU_1")) stop("species 1 is not OTU_1")
if (!identical(fit$infos$model, "binary")) stop("the fit's model type is not binary")

# Species 1's coefficient on column 1 and its parts, kept sweeps by chains.
# With no traits the trait matrix is a column of ones, so the community part is
# row 1 of G; the latent-trait part is A[1, ] times column 1 of C; the
# leftover is the rest (computeBtcoef() in R/jsdmfun.R).
b11 <- jo$B_output[1, 1, , ]
community <- jo$G_output[1, 1, , ]
dims <- c(dim(jo$C_output)[1], niter, MCMCparams$nchain)   # latent traits, sweeps, chains
lattrait <- colSums(array(jo$A_output[1, , , , drop = FALSE], dims) *
                      array(jo$C_output[, 1, , , drop = FALSE], dims))
leftover <- b11 - community - lattrait
q99 <- function(x) unname(stats::quantile(abs(x), 0.99))
reading <- if (is.finite(cp$slab_scale)) {
  if (q99(b11) < LIMIT) "passes" else "fails"
} else {
  if (q99(b11) > LIMIT) "the tail did harm here" else "no harm seen in this data set"
}

conv <- returnGroupedConvergence(fit)
res <- list(
  label = label, slab_scale = cp$slab_scale, listPriors = priors,
  b11 = b11, q99 = q99(b11), max_abs = max(abs(b11)), reading = reading,
  parts = list(community = community, lattrait = lattrait, leftover = leftover),
  parts_q99 = c(community = q99(community), lattrait = q99(lattrait),
                leftover = q99(leftover)),
  parts_max_abs = c(community = max(abs(community)), lattrait = max(abs(lattrait)),
                    leftover = max(abs(leftover))),
  axis1 = jo$lambda_output[1, , ],          # effective axis scale of column 1
  sigma_emb = jo$sigmaemb_output,
  b_species1 = jo$B_output[, 1, , ],        # species 1 on every column
  convergence = conv, verdict = conv$verdict,
  axis_move = fit$infos$axis_move,
  minutes = minutes, power_start = power_start, power_end = power_end,
  started = format(started, tz = "GMT", usetz = TRUE),
  finished = format(finished, tz = "GMT", usetz = TRUE),
  lib = lib, commit = commit, R = as.character(getRversion()),
  settings = list(MCMCparams = MCMCparams, data_seed = 3001L, fit_seed = 3101L,
                  n_lattrait = 2L, n_factors = 2L, limit = LIMIT),
  data = dat)
saveRDS(res, out_file)

cat(sprintf("\n%s: slab_scale %s, |species 1 on column 1|: 99%% point %.3f, largest %.3f -> %s\n",
            label, format(cp$slab_scale), res$q99, res$max_abs, reading))
cat("parts, 99% point of |part|:", sprintf("%s %.3f", names(res$parts_q99), res$parts_q99), "\n")
cat("effective axis scale of column 1: median", sprintf("%.4f", stats::median(res$axis1)),
    ", 99% point", sprintf("%.4f", unname(stats::quantile(res$axis1, 0.99))),
    ", largest", sprintf("%.4f", max(res$axis1)), "\n")
cat("returnGroupedConvergence() verdict:", conv$verdict, "\n")
print(conv$checks[, c("check", "value", "limit", "flagged")])
cat(sprintf("%.1f minutes; power at end: %s\n", minutes, power_end))
cat("saved", out_file, "\n")
