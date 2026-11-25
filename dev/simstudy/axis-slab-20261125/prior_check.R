# Prior check of the cap on the embedding scales (evidence step 1 of
# dev/superpowers/specs/2026-11-25-axis-slab-design.md).
#
# Nothing is fitted. Each draw is one draw from the prior of one species'
# coefficients on K embedding columns, with no measured traits, so the trait
# matrix of the covariate side is the single column of ones (the community row
# of G). The number reported is the standard deviation across sites of the
# embedding's contribution to that species' linear predictor, in logits.
#
# Adapted from dev/simstudy/grouped-prior-checks-20261111/prior_check.R: same
# sites (250), latent traits (8), draws (2,000), seeds (20261111 at 64
# columns, 20261112 at 128) and site by column matrix. Two changes:
#   - The package is loaded from --lib with library(), not devtools::load_all().
#   - The prior is read by the package's read_coef_prior() and its standard
#     deviations come from the package's coef_prior_sds(), which applies the
#     cap. The cap arm passes no slab_scale, so it uses the default (2); the
#     other arm passes slab_scale = Inf, the prior before the cap.
# Within a column count both arms use the same site by column matrix and the
# same random numbers for their prior draws (the RNG state is reset to the
# state after the matrix was drawn), so they differ only by the cap.
#
# Reading rule, fixed in the design before the run: the check passes if, with
# the cap, the median is between 1 and 2 and the 95% point under 8 at both
# column counts. The Inf rows are reported beside it and should give a 95%
# point near 27.
#
# Usage, from the package root:
#   Rscript --vanilla dev/simstudy/axis-slab-20261125/prior_check.R --lib PATH [--out FILE]
# --out saves the draws and the table as an .rds file.

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
out <- arg("out", NA_character_)
library(occJSDM, lib.loc = lib)
if (normalizePath(dirname(system.file(package = "occJSDM"))) != lib)
  stop("occJSDM was not loaded from ", lib)

read_coef_prior <- occJSDM:::read_coef_prior
coef_prior_sds <- occJSDM:::coef_prior_sds
computeBtcoef <- occJSDM:::computeBtcoef

N_SITES <- 250L
N_DRAWS <- 2000L
GT <- 8L
PROBS <- c(0.05, 0.5, 0.95, 0.99)

half_cauchy <- function(n, scale) abs(stats::rcauchy(n, 0, scale))

# The simulator's embedding columns (simulateData() in R/jsdmfun.R), as in the
# earlier script.
embedding_columns <- function(n, K) {
  Fx <- min(8L, K)
  W <- matrix(rnorm(K * Fx), K, Fx) %*% diag(1 / seq_len(Fx), Fx)
  scale(matrix(rnorm(n * Fx), n, Fx) %*% t(W) +
          0.3 * matrix(rnorm(n * K), n, K))[, , drop = FALSE]
}

# The package's prior for K embedding columns and gt latent traits. `priors`
# is what a user would pass as listPriors; empty gives the defaults.
grouped_prior <- function(K, gt, priors = list()) {
  emb <- paste0("X_emb.", seq_len(K))
  columns <- data.frame(column = emb, variable = emb, kind = "numeric")
  read_coef_prior(priors, list(embeddingVars = emb), S = 2L, gt = gt,
                  columns = columns)
}

# Standard deviation across sites of X %*% beta for n_draws prior draws of
# beta, with the scales drawn from their half-Cauchy priors and the parts'
# standard deviations from coef_prior_sds() (capped when the record's
# slab_scale is finite).
prior_contribution_sd <- function(X, cp, gt, n_draws = N_DRAWS) {
  K <- ncol(X)
  vapply(seq_len(n_draws), function(i) {
    lambda <- half_cauchy(K, cp$axis_scale)
    sigma_emb <- half_cauchy(1L, cp$sigma_emb_scale)
    # sigma_b is not used: every column here is an embedding column.
    sds <- coef_prior_sds(cp, sigma_b = 1, sigma_emb = sigma_emb,
                          lambda = lambda, g_b = 1L, p = K)
    G <- matrix(rnorm(K, 0, sds$g[1, ]), 1L, K)            # community row
    C <- matrix(rnorm(gt * K, 0, rep(sds$c, each = gt)), gt, K)
    A <- matrix(rnorm(gt), 1L, gt)                         # one species
    leftover <- matrix(rnorm(K, 0, sds$b), 1L, K)
    beta <- computeBtcoef(G, matrix(1, 1L, 1L), A, C, leftover)
    stats::sd(X %*% t(beta))
  }, numeric(1))
}

settings <- list(list(K = 64L, seed = 20261111L),
                 list(K = 128L, seed = 20261112L))
arms <- list(cap = list(), inf = list(slab_scale = Inf))

cat("Prior check of the cap on the embedding scales\n")
cat("occJSDM ", as.character(utils::packageVersion("occJSDM", lib.loc = lib)),
    " from ", lib, "\n", sep = "")
cat("R ", as.character(getRversion()), ", ", N_DRAWS, " prior draws, ",
    N_SITES, " sites, ", GT, " latent traits\n", sep = "")

rows <- list()
draws <- list()
for (s in settings) {
  set.seed(s$seed)
  X <- embedding_columns(N_SITES, s$K)
  after_X <- .Random.seed
  cat(sprintf("\nK = %d embedding columns, gt = %d, seed %d\n", s$K, GT, s$seed))
  for (a in names(arms)) {
    cp <- grouped_prior(s$K, GT, arms[[a]])
    if (a == "cap" && !identical(cp$slab_scale, 2)) stop("the default cap is not 2")
    if (a == "inf" && !is.infinite(cp$slab_scale)) stop("the Inf arm has a finite cap")
    assign(".Random.seed", after_X, envir = globalenv())
    x <- prior_contribution_sd(X, cp, GT)
    q <- stats::quantile(x, PROBS, names = FALSE)
    cat(sprintf("  slab_scale = %-4s axis_scale = %.5f (%.2f / K)   5%% %7.3f   median %7.3f   95%% %8.3f   99%% %9.3f   max %10.3f\n",
                format(cp$slab_scale), cp$axis_scale, cp$axis_scale * s$K,
                q[1], q[2], q[3], q[4], max(x)))
    rows[[length(rows) + 1L]] <- data.frame(K = s$K, slab_scale = cp$slab_scale,
                                            q05 = q[1], median = q[2], q95 = q[3],
                                            q99 = q[4], max = max(x))
    draws[[sprintf("K%d_%s", s$K, a)]] <- x
  }
}
table <- do.call(rbind, rows)

capped <- table[is.finite(table$slab_scale), ]
passes <- all(capped$median >= 1 & capped$median <= 2 & capped$q95 < 8)
cat("\nStep 1, with the cap (median between 1 and 2 and 95% point under 8 at both column counts): ",
    if (passes) "passes" else "fails", "\n", sep = "")
if (!is.na(out)) {
  saveRDS(list(table = table, draws = draws, passes = passes, lib = lib,
               R = as.character(getRversion())), out)
  cat("saved", out, "\n")
}
