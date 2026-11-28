# Shared by the tests of the two scale-expansion moves of the grouped prior
# (test-axis-scale-move.R and test-sigma-emb-move.R).

# A state of the coefficient block drawn from its prior at the given scales,
# with an arbitrary likelihood state: `rest` stands for every term of the
# linear predictor other than X B, kappa and Omega for the response minus one
# half and the Polya-Gamma weights.
axis_move_state <- function(seed = 1L, n = 25L, S = 7L) {
  set.seed(seed); gb <- 2L; gt <- 2L; p <- 5L
  cp <- list(type = "grouped", is_emb = c(FALSE, FALSE, TRUE, TRUE, TRUE),
             sigma_b_scale = 1, sigma_emb_scale = 1, axis_scale = 0.3,
             c_base = c(1, 1, rep(1 / sqrt(gt), 3)))
  lambda <- c(1, 1, 0.5, 0.05, 0.9); sigma_b <- 0.7; sigma_emb <- 0.4
  sds <- occJSDM:::coef_prior_sds(cp, sigma_b, sigma_emb, lambda, gb, p)
  X <- matrix(rnorm(n * p), n, p)
  Tr_b <- cbind(1, rnorm(S)); A <- matrix(rnorm(S * gt), S, gt)
  G <- matrix(rnorm(gb * p), gb, p) * sds$g
  C <- matrix(rnorm(gt * p), gt, p) * rep(sds$c, each = gt)
  R <- matrix(rnorm(S * p), S, p) * rep(sds$b, each = S)
  B <- t(Tr_b %*% G + A %*% C + R)
  rest <- matrix(rnorm(n * S, sd = 0.5), n, S)
  list(cp = cp, lambda = lambda, sigma_b = sigma_b, sigma_emb = sigma_emb, X = X, Tr_b = Tr_b,
       A = A, G = G, C = C, B = B, rest = rest, eta = rest + X %*% B,
       kappa = matrix(rbinom(n * S, 1, 0.4) - 0.5, n, S), Omega = matrix(rexp(n * S) / 4, n, S))
}

# The log of the density the sampler targets, written out term by term: the
# Normal priors of the leftovers, G and C at the standard deviations
# coef_prior_sds() gives, the half-Cauchy priors of each axis scale and of
# sigma_emb and of sigma_b expressed as densities of the log scale (hence the
# added log of each), and the likelihood given the Polya-Gamma weights. The
# prior of A, which no scale move changes, is left out.
axis_move_log_target <- function(x) {
  S <- ncol(x$B); gt <- nrow(x$C); emb <- which(x$cp$is_emb)
  sds <- occJSDM:::coef_prior_sds(x$cp, x$sigma_b, x$sigma_emb, x$lambda, nrow(x$G), ncol(x$X))
  R <- t(x$B) - x$Tr_b %*% x$G - x$A %*% x$C
  eta <- x$rest + x$X %*% x$B
  sum(dnorm(R, 0, rep(sds$b, each = S), log = TRUE)) +
    sum(dnorm(x$G, 0, sds$g, log = TRUE)) +
    sum(dnorm(x$C, 0, rep(sds$c, each = gt), log = TRUE)) +
    sum(dcauchy(x$lambda[emb], 0, x$cp$axis_scale, log = TRUE) + log(x$lambda[emb])) +
    dcauchy(x$sigma_emb, 0, x$cp$sigma_emb_scale, log = TRUE) + log(x$sigma_emb) +
    dcauchy(x$sigma_b, 0, x$cp$sigma_b_scale, log = TRUE) + log(x$sigma_b) +
    sum(x$kappa * eta - x$Omega * eta^2 / 2)
}

# A short grouped binary fit, with the scale moves on. Used by
# test-sigma-emb-move.R and test-sigma-b-move.R.
sigma_move_fit <- function(sim, MCMCparams = FIXTURE_MCMC) {
  dat <- sim$data_list
  suppressMessages(suppressWarnings(runOccJSDM(dat,
    listParams = list(n_factors = 1L, embeddingVars = dat$embeddingVars),
    occCovariates = grouped_occ_covariates(dat), MCMCparams = MCMCparams)))
}
