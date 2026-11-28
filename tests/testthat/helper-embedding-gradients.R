# An n x K embedding block with centred, uncorrelated columns of standard
# deviation 1, and K x S coefficients that load on the columns `cols` only, so
# the effect X_emb %*% B has singular values `d` and its gradients are those
# columns.
low_rank_parts <- function(n = 40L, K = 6L, S = 12L, cols = 1:2, d = c(30, 15), seed = 1L) {
  set.seed(seed)
  Q <- qr.Q(qr(scale(matrix(rnorm(n * K), n, K), scale = FALSE)))
  V <- qr.Q(qr(matrix(rnorm(S * length(cols)), S)))
  B <- matrix(0, K, S); B[cols, ] <- diag(d, length(cols)) %*% t(V) / sqrt(n - 1)
  list(X_emb = Q * sqrt(n - 1), B = B)
}

# A fit object holding only what the gradient functions read. `B_emb` is a
# list with one K x S matrix per chain: that chain's posterior mean of the
# embedding coefficients. A measured column comes first.
fake_gradient_fit <- function(X_emb, B_emb, niter = 20L, seed = 1L,
                              site_names = paste0("s", seq_len(nrow(X_emb)))) {
  set.seed(seed)
  n <- nrow(X_emb); K <- ncol(X_emb); S <- ncol(B_emb[[1]]); nchain <- length(B_emb)
  colnames(X_emb) <- paste0("emb", seq_len(K))
  B <- array(0, c(K + 1L, S, niter, nchain))
  for (ch in seq_len(nchain)) for (j in seq_len(S)) for (i in seq_len(K + 1L)) {
    noise <- rnorm(niter, sd = 0.01)
    B[i, j, , ch] <- (if (i == 1L) 0.3 else B_emb[[ch]][i - 1L, j]) + noise - mean(noise)
  }
  list(X_psi = cbind(meas1 = as.numeric(scale(rnorm(n))), X_emb),
       infos = list(coef_prior = list(type = "grouped", is_emb = c(FALSE, rep(TRUE, K))),
                    siteNames = site_names, speciesNames = paste0("sp", seq_len(S))),
       results_output = list(jsdm_output = list(B_output = B)))
}
