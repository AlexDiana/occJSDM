// DIAGNOSTIC ONLY (scratchpad). Blockwise elliptical slice sampler: one slice update per block with its own angle.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
static double col_ll(const arma::vec& y, const arma::vec& off, const arma::vec& h) {
  double v = 0;
  for (arma::uword i = 0; i < y.n_elem; ++i) { double eta = off[i] + h[i];
    v += y[i]*eta - (std::max(eta,0.0) + std::log1p(std::exp(-std::abs(eta)))); }
  return v;
}
// [[Rcpp::export]]
Rcpp::List block_ellipse_draws_cpp(arma::mat lower, arma::mat Y, arma::mat offset, arma::mat weights,
    int n_shared, double indep_scale, arma::mat truth_centred, arma::mat initial, int nburn, int niter, int thin) {
  const arma::uword n = Y.n_rows, S = Y.n_cols, blocks = n_shared + S, d = blocks * n;
  const int nkeep = (niter + thin - 1) / thin;
  arma::cube draws(d, nkeep, initial.n_cols), traces(niter, initial.n_cols, 3 * S);
  arma::vec proposals(initial.n_cols, arma::fill::zeros), truth_ss(S);
  for (arma::uword s = 0; s < S; ++s) truth_ss[s] = arma::dot(truth_centred.col(s), truth_centred.col(s));
  for (arma::uword chain = 0; chain < initial.n_cols; ++chain) {
    arma::vec state = initial.col(chain);
    arma::mat h(n, S, arma::fill::zeros);
    for (arma::uword s = 0; s < S; ++s) {
      for (int k = 0; k < n_shared; ++k) h.col(s) += weights(s,k) * state.subvec(k*n,(k+1)*n-1);
      h.col(s) += indep_scale * state.subvec((n_shared+s)*n,(n_shared+s+1)*n-1);
    }
    arma::vec ll(S); for (arma::uword s = 0; s < S; ++s) ll[s] = col_ll(Y.col(s), offset.col(s), h.col(s));
    for (int it = 0; it < nburn + niter; ++it) {
      for (arma::uword b = 0; b < blocks; ++b) {
        arma::vec cur = state.subvec(b*n,(b+1)*n-1), nu(n);
        arma::vec z(n); for (arma::uword j = 0; j < n; ++j) z[j] = R::rnorm(0,1); nu = lower * z;
        // species affected and their per-unit weight on this block
        arma::vec wt(S, arma::fill::zeros);
        if ((int)b < n_shared) wt = weights.col(b); else wt[b - n_shared] = indep_scale;
        double total = arma::sum(ll);
        double threshold = total + std::log(R::runif(0,1));
        double theta = R::runif(0, 2*M_PI), left = theta - 2*M_PI, right = theta; int attempts = 0;
        while (true) {
          if (++attempts > 10000) Rcpp::stop("bracket");
          arma::vec cand = cur*std::cos(theta) + nu*std::sin(theta), delta = cand - cur;
          arma::vec ll_new = ll; arma::mat h_new = h; double tot_new = total;
          for (arma::uword s = 0; s < S; ++s) if (wt[s] != 0) {
            h_new.col(s) = h.col(s) + wt[s]*delta; ll_new[s] = col_ll(Y.col(s), offset.col(s), h_new.col(s));
            tot_new += ll_new[s] - ll[s]; }
          if (tot_new >= threshold) { state.subvec(b*n,(b+1)*n-1) = cand; ll = ll_new; h = h_new; break; }
          if (theta < 0) left = theta; else right = theta; theta = R::runif(left, right);
        }
        proposals[chain] += attempts;
      }
      if (it >= nburn) { int kept = it - nburn;
        for (arma::uword s = 0; s < S; ++s) {
          traces(kept,chain,3*s) = arma::mean(h.col(s));
          traces(kept,chain,3*s+1) = truth_ss[s] > 0 ? arma::dot(h.col(s), truth_centred.col(s))/truth_ss[s] : 0.0;
          double p = 0; for (arma::uword i = 0; i < n; ++i) p += 1.0/(1.0+std::exp(-(offset(i,s)+h(i,s))));
          traces(kept,chain,3*s+2) = p/n; }
        if (kept % thin == 0) draws.slice(chain).col(kept/thin) = state; }
    }
  }
  return Rcpp::List::create(Rcpp::Named("draws")=draws, Rcpp::Named("traces")=traces, Rcpp::Named("mean_proposals")=proposals/((nburn+niter)*blocks));
}
