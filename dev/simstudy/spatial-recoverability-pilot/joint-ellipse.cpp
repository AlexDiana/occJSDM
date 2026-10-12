// Joint elliptical-slice sampler for the recoverability pilot's truth 1 oracle.
// State = [f_1 .. f_K (shared fields), g_1 .. g_S (independent parts)], each block n long.
// eta_s = offset_s + sum_k weights(s,k) f_k + indep_scale * g_s. Independent blocks are
// present only when indep_scale > 0. With K = 1, S = 1, indep_scale = 0 the random-number
// sequence and arithmetic match ellipse.cpp exactly.
// Murray, Adams and MacKay (2010), Algorithm 1, full 2*pi bracket.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

static arma::mat species_effects(const arma::vec& state, const arma::mat& weights, int n_shared,
                                 double indep_scale, arma::uword n) {
  const arma::uword S = weights.n_rows;
  arma::mat h(n, S, arma::fill::zeros);
  for (arma::uword s = 0; s < S; ++s) {
    for (int k = 0; k < n_shared; ++k) h.col(s) += weights(s, k) * state.subvec(k * n, (k + 1) * n - 1);
    if (indep_scale > 0) h.col(s) += indep_scale * state.subvec((n_shared + s) * n, (n_shared + s + 1) * n - 1);
  }
  return h;
}

static double joint_loglik(const arma::vec& state, const arma::mat& Y, const arma::mat& offset,
                           const arma::mat& weights, int n_shared, double indep_scale) {
  const arma::uword n = Y.n_rows, S = Y.n_cols;
  arma::mat h = species_effects(state, weights, n_shared, indep_scale, n);
  double value = 0;
  for (arma::uword s = 0; s < S; ++s)
    for (arma::uword i = 0; i < n; ++i) {
      double eta = offset(i, s) + h(i, s);
      double softplus = std::max(eta, 0.0) + std::log1p(std::exp(-std::abs(eta)));
      value += Y(i, s) * eta - softplus;
    }
  return value;
}

// [[Rcpp::export]]
Rcpp::List joint_ellipse_draws_cpp(arma::mat lower, arma::mat Y, arma::mat offset, arma::mat weights,
                                   int n_shared, double indep_scale, arma::mat truth_centred,
                                   arma::mat initial, int nburn, int niter, int thin) {
  const arma::uword n = Y.n_rows, S = Y.n_cols;
  const int n_indep = indep_scale > 0 ? (int) S : 0;
  const arma::uword blocks = n_shared + n_indep, d = blocks * n;
  if (lower.n_rows != n || lower.n_cols != n || offset.n_rows != n || offset.n_cols != S ||
      weights.n_rows != S || weights.n_cols != (arma::uword) n_shared || truth_centred.n_rows != n ||
      truth_centred.n_cols != S || initial.n_rows != d || initial.n_cols < 2 || nburn < 0 || niter < 4 ||
      thin < 1 || n_shared < 1 || indep_scale < 0 || !lower.is_finite() || !initial.is_finite() ||
      !offset.is_finite() || !weights.is_finite() || !truth_centred.is_finite() ||
      arma::any(arma::vectorise(Y) < 0) || arma::any(arma::vectorise(Y) > 1))
    Rcpp::stop("Invalid joint sampler dimensions or values");
  const int nkeep = (niter + thin - 1) / thin;
  arma::cube draws(d, nkeep, initial.n_cols);
  arma::cube traces(niter, initial.n_cols, 3 * S);
  arma::vec proposals(initial.n_cols, arma::fill::zeros);
  arma::vec truth_ss(S);
  for (arma::uword s = 0; s < S; ++s) truth_ss[s] = arma::dot(truth_centred.col(s), truth_centred.col(s));
  for (arma::uword chain = 0; chain < initial.n_cols; ++chain) {
    arma::vec state = initial.col(chain);
    double loglik = joint_loglik(state, Y, offset, weights, n_shared, indep_scale);
    for (int iteration = 0; iteration < nburn + niter; ++iteration) {
      if (iteration % 100 == 0) Rcpp::checkUserInterrupt();
      arma::vec white(d);
      for (arma::uword j = 0; j < d; ++j) white[j] = R::rnorm(0, 1);
      arma::vec direction(d);
      for (arma::uword b = 0; b < blocks; ++b)
        direction.subvec(b * n, (b + 1) * n - 1) = lower * white.subvec(b * n, (b + 1) * n - 1);
      double threshold = loglik + std::log(R::runif(0, 1));
      double theta = R::runif(0, 2 * M_PI);
      double left = theta - 2 * M_PI, right = theta;
      int attempts = 0;
      while (true) {
        if (++attempts > 10000) Rcpp::stop("Elliptical bracket failed to contract");
        arma::vec candidate = state * std::cos(theta) + direction * std::sin(theta);
        double candidate_loglik = joint_loglik(candidate, Y, offset, weights, n_shared, indep_scale);
        if (candidate_loglik >= threshold) { state = candidate; loglik = candidate_loglik; break; }
        if (theta < 0) left = theta; else right = theta;
        theta = R::runif(left, right);
      }
      proposals[chain] += attempts;
      if (iteration >= nburn) {
        const int kept = iteration - nburn;
        arma::mat h = species_effects(state, weights, n_shared, indep_scale, n);
        for (arma::uword s = 0; s < S; ++s) {
          traces(kept, chain, 3 * s) = arma::mean(h.col(s));
          traces(kept, chain, 3 * s + 1) = truth_ss[s] > 0 ? arma::dot(h.col(s), truth_centred.col(s)) / truth_ss[s] : 0.0;
          double p = 0;
          for (arma::uword i = 0; i < n; ++i) p += 1.0 / (1.0 + std::exp(-(offset(i, s) + h(i, s))));
          traces(kept, chain, 3 * s + 2) = p / n;
        }
        if (kept % thin == 0) draws.slice(chain).col(kept / thin) = state;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("draws") = draws, Rcpp::Named("traces") = traces,
                            Rcpp::Named("mean_proposals") = proposals / (nburn + niter));
}
