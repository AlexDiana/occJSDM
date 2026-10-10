// Blockwise elliptical-slice sampler for the recoverability pilot's truth 1 oracle.
// State = [f_1 .. f_K (shared fields), g_1 .. g_S (independent parts)], each block n long.
// eta_s = offset_s + sum_k weights(s,k) f_k + indep_scale * g_s. Independent blocks are
// present only when indep_scale > 0. Each sweep updates one block at a time with its own
// ellipse angle, conditional on the other blocks, under the same prior factor and
// likelihood as joint-ellipse.cpp, which moves every block with one angle. With K = 1,
// S = 1, indep_scale = 0 the random-number sequence and arithmetic match ellipse.cpp exactly.
// Murray, Adams and MacKay (2010), Algorithm 1, full 2*pi bracket, applied per block.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

static double column_loglik(const arma::vec& y, const arma::vec& offset, const arma::vec& h) {
  double value = 0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    double eta = offset[i] + h[i];
    double softplus = std::max(eta, 0.0) + std::log1p(std::exp(-std::abs(eta)));
    value += y[i] * eta - softplus;
  }
  return value;
}

// Effect of the whole state on species s.
static arma::vec species_effect(const arma::vec& state, const arma::mat& weights, int n_shared,
                                double indep_scale, arma::uword n, arma::uword s) {
  arma::vec h(n, arma::fill::zeros);
  for (int k = 0; k < n_shared; ++k) h += weights(s, k) * state.subvec(k * n, (k + 1) * n - 1);
  if (indep_scale > 0) h += indep_scale * state.subvec((n_shared + s) * n, (n_shared + s + 1) * n - 1);
  return h;
}

// Accepted angle folded into (-pi, pi], so that a move near 2*pi counts as small.
static double folded_angle(double theta) { return theta - 2 * M_PI * std::round(theta / (2 * M_PI)); }

// [[Rcpp::export]]
Rcpp::List block_ellipse_draws_cpp(arma::mat lower, arma::mat Y, arma::mat offset, arma::mat weights,
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
    Rcpp::stop("Invalid blockwise sampler dimensions or values");
  const int nkeep = (niter + thin - 1) / thin;
  arma::cube draws(d, nkeep, initial.n_cols);
  arma::cube traces(niter, initial.n_cols, 3 * S);
  arma::vec proposals(initial.n_cols, arma::fill::zeros);
  arma::mat block_proposals(blocks, initial.n_cols, arma::fill::zeros);
  arma::mat block_angle(blocks, initial.n_cols, arma::fill::zeros);
  arma::vec truth_ss(S);
  for (arma::uword s = 0; s < S; ++s) truth_ss[s] = arma::dot(truth_centred.col(s), truth_centred.col(s));
  for (arma::uword chain = 0; chain < initial.n_cols; ++chain) {
    arma::vec state = initial.col(chain);
    arma::mat h(n, S);
    arma::vec loglik(S);
    for (arma::uword s = 0; s < S; ++s) {
      h.col(s) = species_effect(state, weights, n_shared, indep_scale, n, s);
      loglik[s] = column_loglik(Y.col(s), offset.col(s), h.col(s));
    }
    for (int iteration = 0; iteration < nburn + niter; ++iteration) {
      if (iteration % 100 == 0) Rcpp::checkUserInterrupt();
      for (arma::uword b = 0; b < blocks; ++b) {
        const arma::uword first = b * n, last = (b + 1) * n - 1;
        arma::vec current = state.subvec(first, last);
        arma::vec white(n);
        for (arma::uword j = 0; j < n; ++j) white[j] = R::rnorm(0, 1);
        arma::vec direction = lower * white;
        // Species this block feeds, their weight on it, and their effect from the other blocks.
        std::vector<arma::uword> affected; std::vector<double> weight;
        if ((int) b < n_shared) {
          for (arma::uword s = 0; s < S; ++s) if (weights(s, b) != 0) { affected.push_back(s); weight.push_back(weights(s, b)); }
        } else { affected.push_back(b - n_shared); weight.push_back(indep_scale); }
        const arma::uword m = affected.size();
        arma::mat rest(n, m);
        double total = 0;
        for (arma::uword a = 0; a < m; ++a) {
          rest.col(a) = species_effect(state, weights, n_shared, indep_scale, n, affected[a]) - weight[a] * current;
          total += loglik[affected[a]];
        }
        double threshold = total + std::log(R::runif(0, 1));
        double theta = R::runif(0, 2 * M_PI);
        double left = theta - 2 * M_PI, right = theta;
        int attempts = 0;
        arma::mat candidate_h(n, m);
        arma::vec candidate_loglik(m);
        while (true) {
          if (++attempts > 10000) Rcpp::stop("Elliptical bracket failed to contract");
          arma::vec candidate = current * std::cos(theta) + direction * std::sin(theta);
          double candidate_total = 0;
          for (arma::uword a = 0; a < m; ++a) {
            candidate_h.col(a) = rest.col(a) + weight[a] * candidate;
            candidate_loglik[a] = column_loglik(Y.col(affected[a]), offset.col(affected[a]), candidate_h.col(a));
            candidate_total += candidate_loglik[a];
          }
          if (candidate_total >= threshold) {
            state.subvec(first, last) = candidate;
            for (arma::uword a = 0; a < m; ++a) { h.col(affected[a]) = candidate_h.col(a); loglik[affected[a]] = candidate_loglik[a]; }
            break;
          }
          if (theta < 0) left = theta; else right = theta;
          theta = R::runif(left, right);
        }
        proposals[chain] += attempts;
        if (iteration >= nburn) { block_proposals(b, chain) += attempts; block_angle(b, chain) += std::abs(folded_angle(theta)); }
      }
      if (iteration >= nburn) {
        const int kept = iteration - nburn;
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
  const double block_steps = (double) (nburn + niter) * (double) blocks;
  return Rcpp::List::create(Rcpp::Named("draws") = draws, Rcpp::Named("traces") = traces,
                            Rcpp::Named("mean_proposals") = proposals / block_steps,
                            Rcpp::Named("block_proposals") = block_proposals / niter,
                            Rcpp::Named("block_mean_abs_angle") = block_angle / niter);
}
