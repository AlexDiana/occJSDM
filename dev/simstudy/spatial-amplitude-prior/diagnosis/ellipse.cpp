// Independent conditional diagnostic sampler, not production fitting code.
// Murray, Adams and MacKay (2010), Algorithm 1, full 2*pi bracket.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double binary_loglik(const arma::vec& state, const arma::vec& y,
                    int trials, const arma::vec& offset, bool intercept) {
  double value = 0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    double eta = offset[i] + state[i] + (intercept ? state[y.n_elem] : 0);
    double softplus = std::max(eta, 0.0) + std::log1p(std::exp(-std::abs(eta)));
    value += y[i] * eta - trials * softplus;
  }
  return value;
}

// [[Rcpp::export]]
double ellipse_loglik_cpp(arma::vec state, arma::vec y, int trials,
                         arma::vec offset, bool intercept = false) {
  if (offset.n_elem != y.n_elem || state.n_elem != y.n_elem + intercept ||
      trials < 1 || !state.is_finite() || !offset.is_finite() ||
      !y.is_finite() || arma::any(y < 0) || arma::any(y > trials) ||
      arma::any(y != arma::floor(y)))
    Rcpp::stop("Invalid likelihood dimensions or values");
  return binary_loglik(state, y, trials, offset, intercept);
}

// [[Rcpp::export]]
Rcpp::List ellipse_draws_cpp(arma::mat lower, arma::vec y, int trials,
                           arma::vec offset, bool intercept,
                           arma::mat initial, int nburn, int niter) {
  const arma::uword d = y.n_elem + intercept;
  if (lower.n_rows != d || lower.n_cols != d || initial.n_rows != d ||
      initial.n_cols < 2 || nburn < 0 || niter < 4 ||
      !lower.is_finite() || !initial.is_finite())
    Rcpp::stop("Invalid sampler dimensions or values");
  ellipse_loglik_cpp(initial.col(0), y, trials, offset, intercept);
  arma::cube draws(d, niter, initial.n_cols);
  arma::vec proposals(initial.n_cols, arma::fill::zeros);
  for (arma::uword chain = 0; chain < initial.n_cols; ++chain) {
    arma::vec state = initial.col(chain);
    double loglik = binary_loglik(state, y, trials, offset, intercept);
    for (int iteration = 0; iteration < nburn + niter; ++iteration) {
      if (iteration % 100 == 0) Rcpp::checkUserInterrupt();
      arma::vec white(d);
      for (arma::uword j = 0; j < d; ++j) white[j] = R::rnorm(0, 1);
      arma::vec direction = lower * white;
      double threshold = loglik + std::log(R::runif(0, 1));
      double theta = R::runif(0, 2 * M_PI);
      double left = theta - 2 * M_PI;
      double right = theta;
      int attempts = 0;
      while (true) {
        if (++attempts > 10000) Rcpp::stop("Elliptical bracket failed to contract");
        arma::vec candidate = state * std::cos(theta) + direction * std::sin(theta);
        double candidate_loglik = binary_loglik(candidate, y, trials, offset, intercept);
        if (candidate_loglik >= threshold) {
          state = candidate;
          loglik = candidate_loglik;
          break;
        }
        if (theta < 0) left = theta; else right = theta;
        theta = R::runif(left, right);
      }
      proposals[chain] += attempts;
      if (iteration >= nburn) draws.slice(chain).col(iteration - nburn) = state;
    }
  }
  return Rcpp::List::create(Rcpp::Named("draws") = draws,
    Rcpp::Named("mean_proposals") = proposals / (nburn + niter));
}
