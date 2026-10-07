#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace {
double log_add(double a, double b) {
  if (a == -std::numeric_limits<double>::infinity()) return b;
  if (b == -std::numeric_limits<double>::infinity()) return a;
  return std::max(a, b) + std::log1p(std::exp(-std::abs(a - b)));
}

double log_sigmoid(double x) {
  return -std::max(-x, 0.0) - std::log1p(std::exp(-std::abs(x)));
}
}

// Integrate one complete site's community likelihood, not specieswise
// marginals. The R caller has already summed occupancy/collection states
// into log_g0/log_g1 and supplies standard-normal tensor quadrature.
// [[Rcpp::export]]
arma::vec site_loglik_cpp(const arma::mat& eta,
                         const arma::mat& log_g0,
                         const arma::mat& log_g1,
                         const arma::mat& loadings,
                         const arma::mat& nodes,
                         const arma::vec& log_weights) {
  if (eta.n_rows != log_g0.n_rows || eta.n_cols != log_g0.n_cols ||
      eta.n_rows != log_g1.n_rows || eta.n_cols != log_g1.n_cols ||
      loadings.n_cols != eta.n_cols || nodes.n_cols != loadings.n_rows ||
      nodes.n_rows != log_weights.n_elem || nodes.n_rows == 0) {
    Rcpp::stop("Incompatible site likelihood dimensions.");
  }
  const arma::mat shifts = nodes * loadings;
  arma::vec result(eta.n_rows);
  for (arma::uword i = 0; i < eta.n_rows; ++i) {
    Rcpp::checkUserInterrupt();
    double total = -std::numeric_limits<double>::infinity();
    for (arma::uword k = 0; k < nodes.n_rows; ++k) {
      double conditional = log_weights(k);
      for (arma::uword j = 0; j < eta.n_cols; ++j) {
        const double predictor = eta(i, j) + shifts(k, j);
        conditional += log_add(log_sigmoid(-predictor) + log_g0(i, j),
                               log_sigmoid(predictor) + log_g1(i, j));
      }
      total = log_add(total, conditional);
    }
    result(i) = total;
  }
  return result;
}
