// Historical helpers retired on 7 October 2026 from src/functions.cpp.
// Original source: 0c5b1ff19c4de43be329fbd06730aee24548ca60.
// Preserved verbatim for reference; not loaded or compiled by the package.


/////////////////////////////////////////////////
/////////  W SAMPLER
/////////////////////////////////////////////////

// [[Rcpp::export]]
NumericMatrix sample_w_cpp(const NumericMatrix& logy1,
                           double mu0, double sigma0,
                           double mu1, double sigma1,
                           const NumericMatrix& theta,
                           const NumericVector& theta0,
                           const NumericMatrix& p,
                           const NumericMatrix& q,
                           const IntegerVector& M,
                           const IntegerVector& K,
                           const IntegerVector& sumL,
                           const IntegerVector& sumM,
                           const IntegerVector& sumK,
                           int maxL,
                           const NumericMatrix& z) {

  int S = theta.ncol();
  int N = theta.nrow();
  int n = M.size();

  NumericMatrix w(N, S);

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) {
      for (int m = 0; m < M[i]; m++) {

        // compute log p(w = 1)
        double log_p1 = 0.0;
        for (int l = 0; l < maxL; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;
            if(logy1(idxK, s) == 0){
              log_p1 += log(1 - p(l,s));
            } else {
              log_p1 += log(p(l,s)) + R::dnorm(logy1(idxK, s), mu1, sigma1, true);
            }
          }
        }

        // compute log p(w = 0)
        double log_p0 = 0.0;
        for (int l = 0; l < maxL; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;
            if(logy1(idxK, s) == 0){
              log_p0 += log(1 - q(l,s));
            } else {
              log_p0 += log(q(l,s)) + R::dnorm(logy1(idxK, s), mu0, sigma0, true);
            }
          }
        }

        // conditional on z[i, s]
        if (z(i, s) == 1.0) {
          // log_p1 += log(theta(sumM[i] + m, s));
          // log_p0 += log(1 - theta(sumM[i] + m, s));;//R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
          log_p1 += R::dbinom(1.0, 1.0, theta(sumM[i] + m, s), true);
          log_p0 += R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
        } else {
          log_p1 += //log(theta0[s]);
            R::dbinom(1.0, 1.0, theta0[s], true);
          log_p0 += //log(1 - theta0[s]);
            R::dbinom(0.0, 1.0, theta0[s], true);
        }

        // numerical stability
        double maxlog = std::max(log_p1, log_p0);
        double p1exp = std::exp(log_p1 - maxlog);
        double p0exp = std::exp(log_p0 - maxlog);
        double p_ws1 = p1exp / (p1exp + p0exp);

        // sample w
        w(sumM[i] + m, s) = R::rbinom(1.0, p_ws1);
      }
    }
  }

  return w;
}


// [[Rcpp::export]]
arma::mat sample_betatheta_cpp_parallel_old(const arma::mat& w,
                                        const arma::mat& z,
                                        arma::mat beta_theta, // Passed by value from R, safe to modify locally
                                        const arma::uvec& idx_z,
                                        const arma::mat& X_theta,
                                        const arma::vec& b_betatheta,
                                        const arma::mat& B_betatheta) {

  int S = beta_theta.n_cols;

  arma::uvec idx_z_cpp = idx_z - 1;
  arma::mat z_all = z.rows(idx_z_cpp);

  // Tell OpenMP to parallelize this loop.
  // All variables declared inside the loop become private to each thread.
  // #pragma omp parallel for
  for (int s = 0; s < S; ++s) {

    // Get the s-th column of z_all
    arma::vec z_col = z_all.col(s);

    // Find indices where z_all[, s] == 1
    arma::uvec find_ones = arma::find(z_col == 1);

    if (find_ones.is_empty()) continue;

    // Extract the column subset
    arma::vec w_sub = w.elem(find_ones + s * w.n_rows);
    arma::vec k = w_sub - 0.5;

    arma::vec n = arma::ones<arma::vec>(k.n_elem);
    arma::mat X_thetasubset = X_theta.rows(find_ones);

    // Extract current column of beta_theta
    arma::vec beta_sub = beta_theta.col(s);

    beta_theta.col(s) = sample_beta_nocov_cpp_TS(beta_sub, X_thetasubset, b_betatheta, B_betatheta, n, k);

  }

  return beta_theta;
}
