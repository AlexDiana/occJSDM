#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// old code

double aterm(int n, double x, double t) {
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }
  return (double)exp(f);
}

double exprnd(double mu) {
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);

    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }

  return X;
}

double randinvg(double mu) {
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );

  if(R::runif(0.0,1.0) > mu /(mu+out)) {
    out = mu*mu / out;
  }
  return out;
}

double tinvgauss(double z, double t) {
  double X, u;
  double mu = 1.0/z;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();

      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }
  return X;
}

double samplepg(double z) {
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;

  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
double t = MATH_2_PI;

// Compute p, q and the ratio q / (q + p)
// (derived from scratch; derivation is not in the original paper)
double K = z*z/2.0 + MATH_PI2/8.0;
double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
double logK = (double)std::log(K);
double Kt = K * t;
double w = (double)std::sqrt(MATH_PI_2);

double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
double ratio = 1.0 / (1.0 + p_over_q);

double u, X;

// Main sampling loop; page 130 of the Windle PhD thesis
while(1)
{
  // Step 1: Sample X ? g(x|z)
  u = R::runif(0.0,1.0);
  if(u < ratio) {
    // truncated exponential
    X = t + exprnd(1.0)/K;
  }
  else {
    // truncated Inverse Gaussian
    X = tinvgauss(z, t);
  }

  // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
  int i = 1;
  double Sn = aterm(0, X, t);
  double U = R::runif(0.0,1.0) * Sn;
  int asgn = -1;
  bool even = false;

  while(1)
  {
    Sn = Sn + asgn * aterm(i, X, t);

    // Accept if n is odd
    if(!even && (U <= Sn)) {
      X = X * 0.25;
      return X;
    }

    // Return to step 1 if n is even
    if(even && (U > Sn)) {
      break;
    }

    even = !even;
    asgn = -asgn;
    i++;
  }
}
return X;
}

// [[Rcpp::export]]

double rpg(int n, double z){

  double x = 0;
  for(int i = 0; i < n; i++){
    x += samplepg(z);
  }

  return(x);
}

// [[Rcpp::export]]
arma::mat samplePGvariables(arma::mat &Xbeta){

  int n1 = Xbeta.n_rows;
  int n2 = Xbeta.n_cols;

  arma::mat Omega_mat(n1, n2);

  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){

      Omega_mat(i,j) = rpg(1, Xbeta(i, j));

    }
  }

  return(Omega_mat);
}

arma::vec mvrnormArma(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + arma::chol(sigma) * Y;
}

arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
  int ncols = cholsigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + cholsigma * Y;
}

// [[Rcpp::export]]
arma::mat diagMatrixProd(arma::mat& X, arma::vec& D){

  arma::mat result(X.n_rows, D.size());
  for(int i = 0; i < result.n_rows; i++){
    for(int j = 0; j < result.n_cols; j++){
      result(i, j) = X(i,j) * D(j);
    }
  }

  return(result);
}

// [[Rcpp::export]]
arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b, arma::vec& Omega, arma::vec& k){

  // arma::mat cov_matrix = arma::inv(arma::trans(X) * Omega * X + arma::inv(B));
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  // arma::mat cov_matrix = arma::inv(tXOmega * X + arma::inv(B));
  // arma::vec result = mvrnormArma(cov_matrix * (arma::trans(X) * k + arma::inv(B) * b), cov_matrix);

  arma::mat L = arma::trans(arma::chol(tXOmega * X + arma::inv(B)));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));

  return(result);
}

// [[Rcpp::export]]
arma::vec sample_Omega_cpp(arma::mat& X, arma::vec& beta, arma::vec& n){

  int nsize = n.size();
  arma::vec Omega_vec(nsize);

  for(int i = 0; i < nsize; i++){

    arma::vec b = X.row(i) * beta;
    Omega_vec[i] = rpg(n[i], b[0]);

  }

  return(Omega_vec);
}

// [[Rcpp::export]]
arma::vec sample_beta_nocov_cpp(arma::vec beta, arma::mat X, arma::vec b,
                                arma::mat B, arma::vec n, arma::vec k){

  arma::vec Omega = sample_Omega_cpp(X, beta, n);

  beta = sample_beta_cpp(X, B, b, Omega, k);

  return(beta);
}

double log_L_gamma_cpp(arma::vec gamma, arma::mat X, arma::vec indexes_covariates,
                       arma::vec b, arma::mat B, arma::vec Omega, arma::vec k){


  IntegerVector index_present(indexes_covariates.size());
  int l = 0;
  for(int i = 0; i < indexes_covariates.size(); i++){
    if(gamma[indexes_covariates[i]-1] == 1){
      index_present[l] = i;
      l += 1;
    }
  }

  arma::mat X_gamma(X.n_rows, l);
  arma::vec b_gamma(l);
  arma::mat B_gamma(l, l);

  for(int i = 0; i < l; i++){
    X_gamma.col(i) = X.col(index_present[i]);
    b_gamma[i] = b[index_present[i]];
    for(int j = 0; j < l; j++){
      B_gamma(i,j) = B(index_present[i], index_present[j]);
    }
  }

  arma::mat tX = arma::trans(X_gamma);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat cholXgOmX = arma::chol(tXOmega * X_gamma + arma::inv(B_gamma));

  double firstTerm = (.5) * log(det(arma::inv(B_gamma))) - log(det(cholXgOmX));

  arma::mat tXKbplusBb = arma::trans(X_gamma) * k + arma::inv(B_gamma) * b_gamma;

  arma::vec v = solve(arma::trimatl(arma::trans(cholXgOmX)),tXKbplusBb);
  arma::mat vtv = arma::trans(v) * v;

  arma::mat secondTerm = - .5 * ( (arma::trans(b_gamma) * arma::inv(B_gamma) * b_gamma) - vtv);

  // double firstTerm = (.5) * log(det(arma::inv(B_gamma))) -
    //   (.5) * log(det(arma::trans(X_gamma) * Omega * X_gamma + arma::inv(B_gamma)));

  // arma::mat tXKbplusBb = arma::trans(X_gamma) * k + arma::inv(B_gamma) * b_gamma;

  // arma::mat secondTerm = - .5 * ( (arma::trans(b_gamma) * arma::inv(B_gamma) * b_gamma) -
                                       //     arma::trans(tXKbplusBb) * arma::inv(arma::trans(X_gamma) * Omega * X_gamma +
                                                                                    //     arma::inv(B_gamma)) * tXKbplusBb);

  return(firstTerm + secondTerm(0,0));
}

int sample_int(IntegerVector samples) {
  // Rcpp::IntegerVector pool = Rcpp::seq(1, 10);
  std::random_shuffle(samples.begin(), samples.end());
  return samples[0];
}

arma::vec sample_gamma_cpp(arma::vec gamma, arma::mat X, arma::vec Omega,
                           arma::vec b, arma::mat B, int ncov, arma::vec k,
                           arma::vec indexes_covariates, int fixedIndexes,
                           double d_bar){

  arma::vec gamma_star = gamma;

  double h_ratio = 1;

  if(R::runif(0, 1) < .33333){ // add

    if((sum(gamma) - fixedIndexes) != ncov){

      // Find zero covariates
      int numOfZeroCov = ncov - (sum(gamma) - fixedIndexes);
      IntegerVector zeroCov(numOfZeroCov);
      int i = 0;
      for(int l = fixedIndexes; l < gamma.size(); l++){
        if(gamma[l] == 0) {
          zeroCov[i] = l;
          i += 1;
        }
      }

      int covariate_to_update = sample_int(zeroCov);

      gamma_star[covariate_to_update] = 1;

      h_ratio = (ncov -(sum(gamma) - fixedIndexes)) / (ncov - (sum(gamma) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) );

    }

  } else if(R::runif(0, 1) < .5){ // delete

    if((sum(gamma) - fixedIndexes) != 0){

      // Find non zero covariates
      int numOfNonZeroCov = sum(gamma) - fixedIndexes;
      IntegerVector nonZeroCov(numOfNonZeroCov);
      int i = 0;
      for(int l = fixedIndexes; l < gamma.size(); l++){
        if(gamma[l] == 1) {
          nonZeroCov[i] = l;
          i += 1;
        }
      }

      int covariate_to_update = sample_int(nonZeroCov);

      gamma_star[covariate_to_update] = 0;

      h_ratio = (ncov - (sum(gamma) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma) - fixedIndexes) + 1);

    }

  } else { // swap

    if((sum(gamma) - fixedIndexes) != 0 & (sum(gamma) - fixedIndexes) != ncov){

      // Find zero covariates
      int numOfZeroCov = ncov - (sum(gamma) - fixedIndexes);
      IntegerVector zeroCov(numOfZeroCov);
      int i = 0;
      for(int l = fixedIndexes; l < gamma.size(); l++){
        if(gamma[l] == 0) {
          zeroCov[i] = l;
          i += 1;
        }
      }
      int covariates2_to_swap = sample_int(zeroCov);

      // Find non zero covariates
      int numOfNonZeroCov = sum(gamma) - fixedIndexes;
      IntegerVector nonZeroCov(numOfNonZeroCov);
      i = 0;
      for(int l = fixedIndexes; l < gamma.size(); l++){
        if(gamma[l] == 1) {
          nonZeroCov[i] = l;
          i += 1;
        }
      }

      int covariates1_to_swap = sample_int(nonZeroCov);

      gamma_star[covariates1_to_swap] = 0;
      gamma_star[covariates2_to_swap] = 1;

      h_ratio = 1;

    }

  }

  double L_gamma_star = log_L_gamma_cpp(gamma_star, X, indexes_covariates, b, B, Omega, k);

  double L_gamma = log_L_gamma_cpp(gamma, X, indexes_covariates, b, B, Omega, k);

  h_ratio = h_ratio * exp(L_gamma_star - L_gamma);

  if(R::runif(0, 1) < h_ratio){
    gamma = gamma_star;
  }

  return(gamma);
}

// [[Rcpp::export]]
List sample_gamma_beta_cpp(arma::vec gamma, arma::vec beta, arma::mat X,
                           arma::vec b, arma::mat B, int D, arma::vec n, arma::vec k,
                           arma::vec indexes_covariates, int fixedIndexes, double d_bar){

  // resize beta and X
  IntegerVector index_present(indexes_covariates.size());
  int l = 0;
  for(int i = 0; i < indexes_covariates.size(); i++){
    if(gamma[indexes_covariates[i]-1] == 1){
      index_present[l] = i;
      l += 1;
    }
  }

  arma::mat X_gamma(X.n_rows, l);
  arma::vec beta_gamma(l);

  for(int i = 0; i < l; i++){
    X_gamma.col(i) = X.col(index_present[i]);
    beta_gamma[i] = beta[index_present[i]];
  }

  // sample Omega
  arma::vec Omega = sample_Omega_cpp(X_gamma, beta_gamma, n);

  // sample gamma
  gamma = sample_gamma_cpp(gamma, X, Omega, b, B, D, k, indexes_covariates, fixedIndexes, d_bar);

  // resize X, b and B
  l = 0;
  for(int i = 0; i < indexes_covariates.size(); i++){
    if(gamma[indexes_covariates[i]-1] == 1){
      index_present[l] = i;
      l += 1;
    }
  }

  arma::mat X_gamma2(X.n_rows, l);
  arma::vec b_gamma(l);
  arma::mat B_gamma(l, l);

  for(int i = 0; i < l; i++){
    X_gamma2.col(i) = X.col(index_present[i]);
    b_gamma[i] = b[index_present[i]];
    for(int j = 0; j < l; j++){
      B_gamma(i,j) = B(index_present[i], index_present[j]);
    }
  }

  // sample beta
  arma::vec beta_new = sample_beta_cpp(X_gamma2, B_gamma, b_gamma, Omega, k);

  for(int i = 0; i < l; i++){
    beta[index_present[i]] = beta_new[i];
  }

  return List::create(_["gamma"] = gamma,
                      _["beta"] = beta,
                      _["Omega"] = Omega);
}

// [[Rcpp::export]]
NumericMatrix sample_z_cpp(const NumericMatrix& w,
                           const NumericMatrix& psi,
                           const NumericMatrix& theta,
                           const NumericVector& theta0,
                           const IntegerVector& M,
                           const IntegerVector& sumM) {

  int S = psi.ncol();
  int n = psi.nrow();
  NumericMatrix z(n, S);

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) {

      // compute p_zsequal1
      double log_p1 = 0.0;
      for (int m = 0; m < M[i]; m++) {
        int idx = sumM[i] + m;
        log_p1 += R::dbinom(w(idx, s), 1.0, theta(idx, s), true);
      }
      log_p1 += R::dbinom(1.0, 1.0, psi(i, s), true);

      // compute p_zsequal0
      double log_p0 = 0.0;
      for (int m = 0; m < M[i]; m++) {
        int idx = sumM[i] + m;
        log_p0 += R::dbinom(w(idx, s), 1.0, theta0[s], true);
      }
      log_p0 += R::dbinom(0.0, 1.0, psi(i, s), true);

      // probability
      double maxlog = std::max(log_p1, log_p0);
      double p1 = std::exp(log_p1 - maxlog);
      double p0 = std::exp(log_p0 - maxlog);
      double p_1 = p1 / (p1 + p0);

      // sample z[i, s]
      z(i, s) = R::rbinom(1.0, p_1);
    }
  }
  return z;
}

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
NumericMatrix sample_w_cim_cipp(const NumericMatrix& y,
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

  // Precompute logs for p and q
  NumericMatrix log_p(maxL, S), log_1p(maxL, S);
  NumericMatrix log_q(maxL, S), log_1q(maxL, S);

  for(int s = 0; s < S; s++) {
    for(int l = 0; l < maxL; l++) {
      log_p(l, s) = std::log(p(l, s));
      log_1p(l, s) = std::log(1.0 - p(l, s));
      log_q(l, s) = std::log(q(l, s));
      log_1q(l, s) = std::log(1.0 - q(l, s));
    }
  }

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) {
      for (int m = 0; m < M[i]; m++) {

        // compute log p(w = 1)
        double log_p1 = 0.0;
        for (int l = 0; l < maxL; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;

            log_p1 += y(idxK, s) * log_p(l, s) + (1 - y(idxK, s)) * log_1p(l, s);

            // if(y(idxK, s) == 0){
            //   log_p1 += log(1 - p(l,s));
            // } else {
            //   log_p1 += log(p(l,s));
            // }
          }
        }

        // compute log p(w = 0)
        double log_p0 = 0.0;
        for (int l = 0; l < maxL; l++) {
          int idxL = sumL[sumM[i] + m] + l;
          for (int k = 0; k < K[idxL]; k++) {
            int idxK = sumK[idxL] + k;

            log_p0 += y(idxK, s) * log_q(l,s) + (1 - y(idxK, s)) * log_1q(l,s);

            // if(y(idxK, s) == 0){
            //   log_p0 += log(1 - q(l,s));
            // } else {
            //   log_p0 += log(q(l,s));
            // }
          }
        }

        // conditional on z[i, s]
        if (z(i, s) == 1.0) {
          log_p1 += log(theta(sumM[i] + m, s));
          log_p0 += log(1 - theta(sumM[i] + m, s));;//R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
          // log_p1 += R::dbinom(1.0, 1.0, theta(sumM[i] + m, s), true);
          // log_p0 += R::dbinom(0.0, 1.0, theta(sumM[i] + m, s), true);
        } else {
          log_p1 += log(theta0[s]);
            // R::dbinom(1.0, 1.0, theta0[s], true);
          log_p0 += log(1 - theta0[s]);
            // R::dbinom(0.0, 1.0, theta0[s], true);
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
List sample_pq_cpp(NumericMatrix c_imk, NumericMatrix w, IntegerVector primerIdx,
                   IntegerVector idx_k, int maxL, double a_p, double b_p,
                   double a_q, double b_q) {

  int S = w.ncol();
  int n_k = idx_k.size();

  NumericMatrix p(maxL, S);
  NumericMatrix q(maxL, S);

  // Main loop through columns (S)
  for (int s = 0; s < S; ++s) {

    for (int l = 0; l < maxL; ++l) {

      int w1_primerl_cases_1 = 0;
      int w1_primerl_cases_0 = 0;
      int w0_primerl_cases_1 = 0;
      int w0_primerl_cases_0 = 0;

      for (int i = 0; i < n_k; ++i) {

        int idx_ki = idx_k[i] - 1;

        if(primerIdx[i] == (l+1) & w(idx_ki, s) == 1 & c_imk(i,s) == 1){
          w1_primerl_cases_1 += 1;
        } else if(primerIdx[i] == (l+1) & w(idx_ki, s) == 1 & c_imk(i,s) == 0){
          w1_primerl_cases_0 += 1;
        } else if(primerIdx[i] == (l+1) & w(idx_ki, s) == 0 & c_imk(i,s) == 2){
          w0_primerl_cases_1 += 1;
        } else {
          w0_primerl_cases_0 += 1;
        }

      }

      p(l, s) = R::rbeta(a_p + w1_primerl_cases_1, b_p + w1_primerl_cases_0);
      q(l, s) = R::rbeta(a_q + w0_primerl_cases_1, b_q + w0_primerl_cases_0);
    }
  }

  return List::create(
    _["p"] = p,
    _["q"] = q
  );
}

double dlaplace(double x, double sigma){

  return(- log(sigma) - x / sigma);

}


// // [[Rcpp::export]]
// arma::mat computeGenzLoglik_cpp(const arma::mat& lower, const arma::mat& upper,
//                            const arma::mat& Y, const arma::cube& u,
//                            const arma::mat& L) {
//
//   int n = lower.n_rows;
//   int d = lower.n_cols;
//   int M = u.n_cols;
//
//   arma::mat y_sign = 2 * Y - 1;
//
//   arma::mat logP_all = arma::zeros(n, M);
//
//   for (int i = 0; i < n; ++i) {
//     for (int m = 0; m < M; ++m) {
//
//       arma::vec y = arma::zeros(d);
//       arma::vec logP_m = arma::zeros(d);
//       // arma::cube grad_aistar = arma::zeros(d, d, d);
//       // arma::cube grad_bistar = arma::zeros(d, d, d);
//       // arma::cube grad_y = arma::zeros(d, d, d);
//
//       arma::vec abistar(d);
//
//       for (int idx_d = 0; idx_d < d; ++idx_d) {
//
//         double ai = lower(i, idx_d);
//         double bi = upper(i, idx_d);
//
//         double ai_star, bi_star;
//
//         if (idx_d == 0) {
//           ai_star = ai / L(idx_d, idx_d);
//           bi_star = bi / L(idx_d, idx_d);
//
//         } else {
//
//           double sum_ai = dot(L.submat(idx_d, 0, idx_d, idx_d - 1), y.subvec(0, idx_d - 1));
//           double sum_bi = sum_ai; // same dot product
//
//           ai_star = (ai - sum_ai) / L(idx_d, idx_d);
//           bi_star = (bi - sum_bi) / L(idx_d, idx_d);
//
//         }
//
//         double Y_idxd = (y_sign(i, idx_d) + 1) / 2.0;
//
//         // double u_val = u(i, m, idx_d);
//
//         if (y_sign(i, idx_d) == 1) {
//           y[idx_d] = sim_rnormtrunc_modified2(ai_star, Y_idxd, u(i, m, idx_d));
//         } else {
//           y[idx_d] = sim_rnormtrunc_modified2(bi_star, Y_idxd, u(i, m, idx_d));
//         }
//
//         double logP;
//         arma::mat grad_logP;
//
//         if (Y_idxd == 1) {
//           logP = R::pnorm(ai_star, 0.0, 1.0, false, true);
//         } else {
//           logP = R::pnorm(bi_star, 0.0, 1.0, true, true);
//         }
//
//         logP_m[idx_d] = logP;
//
//       }
//
//       logP_all(i, m) = accu(logP_m);
//     }
//   }
//
//   // double loglik = accu(logP_all);
//
//   return logP_all;
// }



// // [[Rcpp::export]]
// NumericMatrix sample_cimk_cpp(const NumericMatrix& y,
//                               const NumericMatrix& logy1,
//                               NumericMatrix w,
//                               double mu1,
//                               double sigma1,
//                               double mu0,
//                               double sigma0,
//                               const NumericMatrix& p,
//                               const NumericMatrix& q,
//                               IntegerVector idx_k,
//                               IntegerVector primerIdx,
//                               int maxL) {
//
//   int N3 = idx_k.size();
//   int S = p.ncol();
//
//   NumericMatrix c_imk(N3, S);
//
//   for (int s = 0; s < S; s++) {
//     for (int i = 0; i < N3; i++) {
//
//       double term1_loglik = R::dnorm(logy1(i,s), mu1, sigma1, 1);
//       double term2_loglik = 0;
//
//       if(logy1(i,s) == 0){
//         term2_loglik = log(pi0);
//       } else {
//         term2_loglik = R::dnorm(y(i,s), 0, sigma0, 1) - log(.5);
//         // term2_loglik = dlaplace(logy1(i,s), sigma0);
//       }
//
//       int idx_primer = primerIdx[i] - 1;
//       double term1_prior, term2_prior;
//
//       if(w(idx_k[i] - 1,s) == 1){
//         term1_prior = log(p(idx_primer,s));
//         term2_prior = log(1 - p(idx_primer,s));
//         // term2_prior = R::dbinom(0, 1, p(idx_primer,s), 1);
//       } else {
//         term1_prior = log(q(idx_primer,s));
//         term2_prior = log(1 - q(idx_primer,s));
//         // term1_prior = R::dbinom(1, 1, q(idx_primer,s), 1);
//         // term2_prior = R::dbinom(0, 1, q(idx_primer,s), 1);
//       }
//
//       double term12_diff = (term1_loglik + term1_prior) - (term2_loglik + term2_prior);
//
//       // double p_cimk1 = exp(term12_diff) / (exp(term12_diff) + 1);
//       double p_cimk1 = 1 / (exp(-term12_diff) + 1);
//
//       c_imk(i,s) = R::rbinom(1, p_cimk1);
//
//     }
//   }
//
//   return(c_imk);
// }

