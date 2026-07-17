#include <RcppArmadilloExtensions/sample.h>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

inline std::mt19937& get_rng() {
  thread_local std::mt19937 rng(
#ifdef _OPENMP
      12345 + omp_get_thread_num()
#else
  12345
#endif
  );
  return rng;
}

inline double runif() {
  static thread_local std::uniform_real_distribution<double> dist(0.0,1.0);
  return dist(get_rng());
}

// inline double exprnd(double rate) {
//   std::exponential_distribution<double> dist(rate);
//   return dist(get_rng());
// }

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
  // return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
  return -mu * (double)std::log(1.0 - (double)runif());
}

double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);

    // if(R::runif(0.0,1.0) <= gX) {
    if(runif() <= gX) {
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

  // if(R::runif(0.0,1.0) > mu /(mu+out)) {
  if(runif() > mu /(mu+out)) {
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
      // u = R::runif(0.0, 1.0);
      u = runif();
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
    // u = R::runif(0.0,1.0);
    u = runif();
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
    // double U = R::runif(0.0,1.0) * Sn;
    double U = runif() * Sn;
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
double rinvgamma_cpp(double a, double b){
  return 1 / R::rgamma(a, 1 / b);
}


// [[Rcpp::export]]
bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) < x_grid[i + 1]){
        return(true);
      }
    }

  }

  return(false);
}

// [[Rcpp::export]]
bool isPointInBandLeft(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j) {

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) > x_grid[i - 1]){
        return(true);
      }
    }

  }

  return(false);
}

// [[Rcpp::export]]
bool isPointInBandUp(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) > y_grid[j-1]){
        return(true);
      }
    }

  }

  return(false);

}

// [[Rcpp::export]]
bool isPointInBandDown(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) < y_grid[j+1]){
        return(true);
      }
    }

  }

  return(false);

}

// [[Rcpp::export]]
IntegerVector findClosestPoint(arma::mat XY_sp, arma::mat X_tilde){

  IntegerVector closestPoint(XY_sp.n_rows);

  for(int k = 0; k < XY_sp.n_rows; k++){

    double newDistance = 0;
    double minDistance = exp(50);
    int bestIndex = 0;

    for(int i = 0; i < X_tilde.n_rows; i++){
      newDistance = pow(X_tilde(i, 0) - XY_sp(k, 0), 2) + pow(X_tilde(i, 1) - XY_sp(k, 1), 2);

      if(newDistance < minDistance){
        minDistance = newDistance;
        bestIndex = i + 1;
      }
    }

    closestPoint[k] = bestIndex;

  }

  return(closestPoint);
}


// [[Rcpp::export]]
arma::mat dist_matrix(const arma::mat& coords) {

  // coords: n x d matrix (e.g. x,y or x,y,z)
  int n = coords.n_rows;
  int d = coords.n_cols;

  arma::mat D(n, n, arma::fill::zeros);

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {

      double dist_ij = 0.0;

      for (int k = 0; k < d; k++) {
        double diff = coords(i, k) - coords(j, k);
        dist_ij += diff * diff;
      }

      dist_ij = std::sqrt(dist_ij);

      D(i, j) = dist_ij;
      D(j, i) = dist_ij; // symmetry
    }
  }

  return D;
}

// [[Rcpp::export]]
arma::mat gpCovMatrix(const arma::mat& D,
                      double sigma2,
                      double rho) {

  int n = D.n_rows;
  arma::mat Sigma(n, n);

  // Gaussian exponential covariance
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Sigma(i, j) = sigma2 * std::exp(-D(i, j) / rho);
    }
  }

  return Sigma;
}

// GAUSSIAN PROCESS FUNCTIONS

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

// [[Rcpp::export]]
arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());

  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }
  }

  return res;
}

double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
}

// [[Rcpp::export]]
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
  arma::mat res(x1.n_rows, x2.n_rows);

  for(int i = 0; (unsigned)i < x1.n_rows; i++){
    for(int j = 0; (unsigned)j < x2.n_rows; j++){
      res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
    }
  }

  return res;
}

// [[Rcpp::export]]
arma::mat samplePGvariables(arma::mat &Xbeta){

  int n1 = Xbeta.n_rows;
  int n2 = Xbeta.n_cols;

  arma::mat Omega_mat(n1, n2);

#pragma omp parallel for collapse(2)
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      Omega_mat(i,j) = rpg(1, Xbeta(i,j));
    }
  }
  // #pragma omp parallel for
  //  for(int k=0; k<n1*n2; k++){
  //    Omega_mat[k]=rpg(1,Xbeta[k]);
  //  }

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
arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b,
                          arma::vec& Omega, arma::vec& k){

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
arma::vec sample_beta_nocov_cpp(arma::vec beta, arma::mat& X, arma::vec b,
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
arma::mat sample_betatheta_cpp(const arma::mat& w,
                               const arma::mat& z,
                               arma::mat beta_theta,
                               const arma::uvec& idx_z,
                               const arma::mat& X_theta,
                               const arma::vec& b_betatheta,
                               const arma::mat& B_betatheta) { // Pass the R function if it's not natively in C++

  int S = beta_theta.n_cols;

  arma::uvec idx_z_cpp = idx_z - 1;

  arma::mat z_all = z.rows(idx_z_cpp);

  // Loop over columns (S)
  for (int s = 0; s < S; ++s) {

    // Get the s-th column of z_all
    arma::vec z_col = z_all.col(s);

    // Find indices where z_all[, s] == 1
    // Armadillo's find() returns a uvec (unsigned vector of indices)
    arma::uvec find_ones = find(z_col == 1);

    if (find_ones.is_empty()) continue; // Safeguard if no rows equal 1

    // k <- as.vector(t(w[z_all[,s]==1,s])) - .5
    // In C++, w.submat(find_ones, uvec({static_cast<uword>(s)})) extracts the column subset
    arma::vec w_sub = w.elem(find_ones + s * w.n_rows);
    arma::vec k = w_sub - 0.5;

    // n <- rep(1, length(k))
    arma::vec n = arma::ones<arma::vec>(k.n_elem);

    // X_thetasubset <- X_theta[z_all[,s]==1,,drop=F]
    arma::mat X_thetasubset = X_theta.rows(find_ones);

    // Extract current column of beta_theta
    arma::vec beta_sub = beta_theta.col(s);

    beta_theta.col(s) = sample_beta_nocov_cpp(beta_sub, X_thetasubset, b_betatheta, B_betatheta, n, k);

  }

  return beta_theta;
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

arma::vec logistic(arma::vec x) {
  return 1.0 / (1.0 + exp(-x));
}

inline void logistic_inplace(arma::rowvec& x)
{
  x.transform([](double z) {
    return 1.0 / (1.0 + std::exp(-z));
  });
}

double quantile(arma::vec x, double p) {

  const int n = x.n_elem;
  if (p <= 0.0) return x(0);
  if (p >= 1.0) return x(n - 1);

  double h = (n - 1) * p;
  int j = std::floor(h);
  double g = h - j;

  if (j == n - 1)
    return x(n - 1);

  return (1.0 - g) * x(j) + g * x(j + 1);
}

// [[Rcpp::export]]
arma::cube computePsiOutput(
    const arma::mat& X_psi,
    const arma::cube& beta_psi_output,
    const arma::mat& X_ord,
    const arma::cube& beta_ord_output,
    const arma::cube& LL_output,
    const arma::vec& conflevels)
{

  int n = X_psi.n_rows;
  int S = beta_psi_output.n_cols;
  int niter = beta_psi_output.n_slices;

  arma::cube psi_output(3, n, S);

  arma::mat mcmc_output(n, niter);

  for (int j = 0; j < S; j++) {

    Rcout << "Computing species " << j + 1 << " out of " << S << std::endl;

    for (int iter = 0; iter < niter; iter++) {

      arma::vec bpsi = beta_psi_output.slice(iter).col(j);

      arma::vec bord = beta_ord_output.slice(iter) * LL_output.slice(iter).col(j);

      arma::vec linpred = X_psi * bpsi + X_ord * bord;
      mcmc_output.col(iter) = logistic(linpred);

    }

    for (int i = 0; i < n; i++) {

      arma::vec mcmc_output_i_sorted = arma::sort(mcmc_output.row(i).t());

      psi_output(0, i, j) = quantile(mcmc_output_i_sorted, conflevels[0]);
      psi_output(1, i, j) = quantile(mcmc_output_i_sorted, conflevels[1]);
      psi_output(2, i, j) = quantile(mcmc_output_i_sorted, conflevels[2]);

    }
  }

  return psi_output;
}

// JSDM sampling functions

// Sample B in a regression model Y ~ N(XB, sigma^2 I)

// [[Rcpp::export]]
arma::vec sampleBuniv(arma::mat& X, arma::mat& B,
                      arma::vec& b, arma::vec& y,
                      double sigma){

  // arma::mat cov_matrix = arma::inv(arma::trans(X) * Omega * X + arma::inv(B));
  arma::mat tX = arma::trans(X);
  // arma::mat cov_matrix = arma::inv(tXOmega * X + arma::inv(B));
  // arma::vec result = mvrnormArma(cov_matrix * (arma::trans(X) * k + arma::inv(B) * b), cov_matrix);

  arma::mat L = arma::trans(arma::chol(tX * X * (1 / sigma*sigma) + arma::inv(B)));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * y + arma::inv(B) * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));

  return(result);
}

// Sample B in a regression model Y ~ N(XB, Omega)

// [[Rcpp::export]]
arma::vec sampleB(arma::mat& X, arma::mat& B, arma::vec& b,
                  arma::vec& Omega, arma::vec& k){

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
arma::mat sample_U_cpp(arma::mat& k,
                       arma::mat& L,
                       arma::mat& X,
                       arma::mat& B,
                       arma::mat& Omega,
                       std::string model) {

  int d = L.n_rows;
  int S = k.n_cols;
  int n = k.n_rows;

  arma::mat U(n, d, arma::fill::zeros);

  if (d == 0)
    return U;

  // Compute k_new
  arma::mat XB = X * B;
  arma::mat k_new(n, S);

  if (model == "continuous") {

    k_new = (k - XB) % Omega;

  } else if (model == "binary") {

    k_new = k - Omega % XB;
  } else {
    Rcpp::stop("Unknown model");
  }

  arma::mat B_current = arma::eye(d,d);
  arma::vec b_current(d, arma::fill::zeros);

  arma::mat transL = arma::trans(L);

  // Parallel loop
  // #ifdef _OPENMP
  //   omp_set_num_threads(n_threads);
  // #endif
  //
  // #pragma omp parallel for
  for(int i = 0; i < n; i++) {

    arma::vec Omega_i = arma::conv_to<arma::vec>::from(Omega.row(i));
    arma::vec k_new_i = arma::conv_to<arma::vec>::from(k_new.row(i));

    arma::vec result = sampleB(
      transL,
      B_current,
      b_current,
      Omega_i,
      k_new_i
    );

    U.row(i) = arma::conv_to<arma::rowvec>::from(result);
  }


  return U;
}

// SPATIAL APPROXIMATOR FUNCTIONS

// [[Rcpp::export]]
arma::mat XsBs(arma::mat &A,
               arma::mat &B,
               arma::mat &X_s_centers){

  int n = A.n_rows;
  int m = B.n_rows;
  int maxPoints = X_s_centers.n_cols;

  arma::mat AB_out = arma::mat(n, maxPoints);

  for(int i = 0; i < n; i++){
    for(int j = 0; j < maxPoints; j++){
      int colIndex = X_s_centers(i, j) - 1;
      for(int l = 0; l < m; l++){
        AB_out(i, j) += A(i, l) * B(l, colIndex);
      }
    }
  }

  return(AB_out);

}

// [[Rcpp::export]]
arma::mat KsBproduct(arma::mat &Ks,
                     arma::mat &B,
                     arma::mat &X_s_centers){

  int n = Ks.n_rows;
  int S = B.n_cols;
  int maxPoints = X_s_centers.n_cols;

  arma::mat XB = arma::zeros(n, S);

  for(int s = 0; s < S; s++){
    for(int i = 0; i < n; i++){
      for(int j = 0; j < maxPoints; j++){
        XB(i, s) += Ks(i, j) * B(X_s_centers(i, j) - 1, s);
      }
    }
  }

  return(XB);
}

// [[Rcpp::export]]
arma::mat XtOmegaX_SoR(arma::mat X,
                       int X_centers,
                       arma::vec Omega,
                       arma::mat X_s_index,
                       arma::mat &X_s_sor){

  int p = X.n_cols;

  arma::mat XtOmegaX = arma::zeros(X_centers + p,
                                   X_centers + p);

  int maxPoints = X_s_index.n_cols;

  // spatial  covariates times spatial covariates
  for(int l = 0; l < maxPoints; l++){

    for(int l2 = 0; l2 < maxPoints; l2++){

      for (int i = 0; (unsigned)i < Omega.size(); i++){

        XtOmegaX(p + X_s_index(i,l) - 1, p + X_s_index(i,l2) - 1) +=
          X_s_sor(i, l) * Omega[i] * X_s_sor(i, l2);

      }

    }

  }

  // spatial covariates times standard covariates
  for (int i = 0; (unsigned)i < Omega.size(); i++){

    for(int j = 0; j < p; j++){

      for(int l = 0; l < maxPoints; l++){

        XtOmegaX(j, p + X_s_index(i,l) - 1) +=
          X(i, j) * X_s_sor(i, l) * Omega[i];

      }
    }
  }

  for (int i = 1; i <= X_centers; i++) {

    for (int j = 0; j < p; j++){

      XtOmegaX(p + i - 1, j) = XtOmegaX(j, p + i - 1);
    }

  }

  // standard covariates times standard covariates

  for(int i = 0; i < p; i++){
    for (int j = 0; j <= i; j++) {
      for(int l = 0; l < Omega.size(); l++){
        XtOmegaX(i, j) +=
          Omega[l] * X(l, i) * X(l, j);
      }
    }
  }

  for (int i = 0; i < (p- 1); i++) {
    for (int j = i; j < p; j++) {
      XtOmegaX(i, j) = XtOmegaX(j, i);
    }
  }

  return(XtOmegaX);
}

arma::vec XtK_SoR(arma::mat X, arma::mat &X_s_index, arma::mat &X_s_sor,
                  arma::vec &k, int centers){

  int p = X.n_cols;

  arma::vec Xk = arma::zeros(p + centers);

  for(int i = 0; i < p; i++){
    Xk(i) = as_scalar(k.t() * X.col(i));
  }

  for(int i = 0; i < k.size(); i++){
    for(int l = 0; l < X_s_index.n_cols; l++){
      Xk(p + X_s_index(i,l) - 1) += X_s_sor(i, l) * k[i];
    }
  }

  return(Xk);
}


// Sample B in a regression model Y ~ N((X|Xs)(B|Bs), Omega) where Xs is a
// subset of regressor approximator

// [[Rcpp::export]]
arma::vec sampleB_SoR(arma::mat X, arma::mat &invB, arma::vec &b,
                      arma::vec &k, arma::vec Omega,
                      arma::mat &X_s_index,
                      arma::mat &Ks,
                      int X_centers){

  arma::mat XtOmegaX = XtOmegaX_SoR(X, X_centers, Omega, X_s_index, Ks);
  arma::mat tXk = XtK_SoR(X, X_s_index, Ks, k, X_centers);

  arma::mat Lambda_B = XtOmegaX + invB;
  arma::vec mu_B = tXk + invB * b;

  arma::mat L = arma::trans(arma::chol(Lambda_B));
  arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);

  arma::vec z = arma::randn(invB.n_cols);
  arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);

  arma::vec result = v + alpha;

  return(result);
}

// [[Rcpp::export]]
arma::mat spatEffectMeanCpp(arma::cube& Bs_output,
                            arma::mat& Ks,
                            arma::mat& Xs_centers) {

  int dim1 = Ks.n_rows;
  int dim2 = Bs_output.n_cols;
  int niter = Bs_output.n_slices;

  arma::mat spatEffect_mean(dim1, dim2, arma::fill::zeros);

  for (int iter = 0; iter < niter; iter++) {

    arma::mat B = Bs_output.slice(iter);

    spatEffect_mean += (1.0/niter) *
      KsBproduct(Ks, B, Xs_centers);

  }

  return spatEffect_mean;
}
