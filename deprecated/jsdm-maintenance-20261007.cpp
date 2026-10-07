// Historical helpers retired on 7 October 2026 from src/jsdm.cpp.
// Original source: 0c5b1ff19c4de43be329fbd06730aee24548ca60.
// Preserved verbatim for reference; not loaded or compiled by the package.


// [[Rcpp::export]]
bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){

  for(int k = 0; k < X_tilde.n_rows; k++){

    if((X_tilde(k,1) < y_grid[j + 1]) && (X_tilde(k,1) > y_grid[j - 1])){
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

    if((X_tilde(k,1) < y_grid[j + 1]) && (X_tilde(k,1) > y_grid[j - 1])){
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

    if((X_tilde(k,0) < x_grid[i + 1]) && (X_tilde(k,0) > x_grid[i - 1])){
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

    if((X_tilde(k,0) < x_grid[i + 1]) && (X_tilde(k,0) > x_grid[i - 1])){
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


// [[Rcpp::export]]
arma::cube convert_to_correlation(arma::cube L_output_vec, int niter, int S, int d) {
  // Initialize the output 3D array (cube) with dimensions: niter x S x S
  arma::cube Lambda_output(S, S, niter, arma::fill::zeros);

  for (int iter = 0; iter < niter; ++iter) {
    // 1. Extract the slice for the current iteration and treat it as a matrix.
    // In Armadillo, cube slicing gives us an S x d view.
    // Because R stores arrays in column-major order, we transpose or copy carefully
    // to match R's `matrix(..., byrow = TRUE)` behavior.
    arma::mat L_output_current(S, d);

    // Replicating byrow = TRUE from the original R code
    for (int i = 0; i < S; ++i) {
      for (int j = 0; j < d; ++j) {
        L_output_current(i, j) = L_output_vec(iter, i, j);
      }
    }

    // 2. Compute the covariance matrix: L %*% t(L)
    arma::mat cov_mat = L_output_current * L_output_current.t();

    // 3. Convert Covariance to Correlation (Equivalent to cov2cor)
    // R's cov2cor divides by the square root of the diagonal elements: inv(diag(sqrt(cov)))
    arma::vec inv_sqrt_diag = 1.0 / arma::sqrt(cov_mat.diag());
    arma::mat cor_mat = cov_mat % (inv_sqrt_diag * inv_sqrt_diag.t());

    // 4. Store the result in the current slice of our output cube
    Lambda_output.slice(iter) = cor_mat;
  }

  // Note: RcppArmadillo automatically converts an arma::cube back into a 3D R array.
  // The dimensions in R will be [S, S, niter].
  return Lambda_output;
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
