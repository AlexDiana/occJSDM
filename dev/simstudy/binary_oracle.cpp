#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double loglike(const arma::vec& f, const arma::vec& y, const arma::vec& offset,
               const Rcpp::IntegerVector& loc) {
  double out = 0;
  for (unsigned i = 0; i < y.n_elem; ++i) {
    double eta = offset(i) + f(loc[i]);
    out += y(i)*eta - std::max(0.0,eta) - std::log1p(std::exp(-std::abs(eta)));
  }
  return out;
}

// Independent Gaussian-prior elliptical slice sampler. L L' is the specified
// covariance; this does not call any occJSDM sampler or Polya-Gamma code.
// Algorithm: Murray, Adams, MacKay (2010), PMLR 9:541-548.
// [[Rcpp::export]]
Rcpp::List binary_oracle(const arma::mat& L, const arma::mat& y,
                        const arma::mat& offset, const Rcpp::IntegerVector& loc,
                        int burn, int kept, int thin, double initial_scale) {
  const unsigned m=L.n_rows, S=y.n_cols, n=y.n_rows;
  arma::cube draws(m,S,kept,arma::fill::zeros);
  arma::mat qmean(n,S,arma::fill::zeros);
  double evaluations=0;
  for (unsigned s=0; s<S; ++s) {
    arma::vec ys=y.col(s), off=offset.col(s), f(m);
    for (unsigned j=0;j<m;++j) f(j)=R::rnorm(0,1);
    f=initial_scale*L*f;
    double current=loglike(f,ys,off,loc);
    for (int it=0;it<burn+kept*thin;++it) {
      arma::vec nu(m);
      for (unsigned j=0;j<m;++j) nu(j)=R::rnorm(0,1);
      nu=L*nu;
      double target=current+std::log(R::runif(0,1));
      double angle=R::runif(0,2*M_PI), lo=angle-2*M_PI, hi=angle;
      for (int proposals=0;;++proposals) {
        arma::vec proposed=f*std::cos(angle)+nu*std::sin(angle);
        double next=loglike(proposed,ys,off,loc); ++evaluations;
        if (next>target) {f=proposed; current=next; break;}
        if (angle<0) lo=angle; else hi=angle;
        angle=R::runif(lo,hi);
        if (proposals>10000) Rcpp::stop("slice failed");
      }
      if (it>=burn && (it-burn)%thin==0) {
        int k=(it-burn)/thin; draws.slice(k).col(s)=f;
        for (unsigned i=0;i<n;++i)
          qmean(i,s)+=1/(1+std::exp(-off(i)-f(loc[i])))/kept;
      }
      if (it%1000==0) Rcpp::checkUserInterrupt();
    }
  }
  return Rcpp::List::create(Rcpp::_ ["field_draws"]=draws,
                          Rcpp::_ ["probability_mean"]=qmean,
                          Rcpp::_ ["evaluations"]=evaluations);
}
