#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix collapsed_grid(NumericVector b, NumericVector offset, NumericMatrix l0, NumericMatrix l1) {
  if(l0.nrow()!=offset.size() || l1.nrow()!=offset.size() || l0.ncol()!=l1.ncol()) stop("dimension mismatch");
  NumericMatrix ans(b.size(),l0.ncol());
  for(int a=0;a<b.size();++a) for(int c=0;c<l0.ncol();++c) {
    double total=0;
    for(int i=0;i<offset.size();++i) {
      double eta=b[a]+offset[i];
      double soft=std::max(eta,0.0)+std::log1p(std::exp(-std::abs(eta)));
      double x=l0(i,c)-soft,y=l1(i,c)+eta-soft,m=std::max(x,y);
      total+=m+std::log(std::exp(x-m)+std::exp(y-m));
    }
    ans(a,c)=total;
  }
  return ans;
}
