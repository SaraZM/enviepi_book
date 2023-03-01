#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]


arma::mat ComputeDistMat(arma::mat X) {
  
  unsigned int N = X.n_rows;
  unsigned int M = X.n_cols;
  
  arma::mat dist = zeros(N,N);
  
  for (unsigned int i=0; i<(N-1); i++) {
    for(unsigned int j=(i+1); j<N; j++) {
      
      double total = 0;
      for (unsigned int m=0; m<M; m++) {
        total += pow(X(i,m) - X(j,m), 2);
      }
      
      dist(i,j) = sqrt(total);
      dist(j,i) = dist(i,j);
    }
  }
  
  return dist;
}
