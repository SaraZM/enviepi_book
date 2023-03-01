#include <RcppArmadillo.h>
#include <ComputeDistMat.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]

arma::mat prediction_marginal_gp12(arma::colvec y, arma::mat obsMu, arma::mat predMu,
                                   arma::mat obsDistMat, arma::mat pred2ObsDistMat,
                                   arma::vec sigmasq, arma::vec phi,
                                   arma::vec tausq, unsigned int iterprint)
{
  
  unsigned int L = predMu.n_cols;     // number of posterior samples
  unsigned int N = obsMu.n_rows;     // number of data locations
  unsigned int N0 = predMu.n_rows;    // number of prediction locations
  
  arma::mat ypredsamples = zeros(N0,L);
  arma::vec quant = zeros(3);
  arma::vec prob = {0.025,0.5,0.975};
  arma::mat C, Cinv;
  arma::rowvec Weight;
  arma::vec resid;
  arma::rowvec c0, thisXdist;
  double m, v;
  
  for(unsigned int l=0; l<L; l++) { // l index for MCMC samples
    C = sigmasq(l) * exp(-obsDistMat/phi(l)) + tausq(l) * eye(N,N);
    Cinv = inv(C);
    resid = y - obsMu.col(l);
    for(unsigned int i=0; i<N0; i++) { // i index for prediction locations
      thisXdist = pred2ObsDistMat.row(i);
      c0 = sigmasq(l) * exp(-thisXdist/phi(l));
      Weight = c0 * Cinv;
      m = as_scalar(Weight * resid);
      v = as_scalar(sigmasq(l) + tausq(l) - Weight * c0.t());
      ypredsamples(i,l) = predMu(i,l) + m + pow(v, 0.5) * randn();
    }
    // Monitoring the processes
    if(l%iterprint==0){
      Rcout << "Prediction upto the " << l <<"th MCMC sample is completed" << std::endl;
    }
  }
  Rcout << std::endl;
  return ypredsamples;
}
