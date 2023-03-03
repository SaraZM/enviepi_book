#include <RcppArmadillo.h>
#include <ComputeDistMat.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]

arma::mat predSumNNGPMatern12(arma::colvec obsY, arma::mat obsMuY, arma::mat prdMuY,
                      arma::mat obsCoords, arma::mat prdNeiDist,
                      arma::umat prdNeiID, arma::vec sigmasq, arma::vec lscale,
                      arma::vec tausq, unsigned int iterprint)
  {
  unsigned int L = prdMuY.n_cols;     // number of posterior samples
  unsigned int N0 = prdMuY.n_rows;   // number of prediction locations
  //int K = prdNeiID.n_cols;   // number of nearest neighbor
  arma::mat ypredsummary = zeros(N0,5); // initialize prediction summary matrix
  arma::vec ypred = zeros(L); // initialize predictions vector
  arma::vec quant = zeros(3); // initialize quadrants vector
  arma::vec prob = {0.025,0.5,0.975}; // desired output probabilities
  
  // arma::mat thisCoords = zeros(K,2); not sure I need this
  // arma::mat thisDist, thisNeiDist, CNei, thisNeiObsMuY; not sure I need this
 // arma::rowvec Weight;
 // arma::rowvec c0, thisXdist; 
//  arma::urowvec thisNeiID;
//  arma::colvec thisNeiObsY;
  double m;
//  arma::umat NeiID = prdNeiID - 1;  // index in C is equal to (index in R - 1)
  
  for(unsigned int i=0; i<N0; i++) { // i index for the number of prediction locations
    thisNeiID = NeiID.row(i);
    thisCoords = obsCoords.rows(thisNeiID);
    thisNeiDist = ComputeDistMat(thisCoords);
    thisXdist = prdNeiDist.row(i);
    thisNeiObsY = obsY(thisNeiID);
    thisNeiObsMuY = obsMuY.rows(thisNeiID);
    
    for(unsigned int l=0; l<L; l++) { // l index for the number of MCMC samples
      
      c0 = sigmasq(l) * exp(-thisXdist/lscale(l));
      CNei = sigmasq(l) * exp(-thisNeiDist/lscale(l)) + tausq(l)*eye(K,K);
      Weight = c0 * inv(CNei);
      
      m = as_scalar(Weight * (thisNeiObsY - thisNeiObsMuY.col(l)));
     
      ypred(l) = prdMuY(i,l) + m + pow(v, 0.5) * randn();
      
      
      }
    quant = quantile(ypred,prob);
    ypredsummary(i,0) = mean(ypred);
    ypredsummary(i,1) = stddev(ypred);
    ypredsummary(i,2) = quant(0);
    ypredsummary(i,3) = quant(1);
    ypredsummary(i,4) = quant(2);
    
    // Monitoring the processes
    if(i%iterprint==0){
      Rcout << "Prediction upto the " << i <<"th locations is completed" << std::endl;
    }
    
  }
  Rcout << std::endl;
  Rcout << std::endl;
  Rcout << "Output is a matrix with columns mean, sd, q2.5, q50, q97.5, respectively." << std::endl;
  return ypredsummary;
}
