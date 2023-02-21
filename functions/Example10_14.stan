data {
  int<lower=1> n; //  number of villages
  vector[n] N; //  total number of surveys
  vector[n] y; //  number of positive cases
  matrix[n,n] dist_matrix; //  distance matrix
  vector[n] X; // vector with altitudes
}

parameters {
  real<lower=0> phi;
  real<lower=0> tau;
  real<lower=0> sigma;
  vector[n] S; // spatial random effect
  vector[n] p; //  probabilities by village
  real beta;
  real beta0;
}

model {
  matrix[n, n] L;
  matrix[n, n] Sigma;
  vector[n] zeroes; // vector of zeroes

  real sigma_sq = square(sigma);
  real tau_sq = square(tau);

 for(i in 1:(n-1)){
   for(j in (i+1):n){
     Sigma[i,j] = sigma_sq*exp(-dist_matrix[i,j]/phi);
     Sigma[j,i] = Sigma[i,j];
   }
 }
   // diagonal elements covariances
for(i in 1:n){
  Sigma[i,i] = sigma_sq + tau_sq;
}

// sample spatial random effect
  // L = cholesky_decompose(Sigma);
  // zeroes = rep_vector(0, n);
  // S ~ multi_normal_cholesky(zeroes, L);
  // 
  for(i in 1:n) {
    logit(p[i]) = beta0 + X[i]*beta + S[i];
    y[i] ~ binomial(N[i], p[i]);
  }
  // 
  // 
  // beta0 ~ normal(0,10);
  // beta ~ normal(0,10);
  // phi ~ exp(rate = 1/45.51057);
  // tau ~ cauchy(0,1);
  // sigma ~ cauchy(0,1);

}

