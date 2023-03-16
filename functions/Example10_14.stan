data {
  int<lower=1> n; //  number of children
  int<lower=1> k; // number of covariates
  int N; //  total number of villages
  int y[n]; //  tests
  matrix[N,N] dist_matrix; //  distance matrix
  matrix[n, k] X; // matrix with covariates
  int index_village[n];
}

parameters {
  real<lower=0> phi_inv;
  real<lower=0> tau;
  real<lower=0> sigma;
  vector[N] S; // spatial random effect
  vector[k] betas;
  real beta0;
}

transformed parameters {
  // vector[n] p = inv_logit(logit_p);
  real sigma_sq = square(sigma);
  real tau_sq = square(tau);
  real<lower=0> phi = 1/phi_inv;
}

model {
  matrix[N, N] L;
  matrix[N, N] Sigma;
  vector[N] zeroes; // vector of zeroes
  vector[n] logit_p;

 for(i in 1:(N-1)){
   for(j in (i+1):N){
     //Sigma[i,j] = sigma_sq*exp(-dist_matrix[i,j]/phi);
     Sigma[i,j] = sigma_sq*(1 + (sqrt(3)*dist_matrix[i,j])/phi) * exp(-sqrt(3)*dist_matrix[i,j] / phi);
     Sigma[j,i] = Sigma[i,j];
   }
 }
   // diagonal elements covariances
for(i in 1:N){
  Sigma[i,i] = sigma_sq + tau_sq;
}

// sample spatial random effect
  L = cholesky_decompose(Sigma);
  zeroes = rep_vector(0, N);
  S ~ multi_normal_cholesky(zeroes, L);

  for(i in 1:n) {
    logit_p[i] = beta0 + X[i,]*betas + S[index_village[i]];
    y[i] ~ binomial_logit(1, logit_p[i]);
  }

  beta0 ~ normal(0,10);
  betas ~ normal(0,10);
  phi_inv ~ gamma(5, 5);
  tau ~ cauchy(0,1);
  sigma ~ cauchy(0,1);
  

}
