
data {
  int<lower=1> N;
  int<lower=0> p;
  vector[N] y;
  matrix[N,N] dist_matrix;
  matrix[N,p] X; 
}

parameters {
  real<lower=0> phi;
  real<lower=0> tau;
  real<lower=0> sigma;
  vector[p] beta;
  real beta0;

}
transformed parameters{
  real<lower=0> sigma_sq = square(sigma);
  real<lower=0> tau_sq = square(tau);
}
model {
  vector[N] mu;
  matrix[N, N] L;
  matrix[N, N] Sigma;

 for(i in 1:(N-1)){
   for(j in (i+1):N){
     Sigma[i,j] = sigma_sq*exp(-dist_matrix[i,j]/phi);
     Sigma[j,i] = Sigma[i,j];
   }
 }
   // diagonal elements

  for(i in 1:N) {
    Sigma[i,i] = sigma_sq + tau_sq;
    mu[i] = beta0 + X[i,]*beta;
  }
  
  L = cholesky_decompose(Sigma);
  
  beta0 ~ normal(0,10);
  beta ~ normal(0,10);
  phi ~ inv_gamma(5, 5);
  tau ~ cauchy(0,1);
  sigma ~ cauchy(0,1);

  y ~ multi_normal_cholesky(mu, L);
}
