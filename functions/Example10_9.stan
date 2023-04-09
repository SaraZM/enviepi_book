
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
  
  Sigma = sigma_sq * exp(-dist_matrix/ phi) + tau_sq *diag_matrix(rep_vector(1, N));


  for(i in 1:N) {
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
