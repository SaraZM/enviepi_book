// Autoregressive model
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N,p] X;
  vector[N] y;
  real m0;
  real<lower=0> C0;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma; // standard deviation
  real theta0;
  vector[N] theta;
  real<lower=0> W;
}
model {
  // Priors
  beta ~ normal(0, 10);
  theta0 ~ normal(m0, C0);
  theta[1] ~ normal(theta0, W);
  y[1] ~ normal(theta[1] +  X[1,]*beta, sigma);
  
  for (n in 2:N){
    theta[n] ~ normal(theta[n-1], W);
    y[n] ~ normal(theta[n] + X[n,]*beta, sigma);
  }
}
generated quantities{
  vector[N] log_lik;
  vector[N] yfit;
  
  yfit[1] = normal_rng (theta[1] + X[1,]*beta, sigma);
  
  for (n in 2:N){
    log_lik[n] = normal_lpdf(y[n] | theta[n-1] +  X[n,]*beta, sigma);
    yfit[n]  = normal_rng (theta[n-1] + X[n,]*beta, sigma);
  }

}
