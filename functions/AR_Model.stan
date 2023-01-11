// Autoregressive model
data {
  int<lower=0> N;
  vector[N] y;
  real m0;
  real<lower=0> C0;
}
parameters {
  real beta;
  real<lower=0> sigma; // standard deviation
  real theta0;
  vector[N] theta;
  real<lower=0> W;
}
model {
  
  beta ~ uniform(-1, 1);
  theta0 ~ normal(m0, C0);
  theta[1] ~ normal(theta0, W);
  
  for (n in 2:N){
    theta[n] ~ normal(theta[n-1], W);
    y[n] ~ normal(theta[n] + beta * y[n-1], sigma);
  }
}
generated quantities{
  vector[N] log_lik;
  vector[N] yfit;
  
  yfit[1] = y[1];
  
  for (n in 2:N){
    log_lik[n] = normal_lpdf(y[n] | theta[n-1] + beta * y[n-1], sigma);
    yfit[n]  = normal_rng (theta[n-1] + beta * yfit[n-1], sigma);
  }

}
