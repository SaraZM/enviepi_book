// Autoregressive model
data {
  int<lower=0> N;
  int<lower=0> q;
  matrix[N,q] FF;
  matrix[q,q] GG;
  vector[N] y;
  vector[q] m0;
  matrix[q,q] C0;
}
parameters {
  real<lower=0> sigma; // standard deviation
  vector[q] theta0;
  matrix[q, N] theta;
  cov_matrix[q] W;
}
model {
  // Priors
  theta0 ~ multi_normal(m0, C0);
  theta[,1] ~ multi_normal(GG*theta0, W);
  y[1] ~ normal( FF[1,]*theta[,1], sigma);
  
  for (n in 2:N){
    theta[,n] ~ multi_normal(GG*theta[,n-1], W);
    y[n] ~ normal(FF[n,]*theta[,n], sigma);
  }
}
generated quantities{
  vector[N] log_lik;
  vector[N] yfit;

  yfit[1] = normal_rng (FF[1,]*theta[,1], sigma);

  for (n in 2:N){
    log_lik[n] = normal_lpdf(y[n] | FF[n,]*theta[,n-1], sigma);
    yfit[n]  = normal_rng (FF[n,]*theta[,n-1], sigma);
  }

}
