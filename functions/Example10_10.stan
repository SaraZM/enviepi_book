data {
  int<lower=1> N;
  int<lower = 0> N_pred;
  int<lower=0> p;
  vector[N] y;
  matrix[N,N] obs_dist_matrix;
  matrix[N,N_pred] obs_pred_dist_matrix;
  matrix[N_pred, N_pred] pred_dist_matrix;
  matrix[N,p] X; 
  matrix[N_pred,p] X_pred; 
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
  matrix[N, N] L;
  matrix[N, N] Sigma;
  vector[N] mu;

  Sigma = sigma_sq * exp(-obs_dist_matrix/ phi) + tau_sq *diag_matrix(rep_vector(1, N));
   // diagonal elements

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

generated quantities {
  matrix[N_pred, N_pred] L;
  matrix[N, N] Sigma_obs;
  matrix[N_pred, N_pred] Sigma_pred;
  matrix[N_pred, N_pred] Cov_pred;
  matrix[N, N_pred] Sigma_obs_pred;
  vector[N] mu;
  vector[N_pred] mu_unobs;
  vector[N_pred] mu_pred;
  vector[N_pred] y_pred;

 Sigma_obs = sigma_sq * exp(-obs_dist_matrix/ phi) + tau_sq *diag_matrix(rep_vector(1, N));
 Sigma_pred = sigma_sq * exp(-pred_dist_matrix/ phi) + tau_sq *diag_matrix(rep_vector(1, N_pred));
 Sigma_obs_pred =  sigma_sq * exp(-obs_pred_dist_matrix/ phi);

 for (j in 1:N_pred) {
    mu_unobs[j] = beta0 + X_pred[j,]*beta;
  }
  
  for(i in 1:N) {
    mu[i] = beta0 + X[i,]*beta;
  }
  
  mu_pred = mu_unobs + Sigma_obs_pred'*inverse(Sigma_obs)*(y - mu);
  Cov_pred = Sigma_pred - Sigma_obs_pred'*inverse(Sigma_obs)*Sigma_obs_pred;
  L = cholesky_decompose(Cov_pred);
  y_pred  = multi_normal_cholesky_rng(mu_pred, L);
  
 
  
}
