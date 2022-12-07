// Example10_9Code <-  nimbleCode ({
//   # Covariance matrix spatial effect
//   Sigma[1:n, 1:n] <-  (sigma ^ 2) * exp(-obs_dist_mat[1:n, 1:n] / phi)
//   m[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n, 1:n])
//   
//   for (site in 1:n) {
//     mean.site[site] <- beta0 + m[site]
//     y[site] ~ dnorm(mean.site[site], sd = tau)
//   }
//   
//   # Set up the priors for the spatial model
//   sigma ~ T(dt(mu = 0, sigma = 1, df = 1), 0, Inf)
//   tau ~ T(dt(mu = 0, sigma = 1, df = 1), 0, Inf)
//   phi ~ dexp(127.72)
//   
//   # prior for the intercept term
//   beta0 ~ dnorm (0 , 10)
//   
// }) 


data {
  int<lower=1> N;
  vector[N] y;
  matrix[N,N] dist_matrix;
}
transformed data {
  vector[N] mu = rep_vector(0, N);
}
parameters {
  real<lower=0> phi;
  real<lower=0> tau;
  real<lower=0> sigma;
}
model {
  matrix[N, N] L;
  matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  real sq_sigma = square(sigma);

  // diagonal elements
  for (n in 1:N) {
    K[n, n] = K[n, n] + sq_sigma;
  }

  L = cholesky_decompose(K);

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();

  y ~ multi_normal_cholesky(mu, L);
}