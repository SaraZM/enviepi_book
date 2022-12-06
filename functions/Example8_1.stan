data {
  int<lower=0> N;
  real<lower=0> E[N]; // need to indicate that variable is strictly positive
  int<lower=0> Y[N];
  real<lower=0>lambda_a;
  real<lower=0>lambda_b;
}

parameters {
  real<lower=0> theta[N];
  real<lower=0> a;
  real<lower=0> b;
}

transformed parameters{
  real <lower=0> mu[N];
  
  for(i in 1:N){
    mu[i]=E[i]*theta[i];
  }
  
}

model {
  // likelihood function and prior for theta
  for(i in 1:N){
    Y[i] ~ poisson(mu[i]);
    theta[i]~gamma(a,b);
  }
  a~exponential(lambda_a);
  b~exponential(lambda_b);
}

generated quantities {
  vector [N] log_lik;
  int<lower=0> yfit [N];
  
  //computing the log_likelihood for each value of the mean mu and the fitted values
  for(i in 1:N){
    log_lik[i]=poisson_lpmf(Y[i] |mu[i]);
    yfit[i]=poisson_rng(mu[i]);
  }
  
}
