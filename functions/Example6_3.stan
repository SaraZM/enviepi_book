data {
  int<lower=0> N;
  vector[N] E; // expected cases?
  vector[N] X1; // what are these covariates?
  vector[N] X2;
  int Y[N] ;
}

parameters {
  real beta0;
  real beta1;
  real beta2;
}
transformed parameters{
  real base = exp(beta0);
  real RR = exp(beta1);
}

model {
  vector[N] mu;
  
  for(i in 1:N){
    mu[i] = log(E[i])+ beta0 + beta1*X1[i] + beta2*X2[i];
    Y[i] ~ poisson_log(mu[i]);
  }

  beta0 ~ normal(0 , 100);
  beta1 ~ normal(0 , 100);
  beta2 ~ normal(0 , 100);
  
}

