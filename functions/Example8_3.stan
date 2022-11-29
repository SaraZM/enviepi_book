data {
  int<lower=1> N;
  int<lower=1> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            // intercept
  vector[N] s;         // spatial effects
  real<lower=0> sigma_s;        // marginal standard deviation of spatial effects
}

transformed parameters {
  vector[N] b; // latent effect
  b =  sigma_s*s;
}

model {
  y ~ poisson_log(log_E + beta0 + b); 
  // This is the prior for s! (up to proportionality)
  target += -0.5 * dot_self(s[node1] - s[node2]);
  sum(s) ~ normal(0, 0.001 * N); 
  
  beta0 ~ normal(0.0, 10.0);
  sigma_s ~ normal(0.0,1.0);
}

generated quantities {
  vector[N] mu=exp(log_E + beta0 + b);
  vector[N] lik;
  
  for(i in 1:N){
    lik[i] = exp(poisson_lpmf(y[i] | mu[i] ));
  }
}
