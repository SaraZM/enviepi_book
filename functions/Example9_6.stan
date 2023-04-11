data {
  int<lower=1> N;
  int<lower=1> N_edges;
  int<lower=1> p;
  matrix[N,p] X;
  int<lower=1, upper=N> d[N];  // vector of number of neighbours
  int<lower=1, upper=N> node2[N_edges];  // long vector of neighbours
  int<lower=0, upper=N_edges> C[N+1]; // cumulative number of neighbours
  int<lower=0> y[N];              
  vector<lower=0>[N] E;          
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            
  vector[p] beta;
  real<lower=0> sigma;     
  real<lower=0, upper=1> rho; 
  vector[N] s;        
}

transformed parameters {
  vector[N_edges] Ws; // weighted spatial neighbours
  vector[N] mean_s;
  vector[N] sigma_s;
  
  for (j in 1:N_edges)
    Ws[j] = s[node2[j]];
  
  for (i in 1:N){
    mean_s[i] = (rho/(1-rho+rho*d[i]))*sum(Ws[(C[i]+1):(C[i+1])]);
    sigma_s[i] = sigma/sqrt(1-rho+rho*d[i]);
  }
}

model {
  y ~ poisson_log(log_E + beta0 + X*beta + s);
  
  for(j in 1:p){
    beta[j] ~ normal(0.0, 10.0);
  }

  beta0 ~ normal(0.0, 10.0);
  sigma ~ normal(0.0,1.0);
  rho ~ uniform(0.0,1.0);
  for (i in 1:N){
    s[i] ~ normal(mean_s[i], sigma_s[i]);
  }
  // sum(s) ~ normal(0.0, N*0.001); // soft constraint
}

generated quantities {
  vector[N] mu_log = log_E + beta0 + X*beta + s;
  vector[N] lik ;
  
  for(i in 1:N){
    lik[i]= exp(poisson_log_lpmf(y[i] | mu_log[i]));
  }
}

