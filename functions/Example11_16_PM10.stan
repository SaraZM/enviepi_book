functions {
    real uous_dlm_ldensity(array[] real y, array[] real Ft, real G, real V, real W, real m0, real C0, int Tt){
    array[Tt+1] real a;
    array[Tt+1] real R;
    array[Tt] real lldata;
    a[1] = m0;
    R[1] = C0;
    for (i in 1:Tt) {
      real u;
      real Q;
      real A;
      real L;
      u = y[i] - Ft[i] * a[i];
      Q = Ft[i] * R[i] * Ft[i] + V;
      A = G * R[i] * Ft[i] * inv(Q);
      L = G - A * Ft[i];
      lldata[i] = normal_lpdf(u | 0, sqrt(Q));
      a[i+1] = G * a[i] + A * u;
      R[i+1] = G * R[i] * L + W;
    }
  return sum(lldata);
  }
  
  
  array[] real uous_ffbs_rng(array[] real y, array[] real Ft, real G, real V, real W, real m0, real C0, int Tt){
    
    array[Tt] real theta;
    array[Tt] real a;
    array[Tt] real R;
    array[Tt] real m;
    array[Tt] real C;
    
    // Kalman filtering
    
    real mt = m0;
    real Ct = C0;
    
    for(i in 1:Tt){
      
      real ft;
      real Qt;
      real at;
      real Rt;
      real At;
      
      at = G * mt;
      Rt = G * Ct * G + W;
      ft = Ft[i] * at;
      Qt = Ft[i] * Rt * Ft[i] + V;
      At = Rt * Ft[i] * inv(Qt);
      mt = at + At * (y[i] - ft);
      Ct = Rt - At * Qt * At;
      
      //store for backward sampling
      a[i] = at;
      R[i] = Rt;
      m[i] = mt;
      C[i] = Ct;
  }
    // backward sampling
    array[Tt-1] int ind = sort_indices_desc(linspaced_int_array(Tt-1,1,Tt-1));
    theta[Tt] = normal_rng(m[Tt], sqrt(C[Tt]));
    for(i in ind) {
      real Bt;
      real ht;
      real Ht;
      Bt = C[i] * G * inv(R[i+1]);
      ht = m[i] + Bt * (theta[i+1] - a[i+1]);
      Ht = C[i] - Bt * G * C[i];
      theta[i] = normal_rng(ht, sqrt(Ht));
    }
  return theta;
}


  array[] real uous_dlm_one_step_ahead_rng(array[] real y, array[] real Ft, real G, real V, real W, real m0, real C0, int Tt){
    array[Tt] real yfit;
    array[Tt+1] real a;
    array[Tt+1] real R;
    array[Tt] real lldata;
    
    a[1] = m0;
    R[1] = C0;
    
    for (i in 1:Tt) {
      real u;
      real Q;
      real A;
      real L;
      u = y[i] - Ft[i] * a[i];
      Q = Ft[i] * R[i] * Ft[i] + V; //F[i]' * R[i] * F[i] + V;
      A = G * R[i] * Ft[i] * inv(Q);
      L = G - A * Ft[i];
      yfit[i] = normal_rng(Ft[i] * a[i], sqrt(Q)); // univariate
      a[i+1] = G * a[i] + A * u;
      R[i+1] = G * R[i] * L + W;
    }
  return yfit;
  }
  
  
}
  

data{
  int Tt;
  array[Tt] real y;
  array[Tt] real Ft;
  real G;
  real m0;
  real<lower=0> C0;
}

parameters{
  real<lower=0> tau;
  real<lower=0> sqrt_W;
}

model {
  real V = square(tau);
  real W = square(sqrt_W);
  tau ~ std_normal();
  sqrt_W ~ std_normal();
  target += uous_dlm_ldensity(y, Ft, G, V, W, m0, C0, Tt);
}

generated quantities{
  array[Tt] real theta;
  array[Tt] real yfit;

  real V = square(tau);
  real W = square(sqrt_W);
  theta = uous_ffbs_rng(y, Ft, G, V, W, m0, C0, Tt);
  yfit = uous_dlm_one_step_ahead_rng(y, Ft, G, V, W, m0, C0, Tt);

}
