functions {
    real uous_dlm_ldensity(array[] real y, array[] real F, real G, real V, real W, real m0, real C0, int T){
    array[T+1] real a;
    array[T+1] real R;
    array[T] real lldata;
    a[1] = m0;
    R[1] = C0;
    for (i in 1:T) {
      real u;
      real Q;
      real A;
      real L;
      u = y[i] - F[i] * a[i];
      Q = F[i] * R[i] * F[i] + V;
      A = G * R[i] * F[i] * inv(Q);
      L = G - A * F[i];
      lldata[i] = normal_lpdf(u | 0, sqrt(Q));
      a[i+1] = G * a[i] + A * u;
      R[i+1] = G * R[i] * L + W;
    }
  return sum(lldata);
  }
  
  
  array[] real uous_ffbs_rng(array[] real y, array[] real F, real G, real V, real W, real m0, real C0, int T){
    
    array[T] real theta;
    array[T] real a;
    array[T] real R;
    array[T] real m;
    array[T] real C;
    
    // Kalman filtering
    
    real mt = m0;
    real Ct = C0;
    
    for(i in 1:T){
      
      real ft;
      real Qt;
      real at;
      real Rt;
      real At;
      
      at = G * mt;
      Rt = G * Ct * G + W;
      ft = F[i] * at;
      Qt = F[i] * Rt * F[i] + V;
      At = Rt * F[i] * inv(Qt);
      mt = at + At * (y[i] - ft);
      Ct = Rt - At * Qt * At;
      
      //store for backward sampling
      a[i] = at;
      R[i] = Rt;
      m[i] = mt;
      C[i] = Ct;
  }
    // backward sampling
    array[T-1] int ind = sort_indices_desc(linspaced_int_array(T-1,1,T-1));
    theta[T] = normal_rng(m[T], sqrt(C[T]));
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
}
  

data{
  int T;
  array[T] real y;
  array[T] real F;
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
  target += uous_dlm_ldensity(y, F, G, V, W, m0, C0, T);
}

generated quantities{
  array[T] real theta;
  real V = square(tau);
  real W = square(sqrt_W);
  theta = uous_ffbs_rng(y, F, G, V, W, m0, C0, T);
}
