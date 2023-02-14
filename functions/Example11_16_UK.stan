functions {
  
  // FFBS for DLM with univariate observation equation and univariate system equation
  array[] vector uoms_ffbs_rng(array[] real y, array[] vector F, matrix G, real V, matrix W, vector m0, matrix C0, int T, int p){
    
    array[T] vector[p] theta;
    array[T] vector[p] a;
    array[T] matrix[p,p] R;
    array[T] vector[p] m;
    array[T] matrix[p,p] C;
    
    // Kalman filtering
    
    vector[p] mt = m0;
    matrix[p,p] Ct = C0;
    
    for(i in 1:T){
      
      real ft;
      real Qt;
      vector[p] at;
      matrix[p,p] Rt;
      vector[p] At;
      
      at = G * mt;
      Rt = G * Ct * G' + W;
      ft = F[i]' * at;
      Qt = quad_form(Rt, F[i]) + V; //F[i]' * Rt * F[i] + V;
      At = Rt * F[i] * inv(Qt);
      mt = at + At * (y[i] - ft);
      Ct = Rt - At * Qt * At';
      
      //store for backward sampling
      a[i] = at;
      R[i] = Rt;
      m[i] = mt;
      C[i] = Ct;
  }
    // backward sampling
    array[T-1] int ind = sort_indices_desc(linspaced_int_array(T-1,1,T-1));
    theta[T] = multi_normal_rng(m[T], C[T]);
    for(i in ind) {
      matrix[p,p] Bt;
      vector[p] ht;
      matrix[p,p] Ht;
      Bt = C[i] * G' * inverse(R[i+1]);
      ht = m[i] + Bt * (theta[i+1] - a[i+1]);
      Ht = C[i] - Bt * R[i+1] * Bt';
      theta[i] = multi_normal_rng(ht, Ht);
    }
  return theta;
}

  real uoms_dlm_ldensity(array[] real y, array[] vector F, matrix G, real V, matrix W, vector m0, matrix C0, int T, int p){
    array[T+1] vector[p] a;
    array[T+1] matrix[p, p] R;
    array[T] real lldata;
    a[1] = m0;
    R[1] = C0;
    for (i in 1:T) {
      real u;
      real Q;
      real Qinv;
      vector[p] A;
      matrix[p, p] L;
      u = y[i] - F[i]' * a[i];
      Q = quad_form(R[i],F[i]) + V; //F[i]' * R[i] * F[i] + V;
      Qinv = inv(Q); //
      A = G * R[i] * F[i] * Qinv;
      L = G - A * F[i]';
      //lldata[i] = -0.5 * (log(2 * pi()) + log(Q) + Qinv*square(u));
      lldata[i] = normal_lpdf(u | 0, sqrt(Q)); // univariate
      a[i+1] = G * a[i] + A * u;
      R[i+1] = G * R[i] * L' + W;
    }
  return sum(lldata);
  }
  
  array[] real uoms_dlm_one_step_ahead_rng(array[] real y, array[] vector F, matrix G, real V, matrix W, vector m0, matrix C0, int T, int p){
    array[T] real yfit;
    array[T+1] vector[p] a;
    array[T+1] matrix[p, p] R;
    array[T] real lldata;
    a[1] = m0;
    R[1] = C0;
    for (i in 1:T) {
      real u;
      real Q;
      real Qinv;
      vector[p] A;
      matrix[p, p] L;
      u = y[i] - F[i]' * a[i];
      Q = quad_form(R[i],F[i]) + V; //F[i]' * R[i] * F[i] + V;
      Qinv = inv(Q); //
      A = G * R[i] * F[i] * Qinv;
      L = G - A * F[i]';
      yfit[i] = normal_rng(F[i]' * a[i], sqrt(Q)); // univariate
      a[i+1] = G * a[i] + A * u;
      R[i+1] = G * R[i] * L' + W;
    }
  return yfit;
  }
}

data{
  int T;
  int p;
  array[T] real y;
  array[T] vector[p] F;
  matrix[p, p] G;
  vector[p] m0;
  cov_matrix[p] C0;
}

parameters{
  real<lower=0> tau;
  vector<lower=0>[p] sqrt_W_diag;
}

model {
  real V = square(tau);
  matrix[p, p] W = diag_matrix(square(sqrt_W_diag));
  tau ~ std_normal();
  sqrt_W_diag ~ std_normal();
  target += uoms_dlm_ldensity(y, F, G, V, W, m0, C0, T, p);
}

generated quantities{
  array[T] vector[p] theta;
  array[T] real yfit;
  real V = square(tau);
  matrix[p, p] W = diag_matrix(square(sqrt_W_diag));
  theta = uoms_ffbs_rng(y, F, G, V, W, m0, C0, T, p);
  yfit = uoms_dlm_one_step_ahead_rng(y, F, G, V, W, m0, C0, T, p);
}

