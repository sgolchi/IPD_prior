data{
  int j;
  int C;
  int S;
  int N;
  int Nt;
  real<lower = 0> y[N];
  real<lower = 0, upper = 1> w[N];
  real<lower = 0> y_t[Nt];
  int<lower = 0, upper = 1> nu[N];
  int<lower = 0, upper = 1> nu_t[Nt];
  matrix[N, C] X;
  matrix[Nt, C] X_t;
  matrix[N, S] Z;
}
parameters {
  real<lower = 0> alpha; 
  real delta;
  real theta; 
  vector[C] beta;
}
model {
  alpha ~ normal(1, 1000);
  beta ~ normal(0, 1000);
  theta ~ normal(0, 1000);
  delta ~ normal(0, 1000);
    for (n in 1:N) {
      real lambda = X[n]*beta + theta;
      if (nu[n]==1)
          target += w[n] * weibull_lpdf(y[n] | alpha, exp(-(lambda)/alpha));
      else
          target += w[n] * weibull_lccdf(y[n] | alpha, exp(-(lambda)/alpha));
  
  }
  for (n in 1:Nt) {
      real lambda = X_t[n]*beta + theta + delta;
      if (nu_t[n]==1)
          target += weibull_lpdf(y_t[n] | alpha, exp(-(lambda)/alpha));
      else
          target += weibull_lccdf(y_t[n] | alpha, exp(-(lambda)/alpha));
  
  }
}