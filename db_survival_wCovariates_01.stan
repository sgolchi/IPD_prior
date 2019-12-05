data{
  int j;
  int C;
  int S;
  int N;
  int Nt;
  real<lower = 0> y[N];
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
  vector[S] theta; 
  vector[C] beta;
  real mu;
  real<lower = 0> tau;
}
model {
  alpha ~ gamma(1, 1);
  beta ~ normal(0, 100);
  theta ~ normal(mu, tau);
  delta ~ normal(0, 100);
  mu ~ normal(0, 100); 
  tau ~ normal(0, 1000);
    for (n in 1:N) {
      real lambda = X[n]*beta + Z[n]*theta;
      if (nu[n]==1)
          y[n] ~ weibull(alpha, exp(-(lambda)/alpha));
      else
          target += weibull_lccdf(y[n] | alpha, exp(-(lambda)/alpha));
  }
  for (n in 1:Nt) {
      real lambda = X_t[n]*beta + theta[j] + delta;
      if (nu_t[n]==1)
          y_t[n] ~ weibull(alpha, exp(-(lambda)/alpha));
      else
          target += weibull_lccdf(y_t[n] | alpha, exp(-(lambda)/alpha));
  
  }
}