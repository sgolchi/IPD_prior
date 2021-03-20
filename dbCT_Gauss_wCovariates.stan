data{
  int N;
  int S;
  real y[N];
  real x[N];
  real z[N];
  matrix[N, S] Z;
}
parameters {
  real mu;
  real<lower = 0> tau;
  real theta; 
  real beta;
  real<lower = 0> sig;
  vector[S] beta0;
}
model {
  theta ~ normal(0, 1000);
  beta ~ normal(0, 1000);
  beta0 ~ normal(mu, tau);
  mu ~ normal(0, 1000); 
  tau ~ cauchy(0, 1000);
  sig ~ cauchy(0, 1000);
    for (n in 1:N) {
      target += normal_lpdf(y[n] | Z[n]*beta0 + z[n]*theta + beta*x[n], sig);
    }
}
