data{
  int N;
  int Nw;
  real y[N];
  real x[N];
  real z[N];
  real<lower = 0, upper = 1> w[Nw];
}
parameters {
  real theta; 
  real beta;
  real beta0;
  real<lower = 0> sig;
}
model {
  theta ~ normal(0, 1000);
  beta ~ normal(0, 1000);
  beta0 ~ normal(0, 1000);
  sig ~ cauchy(0, 1000);
  for (n in 1:N)  target += w[n] * normal_lpdf(y[n] | beta0 + theta*z[n] + beta*x[n], sig);
    }
