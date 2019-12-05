data{
  int C;
  int N;
  real<lower = 0> y[N];
  int<lower = 0, upper = 1> nu[N];
  real x1[N];
  int x2[N,C-1];
}
parameters {
  real<lower = 0> theta1; 
  real<lower = 0> theta2; 
  real beta1;
  vector[C-1] beta2;
  real mu_x;
  real<lower = 0, upper = 1> p_x[C-1];
  real<lower = 0, upper = 1> p_nu;
  real<lower = 0> sig_x;
  real<lower = 0> sig_y;
}
model {
  theta1 ~ normal(0, 1000);
  beta1 ~ normal(0, 1000);
  beta2 ~ normal(0, 1000);
  theta2 ~ normal(0, 1000);
  mu_x ~ normal(0, 1000); 
  p_x ~ beta(1, 1);
  sig_x ~ cauchy(0, 1000);
  sig_y ~ cauchy(0, 1000);
    for (n in 1:N) {
      real lambda;
      real temp[C-1];
      x1[n] ~ normal(mu_x, sig_x);
      nu[n] ~ bernoulli(p_nu);
      for (j in 1:(C-1)) {
        x2[n,j] ~ bernoulli(p_x[j]);
        temp[j] = beta2[j]*x2[n,j];
      }
      lambda = beta1*x1[n] + sum(temp);
      if (nu[n]==1)
          target += normal_lpdf(y[n] | lambda + theta1, sig_y);
      else
          target += normal_lpdf(y[n] | lambda + theta2, sig_y);
  }
}