
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  sigma ~ exponential(1);
  y ~ normal(alpha + X * beta, sigma);
}
generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  mu = alpha + X * beta;
  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }
}

