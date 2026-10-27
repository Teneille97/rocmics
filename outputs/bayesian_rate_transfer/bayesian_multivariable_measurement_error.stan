
data {
  int<lower=2> N;
  int<lower=1> P;
  vector[N] y;
  vector<lower=0>[N] y_se;
  matrix[N, P] x_obs;
  matrix<lower=0>[N, P] x_se;
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> sigma;
  matrix[N, P] z_x;
}
transformed parameters {
  matrix[N, P] x_true;
  for (n in 1:N)
    for (p in 1:P)
      x_true[n, p] = x_obs[n, p] + x_se[n, p] * z_x[n, p];
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  sigma ~ exponential(1);
  to_vector(z_x) ~ std_normal();
  for (n in 1:N)
    y[n] ~ normal(alpha + x_true[n] * beta,
                  sqrt(square(y_se[n]) + square(sigma)));
}
generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  for (n in 1:N) {
    mu[n] = alpha + x_true[n] * beta;
    y_rep[n] = normal_rng(mu[n], sqrt(square(y_se[n]) + square(sigma)));
  }
}

