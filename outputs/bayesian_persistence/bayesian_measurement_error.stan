
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  vector<lower=0>[N] y_se;
  matrix[N, K] x_obs;
  matrix<lower=0>[N, K] x_se;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
  matrix[N, K] z_x;
}
transformed parameters {
  matrix[N, K] x_true;

  for (n in 1:N) {
    for (k in 1:K) {
      x_true[n, k] = x_obs[n, k] + x_se[n, k] * z_x[n, k];
    }
  }
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  sigma ~ exponential(1);

  to_vector(z_x) ~ std_normal();

  for (n in 1:N) {
    y[n] ~ normal(
      alpha + dot_product(beta, to_vector(x_true[n])),
      sqrt(square(y_se[n]) + square(sigma))
    );
  }
}
generated quantities {
  vector[N] mu;
  vector[N] y_rep;

  for (n in 1:N) {
    mu[n] = alpha + dot_product(beta, to_vector(x_true[n]));
    y_rep[n] = normal_rng(
      mu[n],
      sqrt(square(y_se[n]) + square(sigma))
    );
  }
}

