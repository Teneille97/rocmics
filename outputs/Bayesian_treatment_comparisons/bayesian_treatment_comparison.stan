
data {
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> treatment;
  vector[N] y;
}

parameters {
  vector[K] mu;
  real<lower=0> sigma;
  real<lower=2> nu;
}

model {
  mu ~ normal(0, 1);
  sigma ~ student_t(3, 0, 1);
  nu ~ gamma(2, 0.1);

  for (n in 1:N) {
    y[n] ~ student_t(nu, mu[treatment[n]], sigma);
  }
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;

  for (n in 1:N) {
    y_rep[n] = student_t_rng(nu, mu[treatment[n]], sigma);
    log_lik[n] = student_t_lpdf(y[n] | nu, mu[treatment[n]], sigma);
  }
}

