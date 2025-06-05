data {
  int<lower=1> N;                 // number of observations
  vector[N] rt;                   // response times
  vector[N] t_target;            // target times (fixed)
}

parameters {
  real<lower=0, upper = 1> theta;            // threshold multiplier
  real<lower=0> sigma;            // noise std dev
}

transformed parameters {
  vector[N] mu;
  vector[N] lambda;

  for (i in 1:N) {
    mu[i] = theta * t_target[i];
    lambda[i] = square(mu[i]) / square(sigma);
  }
}

model {
  // Priors
  theta ~ beta(2, 2);
  sigma ~ gamma(1, 1);

  // Custom Inverse Gaussian log-likelihood
  for (i in 1:N) {
    target += -0.5 * log(2 * pi())
              - 1.5 * log(rt[i])
              + 0.5 * log(lambda[i])
              - (lambda[i] * square(rt[i] - mu[i])) / (2 * mu[i]^2 * rt[i]);
  }
}

