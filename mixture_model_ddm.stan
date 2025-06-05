functions {
  real inv_gaussian_lpdf(real x, real mu, real lambda) {
    return -0.5 * log(2 * pi())
           - 1.5 * log(x)
           + 0.5 * log(lambda)
           - (lambda * square(x - mu)) / (2 * mu^2 * x);
  }

  real rt_mixture_lpdf(real x, real mu, real lambda, real t_target, real g) {
    real log_ddm = inv_gaussian_lpdf(x | mu, lambda);
    real log_uniform = (x >= 0 && x <= t_target) ? -log(t_target) : negative_infinity();
    return log_mix(g, log_uniform, log_ddm);
  }
}

data {
  int<lower=1> N;                 // number of observations
  vector[N] rt;                   // response times
  vector[N] t_target;             // target times
}

parameters {
  real<lower=0, upper=1> theta;   // threshold multiplier
  real<lower=0> sigma;            // noise std dev
  real<lower=0, upper=1> g;       // guessing rate
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
  g ~ beta(1, 5);  // favor low guess rates

  // Likelihood
  for (i in 1:N) {
    target += rt_mixture_lpdf(rt[i] | mu[i], lambda[i], t_target[i], g);
  }
}
