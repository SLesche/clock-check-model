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
  int<lower=1> J;                 // number of participants
  int<lower=1, upper=J> subj[N];  // subject index for each observation
  vector[N] rt;                   // response times
  vector[N] t_target;             // target times
}

parameters {
  // Group-level (population) parameters for theta
  real<lower=0, upper=1> mu_theta;
  real<lower=0> sigma_theta;

  // Group-level parameters for sigma (noise std dev)
  real<lower=0> mu_sigma;
  real<lower=0> sigma_sigma;

  // Group-level parameters for guess rate g
  real<lower=0, upper=1> mu_g;
  real<lower=0> sigma_g;

  // Subject-level parameters
  vector<lower=0, upper=1>[J] theta;
  vector<lower=0>[J] sigma;
  vector<lower=0, upper=1>[J] g;
}

transformed parameters {
  vector[N] mu;
  vector[N] lambda;

  for (i in 1:N) {
    mu[i] = theta[subj[i]] * t_target[i];
    lambda[i] = square(mu[i]) / square(sigma[subj[i]]);
  }
}

model {
  // Hyperpriors
  mu_theta ~ beta(2, 2);
  sigma_theta ~ normal(0, 0.5);

  mu_sigma ~ gamma(1, 1);
  sigma_sigma ~ normal(0, 0.5);

  mu_g ~ beta(1, 5);
  sigma_g ~ normal(0, 0.5);

  // Subject-level priors
  theta ~ normal(mu_theta, sigma_theta);
  sigma ~ lognormal(log(mu_sigma), sigma_sigma);
  g ~ normal(mu_g, sigma_g);

  // // Truncate g between 0 and 1 (soft truncation using target intervals)
  // for (j in 1:J) {
  //   target += uniform_lpdf(g[j] | 0, 1);
  // }

  // Likelihood
  for (i in 1:N) {
    target += rt_mixture_lpdf(rt[i] | mu[i], lambda[i], t_target[i], g[subj[i]]);
  }
}
