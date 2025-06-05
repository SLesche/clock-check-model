data {
  int<lower=1> N;                     // number of observations
  int<lower=1> J;                     // number of participants
  int<lower=1, upper=J> subj[N];      // subject index for each observation
  vector[N] rt;                       // response times
  vector[N] t_target;                // target times
}

parameters {
  // Group-level (population) parameters
  real<lower=0, upper=1> mu_theta;
  real<lower=0> sigma_theta;

  real<lower=0> mu_sigma;
  real<lower=0> sigma_sigma;

  // Subject-level parameters
  vector<lower=0, upper=1>[J] theta;
  vector<lower=0>[J] sigma;
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

  mu_sigma ~ normal(0.5, 0.5);
  sigma_sigma ~ normal(0, 0.5);

  // Subject-level priors
  theta ~ normal(mu_theta, sigma_theta);
  sigma ~ lognormal(log(mu_sigma), sigma_sigma);

  // Likelihood
  for (i in 1:N) {
    target += -0.5 * log(2 * pi())
              - 1.5 * log(rt[i])
              + 0.5 * log(lambda[i])
              - (lambda[i] * square(rt[i] - mu[i])) / (2 * mu[i]^2 * rt[i]);
  }
}
