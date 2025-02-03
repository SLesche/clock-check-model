data {
  int<lower=1> N;                 // Number of observations
  int<lower=1> S;                 // Number of subjects
  int<lower=1> subject[N];         // Subject ID for each observation

  vector<lower=0>[N] clock_check_time;  // Observed times
  
  int<lower=1> K_k;
  int<lower=1> K_lambda;
  
  matrix[N, K_k] X_k;               // Covariates for k
  matrix[N, K_lambda] X_lambda;      // Covariates for lambda
}

parameters {
  // Population-level (fixed) effects
  vector[K_k] k_pop;
  vector[K_lambda] lambda_pop;

  // Subject-level (random) effects
  vector[S] k_raw;
  vector[S] lambda_raw;

  // Hyperparameters for subject-level variation
  real<lower=0> sigma_k;
  real<lower=0> sigma_lambda;
}

transformed parameters {
  vector<lower=0>[N] k_pred;
  vector<lower=0>[N] lambda_pred;

  // Compute subject-specific k and lambda
  vector[S] k_subject = k_pop[1] + sigma_k * k_raw;
  vector[S] lambda_subject = lambda_pop[1] + sigma_lambda * lambda_raw;

  // Compute predicted k and lambda for each observation
  for (n in 1:N) {
    k_pred[n] = exp(X_k[n] * k_pop + k_subject[subject[n]]);
    lambda_pred[n] = exp(X_lambda[n] * lambda_pop + lambda_subject[subject[n]]);
  }
}

model {
  // Priors
  k_pop ~ normal(0, 1);
  lambda_pop ~ normal(0, 1);

  k_raw ~ normal(0, 0.1);
  lambda_raw ~ normal(0, 0.1);

  sigma_k ~ normal(0, 1);
  sigma_lambda ~ normal(0, 0.1);

  // Weibull likelihood
  target += weibull_lpdf(clock_check_time | k_pred, lambda_pred);
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N) {
    y_rep[n] = weibull_rng(k_pred[n], lambda_pred[n]);
  }
}
