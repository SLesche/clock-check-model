data {
  int<lower=1> N;                      // Number of observations
  int<lower=1> S;                      // Number of subjects
  int<lower=1, upper=S> subject_id[N];  // Subject indices for each observation
  vector<lower=0>[N] clock_check_time;  // Observation times

  // Level 2 (Subject-Level) Predictors
  int<lower=1> P_alpha;                 // Num of predictors for alpha (subject-level)
  matrix[S, P_alpha] X_alpha;           // Subject-level predictor matrix for alpha

  int<lower=1> P_sigma;                 // Num of predictors for sigma (subject-level)
  matrix[S, P_sigma] X_sigma;           // Subject-level predictor matrix for sigma

  // Level 1 (Observation-Level) Predictors
  int<lower=1> P_alpha_obs;             // Num of predictors for alpha (obs-level)
  matrix[N, P_alpha_obs] X_alpha_obs;   // Observation-level predictor matrix for alpha

  int<lower=1> P_sigma_obs;             // Num of predictors for sigma (obs-level)
  matrix[N, P_sigma_obs] X_sigma_obs;   // Observation-level predictor matrix for sigma
}

parameters {
  // Regression coefficients for alpha
  vector[P_alpha] beta_alpha;           // Subject-level effects
  vector[P_alpha_obs] gamma_alpha;      // Observation-level effects

  // Regression coefficients for sigma
  vector[P_sigma] beta_sigma;           // Subject-level effects
  vector[P_sigma_obs] gamma_sigma;      // Observation-level effects

  // Random effects (subject-level residuals)
  real<lower=0> alpha_sd;               // Std dev for alpha random effects
  real<lower=0> sigma_sd;               // Std dev for sigma random effects
  vector[S] alpha_raw;                  // Subject-specific residuals for alpha
  vector[S] sigma_raw;                  // Subject-specific residuals for sigma
}

transformed parameters {
  vector<lower=0>[N] alpha_obs;  // Final alpha for each observation
  vector<lower=0>[N] sigma_obs;  // Final sigma for each observation

  for (n in 1:N) {
    int s = subject_id[n];  // Get subject index for this observation
    alpha_obs[n] = exp(X_alpha[s] * beta_alpha + X_alpha_obs[n] * gamma_alpha + alpha_sd * alpha_raw[s]);
    sigma_obs[n] = exp(X_sigma[s] * beta_sigma + X_sigma_obs[n] * gamma_sigma + sigma_sd * sigma_raw[s]);
  }
}

model {
  // Priors on regression coefficients
  beta_alpha ~ normal(0, 1);
  gamma_alpha ~ normal(0, 1);
  beta_sigma ~ normal(0, 1);
  gamma_sigma ~ normal(0, 1);

  // Priors on random effects
  alpha_raw ~ normal(0, 1);
  sigma_raw ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    target += weibull_lpdf(clock_check_time[n] | alpha_obs[n], sigma_obs[n]);
  }
}
