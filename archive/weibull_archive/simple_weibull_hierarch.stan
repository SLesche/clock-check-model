data {
  int<lower=1> N;                  // Number of observations
  int<lower=1> S;                  // Number of subjects
  int<lower=1, upper=S> subject_id[N]; // Subject indices for each observation
  vector<lower=0>[N] clock_check_time; // Observation times
}

parameters {
  real<lower=0> alpha_mean;            // Mean shape parameter for flexible hazard Weibull
  real<lower=0> sigma_mean;            // Mean scale parameter for flexible hazard Weibull
  real<lower=0> alpha_sd;              // Standard deviation for alpha across subjects
  real<lower=0> sigma_sd;              // Standard deviation for sigma across subjects

  vector<lower=0>[S] alpha_subject;   // Subject-level shape parameters for flexible hazard Weibull
  vector<lower=0>[S] sigma_subject;   // Subject-level scale parameters for flexible hazard Weibull
}

model {
  // Priors
  alpha_mean ~ normal(1, 1);          // Prior for mean alpha (shape)
  sigma_mean ~ normal(1, 1);          // Prior for mean sigma (scale)
  alpha_sd ~ normal(0, 1);            // Prior for standard deviation of alpha
  sigma_sd ~ normal(0, 1);            // Prior for standard deviation of sigma
  
  // Subject-level priors
  alpha_subject ~ normal(alpha_mean, alpha_sd); // Subject-level alpha
  sigma_subject ~ normal(sigma_mean, sigma_sd); // Subject-level sigma

  // Likelihood
  for (n in 1:N) {
    target += weibull_lpdf(clock_check_time[n] | alpha_subject[subject_id[n]], sigma_subject[subject_id[n]]);
  }
}
