functions {
  real guessing_weibull_lpdf(real t, real b, real k, real g) {
    // Compute hazard and cumulative hazard using element-wise operations
    real hazard = g + (1 - g) * b * k .* t .^ (k - 1);
    real cum_hazard = g * t + (1 - g) * b .* t .^ k;

    // Compute log-likelihood sum
    return log(hazard) - cum_hazard;
  }
}

data {
  int<lower=1> N;                  // Number of observations
  int<lower=1> J;                  // Number of subjects
  int<lower=1, upper=J> subject[N]; // Subject index for each observation
  vector<lower=0>[N] clock_check_time; // Observation times
}

parameters {
  // Group-level parameters (hyperparameters)
  real<lower=0, upper=1> g_group;  // Mean of g across subjects
  real<lower=0> b_group;           // Mean of b across subjects
  real<lower=0> k_group;           // Mean of k across subjects

  // Standard deviations of group-level parameters
  real<lower=0> sigma_g;           // Standard deviation of g across subjects
  real<lower=0> sigma_b;           // Standard deviation of b across subjects
  real<lower=0> sigma_k;           // Standard deviation of k across subjects

  // Individual-level parameters (g, b, k for each subject)
  real<lower=0, upper=1> g[J];     // g for each subject
  real<lower=0> b[J];              // b for each subject
  real<lower=0> k[J];              // k for each subject
}

model {
  // Priors for group-level parameters (hyperparameters)
  g_group ~ beta(2, 2);   // Mean for g, prior
  b_group ~ gamma(1, 1);  // Mean for b, prior
  k_group ~ gamma(1, 1);  // Mean for k, prior

  // Priors for standard deviations (hierarchical model)
  sigma_g ~ gamma(1, 4);
  sigma_b ~ gamma(2, 2);
  sigma_k ~ gamma(2, 2);

  // Priors for individual-level parameters (from normal distributions)
  for (j in 1:J) {
    g[j] ~ normal(g_group, sigma_g);  // g for each individual
    b[j] ~ normal(b_group, sigma_b);  // b for each individual
    k[j] ~ normal(k_group, sigma_k);  // k for each individual
  }
  
  // Likelihood: Compute log-likelihood for each individual
  for (i in 1:N) {
    int subject_id = subject[i];  // Get subject index for each observation
    target += guessing_weibull_lpdf(clock_check_time[i] | b[subject_id], k[subject_id], g[subject_id]);
  }
}
