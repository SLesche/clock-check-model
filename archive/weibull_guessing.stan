functions {
  real guessing_weibull_lpdf(vector t, real b, real k, real g) {
  // Compute hazard and cumulative hazard using element-wise operations
  vector[num_elements(t)] hazard = g + (1 - g) * b * k .* t .^ (k - 1);
  vector[num_elements(t)] cum_hazard = g * t + (1 - g) * b .* t .^ k;

  // Compute log-likelihood sum
  return sum(log(hazard) - cum_hazard);
}

}

data {
  int<lower=1> N;                  // Number of observations
  vector<lower=0>[N] clock_check_time; // Observation times
  // vector<lower=0>[N] known_t_to_target; // Time till target known by individuals
}

parameters {
  real<lower=0, upper=1> g;       // Mixture weight
  real<lower=0> b;                // Scale parameter
  real<lower=0> k;                // Shape parameter (should be positive)
}

model {
  // Priors
  g ~ beta(2, 2);
  b ~ normal(1, 2);
  k ~ normal(1, 2);

  // Likelihood
  target += guessing_weibull_lpdf(clock_check_time | b, k, g);
}
