functions {
  real guessing_weibull_lpdf(vector t, real g, vector b, real k) {
    vector[num_elements(t)] hazard = g + (1 - g) .* b .* k .* t .^ (k - 1);
    vector[num_elements(t)] cum_hazard = g .* t + (1 - g) .* b .* t .^ k;
    
    return sum(log(hazard) - cum_hazard);
  }
}

data {
  int<lower=1> N;                  // Number of observations
  vector<lower=0>[N] clock_check_time;  // Observation times
  vector<lower=0>[N] known_t_to_target; // Predictor for g
}

parameters {
  real<lower=0, upper=1> g;  // Intercept for g regression
  real b_0;  // Scale parameter for Weibull
  real b_1;  // Scale parameter for Weibull
  real<lower=0> k;  // Shape parameter for Weibull
}

transformed parameters {
  vector[N] b;  
  b = b_0 + b_1 * known_t_to_target;
}

model {
  // Priors
  g ~ beta(2, 2);
  b_0 ~ normal(1, 2);
  b_1 ~ normal(1, 2);
  k ~ normal(1, 2);

  // Vectorized likelihood
  target += guessing_weibull_lpdf(clock_check_time | g, b, k);
}
