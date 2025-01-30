functions {
  real guessing_weibull_lpdf(real[] t, real b, real k, real g) {
    real log_prob = 0;
    for (i in 1:num_elements(t)) {
      // Hazard function
      real hazard = g + (1 - g) * b * k * t[i]^(k - 1);
      // Cumulative hazard function
      real cum_hazard = g * t[i] + (1 - g) * b * t[i]^k;
      // Log-likelihood
      log_prob += log(hazard) - cum_hazard;
    }
    return log_prob;
  }
}

data {
  int<lower=1> N;                  // Number of observations
  real<lower=0> clock_check_time[N]; // Observation times
  real<lower=0> known_t_to_target[N]; // Time till target known by individuals
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
