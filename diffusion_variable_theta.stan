functions {
  vector get_predicted_time(vector t_target, real sigma_0, real k, vector alpha) {
    return t_target ./ (1 + alpha * sigma_0 * k); // Vectorized element-wise division
  }
}

data {
  int<lower=0> N; // Number of checking events
  real<lower=0> known_t_to_target[N]; // The time of target
  real<lower=0> observed_time[N]; // The time of the actual clock check
}

parameters {
  real<lower=0> k; // Noise scalar k
  real<lower=0> sigma_0; // Initial noise at t = 0
  real<lower=0, upper=1> base_theta; // Baseline threshold (capped at 0.5)
  real<lower=0> growth_rate; // Growth rate for threshold
  real<lower=0> sigma_err; // Error level
}

transformed parameters {
  vector[N] theta; // Time-dependent thresholds
  vector[N] alpha; // Corresponding alpha values

  for (i in 1:N){
    theta[i] = fmin(0.4999, base_theta * growth_rate * known_t_to_target[i]);
    alpha[i] = inv_Phi(1 - theta[i]);
  }
}

model {
  // Priors
  k ~ normal(1, 0.5);
  sigma_0 ~ gamma(1, 1);
  base_theta ~ beta(2, 2); // Prior for baseline threshold
  growth_rate ~ normal(0.1, 0.05); // Prior for moderate growth rate
  sigma_err ~ gamma(1, 1);

  // Vectorized computation of predicted times
  vector[N] pred_clock_check_time = get_predicted_time(to_vector(known_t_to_target), sigma_0, k, to_vector(alpha));

  // Vectorized likelihood
  observed_time ~ normal(pred_clock_check_time, sigma_err);
}

// generated quantities {
//   vector[N] predicted_time;
//   predicted_time = get_predicted_time(to_vector(known_t_to_target), sigma_0, k, alpha);
// }
