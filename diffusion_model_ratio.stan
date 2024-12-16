functions {
  real get_predicted_ratio(real sigma_0, real k, real alpha) {
    return 1 ./ (1 + alpha * sigma_0 * k); // Vectorized element-wise division
  }
}

data {
  int<lower=0> N; // Number of checking events
  // real<lower=0> known_t_to_target[N]; // The time of target
  real<lower=0, upper=1> observed_ratio[N]; // The ratio of the actual clock check corresponding to target time
}

parameters {
  real<lower=0> k; // Noise scalar k
  real<lower=0> sigma_0; // Initial noise at t = 0
  real<lower=0, upper=1> raw_theta; // Threshold (raw)
  real<lower=0> sigma_err; // Error level
}

transformed parameters {
  real<lower=0, upper=0.5> theta;
  real alpha;
  theta = raw_theta / 2; // Scale the raw threshold to [0, 0.5]
  alpha = inv_Phi(1 - theta); // Precompute alpha
}

model {
  // Priors
  k ~ normal(1, 0.5);
  sigma_0 ~ gamma(1, 1);
  raw_theta ~ beta(1, 1);
  sigma_err ~ gamma(1, 1);

  // Vectorized computation of predicted times
  real pred_check_ratio = get_predicted_ratio(sigma_0, k, alpha);

  // Vectorized likelihood
  observed_ratio ~ normal(pred_check_ratio, sigma_err);
}

// generated quantities {
//   vector[N] predicted_time;
//   predicted_time = get_predicted_time(to_vector(known_t_to_target), sigma_0, k, alpha);
// }
