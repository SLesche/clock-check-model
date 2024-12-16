functions {
  real get_predicted_time(real t_target, real sigma_0, real k, real alpha) {
    return t_target ./ (1 + alpha * sigma_0 * k); // Vectorized element-wise division
  }
}

data {
  int<lower=0> N; // Number of checking events
  real<lower=0, upper = 1> known_t_to_target[N]; // The time of target
  real<lower=0, upper = 1> observed_time[N]; // The time of the actual clock check
}

parameters {
  real<lower=0> k; // Noise scalar k
  real<lower=0> sigma_0; // Initial noise at t = 0
  real<lower=0, upper=1> raw_theta_initial; // Baseline threshold (capped at 0.5)
  real<lower=-1, upper=1> raw_theta; // Growth rate for threshold
  real<lower=0> sigma_err; // Error level
}

transformed parameters {
  real<lower=0, upper=0.5> theta_initial;
  theta_initial = raw_theta_initial/2;
  
  real<lower=0, upper=0.5> theta;
  theta = raw_theta / 2;
}

model {
  // Priors
  k ~ normal(1, 0.5);
  sigma_0 ~ gamma(1, 1);
  raw_theta ~ beta(1, 1); // Prior for baseline threshold
  theta_initial ~  beta(1, 1); // Prior for moderate growth rate
  sigma_err ~ gamma(1, 1);
  

  // Vectorized computation of predicted times
  for (i in 1:N){
    if (known_t_to_target[i] == 1){
      real pred_clock_check_time = get_predicted_time(known_t_to_target[i], sigma_0, k, theta_initial);
      observed_time[i] ~ normal(pred_clock_check_time, sigma_err);
    } else {
      real pred_clock_check_time = get_predicted_time(known_t_to_target[i], sigma_0, k, theta);
      observed_time[i] ~ normal(pred_clock_check_time, sigma_err);
    }
  }
}

// generated quantities {
//   vector[N] predicted_time;
//   predicted_time = get_predicted_time(to_vector(known_t_to_target), sigma_0, k, alpha);
// }
