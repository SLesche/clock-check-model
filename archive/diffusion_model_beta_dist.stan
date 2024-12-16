functions {
  real get_predicted_time(real t_target, real sigma_0, real k, real theta){
    real alpha = inv_Phi(1 - theta); // Inverse CDF of normal for the threshold
    real pred_clock_check_time = t_target / (1 + alpha * sigma_0 * k); // Predicted time for the clock check
  }
}
data {
  int<lower=0> N; // number of checking events
  real<lower=0> known_t_to_target[N]; // the time of target
  real<lower=0> ratio[N]; // the time of the actual clock checke relative to target_time
}

parameters {
  real<lower=0> k;            // Noise scalar k
  real<lower=0> sigma_0;      // Initial noise at t=0
  real<lower=0, upper=1> theta;  // Threshold
  real<lower=0, upper=1> mu_ratio;  // Mean of the ratio (0 < mu_ratio < 1)
  real<lower=0, upper=1> sigma_ratio; // Standard deviation of the ratio
}

transformed parameters {
  real alpha_ratio = mu_ratio^2 * (((1 - mu_ratio) / sigma_ratio^2) - 1 / mu_ratio);
  real beta_ratio = alpha_ratio * ((1 / mu_ratio) - 1);

}

model {
  // Priors
  k ~ normal(1, 0.5);
  sigma_0 ~ gamma(1, 1);
  theta ~ beta(1, 1);
  
  // Priors for mean and standard deviation of the ratio
  mu_ratio ~ beta(1, 1); // Beta distribution for the mean of the ratio
  sigma_ratio ~ uniform(0, 1);  // Standard deviation must be between 0 and 1


  for (i in 1:N) {
    // Likelihood: Beta distribution for the ratio
    ratio[i] ~ beta(alpha_ratio, beta_ratio); // Beta distribution for the ratio between 0 and 1
  }
}
