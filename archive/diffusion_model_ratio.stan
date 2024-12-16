functions {
  real get_predicted_ratio(real sigma_0, real k, real alpha) {
    return 1 ./ (1 + alpha * sigma_0 * k); // Predicted mean ratio
  }
}

data {
  int<lower=0> N; // Number of checking events
  real<lower=0, upper=1> observed_ratio[N]; // The ratio of actual clock check corresponding to target time
}

parameters {
  real<lower=0> k; // Noise scalar k
  real<lower=0> sigma_0; // Initial noise at t = 0
  real<lower=0, upper=1> raw_theta; // Threshold (raw)
  real<lower=0, upper=1> sigma_err; // Precision parameter for Beta distribution
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
  sigma_err ~ uniform(0, 1); // Precision prior for Beta distribution (tweak as needed)

  // Likelihood
  real pred_check_ratio = get_predicted_ratio(sigma_0, k, alpha);
  real alpha_beta = pred_check_ratio^2 * ((1 - pred_check_ratio) / sigma_err^2 - 1 / pred_check_ratio);
  real beta_beta = alpha_beta * (1 / pred_check_ratio - 1);
  
  observed_ratio ~ beta(alpha_beta, beta_beta);
}
