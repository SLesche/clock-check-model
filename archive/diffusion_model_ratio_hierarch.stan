functions {
  real get_predicted_ratio(real sigma_0, real k, real alpha) {
    return 1 ./ (1 + alpha * sigma_0 * k); // Predicted mean ratio
  }
}

data {
  int<lower=0> N; // Number of checking events
  int<lower=1> P; // Number of subjects
  int<lower=1> subject_id[N];
  real<lower=0, upper=1> observed_ratio[N]; // The ratio of actual clock check corresponding to target time
}

parameters {
  real<lower=0> mu_k; // Noise scalar k
  real<lower=0> sigma_k;
  real<lower=0> mu_sigma_0; // Initial noise at t = 0
  real<lower=0> sigma_sigma_0;
  real<lower=0, upper=1> mu_raw_theta; // Threshold (raw)
  real<lower=0, upper=1> sigma_raw_theta;
  real<lower=0, upper=1> mu_sigma_err; // Precision parameter for Beta distribution
  real<lower=0, upper=1> sigma_sigma_err;
  
  real<lower=0> k[P];
  real<lower=0> sigma_0[P];
  real<lower=0, upper=1> raw_theta[P];
  real<lower=0, upper=1> sigma_err[P];
}

transformed parameters {
  real<lower=0, upper=0.5> theta[P];
  real alpha[P];
  
  for (p in 1:P){
    theta[p] = raw_theta[p] / 2;
    alpha[p] = inv_Phi(1 - theta[p]);
  }
}

model {
  // Priors
  mu_k ~ normal(1, 0.5);
  sigma_k ~ gamma(1, 1);
  mu_sigma_0 ~ gamma(1, 1);
  sigma_sigma_0 ~ gamma(1, 1);
  mu_raw_theta ~ beta(2, 3);
  sigma_raw_theta ~ beta(2, 10);
  mu_sigma_err ~ beta(2, 6); 
  sigma_sigma_err ~ beta(2, 10);
  
  k ~ normal(mu_k, sigma_k);
  sigma_0 ~ normal(mu_sigma_0, sigma_sigma_0);
  raw_theta ~ normal(mu_raw_theta, sigma_raw_theta);
  sigma_err ~ normal(mu_sigma_err, sigma_sigma_err);

  for (i in 1:N) {
    real pred_check_ratio = get_predicted_ratio(sigma_0[subject_id[i]], k[subject_id[i]], alpha[subject_id[i]]);
    
    // // Compute alpha and beta parameters for Beta distribution
    // real alpha_beta = pred_check_ratio^2 * ((1 - pred_check_ratio) / sigma_err[subject_id[i]]^2 - 1 / pred_check_ratio);
    // real beta_beta = alpha_beta * (1 / pred_check_ratio - 1);
    // // Apply Beta likelihood for all observations
    // observed_ratio[i] ~ beta(alpha_beta, beta_beta);
    
    observed_ratio[i] ~ normal(pred_check_ratio, sigma_err);
  }
}
