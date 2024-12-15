functions {
  real get_predicted_time(real t_target, real sigma_0, real k, real theta) {
    real alpha = inv_Phi(1 - theta); // Inverse CDF of normal for the threshold
    real pred_clock_check_time = t_target / (1 + alpha * sigma_0 * k); // Predicted time for the clock check
    return(pred_clock_check_time);
  }
}

data {
  int<lower=0> N; // Number of checking events
  int<lower=1> P; // Number of persons
  int<lower=1, upper=P> person_id[N]; // Person ID for each observation
  real<lower=0> known_t_to_target[N]; // Time to target
  real<lower=0> observed_time[N]; // Observed time of the actual clock check
}

parameters {
  // Group-level parameters
  real<lower=0> mu_k; // Group mean for k
  real<lower=0> mu_sigma_0; // Group mean for sigma_0
  real<lower=0, upper=1> mu_theta; // Group mean for theta
  
  real<lower=0> sigma_k; // Standard deviation for k
  real<lower=0> sigma_sigma_0; // Standard deviation for sigma_0
  real<lower=0> sigma_theta; // Standard deviation for theta
  real<lower=0> sigma_err; // Error of prediction
  
  // Person-level parameters
  vector<lower=0>[P] k; // Person-specific k
  vector<lower=0>[P] sigma_0; // Person-specific sigma_0
  vector<lower=0, upper=1>[P] theta; // Person-specific theta
}

model {
  // Priors for group-level parameters
  mu_k ~ normal(1, 0.5);
  mu_sigma_0 ~ gamma(1, 1);
  mu_theta ~ beta(1, 1);
  
  sigma_k ~ gamma(1, 1);
  sigma_sigma_0 ~ gamma(1, 1);
  sigma_theta ~ uniform(0, 1);
  
  sigma_err ~ gamma(1, 1);

  // Priors for person-level parameters
  k ~ normal(mu_k, sigma_k);
  sigma_0 ~ normal(mu_sigma_0, sigma_sigma_0);
  theta ~ normal(mu_theta, sigma_theta);
  
  // Likelihood for observations
  for (i in 1:N) {
    real pred_clock_check_time = get_predicted_time(
      known_t_to_target[i], 
      sigma_0[person_id[i]], 
      k[person_id[i]], 
      theta[person_id[i]]
    );
    
    observed_time[i] ~ normal(pred_clock_check_time, sigma_err);
  }
}
