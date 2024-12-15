data {
  int<lower=0> N; // number of checking events
  real<lower=0> known_t_to_target[N]; // the time of target
  real<lower=0> observed_time[N]; // the time of the actual clock check
}

parameters {
  real<lower=0> k; // Noise scalar k
  real<lower=0> sigma_0; // Initial noise at t=0
  real<lower=0, upper=1> theta; // threshold
}

transformed parameters {
  real alpha = inv_Phi(1 - theta); // Inverse CDF of normal for the threshold
}

model {
  // Priors
  k ~ normal(1, 0.5);
  sigma_0 ~ gamma(1, 1);
  theta ~ beta(1,1);
  
  // Likelihood of the observed response times
  for (i in 1:N) {
    real t = known_t_to_target[i];
    
    // Compute the predicted clock check time using the derived formula
    real pred_clock_check_time = t / (1 + alpha * sigma_0 * k); // Formula for t_step
    
    // Likelihood of the observed response time
    observed_time[i] ~ normal(pred_clock_check_time, sigma_0);
  }
}

