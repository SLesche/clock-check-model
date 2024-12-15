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

model {
  // Priors
  k ~ normal(1, 0.5);
  sigma_0 ~ gamma(1, 1);
  theta ~ beta(1,1);
  
   for (i in 1:N) {
    real t = known_t_to_target[i];
    real pred_clock_check_time = 0;
    
    real t_step = 1;
    real t_step_size = 1;
    
    // Simulate the accumulation process
    while (t_step < t) {
      real noise = sigma_0 * k * t_step;  // Noise increases with time
      
      real prob_included = 1 - normal_cdf(known_t_to_target[i], t_step, noise);
      
      // Check if threshold is crossed
      if (prob_included >= theta) {
        pred_clock_check_time = t_step;
        break;
      }
      t_step = t_step + t_step_size;
    }
    
    // Likelihood of the observed response time
    observed_time[i] ~ normal(pred_clock_check_time, sigma_0);
  }
}
