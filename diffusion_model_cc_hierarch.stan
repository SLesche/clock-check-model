data {
  int<lower=0> N; // number of checking events
  int<lower=1> Nsubj; // number of subjects
  int<lower=1, upper=Nsubj> id[N]; // Subject ID for each trial
  real<lower=0> known_t_to_target[N]; // the time of target
  real<lower=0> observed_time[N]; // the time of the actual clock check
}

parameters {
  real<lower=0> k_mu; // Mean noise scalar k (population level)
  real<lower=0> sigma_k; // Std dev of noise scalar k (population level)
  real<lower=0> sigma_0_mu; // Mean initial noise sigma_0 (population level)
  real<lower=0> sigma_sigma_0; // Std dev of initial noise sigma_0 (population level)
  real<lower=0, upper=1> theta_mu; // Mean threshold (population level)
  real<lower=0, upper=1> sigma_theta; // Std dev of threshold (population level)
  
  real<lower=0> k[Nsubj]; // Subject-specific noise scalar k
  real<lower=0> sigma_0[Nsubj]; // Subject-specific initial noise sigma_0
  real<lower=0, upper=1> theta[Nsubj]; // Subject-specific threshold
  
}

model {
  // Priors for population-level parameters
  k_mu ~ normal(1, 0.5);
  sigma_k ~ normal(0, 1);
  sigma_0_mu ~ gamma(1, 1);
  sigma_sigma_0 ~ gamma(1, 1);
  theta_mu ~ beta(1, 1);
  sigma_theta ~ gamma(1, 1);
  
  // Priors for subject-level parameters
  for (id in 1:Nsubj) {
    k[id] ~ normal(k_mu, sigma_k);
    sigma_0[id] ~ normal(sigma_0_mu, sigma_sigma_0);
    theta[id] ~ normal(theta_mu, sigma_theta);
  }

  // Likelihood for each subject's trials
  for (i in 1:N) {
    int sub_id = id[i];  // get the subject ID for trial i
    real t = known_t_to_target[i];
    real pred_clock_check_time = 0;
    
    real t_step = 1;
    real t_step_size = 1;
    
    // Simulate the accumulation process for each subject
    while (t_step < t) {
      real noise = sigma_0[sub_id] * k[sub_id] * t_step;  // Noise increases with time
      
      real prob_included = 1 - normal_cdf(known_t_to_target[i], t_step, noise);
      
      // Check if threshold is crossed
      if (prob_included >= theta[sub_id]) {
        pred_clock_check_time = t_step;
        break;
      }
      t_step = t_step + t_step_size;
    }
    
    // Likelihood of the observed response time for each subject
    observed_time[i] ~ normal(pred_clock_check_time, sigma_0[sub_id]);
  }
}
