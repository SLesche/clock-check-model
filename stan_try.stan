data {
  int<lower=0> N;                  // Number of clock-checking events
  int<lower=1> Nsubj;                  // Number of participants
  int<lower=1, upper=Nsubj> id[N];     // Participant IDs for each event
  real t_since_last_check[N];          // Time since last check
  real known_t_to_target[N];              // Target times (e.g., 5 minutes for all)
  int<lower=0, upper=1> check[N];                 // Binary (1 = c heck, 0 = no check)
  // real delta;                      // Margin of error for target time
}

parameters {
   // Group-level parameters
  real<lower=0> mu_sigma_0;       // Mean baseline noise
  real<lower=0> sigma_sigma_0;    // SD for baseline noise

  real<lower=0> mu_k;             // Mean drift rate for noise growth
  real<lower=0> sigma_k;          // SD for drift rate

  // Participant-specific parameters
  real<lower=0> sigma_0[Nsubj];     // Baseline noise
  real<lower=0> k[Nsubj];           // Drift rate for noise growth
}

transformed parameters {
  real<lower=0> sigma_est[N];                // Uncertainty over time
  real<lower=0, upper=1> prob_check[N];          // Probability target is within uncertainty range
  

  for (itrial in 1:N) {
    // Estimation noise grows exponentially with time since last check
    sigma_est[itrial] = sigma_0[id[itrial]] * exp(k[id[itrial]] * t_since_last_check[itrial]);
    
    // Probability that the target time is within the estimation noise range
    //prob_check[itrial] = 1 - normal_cdf(known_t_to_target[itrial] - delta, t_since_last_check[itrial], sigma_est[itrial]) + 
    //                       normal_cdf(known_t_to_target[itrial] + delta, t_since_last_check[itrial], sigma_est[itrial]);
    prob_check[itrial] = 1 - normal_cdf(known_t_to_target[itrial], t_since_last_check[itrial], sigma_est[itrial]);
  }

}

model {
  // Priors for shared parameters
  mu_sigma_0 ~ gamma(1, 1);
  sigma_sigma_0 ~ gamma(1, 1);
  
  mu_k ~ gamma(1, 1);
  sigma_k ~ gamma(1, 1);
  
  // Priors for participant-specific parameters
  sigma_0 ~ normal(mu_sigma_0, sigma_sigma_0);
  k ~ normal(mu_k, sigma_k);
  
  // Likelihood of clock-checking behavior
  for (itrial in 1:N) {
    check[itrial] ~ bernoulli(prob_check[itrial]);
  }
}

generated quantities {
  int<lower=0, upper=1> y_rep[N];  // Simulated clock checks

  // Simulate clock-checking behavior based on the predicted probabilities
  for (itrial in 1:N) {
    y_rep[itrial] = bernoulli_rng(prob_check[itrial]);  // Correct use of bernoulli_rng
  }
}
