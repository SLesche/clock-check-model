data {
  int<lower=0> N;                  // Number of clock-checking events
  int<lower=1> Nsubj;                  // Number of participants
  int<lower=1,upper=Nsubj> id[N];     // Participant IDs for each event
  vector[N] t_since_last_check;          // Time since last check
  vector[N] t_to_target;              // Target times (e.g., 5 minutes for all)
  vector[N] check;                 // Binary (1 = c heck, 0 = no check)
  real delta;                      // Margin of error for target time
}

parameters {
  real<lower=0> sigma_0[Nsubj];        // Participant-specific baseline uncertainty
  real<lower=0> k;                // Shared drift rate
  real<lower=0> alpha;            // Shared initial tolerance
  real<lower=0> beta;             // Shared sensitivity to target proximity
}

transformed parameters {
  vector[N] sigma;                // Uncertainty over time
  vector[N] prob_within;          // Probability target is within uncertainty range
  vector[N] threshold;            // Threshold probability for clock checks
  
  for (n in 1:N) {
    int p = id[n];               // Participant ID for this event
    
    // Uncertainty grows with time since last check
    sigma[n] = sqrt(sigma_0[p]^2 + k * t_last_check[n]);
    
    // Probability target is within uncertainty range
    prob_within[n] = normal_cdf(t_target[n], 0, sigma[n]) - 
                     normal_cdf(t_target[n] - delta, 0, sigma[n]);
    
    // Threshold probability for clock check
    threshold[n] = alpha * exp(-beta * (t_target[n] - t_last_check[n]));
  }
}

model {
  // Priors for shared parameters
  k ~ normal(0, 1);
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  
  // Priors for participant-specific parameters
  sigma_0 ~ normal(0, 1);
  
  // Likelihood of clock-checking behavior
  for (n in 1:N) {
    check[n] ~ bernoulli(prob_within[n] < threshold[n]);
  }
}
