data {
  int<lower=1> N;                  // Total number of clock-checking events
  int<lower=1> Nsubj;              // Number of participants
  int<lower=1, upper=Nsubj> id[N]; // Participant IDs for each event
  real<lower=0> check_time[N];     // Times of clock checks
  real<lower=0> block_length[N];
}

parameters {
  real<lower=0> sigma_0[Nsubj];    // Baseline noise
  real<lower=0> k[Nsubj];          // Drift rate for noise growth
  real<lower=0> theta[Nsubj];      // Participant-specific hazard threshold
}

transformed parameters {
  real cumulative_hazard[N];       // Cumulative hazard at each event time
  real hazard_rate[N];             // Hazard rate for each time interval
  real survival_prob[N];           // Survival probability until each event

  // Initialize the cumulative hazard
  cumulative_hazard[1] = 0; // No hazard accumulation before the first event

  for (n in 2:N) {
    // Calculate the time interval since the last clock check
    real time_interval = check_time[n] - check_time[n - 1];

    // Estimation noise grows over time
    real sigma = sigma_0[id[n]] * exp(k[id[n]] * time_interval);

    // Hazard rate (instantaneous probability of a clock check at this time)
    hazard_rate[n] = 1 - normal_cdf(block_length[n], time_interval, sigma);

    // Update the cumulative hazard
    cumulative_hazard[n] = cumulative_hazard[n - 1] + hazard_rate[n] * time_interval;

    // Reset the cumulative hazard to zero if the threshold is exceeded
    if (cumulative_hazard[n] > theta[id[n]]) {
      cumulative_hazard[n] = 0;
    }

    // Calculate survival probability (probability of no event until this time)
    survival_prob[n] = exp(-cumulative_hazard[n]);
  }
}

model {
  // Priors for participant-specific parameters
  sigma_0 ~ normal(1, 1);       // Prior for baseline noise
  k ~ normal(1, 1);             // Prior for drift rate
  theta ~ normal(5, 2);         // Prior for hazard threshold

  // Likelihood: Log-likelihood for clock-check times
  for (n in 2:N) {
    target += log(hazard_rate[n]) + log(survival_prob[n - 1]);
  }
}
