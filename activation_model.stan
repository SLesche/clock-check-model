functions {
  // Define a single function to compute everything in one pass
  vector prob_check_time_optimized(real[] t_steps, real t_target, real k, real g, real threshold, real a, real dstep) {
    int num_steps = num_elements(t_steps);
    vector[num_steps] check_time_probs;
    
    // Initialize the prob_target, prob_action, and cumulative prob
    vector[num_steps] prob_tgt;
    vector[num_steps] prob_act;
    vector[num_steps] cum_prob;

    // Compute prob_target and prob_action in one pass
    prob_tgt[1] = 1 - normal_cdf(t_target, t_steps[1], k * t_steps[1]);
    prob_act[1] = g + (1 - g) / (1 + exp(-a * (prob_tgt[1] - threshold)));
    cum_prob[1] = g; // Cumulative probability at first step
    
    // Loop through time steps to calculate probabilities and cumulative probability
    for (i in 2:num_steps) {
      prob_tgt[i] = 1 - normal_cdf(t_target, t_steps[i], k * t_steps[i]);
      prob_act[i] = g + (1 - g) / (1 + exp(-a * (prob_tgt[i] - threshold)));
      
      cum_prob[i] = cum_prob[i - 1] + prob_act[i - 1] * dstep; // Cumulative probability
    }

    // Compute final check-time probabilities
    check_time_probs[1] = prob_act[1];  // First step, no previous action
    for (i in 2:num_steps) {
      real prob_no_previous_action = 1 - cum_prob[i - 1];
      check_time_probs[i] = prob_act[i] * prob_no_previous_action;
    }

    return check_time_probs;
  }
}

data {
  int<lower=1> N;               // Number of observed data points
  real<lower=0> t_target;       // Target time
  int<lower=1> num_steps;      // Number of steps for checking
  real<lower=0> t_steps[N];     // Observed time steps
}

parameters {
  real<lower=0> k;              // Scaling factor for prob_target
  real<lower=0,upper=1> g;      // Base probability for prob_action
  real<lower=0,upper=1> threshold; // Threshold for prob_action
  real<lower=0> a;              // Sharpness of sigmoid in prob_action
}

transformed parameters {
  real prob_check_time_vals[num_steps];            // Probabilities at each step
  real <lower=0> dstep;
  
  dstep = t_target / num_steps

  // Initialize t_steps (for the time grid)
  for (i in 1:num_steps) {
    t_steps[i] = (i - 1) * dstep;
  }

  // Compute check-time probabilities for all steps
  prob_check_time_vals = prob_check_time_optimized(t_steps, t_target, k, g, threshold, a, dstep);
}

model {
  // Priors (optional, adjust based on domain knowledge)
  k ~ normal(0, 1);
  g ~ beta(2, 2);
  threshold ~ beta(2, 2);
  a ~ normal(0, 1);

  // Likelihood of observed data (for known check times)
  for (n in 1:N) {
    target += log(prob_check_time_vals[n]);
  }
}
