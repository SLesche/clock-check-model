functions {
  real mixture_weibull_lpdf(real[] t, real pi, real sigma1, real alpha2, real sigma2) {
    real log_prob = 0;
    for (i in 1:num_elements(t)) {
      // Density of Weibull component 1 (constant hazard: alpha1 = 1)
      real f1 = exp(weibull_lpdf(t[i] | 1, sigma1));
      
      // Density of Weibull component 2 (flexible hazard: alpha2, sigma2 free)
      real f2 = exp(weibull_lpdf(t[i] | alpha2, sigma2));
      
      // Mixture probability
      log_prob += log(pi * f1 + (1 - pi) * f2);
    }
    return log_prob;
  }
}

data {
  int<lower=1> N;                     // Number of observations
  real<lower=0> clock_check_time[N];  // Observation times
  real<lower=0> known_t_to_target[N]; // Predictor variable
}

parameters {
  // Regression coefficients for mixture weight
  real beta_pi_0;
  real beta_pi_1;

  // Regression coefficients for sigma1 (scale of Weibull 1)
  real beta_sigma1_0;
  real beta_sigma1_1;

  // Regression coefficients for alpha2 (shape of Weibull 2)
  real beta_alpha2_0;
  real beta_alpha2_1;

  // Regression coefficients for sigma2 (scale of Weibull 2)
  real beta_sigma2_0;
  real beta_sigma2_1;
}

model {
  // Priors for regression coefficients
  beta_pi_0 ~ normal(0, 1);
  beta_pi_1 ~ normal(0, 1);
  
  beta_sigma1_0 ~ normal(0, 1);
  beta_sigma1_1 ~ normal(0, 1);
  
  beta_alpha2_0 ~ normal(0, 1);
  beta_alpha2_1 ~ normal(0, 1);
  
  beta_sigma2_0 ~ normal(0, 1);
  beta_sigma2_1 ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    // Compute parameter values as functions of known_t_to_target[n]
    real pi = inv_logit(beta_pi_0 + beta_pi_1 * known_t_to_target[n]);
    real sigma1 = exp(beta_sigma1_0 + beta_sigma1_1 * known_t_to_target[n]);
    real alpha2 = exp(beta_alpha2_0 + beta_alpha2_1 * known_t_to_target[n]);
    real sigma2 = exp(beta_sigma2_0 + beta_sigma2_1 * known_t_to_target[n]);
    
    // Add to the target log-probability
    target += mixture_weibull_lpdf({clock_check_time[n]} | pi, sigma1, alpha2, sigma2);
  }
}
